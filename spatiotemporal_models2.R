
## libraries
library(sdmTMB)
#library(sdmTMBextra)
library(fmesher)
library(INLA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(sf)
library(future)
library(furrr)


########################
##
## Visualise survey data
##
########################

## Processed Survey Data
full_data <- readRDS("data/Allspecies_15_24_Datras.Rds")

## Add UTM cols
full_data <- add_utm_columns(full_data, ll_names = c("lon", "lat"), 
                             ll_crs = 4326, units = "km")

# plot survey data
ggplot(full_data[full_data$biomass_g > 0,], aes(X, Y, col = log(biomass_g))) +
  geom_point() +
  coord_fixed() +
  scale_colour_viridis_c() +
  facet_wrap(~ cmn_nm) +
  labs(colour = "log(Biomass (g))", title = "Observed survey tows")


#############
## Covariates
#############

# log transform
full_data$logHaulDur    <- log(full_data$HaulDur)
full_data$logBathymetry <- log(full_data$Bathymetry)

covariates <- c("logBathymetry", "dist2coast",  
                "thetao_surface.rcp45", "current_surface.rcp45",
                "current_bottom.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45")
## rescale
full_data[,covariates] <- apply(full_data[,covariates],2,function(x){(x-mean(x))/sd(x)}) 



##############
## Create Mesh
##############

data    <- split(full_data, full_data$cmn_nm)
sp_mesh <- NULL

for (i in 1:length(data)) {
  
  print(data[[i]]$cmn_nm[1])
  
  ## define mesh ##
  non_convex_bdry <- INLA::inla.nonconvex.hull(cbind(data[[i]]$X, data[[i]]$Y), 
                                               -0.03, resolution = c(200, 50))
  m <- fmesher::fm_mesh_2d_inla(
    loc = cbind(data[[i]]$X, data[[i]]$Y), # coordinates4
    max.edge = c(120), # max triangle edge length; inner and outer meshes
    cutoff = 60, # minimum triangle edge length
    boundary = non_convex_bdry
  )
  m <- make_mesh(data[[i]], c("X", "Y"), mesh = m)
  
  ## save model
  nm <- paste0(data[[i]]$cmn_nm[1])
  sp_mesh[[nm]] <- m
}

# order meshes by species data
sp_mesh <- sp_mesh[names(data)]


################
## Model Fitting
################


## fit modeldata## fit models for each species
SDM_fitting <- function(dat = dat, mesh = mesh) {
  
  start_time <- Sys.time()
  sp <- dat$cmn_nm[1]
  all_models <- list()
  
  print(sp)
  
  ## full models ##
  
  # delta lognormal
  all_models$dln_full <- sdmTMB(biomass_g ~  logBathymetry + Quarter + thetao_surface.rcp45 + current_surface.rcp45 +
                       current_bottom.rcp45 + ph_surface.rcp45 + ph_bottom.rcp45 +      
                       totoc_bottom.rcp45 + zooc_surface.rcp45 + Gravel.percentage +
                       dist2coast,
                     data = dat,
                     mesh = mesh,
                     family = delta_lognormal(),
                     spatiotemporal = "rw",
                     offset = "logHaulDur",
                     time = "Year",
                     spatial = "on"
  )
  
  # delta gamma 
  all_models$dg_full <- update(all_models$dln_full, family = delta_gamma()) 
  
  # tweedie 
  all_models$tw_full <- update(all_models$dln_full, family = tweedie()) 
  
  
  ## reduced models ## 
  
  for (m in 1:length(all_models)) {
    
    # check if hurdle model
    hurdle_mod <- length(all_models[[m]]$formula) == 2
    
    # model parameters
    mod1 <- tidy(all_models[[m]], model=1) %>%
      mutate(crosses_zero = (conf.low <= 0 & conf.high >= 0))
    
    if (hurdle_mod) {
      mod2 <- tidy(all_models[[m]], model=2) %>%
        mutate(crosses_zero = (conf.low <= 0 & conf.high >= 0))
    } else { 
      mod2 <- NULL 
    }
    
    i <- 0
    
    while (any(mod1$crosses_zero[-1] == TRUE) ||
           hurdle_mod && any(mod2$crosses_zero[-1] == TRUE)) {
      
      i <- i + 1
      
      # update model formula
      covs1 <- mod1[mod1$crosses_zero == FALSE,]$term[-1]
      f     <- paste("biomass_g ~", paste(covs1, collapse = " + "))
      
      ## hurdle models ##
      
      if (hurdle_mod) {
        # update formula    
        covs2    <- mod2[mod2$crosses_zero == FALSE,]$term[-1]
        f2 <- paste("biomass_g ~", paste(covs2, collapse = " + "))
        f <- list(f, f2)
        f <- gsub("Quarter3", "Quarter", f)
        f <- lapply(f, as.formula)
        
        # re-run reduced model
        mod_red <- sdmTMB(formula = f,
                          data = dat,
                          mesh = mesh,
                          family = all_models[[m]]$family,
                          spatiotemporal = "rw",
                          offset = "logHaulDur",
                          time = "Year",
                          spatial = "on"
        )
      } else {  
        
        ## tweedie models ##
        
        f <- gsub("Quarter3", "Quarter", f)
        f <- as.formula(f)
        
        # re-run reduced model
        mod_red <- update(all_models[[m]], formula = f)
      }
      
      ## save model
      nm <- paste0(gsub("full", "red_", names(all_models)[m]), i)
      all_models[[nm]] <- mod_red
      print(nm)
      
      # update parameter estimates
      mod1 <- tidy(mod_red, model=1) %>%
        mutate(crosses_zero = (conf.low <= 0 & conf.high >= 0))
      
      if (hurdle_mod) {
        mod2 <- tidy(mod_red, model=2) %>%
          mutate(crosses_zero = (conf.low <= 0 & conf.high >= 0))
      }
    }
  }
  
  # clear memory
  #rm(mod_red); rm(mod1); rm(mod2); rm(covs1); rm(covs2); rm(hurdle_mod)
  
  # model checks
  cAIC   <- sapply(all_models, cAIC)
  sanity <- sapply(all_models, sanity)
  
  # best models per family
  best_models <- sapply( c("tw", "dln", "dg"), function(x) {
      y <- cAIC[grep(x, names(cAIC))]
      names(y)[which.min(y)]
    }
  )
  
  names(best_models) <- NULL
  
  
  ## spatially-varying and random intercept models ## 
  
  for (i in 1:length(best_models)) {
    
    # check if hurdle model
    mod <- all_models[[best_models[i]]]
    hurdle_mod <- length(mod$formula) == 2
    
    ## spatially-varying Quarter
    f <- gsub("+ Quarter", "", mod$formula, fixed = TRUE)
    if (hurdle_mod) { f <- lapply(f, as.formula) } else { f <- as.formula(f) }
    
    # re-run model
    mod_svQ <- sdmTMB(formula = f,
                      data = dat,
                      mesh = mesh,
                      family = mod$family,
                      spatiotemporal = "rw",
                      offset = "logHaulDur",
                      time = "Year",
                      spatial_varying = ~ Quarter,
                      spatial = "on"
    )
    
    # save model to list
    nm <- paste0(best_models[i], "_svQ")
    all_models[[nm]] <- mod_svQ
    print(nm)
    
    ## Ship random intercept
    if (hurdle_mod) { 
      f <- gsub("biomass_g ~", "biomass_g ~ (1 | Ship) +", f, fixed = TRUE)
      f <- lapply(f, as.formula) 
    } else { f <- update(f, . ~ . + (1 | Ship)) }
    
    # re-run model
    mod_svQ_Sri <- sdmTMB(formula = f,
                          data = dat,
                          mesh = mesh,
                          family = mod$family,
                          spatiotemporal = "rw",
                          offset = "logHaulDur",
                          time = "Year",
                          spatial_varying = ~ Quarter,
                          spatial = "on"
    )
    
    # save model to list
    nm <- paste0(best_models[i], "_svQ_Sri")
    all_models[[nm]] <- mod_svQ_Sri
    print(nm)
    
  }
  
  #rm(mod_svQ); rm(mod_svQ_Sri); rm(nm); rm(f); rm(hurdle_mod)
  
  
  # model results - THIS WOULD NEED CHANGING TO RUN FOR ALL SPECIES
  model    <- names(all_models)
  cAIC     <- sapply(all_models, cAIC)
  sanity   <- sapply(all_models, sanity)
  sanity   <- t(as.data.frame(sanity)[9,])
  formulas <- lapply(all_models, function(x) x[["formula"]])
  
  
  model_results <- data.frame(model    = rep(model, each = 2)[-c(2,4,6,12,14,16)],#add 18 for plaice
                              cAIC     = rep(cAIC, each = 2)[-c(2,4,6,12,14,16)],
                              sanity   = rep(unlist(sanity), each = 2)[-c(2,4,6,12,14,16)],
                              formulas = as.character(unlist(formulas))[-c(2,4)])
  model_results$sp <- sp
  
  # Create file name: file_1.rds, file_2.rds, ...
  file_name1 <- paste0("res/", paste0(gsub(" ", "", sp), "_models.RData"))
  file_name2 <- paste0("res/", paste0(gsub(" ", "", sp), "_model_results.RData"))
  
  # save models and results
  save(all_models, file = file_name1)
  save(model_results, file = file_name2)
  
  end_time <- Sys.time()
  
  return(list(all_models = all_models,
              model_results = model_results,
              time = end_time - start_time))
}



# Set parallel backend
plan(multisession, workers = 2)
  
# Run in parallel
start <- Sys.time()
results <- future_map2(data[c(2)], sp_mesh[c(2)], SDM_fitting) # removed n.pout
end <- Sys.time()
end - start

  
  

#RhpcBLASctl::blas_set_num_threads(1)

