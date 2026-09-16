setwd("C:/Users/Public/010 Multispecies Modelling")

## libraries
library(sdmTMB)
library(fmesher)
library(INLA)
library(dplyr)
library(ggplot2)
library(tidyr)
library(sf)
library(visreg)


########################
##
## Visualise survey data
##
########################

## Processed Survey Data
rdat <- readRDS("data/Cod_15_24_Datras.Rds")

## Add UTM cols
rdat <- add_utm_columns(rdat, ll_names = c("lon", "lat"), 
                       ll_crs = 4326, units = "km") # 32631

uk <- map_data(map = "world", region = "UK")

# plot survey data
ggplot() +
  geom_point(data = rdat, aes(lon, lat, col = biomass_g)) +
  geom_polygon(data = uk, aes(x = long, y = lat, group = group)) +
  coord_fixed() +
  scale_colour_viridis_c(trans = "sqrt") +
  facet_wrap(~ Year + Quarter) +
  labs(colour = "Biomass (kg)", title = "Observed survey tows")



#############
## Covariates
#############

# log haul duration
rdat$logHaulDur <- log(rdat$HaulDur)

rdat$logBathymetry <- log(rdat$Bathymetry)

dat <- rdat

covariates <- c("Bathymetry", "dist2coast",  
                "thetao_surface.rcp45", "current_surface.rcp45",
                "current_bottom.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45","logBathymetry")
## rescale
dat[,covariates] <- apply(dat[,covariates],2,function(x){(x-mean(x))/sd(x)}) 


##############
## Define Mesh
##############

# define a mesh directly with fmesher (formerly INLA):
non_convex_bdry <- INLA::inla.nonconvex.hull(cbind(dat$X, dat$Y), 
                                             -0.025, resolution = c(200, 50))
inla_mesh <- fmesher::fm_mesh_2d_inla(
  loc = cbind(dat$X, dat$Y), # coordinates
  max.edge = c(150), # max triangle edge length; inner and outer meshes
  offset = c(0),  # inner and outer border width
  cutoff = 60, # minimum triangle edge length
  boundary = non_convex_bdry
)
mesh <- make_mesh(dat, c("X", "Y"), mesh = inla_mesh)

plot(mesh)
mesh$mesh$n # 133




################
## Model Fitting
################

## full model ##
dln_cod_full <- sdmTMB(biomass_g ~  logBathymetry + Quarter + thetao_surface.rcp45 + current_surface.rcp45 +
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

sanity(dln_cod_full)
tidy(dln_cod_full, model=1)
tidy(dln_cod_full, model=2)

# store results
model   <- "dln_cod_full"
cAIC    <- cAIC(dln_cod_full)
formula <- as.character(dln_cod_full$formula[1])



## reduced model ##
dln_cod_red <- update(dln_cod_full, formula = biomass_g ~ logBathymetry + Quarter + 
                        thetao_surface.rcp45 + current_surface.rcp45 +
                        ph_surface.rcp45 +    ph_bottom.rcp45 + Gravel.percentage +
                        totoc_bottom.rcp45 + zooc_surface.rcp45)

sanity(dln_cod_red)
tidy(dln_cod_red, model=1)
tidy(dln_cod_red, model=2)

# store results
model   <- c(model, "dln_cod_red")
cAIC    <- c(cAIC, cAIC(dln_cod_red))
formula <- c(formula, as.character(dln_cod_red$formula[1]))



## reduced model 1 ##
dln_cod_red1 <- sdmTMB(list(biomass_g ~  logBathymetry + Quarter + thetao_surface.rcp45 + current_surface.rcp45 +
                                 ph_surface.rcp45 + ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + Gravel.percentage,
                            biomass_g ~  logBathymetry + Quarter),
                          data = dat,
                          mesh = mesh,
                          family = delta_lognormal(),
                          spatiotemporal = "rw",
                          offset = "logHaulDur",
                          time = "Year",
                          spatial = "on"
)

sanity(dln_cod_red1)
tidy(dln_cod_red1, model=1)
tidy(dln_cod_red1, model=2)

# store results
model   <- c(model, rep("dln_cod_red1", 2))
cAIC    <- c(cAIC, rep(cAIC(dln_cod_red1), 2))
formula <- c(formula, as.character(dln_cod_red1$formula))





## Tweedie Models ##################


## full model ##
tw_cod_full <- update(dln_cod_full, family = tweedie())

sanity(tw_cod_full)
tidy(tw_cod_full)

# store results
model   <- c(model, "tw_cod_full")
cAIC    <- c(cAIC, cAIC(tw_cod_full))
formula <- c(formula, as.character(tw_cod_full$formula))


## reduced model ##
tw_cod_red <- update(tw_cod_full, formula = biomass_g ~ logBathymetry + Quarter + thetao_surface.rcp45 + current_surface.rcp45 +
                       ph_surface.rcp45 + ph_bottom.rcp45 + totoc_bottom.rcp45)

sanity(tw_cod_red)
tidy(tw_cod_red)

# store results
model   <- c(model, "tw_cod_red")
cAIC    <- c(cAIC, cAIC(tw_cod_red))
formula <- c(formula, as.character(tw_cod_red$formula))



## Delta-Gamma ##################

## full model ##
dg_cod_full <- update(dln_cod_full, family = delta_gamma())

sanity(dg_cod_full)
tidy(dg_cod_full, model=1)
tidy(dg_cod_full, model=2)

# store results
model   <- c(model, rep("dg_cod_full", 2))
cAIC    <- c(cAIC, rep(cAIC(dg_cod_full), 2))
formula <- c(formula, as.character(dg_cod_full$formula))



## reduced model ##
dg_cod_red <- sdmTMB(list(biomass_g ~ logBathymetry + Quarter + thetao_surface.rcp45 + current_surface.rcp45 +
                              ph_surface.rcp45 + ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + Gravel.percentage,
                            biomass_g ~  logBathymetry + Quarter),
                       data = dat,
                       mesh = mesh,
                       family = delta_gamma(),
                       spatiotemporal = "rw",
                       offset = "logHaulDur",
                       time = "Year",
                       spatial = "on"
)

sanity(dg_cod_red)
tidy(dg_cod_full, model=1)
tidy(dg_cod_full, model=2)

# store results
model   <- c(model, rep("dg_cod_red", 2))
cAIC    <- c(cAIC, rep(cAIC(dg_cod_red), 2))
formula <- c(formula, as.character(dg_cod_red$formula))



## Model Comparison ##################

models <- data.frame(model = model,
                     cAIC = cAIC,
                     formula = formula)

View(models) # Best Model: dg_cod_red




## Quarter and Ship Variables ##################


## spatially-varying Quarter ##
dg_cod_red_svQ <- sdmTMB(list(biomass_g ~ logBathymetry + thetao_surface.rcp45 + current_surface.rcp45 + ph_surface.rcp45 + 
                                ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + Gravel.percentage,
                              biomass_g ~ logBathymetry),
                         data = dat,
                         mesh = mesh,
                         family = delta_gamma(),
                         spatial_varying = ~ Quarter,
                         spatiotemporal = "rw",
                         offset = "logHaulDur",
                         time = "Year",
                         spatial = "on"
)
  
sanity(dg_cod_red_svQ)
tidy(dg_cod_red_svQ, model=1)
tidy(dg_cod_red_svQ, model=2)

# store results
model   <- c(model, rep("dg_cod_red_svQ", 2))
cAIC    <- c(cAIC, rep(cAIC(dg_cod_red_svQ), 2))
formula <- c(formula, as.character(dg_cod_red_svQ$formula))

extra_years <- as.numeric(c(2025:2045))

## ship random-intercept ##
dg_cod_red_svQ_Sri <- sdmTMB(list(biomass_g ~ logBathymetry + thetao_surface.rcp45 + current_surface.rcp45 + ph_surface.rcp45 + 
                                    ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + (1 | Ship),
                                  biomass_g ~  logBathymetry + (1 | Ship)),
                             data = dat,
                             mesh = mesh,
                             time = "Year",
                             spatiotemporal = "rw",
                             family = delta_gamma(),
                             spatial_varying = ~ Quarter,
                             offset = "logHaulDur",
                             spatial = "on",
)

test <- update(dg_cod_red_svQ_Sri, extra_time = c(2030, 2040, 2050, 2060))


sanity(dg_cod_red_svQ_Sri)
tidy(dg_cod_red_svQ_Sri, model=1)
tidy(dg_cod_red_svQ_Sri, model=2) #49909.75

# store results
model   <- c(model, rep("dg_cod_red_svQ_Sri", 2))
cAIC    <- c(cAIC, rep(cAIC(dg_cod_red_svQ_Sri), 2))
formula <- c(formula, as.character(dg_cod_red_svQ_Sri$formula))




## update distribution - dln #########
dln_cod_red_svQ_Sri <- update(dg_cod_red_svQ_Sri, family = delta_lognormal())

sanity(dln_cod_red_svQ_Sri)
tidy(dln_cod_red_svQ_Sri, model=1)
tidy(dln_cod_red_svQ_Sri, model=2)

# store results
model   <- c(model, rep("dln_cod_red_svQ_Sri", 2))
cAIC    <- c(cAIC, rep(cAIC(dln_cod_red_svQ_Sri), 2))
formula <- c(formula, as.character(dln_cod_red_svQ_Sri$formula))




## update distribution - tw ##########
tw_cod_red_svQ_Sri <- update(tw_cod_red, formula = biomass_g ~ logBathymetry + thetao_surface.rcp45 + current_surface.rcp45 +
                               ph_surface.rcp45 + ph_bottom.rcp45 + totoc_bottom.rcp45 + (1| Ship), spatial_varying = ~ Quarter)

sanity(tw_cod_red_svQ_Sri)
tidy(tw_cod_red_svQ_Sri)

# store results
model   <- c(model, "tw_cod_red_svQ_Sri")
cAIC    <- c(cAIC, cAIC(tw_cod_red_svQ_Sri))
formula <- c(formula, as.character(tw_cod_red_svQ_Sri$formula))





## Model Comparison ##################

models <- data.frame(model = model,
                     cAIC = cAIC,
                     formula = formula)

View(models) # Best Model: dg_cod_red





## Diagnostic plots ##################

model_diagnostics <- list(dg_cod_red_svQ_Sri, dln_cod_red_svQ_Sri, tw_cod_red_svQ_Sri) # best three models
model_diagnostics <- c("dg_red_svQ_Sri", "dln_red_svQ_Sri", "tw_red_svQ_Sri") # best three models


model_diagnostics <- list(dg_cod_red_svQ_Sri)

dharma_r <- list()

# simulate dharmr residuals
for (i in 1:length(model_diagnostics)) {
  
  sim <- simulate(model_diagnostics[[i]], nsim = 500, seed = 123, type = "mle-mvn")
  dharma_r <- c(dharma_r, list(dharma_residuals(sim, model_diagnostics[[i]], return_DHARMa = TRUE)))
  
}


# plot residuals
plot(dharma_r[[1]])

DHARMa::testResiduals(dharma_r[[1]]) ### fails first two tests -- is this alright??
DHARMa::testZeroInflation(dharma_r[[1]]) ### seems okay
DHARMa::testSpatialAutocorrelation(dharma_r[[1]], x = dat$X, y = dat$Y) ## didn't run

DHARMa::plotResiduals(dharma_r[[1]]) ## seems okay
DHARMa:: testTemporalAutocorrelation(dharma_r, time = dat$Year) ## didn't run





## Predictions ##################

# best model: dg_cod_red_svQ_Sri
csq <- p_grid <- readRDS(file = "cod_test/prediction_grid.Rds")

p_grid <- p_grid[,-c(9,10)] %>% 
  mutate(logBathymetry = log(Bathymetry),
         Quarter = as.factor("1"),
         Ship = as.factor("58G2"),
         logHaulDur = log(30)) %>%
  select("lon", "lat", "csquare_area", "id", ,"c_square", "stat_rec", 
         "ices_area",  "ecoregion", "res", "source", "X", "Y", "within_mesh",
         "Gravel.percentage", "Year", "logBathymetry", "Quarter", "Ship",
         "current_surface.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45", 
         "totoc_bottom.rcp45", "zooc_surface.rcp45", "thetao_surface.rcp45",
         "logHaulDur", "Bathymetry")

# remove NAs
p_grid <- p_grid[complete.cases(p_grid),]

plot(p_grid)
  


## rescale covariates
covariates <- c("thetao_surface.rcp45", "current_surface.rcp45",
                "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45","logBathymetry")
## rescale
p_grid[,covariates] <- sapply(covariates, function(x){(p_grid[,x] - mean(rdat[,x]))/sd(rdat[,x])})


# ## test with distance to observations
# test <- dat[,c(47,48)]
# 
# chull(test$X,test$Y)
# 
# test2 <- sapply(1:nrow(csq),function(n){min(sqrt((dat$X - csq$X[n])^2 + (dat$Y - csq$Y[n])^2))})
# csq$dist <- test2
# 
# p_grid <- left_join(p_grid, csq[,c("c_square","Year","dist")], by = c("c_square", "Year"))

# model predictions
p1 <- predict(dg_cod_red_svQ_Sri, newdata = p_grid, type = "response", 
             offset = p_grid$logHaulDur) #re.form = NA

#[p_grid$dist <= 50,]$


range(p1$est)
# p$est_biomass <- (1/(1+exp(-p$est1))) * exp(p$est2)   



dat2 <- rdat
dat2$Ship <- as.factor("58G2")
dat2$Quarter <- as.factor("1")
dat2[,covariates] <- sapply(covariates, function(x){(dat2[,x] - mean(rdat[,x]))/sd(rdat[,x])})
dat2$logHaulDur <- log(30)

names(dat2)
unique(p_grid$Year)
unique(dat2$Year)

# model predictions
p2 <- predict(dg_cod_red_svQ_Sri, newdata = dat2, type = "response", 
              offset = dat2$logHaulDur) #re.form = NA

## should be:
# > nrow(p_grid)
# [1] 315739
# > nrow(csq)
# [1] 366410


i <- 1
par(mfrow=c(2,1))
hist(p1[,covariates[i]])
hist(p2[,covariates[i]])


res <- left_join(csq[,c("lon", "lat", "Year")], p1[,c("lon","lat","Year","est","Quarter")], by = c("lon","lat","Year"))
res$Quarter <- 1

# save biomass estimates
write.csv(res,file="res/cod_biomass.csv", na = "NA")





range(rdat$Bathymetry)

summary(res)

## plot densities ##
ggplot() +
  geom_point(data = dat, aes(x = X, y = Y, colour = Ship)) +
  facet_wrap(~ Ship) +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")


ggplot() +
  geom_point(data = res, aes(x = lon, y = lat, colour = est)) +
  geom_polygon(data = uk, aes(x = long, y = lat, group = group)) +
  #geom_point(data = p_grid[p_grid$within_mesh == TRUE,], aes(x = X, y = Y), color = "grey") +
  #geom_point(data = dat, aes(x = X, y = Y), colour = "pink") +
  facet_wrap(~Year) +
  scale_colour_viridis_c() +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")





# bathy plot
ggplot() +
  geom_point(data = p, aes(x = X, y = Y, colour = Bathymetry)) +
  geom_point(data = csq[csq$dist <= 20,], aes(x = X, y = Y)) +
  scale_colour_viridis_c() +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")


# Extract the spatially varying slope for depth_scaled
# This is stored in the random effects
svQ <- predict(dg_cod_red_svQ_Sri, newdata = p_grid, re_form = ~0, what = "spatial_varying")

# Combine with coordinates
plot_data <- cbind(p_grid, slope = svQ$Quarter)

# Plot using ggplot2
ggplot(plot_data, aes(X, Y, fill = slope)) +
  geom_raster() +
  coord_fixed() +
  scale_fill_viridis_c(option = "plasma") +
  labs(
    title = "Spatially Varying Effect of Quarter",
    fill = "Slope"
  ) +
  theme_minimal()














## covariate effects ## (not working yet)
visreg_delta(dg_cod_red_svQ_Sri, xvar = "logBathymetry", model = 1)


visreg_delta(d_norm_cod_red,xvar="current_surface.rcp45",model=1)

test_cov <- c("logBathymetry","thetao_surface.rcp45","current_surface.rcp45",
              "current_bottom.rcp45","ph_surface.rcp45","ph_bottom.rcp45",      
              "totoc_bottom.rcp45", "zooc_surface.rcp45")

cov(dat[,test_cov])
