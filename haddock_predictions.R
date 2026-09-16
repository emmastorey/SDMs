
## herring model diagnostics and predictions ###########


## survey data
full_data <- readRDS("data/Allspecies_15_24_Datras.Rds")
full_data <- add_utm_columns(full_data, ll_names = c("lon", "lat"), 
                             ll_crs = 4326, units = "km")
dat <- full_data[full_data$cmn_nm == "haddock",]

# log transform covariates
dat$logHaulDur    <- log(dat$HaulDur)
dat$logBathymetry <- log(dat$Bathymetry) 




## model diagnostics ###########


# model results
load("C:/Users/Public/010 Multispecies Modelling/res/haddock_model_results.RData")
load("C:/Users/Public/010 Multispecies Modelling/res/haddock_models.RData")

View(model_results)

model_diagnostics <- all_models[c("dg_full_svQ_Sri","dln_full_svQ_Sri", "tw_full_svQ_Sri")] # best three models

dharma_r <- list()

# simulate dharmr residuals
for (i in 1:length(model_diagnostics)) {
  
  sim <- simulate(model_diagnostics[[i]], nsim = 500, seed = 123, type = "mle-mvn")
  dharma_r <- c(dharma_r, list(dharma_residuals(sim, model_diagnostics[[i]], return_DHARMa = TRUE)))
  
}


# plot residuals

i <- 3

plot(dharma_r[[i]])

DHARMa::testResiduals(dharma_r[[i]]) ### fails first two tests -- is this alright??
DHARMa::testZeroInflation(dharma_r[[i]]) ### seems okay
DHARMa::testSpatialAutocorrelation(dharma_r[[i]], x = dat$X, y = dat$Y) ## didn't run

DHARMa::plotResiduals(dharma_r[[i]]) ## seems okay
DHARMa:: testTemporalAutocorrelation(dharma_r, time = dat$Year) ## didn't run




## Predictions ##################

# best_model <- all_models$dg_full_svQ_Sri

# Prediction Grid
csq    <- readRDS(file = "cod_test/prediction_grid.Rds")[,-c(1:2)]
p_grid <- csq

p_grid <- p_grid %>% 
  mutate(logBathymetry = log(Bathymetry),
         Quarter = as.factor("1"),
         Ship = as.factor("58G2"),
         logHaulDur = log(30)) %>%
  select("lon", "lat","c_square", "stat_rec", 
         "ices_area",  "ecoregion", "res", "source", "X", "Y",
         "Year", "logBathymetry", "Quarter", "Ship", "ph_surface.rcp45", 
         "ph_bottom.rcp45", "totoc_bottom.rcp45", "zooc_surface.rcp45", 
         "thetao_surface.rcp45", "dist2coast", "logHaulDur", "Bathymetry")

summary(p_grid)

# remove NAs
p_grid <- p_grid[complete.cases(p_grid),]


## rescale covariates
covariates <- c("logBathymetry", "dist2coast",  
                "thetao_surface.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45")
## rescale
p_grid[,covariates] <- sapply(covariates, function(x){(p_grid[,x] - mean(dat[,x]))/sd(dat[,x])})

# model predictions
p1 <- predict(best_model, newdata = p_grid, type = "response", offset = p_grid$logHaulDur) #re.form = NA

range(p1$est)


# add to full grid
res <- csq %>%
  left_join(p1 %>% select(c_square, Year, est), by = c("c_square", "Year"))
res <- res[!duplicated(res), ]
res$Quarter <- 1
res <- res[res$Bathymetry <= 260,]
res <- res %>%
  select(c("Year", "c_square", "lat", "lon", "est", "Quarter",)) %>%
  filter(!if_all(everything(), is.na))


# view results
ggplot() +
  geom_point(data = res, aes(x = lon, y = lat, colour = est)) +
  facet_wrap(~Year) +
  scale_colour_viridis_c() +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")


# save biomass estimates
write.csv(res, file="res/haddock_biomass.csv", na = "NA", row.names = FALSE)


nrow(res)





### rerun model ################



## covariates ##
dat$logHaulDur    <- log(dat$HaulDur)
dat$logBathymetry <- log(dat$Bathymetry)

covariates <- c("logBathymetry", "dist2coast",  
                "thetao_surface.rcp45", "current_surface.rcp45",
                "current_bottom.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45")
## rescale
dat[,covariates] <- apply(dat[,covariates],2,function(x){(x-mean(x))/sd(x)}) 


## mesh ##
load("res/allspp_mesh.RDS")
mesh <- sp_mesh$haddock


## model ##

best_model # 59703.18

dg_red_svQ_Sri <- sdmTMB(formula = biomass_g ~ (1 | Ship) + logBathymetry + thetao_surface.rcp45 + 
                           current_surface.rcp45 + current_bottom.rcp45 + ph_surface.rcp45 + 
                           ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + 
                           dist2coast,
                         data = dat,
                         mesh = mesh,
                         family = delta_gamma(),
                         spatiotemporal = "rw",
                         offset = "logHaulDur",
                         time = "Year",
                         spatial_varying = ~ Quarter,
                         spatial = "on"
                         )

sanity(dg_red_svQ_Sri)
tidy(dg_red_svQ_Sri)
cAIC(dg_red_svQ_Sri) # 59698.82

# remove current_bottom
dg_red_1_svQ_Sri <- update(dg_red_svQ_Sri, formula = biomass_g ~ (1 | Ship) + logBathymetry + 
                             thetao_surface.rcp45 + current_surface.rcp45 + ph_surface.rcp45 + 
                             ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + 
                             dist2coast)

sanity(dg_red_1_svQ_Sri)
tidy(dg_red_1_svQ_Sri)
cAIC(dg_red_1_svQ_Sri) # 59697.99

# remove current_surface (inc current_bottom)
dg_red_2_svQ_Sri <- update(dg_red_svQ_Sri, formula = biomass_g ~ (1 | Ship) + logBathymetry + 
                             thetao_surface.rcp45 + current_bottom.rcp45 + ph_surface.rcp45 + 
                             ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + 
                             dist2coast)

sanity(dg_red_2_svQ_Sri)
tidy(dg_red_2_svQ_Sri, model = 1)
tidy(dg_red_2_svQ_Sri, model = 2)
cAIC(dg_red_2_svQ_Sri) # 59696.96

# remove current_surface and current_bottom
dg_red_3_svQ_Sri <- update(dg_red_svQ_Sri, formula = biomass_g ~ (1 | Ship) + logBathymetry + 
                             thetao_surface.rcp45 + ph_surface.rcp45 + ph_bottom.rcp45 + 
                             totoc_bottom.rcp45 + zooc_surface.rcp45 + dist2coast)

sanity(dg_red_3_svQ_Sri)
tidy(dg_red_3_svQ_Sri, model = 1)
tidy(dg_red_3_svQ_Sri, model = 2)
cAIC(dg_red_3_svQ_Sri) # 59695.85



# adjust gamma model covariates
dg_red_4_svQ_Sri <- sdmTMB(formula = list(biomass_g ~ (1 | Ship) + logBathymetry + thetao_surface.rcp45 + 
                           current_surface.rcp45 + current_bottom.rcp45 + ph_surface.rcp45 + 
                           ph_bottom.rcp45 + totoc_bottom.rcp45 + zooc_surface.rcp45 + 
                           dist2coast,
                           biomass_g ~ (1 | Ship) + logBathymetry + dist2coast),
                         data = dat,
                         mesh = mesh,
                         family = delta_gamma(),
                         spatiotemporal = "rw",
                         offset = "logHaulDur",
                         time = "Year",
                         spatial_varying = ~ Quarter,
                         spatial = "on"
)

sanity(dg_red_4_svQ_Sri)
tidy(dg_red_4_svQ_Sri, model = 2)
cAIC(dg_red_4_svQ_Sri) #59701.22 59699.2

best_model <- dg_red_3_svQ_Sri
