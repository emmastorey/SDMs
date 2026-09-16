
## whiting model diagnostics and predictions ###########


## survey data
full_data <- readRDS("data/Allspecies_15_24_Datras.Rds")
full_data <- add_utm_columns(full_data, ll_names = c("lon", "lat"), 
                             ll_crs = 4326, units = "km")
dat <- full_data[full_data$cmn_nm == "whiting",]

# log transform covariates
dat$logHaulDur    <- log(dat$HaulDur)
dat$logBathymetry <- log(dat$Bathymetry) 




## model diagnostics ###########



best_model <- all_models$dg_red_2_svQ_Sri

# Prediction Grid
csq    <- readRDS(file = "predictions/prediction_grid.Rds")[,-c(1:2)]
p_grid <- csq

p_grid <- p_grid %>% 
  mutate(logBathymetry = log(Bathymetry),
         Quarter = as.factor("1"),
         Ship = as.factor("58G2"),
         logHaulDur = log(30)) %>%
  select("lon", "lat","c_square", "stat_rec", 
         "ices_area",  "ecoregion", "res", "source", "X", "Y",
         "Year", "logBathymetry", "Quarter", "Ship", "current_surface.rcp45", 
         "dist2coast", "logHaulDur", "Bathymetry")

summary(p_grid)

# remove NAs
p_grid <- p_grid[complete.cases(p_grid),]


## rescale covariates
covariates <- c("logBathymetry", "dist2coast","current_surface.rcp45")

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
write.csv(res, file="res/whiting_biomass.csv", na = "NA", row.names = FALSE)







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
mesh <- sp_mesh$whiting


## model ##

best_model # 89601.64

dln_red_2_svQ_Sri <- sdmTMB(formula = list(biomass_g ~ (1 | Ship) + logBathymetry + dist2coast,
                                        biomass_g ~ (1 | Ship) + logBathymetry + current_surface.rcp45),
                         data = dat,
                         mesh = mesh,
                         family = delta_lognormal(),
                         spatiotemporal = "rw",
                         offset = "logHaulDur",
                         time = "Year",
                         spatial_varying = ~ Quarter,
                         spatial = "on"
                         )

test <- update(dg_cod_red_svQ_Sri, extra_time = c(2030, 2040, 2050, 2060))



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
