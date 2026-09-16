



## Predictions ##################

# best model: dg_cod_red_svQ_Sri
csq    <- readRDS(file = "predictions/prediction_grid.Rds")[,-c(1:2)]
p_grid <- csq


# define a mesh directly with fmesher (formerly INLA):
bdry <- INLA::inla.nonconvex.hull(cbind(dat$X, dat$Y), 
                                  -0.025, resolution = c(200, 50))
# Convert to data frame
bdry <- data.frame(
  X = bdry$loc[, 1],
  Y = bdry$loc[, 2]
)

# Bathymetry and mesh boundary
# ggplot() +
#   geom_point(data = p_grid, aes(x = X, y = Y, colour = Bathymetry)) +
#   geom_polygon(data = bdry, aes(x = X, y = Y), fill = "skyblue", alpha = 0.1, color = "blue") +
#   scale_colour_viridis_c() +
#   theme_light() +
#   labs(x = "Longitude", y = "Latitude")




## predictions for cod
p_grid <- p_grid %>% 
  mutate(logBathymetry = log(Bathymetry),
         Quarter = as.factor("1"),
         Ship = as.factor("58G2"),
         logHaulDur = log(30)) %>%
  select("lon", "lat","c_square", "stat_rec", 
         "ices_area",  "ecoregion", "res", "source", "X", "Y", 
         "Year", "logBathymetry", "Quarter", "Ship",
         "current_surface.rcp45", "ph_surface.rcp45", "ph_bottom.rcp45", 
         "totoc_bottom.rcp45", "zooc_surface.rcp45", "thetao_surface.rcp45",
         "logHaulDur", "Bathymetry")

# test <- p_grid[apply(is.na(p_grid), 1, any), ]
# plot(test$lon, test$lat)
# plot(p_grid$lon, p_grid$lat)
# summary(test)

# remove NAs
p_grid <- p_grid[complete.cases(p_grid),]



## rescale covariates
covariates <- c("thetao_surface.rcp45", "current_surface.rcp45",
                "ph_surface.rcp45", "ph_bottom.rcp45",      
                "totoc_bottom.rcp45", "zooc_surface.rcp45","logBathymetry")
## rescale
p_grid[,covariates] <- sapply(covariates, function(x){(p_grid[,x] - mean(rdat[,x]))/sd(rdat[,x])})


# model predictions
p1 <- predict(test, newdata = p_grid, type = "response", #dg_cod_red_svQ_Sri
              offset = p_grid$logHaulDur) #re.form = NA

range(p1$est)

# results
res <- csq %>%
  left_join(p1 %>% select(c_square, Year, est), by = c("c_square", "Year"))
res <- res[!duplicated(res), ]
res$Quarter <- 1
res <- res[res$Bathymetry <= 260,]
res <- res %>%
  select(c("Year", "c_square", "lat", "lon", "est", "Quarter")) %>%
  filter(!if_all(everything(), is.na))


# plot
ggplot() +
  geom_point(data = res, aes(x = lon, y = lat, colour = est)) +
  #geom_polygon(data = uk, aes(x = long, y = lat, group = group)) +
  #geom_point(data = p_grid[p_grid$within_mesh == TRUE,], aes(x = X, y = Y), color = "grey") +
  #geom_point(data = dat, aes(x = X, y = Y), colour = "pink") +
  facet_wrap(~Year) +
  scale_colour_viridis_c() +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")

tidy(test)


# save biomass estimates
write.csv(res,file="res/cod_biomass.csv", na = "NA", row.names = FALSE)

nrow(res)









# checks
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





summary(res)

## plot densities ##
ggplot() +
  geom_point(data = dat, aes(x = X, y = Y, colour = Ship)) +
  facet_wrap(~ Ship) +
  theme_light() +
  labs(x = "Longitude", y = "Latitude")






# bathy plot
ggplot() +
  geom_point(data = csq, aes(x = X, y = Y, colour = Bathymetry)) +
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









