
setwd("C:/Users/Public/010 Multispecies Modelling")

## libraries
library(sf)
## libraries
library(sp)
library(dplyr)
library(spatialEco)
library(terra)
library(sdmTMB)
library(ggplot2)
library(inlabru)

## prediction grid #####

years <- c(2015:2024, 2030, 2040, 2050, 2060)

# load c-squares for prediction
load("North_Sea_csquare_grid_no_land_final.RData")
pred_csq <- North_Sea_csquare_no_land_updated

# change polygons to centroids
pred_csq$geometry <- st_centroid(st_geometry(pred_csq))

## Add UTM cols
pred_csq <- add_utm_columns(pred_csq, ll_names = c("lon", "lat"), units = "km") # 32631



## Covariate data for Prediction Grid #####


## Static variables (spatial, not temporal)
load(file.path("data/new_gridded_env_data_10km_1980to2099_CERES_Copernicus_and_other_sources", 
               "Static_vars.RData"))

static_env <- rast(env_data_static)
static_env <- terra::project(static_env, "EPSG:4326")
pred_csq   <- vect(pred_csq) # "EPSG:4326"

extracted_vals <- terra::extract(static_env, pred_csq)

pred_csq <- cbind(pred_csq, extracted_vals[,c(5:12)])
pred_csq$Bathymetry <- -pred_csq$Bathymetry


## replicate prediction grid over years
grid_yrs <- replicate_df(pred_csq, "Year", years)



## Extract covariate data for temporal variables

# storage
all_yrs_extracted_vals <- NULL

# for each year
for(i in 1: length(years)) {
  
  print(i)
  
  # load env data
  temp_env_files <- list.files(path = "data//new_gridded_env_data_10km_1980to2099_CERES_Copernicus_and_other_sources",
                               pattern = paste(years[i]), 
                               full.names = TRUE)
  load(temp_env_files)
  
  # transform spatial grid to raster and project to WGS84
  temp_env <- rast(future_env_data)
  temp_env <- terra::project(temp_env, "EPSG:4326")
  
  # extract values at each prediction data point
  extracted_vals <- terra::extract(temp_env, 
                            geom(grid_yrs[grid_yrs$Year == paste(years[i])])[,c("x", "y")])
  
  # keep only surface data for RCP4.5 and remove standard deviations 
  extracted_vals <- extracted_vals %>% select(contains("rcp45") & 
                                                #contains("surface") &
                                                !contains("sd"))
  
  # store values for each year 
  all_yrs_extracted_vals <- rbind(all_yrs_extracted_vals, extracted_vals)
  
}

# combine columns for extracted variables with survey data
values(grid_yrs) <- cbind(values(grid_yrs), all_yrs_extracted_vals)

# convert to dataframe
grid_yrs <- cbind(geom(grid_yrs)[,c(3,4)], values(grid_yrs))

# rename x/y to lon/lat
names(grid_yrs)[1:2] <- c("lon", "lat")

# save data
saveRDS(grid_yrs, file = "predictions/prediction_grid.Rds")





## boundary comparison ########

# load c-squares for prediction
load("North_Sea_csquare_grid_no_land_final.RData")

# c-square centroids
locs <- st_centroid(North_Sea_csquare_no_land_updated)

## Add UTM cols
locs <- add_utm_columns(locs, ll_names = c("lon", "lat"), 
                        units = "km") # 32631

# Convert to data frame
bdry <- data.frame(
  X = non_convex_bdry$loc[, 1],
  Y = non_convex_bdry$loc[, 2]
)



  
  geom_point(data = locs, aes(X, Y)) +
  geom_polygon(
    data = bdry,
    aes(x = X, y = Y),
    fill = "skyblue", alpha = 0.4, color = "blue"
  ) +
  geom_point(
    data = dat,
    aes(x = X, y = Y),
    fill = "pink", alpha = 0.4, color = "pink"
  ) 
