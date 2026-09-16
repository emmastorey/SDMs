
########################################
##
## Pre-processing of NS-IBTS survey data
##
########################################

library(dplyr)
library(ggplot2)
library(icesVocab)
library(icesDatras)
library(terra)
library(tidyterra)



## LOAD DATA ########

yrs <- 2015:2024

raw <- readRDS("data/HH_HL_merged_NS-IBTS.rds")
raw <- raw[raw$Year %in% yrs,]



## PRESENCE DATA ########


# Species data frame
species <- data.frame(latin_nm = c("Merlangius merlangus",	"Clupea harengus", 
                                    "Pleuronectes platessa", "Sprattus sprattus",
                                    "Trisopterus esmarki", "Melanogrammus aeglefinus", 
                                    "Solea solea",	"Pollachius virens",
                                    "Eutrigla gurnardus", "Limanda limanda"),
                      cmn_nm   = c("whiting", "herring", "plaice", "sprat", 
                                    "norway pout", "haddock", "common sole", 
                                    "saithe", "grey gurnard", "dab"),
                      a        = c(0.006, 0.002, 0.007, 0.007, 0.009, 0.005,
                                   0.008, 0.007, 0.004, 0.01),
                      b        = c(3.08, 3.429, 3.101, 3.014, 2.941, 3.16,
                                   3.019, 3.075, 3.198, 2.986)) # need to include sandeel

species$aphia <- as.numeric(findAphia(species$latin_nm, latin = TRUE, full = TRUE)[,3])


## filter by species ##
dat <- raw[raw$SpecCode %in% species$aphia,] 
names(dat)[2] <- "aphia"

# add latin & common names
dat <- left_join(dat, species, join_by(aphia))


## convert lengths to Biomass using w = al^b ##
# see supplementary for https://besjournals.onlinelibrary.wiley.com/doi/10.1111/1365-2664.12238
dat$biomass_g <- dat$a * dat$LngtClass^dat$b * dat$HLNoAtLngt 

# remove lengths, number at length and subfactor and summarise by haul id etc.
dat <- dat %>%
  select(-c("LngtClass", "HLNoAtLngt", "SubFactor", "a", "b")) %>%
  group_by(across(-"biomass_g")) %>%
  summarise(biomass_g = sum(biomass_g))

{dat %>%
  group_by(cmn_nm) %>%
  summarise(unique_count = n_distinct(HaulID)) %>%
  summarise(sum_unique_count = sum(unique_count),
            check = (sum_unique_count == nrow(dat)))} # should be true (31790)



## ABSENCE DATA ########


absences <- NULL

# for each species 
for (i in 1:nrow(species)) {
  
  spp <- species$cmn_nm[i]
  print(spp)
  
  # find hauls where species hasn't been recorded
  spp_abs <- raw[raw$HaulID %in% setdiff(unique(raw$HaulID), unique(dat[dat$cmn_nm == spp,]$HaulID)),]
  spp_abs$biomass_g <- 0 
  
  # add columns for species id
  spp_abs <- spp_abs %>%
    mutate(aphia     = dat[dat$cmn_nm == spp,]$aphia[1],
           cmn_nm    = dat[dat$cmn_nm == spp,]$cmn_nm[1],
           latin_nm  = dat[dat$cmn_nm == spp,]$latin_nm[1],
           biomass_g = 0)
  
  absences <- rbind(spp_abs, absences)
  rm(spp_abs)
}


# remove lengths, number at length and subfactor and summarise by haul id etc.
absences <- absences %>%
  select(-c("LngtClass", "HLNoAtLngt", "SubFactor", "SpecCode")) %>%
  group_by(across(-"biomass_g")) %>%
  summarise(biomass_g = sum(biomass_g))

# combine presence and absence data
dat <- rbind(dat, absences)


# add column for ship id
dat$Ship <- unlist(lapply(strsplit(dat$HaulID, ":", fixed = TRUE), function(x) x[5]))
dat$Ship <- as.factor(dat$Ship)

# Quarter as factor variable
dat$Quarter <- as.factor(dat$Quarter)


{
length(unique(dat$HaulID))
dat %>%
  group_by(cmn_nm) %>%
  summarise(Hauls = n_distinct(HaulID))} # check HaulID by Species (5355)



## SURVEY PLOTS ########


# convert to vector data
dat <- vect(dat, geom=c("ShootLong", "ShootLat"), crs="EPSG:4326")
world <- map_data("world")

# plots
p <- list()

for (i in 1:nrow(species)) {
  
  # species
  spp <- species$cmn_nm[i]

  # plot survey locations across years
  p[[i]] <- ggplot() + 
    geom_sf() + 
    geom_polygon(data = world, mapping = aes(x = long, y = lat, group = group), fill = "dark grey", color = "black") +
    geom_spatvector(data = dat[dat$cmn_nm == spp & dat$biomass_g != 0,], aes(colour = (biomass_g))) +
    coord_sf(xlim = ext(dat)[1:2] *  c(1.1,0.9),
             ylim = ext(dat)[3:4] * c(0.99,1.01)) +
    theme_bw() +
    theme(axis.line = element_line(color='black'),
          plot.background = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank()) +
    labs(x = "Longitude", y = "Latitude")+ 
    facet_wrap(~Year) + 
    theme(plot.title = element_text(size=10,face="bold"))
  
}

spp <- 5
species$cmn_nm[spp]
p[[spp]]




## EXTRACT COVARIATES DATA ########


## static variables ##

load(file.path("data/new_gridded_env_data_10km_1980to2099_CERES_Copernicus_and_other_sources", 
               "Static_vars.RData"))

static_env     <- terra::rast(env_data_static)
static_env     <- terra::project(static_env, "EPSG:4326")
extracted_vals <- terra::extract(static_env, geom(dat)[,3:4])

dat <- cbind(dat, extracted_vals[,c(4:8)])

dat$Bathymetry <- -dat$Bathymetry # change to positive


## temporal variables ##


# storage
all_yrs_extracted_vals <- NULL

# for each year
for(i in 1: length(yrs)) {
  
  # load env data
  temp_env_files <- list.files(path = "data//new_gridded_env_data_10km_1980to2099_CERES_Copernicus_and_other_sources",
                               pattern = paste(yrs[i]), 
                               full.names = TRUE)
  load(temp_env_files)
  
  # transform spatial grid to raster and project to WGS84
  temp_env <- terra::rast(future_env_data)
  temp_env <- terra::project(temp_env, "EPSG:4326")
  
  # extract values at each survey data point
  extracted_vals <- terra::extract(temp_env, 
                            geom(dat[dat$Year == paste(yrs[i])])[,c("x", "y")])
  
  # keep only surface data for RCP4.5 and remove standard deviations 
  extracted_vals <- extracted_vals %>% select(contains("rcp45") & 
                                                #contains("surface") &
                                                !contains("sd"))
  
  # store values for each year 
  all_yrs_extracted_vals <- rbind(all_yrs_extracted_vals, extracted_vals)
  
}

# combine columns for extracted variables with survey data
dat <- cbind(dat, all_yrs_extracted_vals)

# convert to dataframe
dat <- cbind(geom(dat)[,c(3,4)], values(dat))

# remove rows with NAs
dat <- dat[complete.cases(dat), ]

# rename x/y to lon/lat
names(dat)[1:2] <- c("lon", "lat")



## SAVE DATA ##########

saveRDS(dat, file = "data/Allspecies_15_24_Datras.Rds")
