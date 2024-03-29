
require(mapdata); require(maptools)
library(rgdal)
#library(ggmap)
library(sp)
library(dplyr)
library(raster) 
#library(reshape)
#library(reshape2)
library(rgeos)
library(rgdal)
require(ggspatial)
library(sf)
library(tidyverse)
library(stringr)

#library(tmap)

# install.packages("ncdf4")
library(ncdf4)


## GLOBALS
#rastname_base = '/Volumes/Padlock/covid/rasts_POC_1/'
rastname_base = '~/covid/output/rasts_POC_any_2019/'

## Air surface
county_ref = read.csv('~/covid/data/california_county_dict.csv')
epa_raw19 = read.csv('~/covid/data/EPA_2019_pm25_CA.csv')
epa_raw20 = read.csv('~/covid/data/EPA_2020_pm25_CA.csv')
epa_raw21 = read.csv('~/covid/data/EPA_2021_pm25_CA.csv')
# filter by POC (comment out to use all sensors)
#epa_raw <- epa_raw %>% filter(POC == 1)
epa_raw = rbind(epa_raw19, epa_raw20)
epa_raw = rbind(epa_raw, epa_raw21)
# How many counties do not have an EPA sensor?
epa_counties <- unique(epa_raw$COUNTY_CODE)
#print(nrow(county_ref) - length(epa_counties))

# Which counties do not have a sensor?
missing_fips <- setdiff(county_ref$county_fips, epa_counties)
counties_no_sensors = county_ref[county_ref$county_fips %in% missing_fips, ]
#print(counties_no_sensors$county_name)

# Distance matrix of county centroids to sensors
distance_matrix <- read.csv('~/covid/data/distance_matrix.csv')
colnames(distance_matrix) <- c('county_fips', 'sensor_id', 'distance')


# Which counties have multiple sensors?
sensor_locations <- epa_raw %>% 
  dplyr::select(Site.ID, Site.Name, CBSA_CODE, CBSA_NAME, STATE_CODE, STATE, COUNTY_CODE, COUNTY, SITE_LATITUDE, SITE_LONGITUDE) %>%
  distinct()

sensors_by_county <- sensor_locations %>%
  group_by(COUNTY_CODE) %>% dplyr::summarise(num_sensors = n())

# add number of EPA sensors to the county reference df
county_ref = merge(county_ref, sensors_by_county, by.x = 'county_fips', by.y = 'COUNTY_CODE', all.x = TRUE)
county_ref[is.na(county_ref)] <- 0

pm25 <- raster("~/covid/geodata/V4NA03_PM25_NA_201801_201812-RH35.nc")
pmdat <- projectRaster(pm25, crs = "+proj=longlat +datum=WGS84")



# determine county multipliers
epa_rast <- raster::extract(pmdat, SpatialPointsDataFrame(sensor_locations[c("SITE_LONGITUDE", 'SITE_LATITUDE')], 
                                                          data = sensor_locations['Site.ID']), df=TRUE)
epa_rast <- cbind(epa_rast, sensor_locations$Site.ID)
colnames(epa_rast) <- c('ID', 'PM2.5', 'Site.ID') #rename columns

epa <- merge(epa_raw, epa_rast, by='Site.ID', all.x = TRUE)

epa <- epa %>% mutate(multiplier = Daily.Mean.PM2.5.Concentration/PM2.5)

# if any Site.IDs are duplicated, select the lowest POC and discard the others 
#  delete this if we are only using POC == 1
epa_rmv <- epa %>% group_by(Site.ID, Date) %>% 
  dplyr::summarize(dup_count = n(), keep_POC = min(POC))

# remove rows where Site.ID = '' and POC = '')
epa  <- merge(epa, epa_rmv, by = c('Site.ID', 'Date'))
epa <- epa %>% filter(POC == keep_POC) %>% dplyr::select(-c('dup_count', 'keep_POC'))

# compute the average multiplier by county
county_mult <- epa %>% dplyr::group_by(COUNTY_CODE, Date) %>%
  dplyr::summarise(
    mean_mult = mean(multiplier), stdev_mult = sd(multiplier), num_measures = n())

# Read in a map of California counties, show number of sensors per county
ca_county = rgdal::readOGR('~/covid/geodata/cb_2018_us_county_5m')
ca_county <- subset(ca_county, STATEFP == '06')
ca_county$county_int = as.integer(ca_county$COUNTYFP)
ca_county <- merge(ca_county, county_ref, by.x = 'county_int', by.y = 'county_fips')
pm25_ca <- mask(pm25, ca_county)

# Multiply each county.

# basic raster with extent
r<- raster(ncol = ncol(pm25_ca), nrow = nrow(pm25_ca))
extent(r) <- extent(ca_county)

#   for each day, rasterize the shapefile; the value of each county is its multiplier
#   multiply the rasters
#   store each day as an output
datelist <- seq(as.Date('2019-11-15'), as.Date('2019-12-31'), by='days')
for (i in 1:length(datelist)) {
  #print(paste(i, 'start', Sys.time()))
  rastname = paste(rastname_base, datelist[i], '.grd', sep = '')
  datemult <- county_mult %>% filter(as.Date(Date, '%m/%d/%Y') == datelist[i])
  if (!file.exists(rastname)) { 
  if (nrow(datemult) < 58) {
    # deal with counties without a sensor
    tmp_epa <- epa %>% filter(as.Date(Date, '%m/%d/%Y') == datelist[i])  # epa sensors available for the current date
    tmp_fips <- setdiff(county_ref$county_fips, datemult$COUNTY_CODE)  # what counties are we missing?
    
    tmp_dist <- distance_matrix %>% filter(county_fips %in% tmp_fips, sensor_id %in% tmp_epa$Site.ID) %>% 
      group_by(county_fips) %>%
      dplyr::summarise(nearest_epa = sensor_id[which.min(distance)])  # find the nearest sensor available
    
    if (nrow(tmp_dist) + nrow(datemult) != 58) {
      print(paste('No sensors for the counties:', setdiff(tmp_fips, tmp_dist$county_fips)))
    }
    
    tmp_add = merge(tmp_epa, tmp_dist, by.x = 'Site.ID' , by.y = 'nearest_epa')
    
    # reassign county, create and select necessary columns
    tmp_add <- tmp_add %>% mutate(COUNTY_CODE = county_fips, 
                                  mean_mult = multiplier, stdev_mult = 0, num_measures = 1) %>% 
      dplyr::select(COUNTY_CODE, Date, mean_mult, stdev_mult, num_measures)
    
    datemult <- rbind(datemult, tmp_add)  # add on values for the missing counties
  }
  
  county_full <- merge(ca_county, datemult, by.x = 'county_int', by.y = 'COUNTY_CODE')
  
  # make a raster of the county file, value set to county multiplier
  rast_mult <- rasterize(county_full, r, 'mean_mult')
  
  rs_rast_mult <- resample(rast_mult, pm25_ca)  # TODO: why is this so SLOW??
  
  # multiply the two rasters and save the output
  rast_date = rs_rast_mult * pm25_ca

  writeRaster(rast_date, filename = rastname)
}
}

















