
require(mapdata); require(maptools)
library(rgdal)
library(ggmap)
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

library(tmap)

# install.packages("ncdf4")
library(ncdf4)


## GLOBALS
#rastname_base = '/Volumes/Padlock/covid/rasts_POC_1/'
rastname_base = '~/covd/output/rasts_POC_1/'

## FUNCTIONS #######

sensor_locate <- function(sensors) {
  locations <- sensors %>% 
    dplyr::select(Site.ID, Site.Name, POC, CBSA_CODE, CBSA_NAME, STATE_CODE, STATE, COUNTY_CODE, COUNTY, SITE_LATITUDE, SITE_LONGITUDE) %>%
    distinct()
  return(locations)
}

sensor_county <- function(sensors_loc) {
  sensors_by_county <- sensors_loc %>%
    group_by(COUNTY_CODE) %>% dplyr::summarise(num_sensors = n())
}


## Air surface
county_ref = read.csv('~/covid/data/california_county_dict.csv')
epa_raw = read.csv('~/covid/data/EPA_2020_pm25_CA.csv')

# Distance matrix of county centroids to sensors
distance_matrix <- read.csv('~/covid/data/distance_matrix.csv')
colnames(distance_matrix) <- c('county_fips', 'sensor_id', 'distance')

# filter by POC (comment out to use all sensors)
poc1 <- epa_raw %>% filter(POC == 1)

sensor_locations <- sensor_locate(epa_raw)
sensors_by_county <- sensor_county(sensor_locations)

poc1_locations <- sensor_locate(poc1)
poc1_by_county <- sensor_county(poc1_locations)

# add number of EPA sensors to the county reference df
county_ref = merge(county_ref, poc1_by_county, by.x = 'county_fips', by.y = 'COUNTY_CODE', all.x = TRUE)


county_ref = merge(county_ref, sensors_by_county, by.x = 'county_fips', by.y = 'COUNTY_CODE', all.x = TRUE)


county_ref[is.na(county_ref)] <- 0
colnames(county_ref) = c('county_fips', 'county_name', 'state', 'POC_1', 'POC_any')
write.csv(county_ref, file = '~/covid/county_sensor_counts.csv', row.names = FALSE)

write.csv(sensor_locations, file = '~/covid/sensor_locations_POC_any.csv', row.names = FALSE)
write.csv(poc1_locations, file = '~/covid/sensor_locations_POC_1.csv', row.names =FALSE)

minpoc <- sensor_locations %>% group_by(COUNTY_CODE, COUNTY) %>%
  summarise(minimum_POC = min(POC)) %>% filter(minimum_POC >1)
write.csv(minpoc, file = '~/covid/minimum_POC.csv', row.names = FALSE)

# In counties with multiple sensors, how well do they agree?

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
datelist <- seq(as.Date('2020-01-01'), as.Date('2020-12-31'), by='days')
for (i in 1:length(datelist)) {
  #print(paste(i, 'start', Sys.time()))
  datemult <- county_mult %>% filter(as.Date(Date, '%m/%d/%Y') == datelist[i])
  
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

  rastname = paste(rastname_base, datelist[i], '.grd', sep = '')
  print(rastname)
  writeRaster(rast_date, filename = rastname)

}

















