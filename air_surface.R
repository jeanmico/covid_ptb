require(RColorBrewer); require(ggplot2)
require(mapdata); require(maptools)
library(raster); library(zipcode) 
library(choroplethr)
library(rgdal)
library(ggmap)
library(sp)
require("plyr"); require(dplyr)
library(tgp)
library(mgcv)
library(gstat)
library(automap)
library(raster) 
library(dismo)
library(reshape)
library(reshape2)
library(leaflet); library(rgeos)
library(leaflet.extras)
library(rgdal)
require(ggspatial)
library(mapview); library(webshot)
library(rasterVis)
library(sf)
library(tidyverse)
library(stringr)
library(viridis)
library(scico)
library(patchwork)

# install.packages("ncdf4")
library(ncdf4)

setwd("~/covid/")

#pm25  <- raster("~/covid/V4NA03_PM25_NA_201801_201812-RH35-NoNegs.asc/V4NA03_PM25_NA_201801_201812-RH35-NoNegs.asc")

pm25 <- raster("V4NA03_PM25_NA_201801_201812-RH35.nc")
pmdat <- projectRaster(pm25, crs = "+proj=longlat +datum=WGS84")

## just pulling out a kind of random box somewhat covering NYC, in your case it'll be the box around MA
bay_area <- as(extent(-123, -120.93, 36.6, 38.6 ), 'SpatialPolygons')

pm_bayarea <- crop(pmdat, bay_area)
plot(pm_bayarea) ## see that it still gives full coverage! so no need for idw!


## here i am just creating a dummy dataset with my home and office addresses, but in your case it'd be the data frame with the ALS cases' addresses
mak.addr <- data.frame(rbind(c(-74.0108925, 40.7190305), c(-73.9433492, 40.8424454)))
mak.addr$id <- c("home", "office")
names(mak.addr)[1:2] <- c("long", "lat")

mm1 <- raster::extract(pm_bayarea, SpatialPoints(mak.addr[c("long", "lat")]), df=T)

mm1

# read in the epa data for 2020
epa_raw = read.csv('EPA_2020_pm25_CA.csv')


# what is county coverage?
num_counties = length(unique(epa_raw$COUNTY_CODE))
# CA has 58 total counties; we have data for 51

# what counties are we missing?

epa_location = epa_raw %>% dplyr::select(Site.ID, SITE_LATITUDE, SITE_LONGITUDE) %>% distinct()

epa_rast <- raster::extract(pmdat, SpatialPointsDataFrame(epa_location[c("SITE_LONGITUDE", 'SITE_LATITUDE')], 
                                                          data = epa_location['Site.ID']), df=TRUE)
epa_rast <- cbind(epa_rast, epa_location$Site.ID)
colnames(epa_rast) <- c('ID', 'PM2.5', 'Site.ID')

a = SpatialPointsDataFrame(epa_location[c("SITE_LONGITUDE", 'SITE_LATITUDE')], 
                           data = epa_location['Site.ID'])
a@data
b = raster::extract(pmdat, a)

epa <- merge(epa_raw, epa_rast, by='Site.ID', all.x = TRUE)

epa <- epa %>% mutate(multiplier = Daily.Mean.PM2.5.Concentration/PM2.5)

county_mult <- epa %>% dplyr::group_by(COUNTY_CODE, Date) %>%
  dplyr::summarise(
    mean_mult = mean(multiplier), stdev_mult = sd(multiplier), num_measures = n())

#histograms of standard deviations (for multiple sensors within county)
sd_hist <- ggplot(filter(county_mult, num_measures>1), aes(x=stdev_mult)) +
  geom_histogram()
sd_hist

n_hist <- ggplot(county_mult, aes(x=num_measures)) + 
  geom_histogram()
n_hist
range(county_mult$num_measures)

b <- county_mult %>% filter(num_measures == 1)

range(county_mult$sd_mult, na.rm = TRUE)
range(county_mult$mean_mult, na.rm = TRUE)
