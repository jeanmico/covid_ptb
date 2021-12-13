library(tmap)
library(ggplot2)
library(tidyr)
library(ggmap)
library(rgdal)
library(gdalUtils)
library(rgeos)
library(maptools)
library(dplyr)
library(tmap)
library(gstat)
library(raster)
library(sp)
library(sf)
library(cartography)


ca_counties  <- readOGR('~/JelliffeWitte/ref_geo/cb_2018_us_county_5m/cb_2018_us_county_5m.shp')
plot(ca_counties)

myrates = read.csv('/Volumes/Padlock/covid/output/20211129_rates_to_map.csv')

ca_counties <- merge(ca_counties, myrates, by.x = 'County_num', by.y = 'countynum', all.x = T)
pm_tiles = quantile(ca_counties$pmmean, probs = seq(0,1,.2))
pm_colors = c("#ffffb2", "#fecc5c", "#fd8d3c", "#f03b20", "#bd0026")

pt_tiles = quantile(ca_counties$ptb_rate, probs = seq(0,1,.2), na.rm = T)
pt_colors = c()



ca_counties_full = ca_counties[ca_counties$birthcount > 20,]
ca_counties_na = ca_counties[ca_counties$birthcount <=20, ]
ca_na = st_as_sf(ca_counties_na)


ca_na_sfc = ca_na %>% hatchedLayer(mode = 'sfc', pattern = 'left2right', density = 1)

ca_na_patterns = st_sf(geometry = ca_na_sfc)

pm_map = 
  tm_shape(ca_na_patterns, bbox = ca_counties) + 
    tm_lines(col = "grey50") + 
  tm_shape(ca_counties_full) + 
  tm_fill("pmmean", breaks = pm_tiles,  
          legend.show = T, 
          title = expression('Mean PM'[2.5]), 
          legend.format = c(digits = 1)) +  
  tm_layout(legend.position = c('right', 'top'),
            legend.title.size = 4,
            legend.text.size = 1.5) + 
  tm_borders() #+ 
  #  tm_legend()
pm_map
tmap_save(pm_map, file = '/Volumes/Padlock/covid/figures/pm_rate_map.png')

ptb_map = 
  tm_shape(ca_na_patterns, bbox = ca_counties) + 
  tm_lines(col = "grey50") + 
  tm_shape(ca_counties_full) + 
  tm_fill("ptb_rate", breaks = pt_tiles,  
          legend.show = T, 
          title = expression('PTB rate (%)'), 
          legend.format = c(digits = 1)) +  
  tm_layout(legend.position = c('right', 'top'),
            legend.title.size = 4,
            legend.text.size = 1.5) + 
  tm_borders() 
ptb_map

tmap_save(ptb_map, file = '/Volumes/Padlock/covid/figures/ptb_rate_map.png')
