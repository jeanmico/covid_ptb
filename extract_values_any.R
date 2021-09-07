library(raster)
library(dplyr)
library(rgdal)
lookup_raw = read.csv('~/covid/data/lookup_points.csv')

filenames = list.files('~/covid/output/rasts_POC_any_2019', pattern = ".grd", full.names = TRUE)
outbase = paste('~/covid/output/values_extracted_any/')

for (dayfile in filenames) {
  tmpdate = strsplit(basename(dayfile), '\\.')[[1]][1]
  print(paste(tmpdate, 'start', Sys.time()))
  #use extract raster and make a dataframe
  day_rast = raster(dayfile)
  tmpdf <- raster::extract(day_rast, lookup_raw, df=TRUE)
  daydf <- cbind(lookup_raw, tmpdf) 
  daydf <- daydf %>% mutate(date = as.Date(tmpdate)) %>% dplyr::select(-ID)
  
  outfile = paste(outbase, tmpdate, '.csv', sep = '')
  write.csv(daydf, file = outfile, row.names = FALSE)
}
