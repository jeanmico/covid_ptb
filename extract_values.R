lookup_raw = read.csv('/Volumes/Padlock/covid/lookup_points.csv')

filenames = list.files('/Volumes/Padlock/covid/rasts_POC_1', pattern = ".grd", full.names = TRUE)
outbase = paste('/Volumes/Padlock/covid/values_extracted1/')

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
