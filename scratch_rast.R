exposure_plot <- ggplot(data = births, aes(x = pm25mean)) + 
  geom_histogram(binwidth = .5) +
  theme_minimal() +
  xlab('Mean PM2.5 Exposure')
exposure_plot
ggsave(exposure_plot, file = image_name('mean_pm25_hist_fullsample'), dpi = 300)


myrast = raster('/Volumes/Padlock/covid/rasts_POC_any/2020-01-01.grd')
myrast_cr <- crop(myrast, bound_box)


bound_box <<-  as(raster::extent( -125, -114, 32, 43), "SpatialPolygons")
proj4string(bound_box) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"

rastplt <- tm_shape(myrast, bbox = bound_box) + 
  tm_raster(palette = "YlOrBr", title = '2020-01-01 PM2.5', style='fixed', breaks=c(0,5,10,15,20,25,30,35,40)) +
  tm_legend(legend.outside=TRUE, legend.show=TRUE, legend.text.size=1.2)
rastplt

tmap_save(rastplt + tm_layout(outer.margins = c(0,0,0,0)), '/Volumes/Padlock/covid/images/raster_20200101.png')
