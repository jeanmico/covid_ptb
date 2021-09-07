# investigating low air pollution values
library(data.table)
library(dplyr)

fp_low <<- '/Volumes/Padlock/covid/lowvalues/'

image_low <- function(fname) {
  return(paste(fp_low, Sys.Date(), '_', fname, '.png', sep = ''))
}

output_low <- function(fname) {
  return(paste(fp_low, Sys.Date(), '_', fname, '.csv', sep = ''))
}


low_pm = births[births$pm25mean>0, ]

range(low_pm$pm25mean)
which.min(low_pm$pm25mean)

meanexp_hist <- ggplot(low_pm, aes(x = pm25mean)) + 
  geom_histogram() + 
  ggtitle('Mean exposure during pregnancy') + 
  theme_bw()
meanexp_hist
ggsave(meanexp_hist, file = image_low('mean_exp_hist'))


lowexp_hist = ggplot(under5, aes(x = pm25mean)) + 
  geom_histogram(binwidth = .5) + 
  ggtitle('Mean exposure - detail of low exposure') + 
  theme_bw()
lowexp_hist
ggsave(lowexp_hist, file = image_low('low_exp_hist'))

under5write = under5 %>% dplyr::select(longitude, latitude, VS_unique, census_tract, pm25mean)

under5_rm <- under5 %>% filter(na_days <20)
write.csv(under5write, file = output_low('pregnancies_under5'), row.names = FALSE)


# missing pm
a = births[births$na_days > 20, ]
hist(births$na_days)
missingdata <- a %>% dplyr::select(longitude, latitude, VS_unique, census_tract, pm25mean, na_days)
write.csv(missingdata, file = '/Volumes/Padlock/covid/many_nas.csv')

timeline_plot(under5)

b = read.csv('/Volumes/Padlock/covid/data/debug_low_values.csv')
hist(b$layer)
ndays = time_length(as.Date(max(b$date)) - as.Date(min(b$date)), unit = "days")
max(b$date)
nrow(b[b$layer > 10, ])


# san mateo's county sensor id = 60811001
sm_epa = epa[epa$Site.ID == 60811001, ]
hist(sm_epa$Daily.Mean.PM2.5.Concentration)
sm_epahigh = sm_epa[sm_epa$Daily.Mean.PM2.5.Concentration > 10, ]

# under5 counties
table(under5$county)
under5_county = data.frame(table(under5$county))
colnames(under5_county) = c('county', 'under5_count')

total_counties = data.frame(table(births$county))
colnames(total_counties) = c('county', 'total_count')

under5_county<- merge(under5_county, total_counties, all = TRUE)
under5_county <- under5_county %>% mutate(ratio = under5_count/total_count, percent = ratio*100)
under5_county[is.na(under5_county)] <- 0
hist(under5_county$ratio)
write.csv(under5_county, file = output_low('mean_exp_by_county'), row.names = FALSE)


epa_low <- epa %>% mutate(day = as.Date(Date, format = '%m/%d/%Y'))
epa_avg <- epa_low %>% group_by(COUNTY_CODE, day) %>% dplyr::summarize(county_avg = mean(Daily.Mean.PM2.5.Concentration))
epa_low_avg  <- epa_avg %>% filter(county_avg<5) %>% group_by(COUNTY_CODE) %>% dplyr::summarize(lowdays = n())
epa_low_avg <- epa_low_avg %>% mutate(lowdays_pct = lowdays*100/951, 
                                      county_num = paste('6', str_pad(COUNTY_CODE, 3, side = 'left', pad = '0'), sep = ''))

write.csv(epa_low_avg, file = output_low('lowpm25_epa_values'), row.names = FALSE)
write.csv(epa_avg, file = output_low('epa_avg'), row.names = FALSE)

# duration of EPA file
time_length(max(epa_low$day) - min(epa_low$day), unit = "days")

county_epa_avg <- epa_low %>% group_by(COUNTY_CODE) %>% dplyr::summarise(county_fullavg = mean(Daily.Mean.PM2.5.Concentration))
county_epa_avg <- county_epa_avg %>% mutate(county_fips = paste('06', str_pad(COUNTY_CODE, 3, side = 'left', pad = '0'), sep = ''))
write.csv(county_epa_avg, file = output_low('epa_county_average'), row.names = FALSE)
county_ind <- merge(county_epa_avg, under5_county, by.x = "county_fips", by.y = 'county')

county_ind <- county_ind %>% mutate(less100 = case_when(total_count < 100 ~ TRUE, TRUE ~ FALSE))

county_ind_scatter = ggplot(county_ind, aes(x = county_fullavg, y = percent, colour = less100)) + 
  geom_point() + 
  xlab('County average PM2.5 (EPA)') +
  ylab('% pregnancies <5ug/m3') +
  ggtitle('% pregnancies with low PM2.5 by county EPA average') +
  labs(colour = "Less than 100 births") +
  theme_minimal()
county_ind_scatter
ggsave(county_ind_scatter, file = image_low('county_ind_scatter'), height = 5, width = 7)
                            