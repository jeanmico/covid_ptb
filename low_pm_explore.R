# investigating low air pollution values
library(data.table)

low_pm = births
hist(low_pm$pm25mean)
range(low_pm$pm25mean)
which.min(low_pm$pm25mean)

under5 = low_pm[low_pm$pm25mean <= 5,]
hist(under5$pm25mean)

under5write = under5 %>% dplyr::select(longitude, latitude, VS_unique, census_tract, pm25mean)

under5_rm <- under5 %>% filter(na_days <20)
write.csv(under5write, file = '/Volumes/Padlock/covid/lowpm25_explore.csv', row.names = FALSE)


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
write.csv(under5_county, file = '/Volumes/Padlock/covid/lowpm25_counties.csv', row.names = FALSE)


# how many days did each county's EPA monitor register under 5?
epa_low = epa %>% filter(Daily.Mean.PM2.5.Concentration<5) %>% group_by(COUNTY_CODE) %>% dplyr::summarize(n())
View(epa_low)


epa_low <- epa %>% mutate(day = as.Date(Date, format = '%m/%d/%Y'))
epa_avg <- epa_low %>% group_by(COUNTY_CODE, day) %>% dplyr::summarize(county_avg = mean(Daily.Mean.PM2.5.Concentration))
epa_low_avg  <- epa_avg %>% filter(county_avg<5) %>% group_by(COUNTY_CODE) %>% dplyr::summarize(lowdays = n())
epa_low_avg <- epa_low_avg %>% mutate(lowdays_pct = lowdays*100/365, 
                                      county_num = paste('6', str_pad(COUNTY_CODE, 3, side = 'left', pad = '0'), sep = ''))

write.csv(epa_low_avg, file = '/Volumes/Padlock/covid/lowpm25_epa_values.csv', row.names = FALSE)

# duration of EPA file
time_length(max(epa_low$day) - min(epa_low$day), unit = "days")
