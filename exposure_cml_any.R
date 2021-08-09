# this R script exists to calculate exposure for every pregnancy in our dataset.
# it's called from exposure_cml.sh
# exposure is cumulative

# read in the pregnancies (ID, longitude, latitude, start, end)
# spread (pivot_longer) over the full date range (by day) of the pregnancy
# read in the air data into one enormous datatable
# do a left join (or datatable equivalent)
# summarize by summing over each pregnancy

library(lubridate)
library(dplyr)
library(data.table)
files = list.files("/c4/home/jcostello/covid/output/values_extracted_any/", full.names = TRUE)
air_vals <- do.call(rbind, lapply(files, fread, header=TRUE))
print(paste('air data read', Sys.time()))
setnames(air_vals, new = c('longitude', 'latitude', 'layer', 'date'))

# now we have a data table (air values) with a few columns: longitude, latitude, layer, date

# read in the births data, simplify it to the columns we need (ID, longitude, latitude, start, end)
# spread/pivot the births data over the date range (by day)
births <- read.csv('/c4/home/jcostello/covid/data/births_exposure_calc.csv')
#births_mini <- sample_n(births, 500)
print(paste('births data read', Sys.time()))
births_exp <- births %>% 
   rowwise() %>%
   do(data.frame(VS_unique = .$VS_unique, longitude = .$longitude, latitude = .$latitude, start_date = .$start_date, baby_DOB =.$baby_DOB, day_date = seq(ymd(.$start_date), ymd(.$baby_DOB), by = 'day')))
print(paste('births data read and spread', Sys.time()))
setDT(births_exp)

# do a simple join
setkey(births_exp, longitude, latitude, day_date)
setkey(air_vals, longitude, latitude, date)
births_exp[air_vals, pm25 := i.layer]
write.csv(births_exp, file = 'exposure_all.csv')
print(paste('join performed', Sys.time()))

# summarize over the pregnancy
births_exp[is.na(pm25), pm25:=0]
#births_summed = births_exp[, sum(pm25), by=VS_unique]

# also add a count of days above a PM2.5 threshold of 100
births_counts = births_exp[, sum(pm25>=100), by=VS_unique]
births_full = births_exp[, sum(pm25==0), by=VS_unique]
print(paste('pregnancies summed', Sys.time()))
# save the output
#write.csv(births_summed, file  = 'exposure_sum_all.csv')
#write.csv(births_counts, file = 'exposure_counts_all.csv')
#write.csv(births_full, file = 'exposure_na_all.csv')
print(paste('output written', Sys.time()))

# write example low values for debug
debug_ex = air_vals %>% filter(longitude == -122.458707 & latitude == 37.686516)
write.csv(debug_ex, file = 'debug_low_values.csv')
