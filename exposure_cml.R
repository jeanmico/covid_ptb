# this R script exists to calculate exposure for every pregnancy in our dataset.
# exposure is cumulative
# for now, filter out pregnancies which last into 2021 since we don't have that data yet

# read in the pregnancies (ID, longitude, latitude, start, end)
# spread (pivot_longer) over the full date range (by day) of the pregnancy
# read in the air data into one enormous datatable
# do a left join (or datatable equivalent)
# summarize by summing over each pregnancy

# figure out how to test this on a subset first; it's going to be horrible
library(lubridate)
library(dplyr)
library(data.table)
print(Sys.time())
files = list.files("/c4/home/jcostello/covid/output/values_extracted1/", full.names = TRUE)
#air_vals <- lapply(files[1:3], fread)  # this works, but each file is a separate element of dt
#air_vals <- rbind(lapply(files[1:5], fread, header=TRUE))
air_vals <- do.call(rbind, lapply(files[1:3], fread, header=TRUE))
print(paste('air data read', Sys.time()))
setnames(air_vals, new = c('longitude', 'latitude', 'layer', 'date'))

# now we have a data table (air values) with a few columns: longitude, latitude, layer, date

# read in the births data, simplify it to the columns we need (ID, longitude, latitude, start, end)
# spread/pivot the births data over the date range (by day)
births <- read.csv('/c4/home/jcostello/covid/data/births_exposure_calc.csv')
births_mini <- sample_n(births, 500)
print(births_mini %>% filter(is.na(start_date)|is.na(baby_DOB)))
print(nrow(births_mini))
births_exp <- births_mini %>% 
   rowwise() %>%
   do(data.frame(VS_unique = .$VS_unique, longitude = .$longitude, latitude = .$latitude, start_date = .$start_date, baby_DOB =.$baby_DOB, day_date = seq(ymd(.$start_date), ymd(.$baby_DOB), by = 'day')))
print(paste('births data read and spread', Sys.time()))
setDT(births_exp)

# do a simple join
setkey(births_exp, longitude, latitude, day_date)
setkey(air_vals, longitude, latitude, date)
births_exp[air_vals, pm25 := i.layer]
write.csv(births_exp, file = 'testing.csv')
print(paste('join performed', Sys.time()))

# summarize over the pregnancy
births_exp[is.na(pm25), pm25:=0]
births_summed = births_exp[, sum(pm25), by=VS_unique]

print(paste('pregnancies summed', Sys.time()))
# save the output
write.csv(births_summed, file  = 'testing_sum.csv')
print(paste('output written', Sys.time()))
