library(tidyverse)
library(lubridate)
library(tableone)
library(epitools)

setwd('/Volumes/Padlock/covid/')

births_raw = read.csv('Deb_COVID_set_022621.csv')


# process the data, factorize where needed, determine sample sizes

print(colnames(births_raw))
table(births_raw$any_mat_COVID)
table(births_raw$normal_BMI)
unique(births_raw$normal_BMI)

decodes = read.csv('colnames_raw.csv')


dropcols <- decodes %>% filter(keep == 0) %>% dplyr::select(rawname) # drop unneeded columns
births_tmp <- births_raw %>% dplyr::select(-(dropcols$rawname))

keepcols <- decodes%>% filter(keep == 1)
# create a named vector to use with dplyr::rename
name_vec <- keepcols$rawname
names(name_vec) <- keepcols$varname

# rename variables to make things clearer
births_tmp <- births_tmp %>% dplyr::rename(name_vec[name_vec %in% names(births_tmp)])

births_tmp <- births_tmp %>% mutate(female = case_when(male == 0 ~1, TRUE~0))

births <- births_tmp

# condense multiple binary columns into a single factor column
factor_defs = read.csv('factor_defs.csv')
factor_comps = read.csv('factor_components.csv')

fact_cols = factor_defs$factname

births_tmp <- births

for (myfact in fact_cols) {
  fact_default <- filter(factor_defs, factname == myfact)$default
  print(myfact)
  
  births_tmp[myfact] = fact_default
  
  components <- filter(factor_comps, factname == myfact)$varname
  
  for (colname in components)({
    coldata=births_tmp[,colname]
    births_tmp[[myfact]][which(coldata==1)]=colname
  })
  
  mylevels <- c(fact_default, components)
  mylevels = mylevels[!duplicated(mylevels)]
  
  # convert the column to a factor and set the default
  births_tmp[[myfact]] <- factor(births_tmp[[myfact]], levels = mylevels)
}

births <- births_tmp

# generate table one 

mytab <- CreateTableOne(vars = c(fact_cols), data = births, factorVars = fact_cols)

mytab

# exploring sample and sample sizes


# determine start date of pregnancy
births <- births %>% mutate(baby_DOB = mdy(baby_DOB)) %>%
  mutate(start_date = ymd(baby_DOB) - weeks(gest_weeks))

births <- births %>% mutate(tri2_start = start_date + make_difftime(weeks = 12), tri3_start = tri2_start + make_difftime(weeks = 12))

# LA county closed on Mar 16 2020
covid_start_date = ymd('2020-04-16')

births <- births %>% mutate(covid_start_in = case_when(
    covid_start_date < start_date ~ 'tri_0',
    covid_start_date < tri2_start ~ 'tri_1',
    covid_start_date < tri3_start ~ 'tri_2',
    covid_start_date < baby_DOB ~ 'tri_3',
  TRUE ~ 'post_partum'))

births_dateplt <- births %>% arrange(desc(start_date)) %>% mutate(date_order = row_number())

timeplt <- ggplot(births_dateplt, aes(as.Date(start_date), date_order, color = covid_start_in)) +
  geom_crossbar(aes(xmin = as.Date(start_date), xmax = as.Date(baby_DOB)), width = 0.2) +
  scale_color_manual(values=c( '#a2191f', '#606060', '#808080', '#E0E0E0')) +
  scale_x_date() + 
  xlab('Date') + 
  ylab('Pregnancy') +
  geom_vline(xintercept = covid_start_date) + 
  #guides(colour = guide_legend(override.aes = list(linetype=c(1,1,1,1), size = 2.2))) + # can't seem to override legend symbol
  theme(axis.ticks.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.background=element_blank(), 
        panel.background=element_rect(fill='white'), 
        panel.grid = element_line(color='grey', size=.1)) 


timeplt

ggsave(timeplt, file = 'timeplt.png')


births_postcovid = births %>% filter(covid_start_in == 'tri_0')

births_postcovid$preterm[is.na(births_postcovid$preterm)] <- 0
a = epitable(births_postcovid$any_mat_covid, births_postcovid$preterm)
epitab(a, method = 'oddsratio')

mytab <- CreateTableOne(vars = c(fact_cols, 'wic'), data = births, factorVars = fact_cols)

mytab

tbl1print = print(mytab)
write.csv(tbl1print, file = 'tbl1print.csv')

points_file <- births %>% dplyr::select(longitude, latitude) %>% distinct()
write.csv(points_file, file = '/Volumes/Padlock/covid/lookup_points.csv', row.names = FALSE)



# read in raster value extracts
filenames = list.files('/Volumes/Padlock/covid/values_month_summed', full.names = TRUE)
extract_vals = do.call('cbind', lapply(filenames, read.csv, header = TRUE))
colnames(extract_vals) = c('m1', 'm2', 'm3', 'm4', 'm5', 'm6', 'm7', 'm8', 'm9', 'm10', 'm11', 'm12')
extract_vals = cbind(points_file, extract_vals)

# compute CUMULATIVE exposure for each pregnancy
# calculate complete months
exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB)

# code months as no pregnant days (0), some pregnant days, all pregnant days
#exposure_calc <- exposure_calc %>% mutate(month_cncp = paste('m', month(start_date), sep = ''),
 #                        month_deliv = paste('m', month(baby_DOB), sep = ''),
  #                       full_cncp = case_when(day(start_date) == 1 ~ TRUE, TRUE ~ FALSE),
   #                      full_deliv = case_when(baby_DOB == ceiling_date(baby_DOB, 'month') - days(1) ~ TRUE, TRUE ~ FALSE))

earliest_date = ymd('2020-01-01')
latest_date = ymd('2021-01-31')

#date_cols = seq.Date(earliest_date, latest_date, by = 'month')
#date_cols = format(date_cols, '%Y-%m')
#date_cols = paste('X', date_cols, sep = '')

 #exposure_calc[,date_cols] = NA
 
 
 #mini_exp = sample_n(exposure_calc, 20000)
 
 # THIS IS SLOWWWWW
# test <- exposure_calc %>% 
#   rowwise() %>%
#   do(data.frame(VS_unique = .$VS_unique, longitude = .$longitude, latitude = .$latitude, start_date = .$start_date, baby_DOB =.$baby_DOB,
#                 month = seq(earliest_date, latest_date, by = 'month')))

exposure_calc <- exposure_calc %>% filter(!is.na(start_date) & !is.na(baby_DOB) & !is.na(latitude) & !is.na(longitude))
#mini_exp <- mini_exp %>% filter(!is.na(start_date))
test1 <- exposure_calc %>% 
   rowwise() %>%
   do(data.frame(VS_unique = .$VS_unique, longitude = .$longitude, latitude = .$latitude, start_date = .$start_date, baby_DOB =.$baby_DOB,
                 month = seq(floor_date(.$start_date, "month"), ceiling_date(.$baby_DOB, 'month') - days(1), by = 'month')))

write.csv(test1, file = '/Volumes/Padlock/covid/20210602exposure_calc.csv', row.names = FALSE)

exposure_calc = test1

test <- mini_exp

extract_values <- extract_vals %>% pivot_longer(cols = starts_with("m"),
                                                names_to = "month",
                                                names_prefix = "m",
                                                values_to = "values")

extract_values <- extract_values %>% mutate(year_month = ym(paste('2020', month , sep = '')))

#mini_exp = exposure_calc[sample(nrow(exposure_calc), 500), ]
exposure_calc <- exposure_calc %>% dplyr::rename(year_month = month)

t <- left_join(exposure_calc, extract_values)
