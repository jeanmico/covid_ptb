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

ggsave(timeplt, file = 'timeplt.png')


tbl1print = print(mytab)
write.csv(tbl1print, file = 'tbl1print.csv')

points_file <- births %>% dplyr::select(longitude, latitude) %>% distinct()
write.csv(points_file, file = '/Volumes/Padlock/covid/lookup_points.csv', row.names = FALSE)



exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB)
exposure_calc <- exposure_calc %>% filter(!is.na(start_date) & !is.na(baby_DOB) & !is.na(latitude) & !is.na(longitude))
write.csv(exposure_calc, file = '/Volumes/Padlock/covid/data/births_exposure_calc.csv', row.names = FALSE)


# read and compare exposures with POC == 1 and any POC
poc1 = read.csv('/Volumes/Padlock/covid/data/cumulative_exposure_POC1_pass1.csv')
#cut out extreme and missing values for initial comparison; follow up with these later
poc1 <- poc1 %>% mutate(poc_type = 'POC_1') %>% filter(V1 > 200 & V1 < 7500) %>% rename(pm25 = V1)
pocany = read.csv('/Volumes/Padlock/covid/data/cumulative_exposure_POCany_pass1.csv')
pocany <- pocany %>% mutate(poc_type = "POC_any") %>% filter(V1 > 200 & V1 < 7500) %>% rename(pm25 = V1)

poc = rbind(poc1, pocany)

pochist <- ggplot() + 
  geom_density(data = poc, aes(x = pm25, group = poc_type, fill = poc_type, color = poc_type), alpha = .5) +
  scale_color_manual(values = c('#009d9a', '#b28600')) + 
  scale_fill_manual(values = c('#009d9a', '#b28600')) +
  theme_minimal() +
  xlab("Cumulative PM2.5 Exposure")
pochist

ggsave(pochist, file = '/Volumes/Padlock/covid/images/20210609_poc1_vs_any_hist.png', dpi = 300)

poc <- poc %>% dplyr::select(-X)
pocw <- poc %>% pivot_wider(id_cols = VS_unique, names_from = poc_type, values_from = pm25) %>%
  mutate(pmdiff_abs = abs(POC_1 - POC_any))
print(mean(pocw$pmdiff_abs, na.rm = TRUE))
print(sd(pocw$pmdiff_abs, na.rm = TRUE))
print(max(abs(pocw$pmdiff_abs), na.rm = TRUE))

write.csv(pocw, file= '/Volumes/Padlock/covid/data/20210609_exposures_by_POC.csv', row.names = FALSE)


