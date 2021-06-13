## LIBRARIES #################

library(tidyverse)
library(lubridate)
library(tableone)
library(epitools)
library(withr)
library(survival)
library(caTools)

## CONSTANTS AND GLOBALS ###########

# LA county closed on Mar 16 2020
covid_start_date <<- ymd('2020-04-16')

# filepaths
fp_base <<- '/Volumes/Padlock/covid/'
fp_data <<- paste(fp_base, 'data/', sep = '')
fp_out <<- paste(fp_base, 'output/', sep = '')
fp_img <<- paste(fp_base, 'images/', sep = '')


births_fp = paste(fp_data, 'Deb_COVID_set_022621.csv', sep = '')
decodes_fp = paste(fp_data, 'colnames_raw.csv', sep = '')
acs_fp = paste(fp_data, 'acs_data_5yr_2019.csv', sep = '')

## FUNCTIONS #############

# functions for consistent naming; prevent overwriting with Sys.Date()
image_name <- function(fname) {
  return(paste(fp_img, Sys.Date(), '_', fname, '.png', sep = ''))
}

output_name <- function(fname) {
  return(paste(fp_out, Sys.Date(), '_', fname, '.csv', sep = ''))
}

timeline_plot <- function(df) {
  births_dateplt <- df %>% arrange(desc(start_date)) %>% mutate(date_order = row_number())
  
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
  
  ggsave(timeplt, file = image_name('timeline_plot'))
}


## START #################

### process files #####
births_raw = read.csv(births_fp)
decodes = read.csv(decodes_fp)
acs = read.csv(acs_fp)

# process the data, factorize where needed, determine sample sizes
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

### create factor columns ############
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
births[is.na(births)] <- 0

### integrate ACS data ########
# add string variable for census tract id
#births <- births %>% mutate(census_tract = str_sub(with_options(c(scipen=999), 
#                                                                str_pad(census_block, 15, "left", pad = "0")), 
#                                                   1,  11))

# add string for ACS census id and pivot wider
#acs <- acs %>% mutate(census_tract = with_options(c(scipen=999), str_pad(GEOID, 11, "left", pad = "0"))) %>% dplyr::select(-c(variable, codename))
#acs <- acs %>% pivot_wider(names_from = rename, names_sep = '_', values_from = c(estimate, moe))

#births_tmp <- births
#births_tmp <- merge(births, acs, by = 'census_tract')

### determine pregnancy start date ##########
births <- births %>% mutate(baby_DOB = mdy(baby_DOB)) %>%
  mutate(start_date = ymd(baby_DOB) - weeks(gest_weeks))

births <- births %>% mutate(tri2_start = start_date + make_difftime(weeks = 12), tri3_start = tri2_start + make_difftime(weeks = 12))

births <- births %>% mutate(covid_start_in = case_when(
  covid_start_date < start_date ~ 'tri_0',
  covid_start_date < tri2_start ~ 'tri_1',
  covid_start_date < tri3_start ~ 'tri_2',
  covid_start_date < baby_DOB ~ 'tri_3',
  TRUE ~ 'post_partum'))


### write a coordinates file
exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB)
write.csv(exposure_calc, file = '/Volumes/Padlock/covid/data/births_exposure_calc.csv', row.names = FALSE)


### table 1 ###############

mytab <- CreateTableOne(vars = c(fact_cols), data = births, factorVars = fact_cols)
tbl1print = print(mytab)
write.csv(tbl1print, file = output_name('tbl1print'))

## TIMELINE PLOT #############

timeline_plot(births)


## ASSIGN EXPOSURE ########
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

set.seed(100)
library(stringi)
newids = stri_rand_strings(nrow(pocw), 8)
pocw <- pocw %>% mutate(new_id = newids)
pocw_scrambled = pocw %>% dplyr::select(-VS_unique)
write.csv(pocw_scrambled, file = '~/covid/data/exposures_by_POC.csv', row.names = FALSE)


births_exposure <- merge(births, pocany, by = "VS_unique")
births_exposure <- births_exposure %>% 
  filter(pm25 > 0 & start_date >= ymd('2020-01-01') & baby_DOB <= ymd('2020-12-31')) %>%
  mutate(pm25mean = pm25/time_length(baby_DOB - start_date, unit = "days"))

exposure_plot <- ggplot(data = births_exposure, aes(x = pm25mean, y = ..density..)) + 
  geom_density() +
  theme_minimal() +
  xlab('Mean PM2.5 Exposure')
exposure_plot
ggsave(exposure_plot, file = '/Volumes/Padlock/covid/images/20200610_meanpm25_density.png', dpi = 300)


## ANALYSIS #######
## STEP 1: test the effects of air pollution on the outcome of PTB

# start with a simple table

s1_train = sample_frac(births_exposure, .75)
s1_id <- as.numeric(rownames(s1_train))
s1_test <- births_exposure[-s1_id]

s1_logit <- glm(preterm ~ pm25mean, family = binomial('logit'), data = s1_train)
s1_logit_full <- glm(preterm ~ race + female + pm25mean, family = binomial('logit'), data = s1_train)





