# This is the main code for the COVID project.
#   exposure is computed elsewhere.
# This code
#   reads data files for births, exposure, ACS data
#   processes data
#   performs analyses
#   generates tables, plots

## LIBRARIES #################

library(tidyverse)
library(lubridate)
library(tableone)
library(epitools)
library(withr)
library(survival)
library(caTools)
library(Hmisc)
library(survminer)
library(ggcorrplot)
library(cowplot)
library(janitor)
library(broom)
library(ggplot2)
library(ggstance)
library(ggfortify)
library(splines)

## CONSTANTS AND GLOBALS ###########
set.seed(123)

# LA county closed on Mar 16 2020
covid_start_date <<- ymd('2020-03-16')

# filepaths
fp_base <<- '/Volumes/Padlock/covid/'
fp_data <<- paste(fp_base, 'data/', sep = '')
fp_out <<- paste(fp_base, 'output/', sep = '')
fp_img <<- paste(fp_base, 'images/', sep = '')

# data filenames
births_fp = paste(fp_data, 'Deb_COVID_set_022621.csv', sep = '')
census_tracts_fp = paste(fp_data, 'fix_census_062421.csv', sep = '')
decodes_fp = paste(fp_data, 'colnames_raw.csv', sep = '')
acs_fp = paste(fp_data, 'acs_data_5yr_2019.csv', sep = '')
factor_defs_fp = paste(fp_data, 'factor_defs.csv', sep = '')
factor_comps_fp = paste(fp_data, 'factor_components.csv', sep = '')
coldict_fp = paste(fp_data, 'columns_dict.csv', sep = '')
pocany_fp = paste(fp_data, 'cumulative_exposure_POCany_pass2.csv', sep = '')
badair_days_fp = paste(fp_data, 'exposure_counts_all.csv', sep = '')
pm25_na_fp = paste(fp_data, 'exposure_na_all.csv', sep = '')
## FUNCTIONS #############

# functions for consistently naming output files; prevent overwriting with Sys.Date()
image_name <- function(fname) {
  return(paste(fp_img, Sys.Date(), '_', fname, '.png', sep = ''))
}

output_name <- function(fname) {
  return(paste(fp_out, Sys.Date(), '_', fname, '.csv', sep = ''))
}

factorize_data <- function(factors, df){
  fact_cols = factors$factname
  
  for (myfact in fact_cols) {
    fact_default <- filter(factors, factname == myfact)$default
    print(myfact)
    
    df[myfact] = fact_default
    
    components <- filter(factor_comps, factname == myfact)$varname
    
    for (colname in components)({
      coldata=df[,colname]
      df[[myfact]][which(coldata==1)]=colname
    })
    
    mylevels <- c(fact_default, components)
    mylevels = mylevels[!duplicated(mylevels)]
    
    # convert the column to a factor and set the default
    df[[myfact]] <- factor(df[[myfact]], levels = mylevels)
  }
  return(df)
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


cox_coxph <- function(df) {
  coxmodel = coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
          ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
          strata(race, insurance, pn_care), data = df)
  return(coxmodel)
  }

hazard_ratios <- function(cox_model, plotname, myvar) {
  cox_termplot = termplot(cox_model, terms = NULL, se=T, plot = F)
  
  cox_df = data.frame(cox_termplot[[myvar]])
  cox_df <- cox_df %>% mutate(hr = exp(y),
                            lci = exp(y - 196*se),
                            uci = exp(y + 1.96*se))
  
  hr_plot <- ggplot(cox_df, aes(x = x, y = hr)) + 
    geom_line(size = 1.5) + 
    geom_line(data=cox_df, aes(x = x, y = lci), color = '#009d9a') + 
    geom_line(data=cox_df, aes(x = x, y = uci), color = '#009d9a') +
    geom_hline(yintercept=1, linetype = 'longdash') +
    xlab("Mean PM2.5") + 
    ylab(paste("Hazard ratio for",  myvar, "using psplines")) +
    theme_bw()
  
  hr_plot_log <- ggplot(cox_df, aes(x = x, y = hr)) + 
    geom_line(size = 1.5) + 
    geom_line(data=cox_df, aes(x = x, y = lci), color = '#009d9a') + 
    geom_line(data=cox_df, aes(x = x, y = uci), color = '#009d9a') +
    geom_hline(yintercept=1, linetype = 'longdash') +
    xlab("Mean PM2.5") + 
    ylab(paste("Log Hazard ratio for",  myvar, "using psplines")) +
    scale_y_log10() +
    theme_bw()
  
  hr_grid <- plot_grid(hr_plot, hr_plot_log)
    
  ggsave(hr_plot, file = image_name(paste('hazard_ratio_', plotname, '_', myvar, sep = '')), dpi = 300)
  
  return(cox_df)
  }

hazard_forest <- function(hr_df, plotname, filename){
  
  renames = read.csv(paste(fp_data, 'hazard_renames.csv', sep = ''))
  
  hr_df <- merge(hr_df, renames)
  hr_fp <- ggplot(hr_df, aes(y = rename, x= var_exp, xmin = lower, xmax = upper)) + 
    geom_linerangeh(size = .5) + 
    geom_point(size = 1) +
    scale_x_log10() +
    xlab('exp(coeff) (log scale)') + 
    ggtitle(plotname) + 
    theme_bw()
  ggsave(hr_fp, file = image_name(filename), dpi = 300)
  
  return(hr_fp)
  }

cox_cph <- function(df) {}

## START #################

### read data files #####
births_raw = read.csv(births_fp)  # births data
census_tracts = read.csv(census_tracts_fp)  # corrected census tracts
acs = read.csv(acs_fp)  # data from American Community Survey 2015-19
pocany = read.csv(pocany_fp)  # PM2.5 exposure data
badair_days = read.csv(badair_days_fp)  # how many days was PM2.5 > 100ug/m3
pm25_na = read.csv(pm25_na_fp)

### read data dictionaries #####
factor_defs = read.csv(factor_defs_fp)  # factor definitions (names, missing, default)
factor_comps = read.csv(factor_comps_fp)  # factor components
decodes = read.csv(decodes_fp)  # column renames
coldict = read.csv(coldict_fp)  # column dictionary

### process files ####
births_raw <- births_raw %>% dplyr::select(-PGB_CENSUS_BLOCK_193)  # remove truncated census tract variable
births_raw <- merge(births_raw, census_tracts)  # add correct census tract

# process the data, factorize where needed, determine sample sizes
dropcols <- decodes %>% filter(keep == 0) %>% dplyr::select(rawname) # drop unneeded columns
births <- births_raw %>% dplyr::select(-(dropcols$rawname))
keepcols <- decodes%>% filter(keep == 1)

# rename variables for readability
name_vec <- keepcols$rawname  # create a named vector to use with dplyr::rename
names(name_vec) <- keepcols$varname
births <- births %>% dplyr::rename(name_vec[name_vec %in% names(births)])

births <- births %>% mutate(female = case_when(male == 0 ~1, TRUE~0))  # add missing sex variable


### create factor columns ############
# condense multiple binary columns into a single factor column
births_tmp <- factorize_data(factor_defs, births)

# collapse any factors we wish to combine
births_tmp$insurance = fct_collapse(births_tmp$insurance, all_other_pay = c("all_other_pay", "self_pay"))
births <- births_tmp

births[is.na(births)] <- 0  # this is how the data are coded 

births <- births %>% filter(singleton == 1)  # select only singleton births

### integrate ACS data ########
# add string variable for census tract id
births <- births %>% mutate(census_tract = str_sub(with_options(c(scipen=999), 
                                                                str_pad(census_block, 15, "left", pad = "0")), 
                                                  1,  11))

# add string for ACS census id and pivot wider
acs <- acs %>% mutate(census_tract = with_options(c(scipen=999), str_pad(GEOID, 11, "left", pad = "0"))) %>% dplyr::select(-c(variable, codename))
acs <- acs %>% pivot_wider(names_from = rename, names_sep = '_', values_from = c(estimate, moe))

births <- births %>% mutate(county = str_sub(with_options(c(scipen=999), 
                                                          str_pad(census_block, 15, "left", pad = "0")), 
                                             1,  5))

births <- merge(births, acs, by = 'census_tract')


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

### integrate exposure data #####
colnames(pm25_na) = c('X', 'VS_unique', 'na_days')
births <- merge(births, pm25_na, by = 'VS_unique')

pocany <- pocany %>% mutate(poc_type = "POC_any") %>% rename(pm25 = V1)
births <- merge(births, pocany, by = "VS_unique")
births <- births %>% 
  mutate(pm25mean = pm25/(time_length(baby_DOB - start_date, unit = "days") - na_days))

# how many bad air days did eavh person experience?
colnames(badair_days) = c("X", "VS_unique", "badair_count")
births <- merge(births, badair_days, by = "VS_unique")
births <- births %>% mutate(badair_binary = case_when(badair_count > 0 ~ 1, TRUE ~ 0))



# Start date of pregnancy
births <- births %>% mutate(year_mth_start = str_sub(start_date, 1, 7))
births <- births %>% mutate(year_mth_fact = as.factor(year_mth_start))

# to use the start year and month of the pregnancy, use a spline
datelist = sort(unique(births$year_mth_start))
births <- births %>% mutate(year_mth_ord = match(year_mth_start, datelist))

# if we choose to right-censor: calculate the time at event; if not preterm, use 37wks = 259 days
births <- births %>% mutate(event_time = case_when(
  gest_weeks > 37 ~ 37,
  TRUE ~ gest_weeks
))


### FILTER BY MONTH??
births <- births %>% filter(year_mth_ord <= 9) # 9 was chosen because it's the latest month with more than 2000 births

## END OF VARIABLE CREATION AND SAMPLE SELECTION ###

### write a coordinates file
exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB)
write.csv(exposure_calc, file = output_name('exposure_coordinates'), row.names = FALSE)


# write a file of the processed columns
write.csv(colnames(births), "/Volumes/Padlock/covid/data/columns_processed.csv", row.names = FALSE)




### TABLE 1 ###############

mytab <- CreateTableOne(vars = c(fact_cols), data = births, factorVars = fact_cols)
tbl1print = print(mytab)
write.csv(tbl1print, file = output_name('tbl1print'))

## PLOTS PT1 #############

timeline_plot(births)

exposure_plot <- ggplot(data = births, aes(x = pm25mean, y = ..density..)) + 
  geom_density() +
  theme_minimal() +
  xlab('Mean PM2.5 Exposure')
exposure_plot
ggsave(exposure_plot, file = image_name('mean_pm25_density_fullsample'), dpi = 300)


# compute mean air exposure by county and write output to generate a figure
pm25_county_mean <- births %>% group_by(county) %>%
  summarise(pm25county_mean = mean(pm25mean),
            pm25county_min = min(pm25mean),
            pm25county_max = max(pm25mean),
            pm25county_sd = sd(pm25mean),
            county_births_count = n())
write.csv(pm25_county_mean, file = output_name('county_summary'), row.names = FALSE)

# histogram of PM2.5 values to identify extremes
pmhist_breaks = read.csv('/Volumes/Padlock/covid/output/hist_breaks.csv')
pmhist_counts = read.csv('/Volumes/Padlock/covid/output/hist_counts.csv')


pm_freq = data.frame(pmhist_breaks$x[1:nrow(pmhist_breaks) -1], pmhist_counts$x)
colnames(pm_freq) = c('pm_val', 'counts')

pmhist = ggplot(pm_freq, aes(pm_val, counts)) + geom_col(width = 10) + 
  theme_bw() + xlab('PM2.5 Concentration (ug/m3)') + 
  ylab('log(count)')
pmhist
ggsave(pmhist, file = image_name('PM25_histogram'))


# ACS histograms
pov_hist <- ggplot(births, aes(x = estimate_poverty)) + geom_histogram() + 
  xlab("Poverty rate") +
  theme_bw()
pov_hist
ggsave(pov_hist, file = image_name('poverty_hist'))

house_hist <- ggplot(births, aes(x = estimate_household_income)) + geom_histogram() + 
  xlab("Household income") +
  theme_bw()
house_hist
ggsave(house_hist, file = image_name('housing_hist'))

badair_days_hist <- ggplot(births %>% filter(badair_count > 0), aes(x = badair_count)) + 
  geom_histogram(binwidth = 1) +
  xlab('Days where mean PM2.5 > 100ug/m3') +
  theme_bw()
badair_days_hist
ggsave(badair_days_hist, file = image_name('badair_days_hist'))


births_full = births # we keep this df with all columns for reference


# correlation matrix for ACS variables
births_acs = births %>% dplyr::select(starts_with('estimate'))
data.frame(colSums(is.na(births_acs)))


mycorr <- round(cor(births_acs, use = 'complete.obs'), 2)
write.csv(mycorr, file = output_name('acs_correlation'), row.names = FALSE)

png(filename = image_name('acs_correlation'), width = 600, height = 900)
ggcorrplot(mycorr, type = "lower", colors = c('#b2182b', '#f7f7f7', '#2166ac'), lab = TRUE)
dev.off()



# frequency tables to check sample sizes

frq_pn <- births %>% tabyl(race, pn_care, insurance)

frq_bmi <- births %>% tabyl(race, bmi, insurance)

## ANALYSIS #######

# center and scale continuous variables
births <- births %>% mutate(poverty_scale = scale(estimate_poverty),
                            pm25_scale = pm25mean/5 - mean(pm25mean))

# split test and train
births_train = sample_frac(births, .75)
births_id <- as.numeric(rownames(births_train))
births_test <- births[-births_id]


# basic log-linear model
baselog = glm(preterm ~ pm25_scale, data = births, family = poisson)
exp(coef(baselog))

countlog = glm(preterm ~ badair_count, data = births, family = poisson)
exp(coef(countlog))

partial_log = glm(preterm ~ pm25_scale + ipi + year_mth_ord + mom_age + inf_sex + estimate_poverty + race + insurance + pn_care, data = births, family = poisson)
exp(coef(partial_log))
summary(partial_log)

partialcount_log = glm(preterm ~ badair_count + ipi + year_mth_ord + mom_age + inf_sex + estimate_poverty + race + insurance + pn_care, data = births, family = poisson)
exp(coef(partialcount_log))
summary(partialcount_log)

nointeraction_log = glm(preterm ~ pm25_scale + any_mat_covid, data = births, family=poisson)
exp(coef(nointeraction_log))
summary(nointeraction_log)

baseinteraction_log = glm(preterm ~ pm25_scale + any_mat_covid + pm25_scale*any_mat_covid, data = births, family=poisson)
exp(coef(baseinteraction_log))
summary(baseinteraction_log)

baseinter_quad_log = glm(preterm ~ pm25_scale + pm25_scale^2 + any_mat_covid + pm25_scale^2*any_mat_covid, data = births, family=poisson)
exp(coef(baseinter_quad_log))
summary(baseinter_quad_log)

# Basic Cox Proportional Hazards model, not adjusting for COVID status

### MODEL SELECTION #####
# build models
cox1 <- coxph(Surv(gest_weeks, preterm) ~  pm25mean + year_mth_ord +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox2 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + year_mth_ord +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox3 <- coxph(Surv(gest_weeks, preterm) ~  pm25mean + pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox4 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox5 <- coxph(Surv(gest_weeks, preterm) ~  pm25mean + pm25mean*pm25mean + pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox6 <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + year_mth_ord + 
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox7 <- coxph(Surv(gest_weeks, preterm) ~  pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox8 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty  +
                strata(race, insurance, pn_care), data = births_train)

cox9 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
               ipi + mom_age  + inf_sex  + estimate_poverty  + badair_binary + 
               strata(race, insurance, pn_care), data = births_train)

cox10 <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + year_mth_ord + 
                 ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                 strata(race, insurance, pn_care), data = births_train)

cox11 <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + year_mth_ord + 
                 ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary +
                 strata(race, insurance, pn_care), data = births_train)

cox12 <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + pspline(year_mth_ord, df = 4) + 
                 ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary +
                 strata(race, insurance, pn_care), data = births_train)

cox13 <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + pspline(year_mth_ord, df = 4) + 
                 ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                 strata(race, insurance, pn_care), data = births_train)


coxinter <- coxph(Surv(gest_weeks, preterm) ~  pm25mean*any_mat_covid + pm25mean + any_mat_covid + year_mth_ord +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

coxinter1 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4)*any_mat_covid + pspline(pm25mean, df = 4) + any_mat_covid + 
                    pspline(year_mth_ord) +
                    ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                    strata(race, insurance, pn_care), data = births_train)


model_list = list(cox1, cox2, cox3, cox4, cox5, cox6, cox7, cox8, cox9, cox10, cox11, cox12, cox13)

spline_termplot <- termplot(cox12, term = 1, se=T, plot = T)


# compare models
aic_cox = lapply(model_list, AIC)

min_aic_cox <- which.min(lapply(model_list, AIC))

# comparing PM25 to PM25 spline
exp((aic_cox[[4]] - aic_cox[[3]])/2)

# comparing linear + quadratic to spline
exp((aic_cox[[4]] - aic_cox[[5]])/2)

log_like = lapply(model_list, logLik)

# looking at the models in further detail
cox_base <- cox9

base_summ <- summary(cox_base)
base_coeff <- round(data.frame(base_summ$coefficients), 2)
base_ci <- round(data.frame(base_summ$conf.int), 2)
base_ci <- base_ci %>% rownames_to_column()
colnames(base_ci) <- c('varname', 'var_exp', 'var_nexp', 'lower', 'upper')

write.csv(base_coeff, file = output_name('base_model_coeff'), row.names = TRUE)
write.csv(base_ci, file = output_name('base_model_ci'), row.names  = TRUE)

hazard_ratios(cox_base, 'all', 'pm25mean')
hazard_ratios(cox_base, 'all', 'year_mth_ord')

png(filename = image_name('partial_pm25'))
partials_plots = termplot(cox_base, terms = 1, se=T, plot = T)
dev.off()

png(filename = image_name('partial_start'))
partials_plots = termplot(cox_base, terms = 2, se=T, plot = T)
dev.off()

hr1 = smoothHR(data = births_train, coxfit = cox_base)
plot(hr1, predictor = 'pm25mean')

partials = termplot

sr <- cox.zph(cox_base)
write.csv(round(sr$table, 3), file = output_name('full_model_sr'))



births_train_neg = births_train %>% filter(any_mat_covid == 0)
cox_neg = coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
                  ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                  strata(race, insurance, pn_care), data = births_train)

births_train_pos <- births_train %>% filter(any_mat_covid == 1)
cox_pos = coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
                  ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                  strata(race, insurance, pn_care), data = births_train_pos)

neg_summ <- summary(cox_neg)
neg_coeff <- round(data.frame(neg_summ$coefficients), 2)
neg_ci <- round(data.frame(neg_summ$conf.int), 2)
neg_ci <- neg_ci %>% rownames_to_column()
colnames(neg_ci) <- c('varname', 'var_exp', 'var_nexp', 'lower', 'upper')

write.csv(neg_coeff, file = output_name('neg_model_coeff'), row.names = TRUE)
write.csv(neg_ci, file = output_name('neg_model_ci'), row.names  = TRUE)

hazard_ratios(cox_neg, 'neg', 'pm25mean')
hazard_ratios(cox_neg, 'neg', 'year_mth_ord')

sr <- cox.zph(cox_neg)
write.csv(round(sr$table, 3), file = output_name('neg_model_sr'))


pos_summ <- summary(cox_pos)
pos_coeff <- round(data.frame(pos_summ$coefficients), 2)
pos_ci <- round(data.frame(pos_summ$conf.int), 2)
pos_ci <- pos_ci %>% rownames_to_column()
colnames(pos_ci) <- c('varname', 'var_exp', 'var_nexp', 'lower', 'upper')

write.csv(pos_coeff, file = output_name('pos_model_coeff'), row.names = TRUE)
write.csv(pos_ci, file = output_name('pos_model_ci'), row.names  = TRUE)

hazard_ratios(cox_pos, 'pos', 'pm25mean')
hazard_ratios(cox_pos, 'pos', 'year_mth_ord')

sr <- cox.zph(cox_pos)
write.csv(round(sr$table, 3), file = output_name('pos_model_sr'))

base_fit <- survfit(cox_base)

hazard_forest(base_ci, 'all', 'all')
hazard_forest(neg_ci, 'negative', 'negative')
hazard_forest(pos_ci, 'positive', 'positive')

autoplot(base_fit)

a = ggsurvplot(base_fit, data = births_train)
a
library(smoothHR)
pm25_smhr = smoothHR(data=births_train, coxfit = cox_base)
plot(pm25_smhr, predictor = 'pm25mean', prob = 0, conf.level = .95, round.x=2)

# Plot the Schoenfeld Residuals
srplt = ggcoxzph(sr, df = 3)  # df >4 leads to an NA/NaN/Inf error

pm25_tp = termplot(cox_base, term = 1, se = T, plot = T)

png(file = image_name('schoenfeld_test'), width = 1000, height = 900)
srplt
dev.off()

srplt
sr

cox_ps = coxph(Surv(gest_weeks, preterm) ~ pspline(pm25mean, df = 4), data = births)
ps_termplot = termplot(cox_base, term = 1, se=T, plot = F)

ps_df = data.frame(ps_termplot$pm25mean)
ps_df <- ps_df %>% mutate(hr = exp(y),
                          lci = exp(y - 196*se),
                          uci = exp(y + 1.96*se))

ps_plot <- ggplot(ps_df, aes(x = x, y = hr)) + 
  geom_line(size = 1.5) + 
  geom_line(data=ps_df, aes(x = x, y = lci), color = '#009d9a', linetype = 'dashed') + 
  geom_line(data=ps_df, aes(x = x, y = uci), color = '#009d9a', linetype = 'dashed') +
  geom_hline(yintercept=1, linetype = 'longdash') +
  xlab("Mean PM2.5") + 
  ylab("Hazard ratio for PM2.5 using psplines") +
  theme_bw()
ps_plot
ggsave(ps_plot, file= image_name('hazard_psplines_plot_pm25_fullmodel'))
  

library(rms)
surv_upd <- Surv(time = births$gest_weeks, event = births$preterm)
cox_upd = cph(surv_upd ~ ipi  + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strat(insurance) + strat(race) + strat(pn_care), 
              , data = births 
              , x = TRUE, y = TRUE, surv  = TRUE)

ddist <- datadist(births)
options(datadist = "ddist")
summary(cox_upd)
plot(summary(cox_upd))

cox.zph(cox_upd)

plot(cox.zph(cox_upd, "identity"))


cox_full = cph(surv_upd ~ ipi  + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strat(insurance) + strat(race) + strat(pn_care) + rcs(pm25mean, 4) + rcs(year_mth_ord, 4), 
              , data = births 
              , x = TRUE, y = TRUE, surv  = TRUE)
summary(cox_full)
plot(summary(cox_full))
cox.zph(cox_full)

residuals(cox_full)

ggplot(Predict(cox_full, pm25mean))
ggplot(Predict(cox_full, year_mth_ord))

AIC(cox_full) # the minimum AIC is with both pm25mean and year_mth_ord as splines

pm_linear_test = ggcoxfunctional(Surv(gest_weeks, preterm) ~  pm25mean, data = births) 
ggsave(pm_linear_test$pm25mean, file = image_name('PM25_martingale'))

start_linear_test = ggcoxfunctional(Surv(gest_weeks, preterm) ~  year_mth_ord, data = births) 
ggsave(start_linear_test$year_mth_ord, file = image_name('start_martingale'))

## separate by covid status ############
births_nocov <- births %>% filter(any_mat_covid == 0)
births_cov <- births %>% filter(any_mat_covid == 1)

cox_cov = cox_coxph(births_cov)
cox_nocov = cox_coxph(births_nocov)

hazard_ratios(cox_cov, 'covid')
hazard_ratios(cox_nocov, 'nocovid')

cov_sr <- cox.zph(cox_cov)
nocov_sr <- cox.zph(cox_nocov)

ggplot(Predict(cox_nocov, pm25mean))



# plot errors...are there any categories for whom no PTB events occurred?
dr2 <- ggcoxdiagnostics(cox_base, type = 'deviance', linear.predictions = FALSE, se = FALSE, ggtheme = theme_bw())
#dr2 this does not work just processes forever

contingency_pm25 <- epitable(births$pm25_tertile, births$preterm)
pm25_rr = epitab(contingency_pm25, method = 'riskratio')

contingency_covid <- epitable(births$any_mat_covid, births$preterm)
covid_rr <- epitab(contingency_covid, method = 'riskratio')

contingency_pm25bin <- epitable(births$pm25binary, births$preterm)
pm25_rr = epitab(contingency_pm25bin, method = 'riskratio')
