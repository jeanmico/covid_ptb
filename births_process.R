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
library(ggExtra)
library(splines)

## CONSTANTS AND GLOBALS ###########
set.seed(123)

# LA county closed on Mar 16 2020
covid_start_date <<- ymd('2020-03-16')

# filepaths
fp_base <<- '/Volumes/Padlock/covid/'
fp_data <<- paste(fp_base, 'data/', sep = '')
fp_exp <<- paste(fp_base, 'exposure/', sep = '')
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

county_pop_fp = paste(fp_data, 'county_pop_census_2019.csv', sep = '')
county_area_fp = paste(fp_data, 'county_area.csv', sep = '')
county_dict_fp = paste(fp_data, 'california_county_dict.csv', sep = '')

pocany_fp = paste(fp_data, 'exposure_sum_final.csv', sep = '')
badair_days_fp = paste(fp_data, 'exposure_counts_final.csv', sep = '')
pm25_na_fp = paste(fp_data, 'exposure_na_final.csv', sep = '')
## FUNCTIONS #############

# functions for consistently naming output files; prevent overwriting with Sys.Date()
image_name <- function(fname) {
  return(paste(fp_img, Sys.Date(), '_', fname, '.png', sep = ''))
}

output_name <- function(fname) {
  return(paste(fp_out, Sys.Date(), '_', fname, '.csv', sep = ''))
}

births_combine <- function(olddf, newdf) {
  oldcol = colnames(olddf)
  newcol = colnames(newdf)
  collist = list()
  for (i in oldcol) {
    if (!(i %in% newcol))  {collist = c(collist, i)}
  }
  
  print(collist)
  
  colselect = oldcol[!(oldcol %in% collist)]
  
  tmpnew = newdf %>% dplyr::select(all_of(colselect))
  tmpold <- olddf %>% dplyr::select(all_of(colselect))
  
  tmpnew = tmpnew %>% filter(!(VS_unique %in% tmpold$VS_unique))
  
  tmp = rbind(tmpold, tmpnew)
  #[is.na(tmp)] <- 0 
  return(distinct(tmp))
  
  
  return(df) 
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
    scale_color_manual(values=c( '#003061', '#8557a8', '#2d7166', '#963c4b')) +
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

exposure_add <- function(df) {
  df = births
  explist = c('exposure')
  
  for (i in explist) {
    suffix = gsub('exposure', '', i)

    na = paste('na_days', suffix, sep = '')
    bad = paste('badair_count', suffix, sep = '')
    badbin = paste('badair_binary', suffix, sep = '')
    pm = paste('pm25', suffix, sep = '')
    pmmean = paste('pm25mean', suffix, sep = '')
    
    pocany = read.csv(paste(fp_exp, i, '_sum_final.csv', sep = ''))
    badair_days = read.csv(paste(fp_exp, i, '_counts_final.csv', sep = ''))
    pm25_na = read.csv(paste(fp_exp, i, '_na_final.csv', sep = ''))
    
    print(nrow(pocany))
    print(head(pocany))

    # integrate exposure data #####
    colnames(pm25_na) = c('VS_unique', na)
    df <- merge(df, pm25_na, by = 'VS_unique')

    colnames(pocany) = c('VS_unique', pm)
    df <- merge(df, pocany, by = "VS_unique")

    #df <- df %>% 
    #  mutate(!!pmmean := !!pm/(time_length(baby_DOB - start_date, unit = "days") - !!na))

    df[[pmmean]] = df[[pm]]/(time_length(df[["baby_DOB"]] - df[["start_date"]], unit = "days") - df[[na]])

    # how many bad air days did eavh person experience?
    colnames(badair_days) = c('VS_unique', bad)
    df <- merge(df, badair_days, by = "VS_unique")

    df[[badbin]] =  ifelse(df[[bad]] == 0, 0, 1)
      #df[[bad]]/df[[bad]]  # map to 0 or 1
    #df <- df %>% mutate(!!badbin := case_when(!!bad > 0 ~ 1, TRUE ~ 0))
    
    }
  return(df)
  }

rsq <- function(x, y) {
  return(cor(x,y)^2)
}

hazard_ratios <- function(cox_model, plotname, myvar, bdf) {
  cox_termplot = termplot(cox_model, terms = NULL, se=T, plot = F)

  cox_df = data.frame(cox_termplot[[myvar]])
  cox_df <- cox_df %>% mutate(hr = exp(y),
                            lci = exp(y - 1.96*se),
                            uci = exp(y + 1.96*se))
  
  hr_plot <- ggplot(cox_df, aes(x = x, y = hr)) + 
    geom_line(size = 1.5) + 
    geom_line(data=cox_df, aes(x = x, y = lci), color = '#009d9a') + 
    geom_line(data=cox_df, aes(x = x, y = uci), color = '#009d9a') +
    geom_hline(yintercept=1, linetype = 'longdash') +
    geom_point(data = bdf, aes(x = pm25mean, y = 0), alpha = 0) + 
    xlab("Mean PM2.5") + 
    ylab(paste("Hazard ratio for",  myvar, "using a natural spline")) +
    theme_bw()
    
  ggsave(hr_plot, file = image_name(paste('hazard_ratio_', plotname, '_', myvar, sep = '')), dpi = 300)
  print(image_name(paste('hazard_ratio_', plotname, '_', myvar, sep = '')))
  
  hr_margin = ggMarginal(hr_plot, type = 'boxplot', margins = 'x')
  ggsave(hr_margin, file = image_name(paste('hazard_ratio_margin_', plotname, '_', myvar, sep = '')), dpi = 300) 
  
  return(hr_plot)
  }

hazard_ratios_compare <- function(cox_model1, cox_model2, plotname, myvar, lbls) {
  cox_df = data.frame(matrix(nrow = 0, ncol = 3))
  colnames(cox_df) = c('x', 'y', 'se')
  mdls = list(cox_model1, cox_model2)
  for (i in 1:length(mdls)) {
    cox_termplot = termplot(mdls[[i]], terms = NULL, se=T, plot = F)
    cox_data = data.frame(cox_termplot[[myvar]])
    cox_data$mdl = lbls[[i]]
    
    cox_df = rbind(cox_df, cox_data)
  }
  
  
  #return(cox_df)
  cox_df$mdl = as.factor(cox_df$mdl)
  print(levels(cox_df$mdl))
  cox_df <- cox_df %>% mutate(hr = exp(y),
                              lci = exp(y - 1.96*se),
                              uci = exp(y + 1.96*se))
  
  mycolors = c('#5e3c99', '#e66101')
  
  hr_plot <- ggplot(cox_df, aes(x = x, y = hr, color = mdl)) + 
    geom_line(size = 1.5, alpha = .75) + 
    geom_line(data=cox_df, aes(x = x, y = lci, color = mdl), alpha = .75) + 
    geom_line(data=cox_df, aes(x = x, y = uci, color = mdl), alpha = .75) +
    geom_hline(yintercept=1, linetype = 'longdash') +
    xlab(expression('Mean PM'[2.5])) + 
    ylab(paste("Hazard ratio for preterm birth")) +
    theme_bw() + 
    theme(text = element_text(size = 20)) +
    #theme_bw() + 
    labs(color= "COVID Status")  +
    scale_color_manual(values = c("COVID-" = "#5e3c99", "COVID+" = "#e66101"))
  
    #facet_grid(cols = vars(model))
  
  hr_plot
  ggsave(hr_plot, file = image_name(paste('hazard_ratio_compare', plotname, '_', myvar, sep = '')), dpi = 300)
  print(image_name(paste('hazard_ratio_', plotname, '_', myvar, sep = '')))
  return(hr_plot)
}


hazard_forest <- function(hr_df, plotname, fname){
  
  renames = read.csv(paste(fp_data, 'hazard_renames.csv', sep = ''))
  
  hr_df <- merge(hr_df, renames)
  hr_fp <- ggplot(hr_df, aes(y = rename, x= var_exp, xmin = lower, xmax = upper)) + 
    geom_linerangeh(size = .5) + 
    geom_point(size = 1) +
    scale_x_log10() +
    xlab('exp(coeff) (log scale)') + 
    ggtitle(plotname) + 
    theme_bw()
  ggsave(hr_fp, file = image_name(fname), dpi = 300)
  print(image_name(fname))
  }

# termplot function to check the PM2.5 distribution
termp <- function(df, fnamestr) {
  
  
  return(termplot(df, term = 1, se = T, plot = T))
  
}

# Cox model function
#  takes arguments as strings, builds literal function call
cox_cph <- function(lbl, df, varlist, stratalist, splinelist) {
  myvars = paste(varlist, collapse = ' + ')
  stratas = paste(stratalist, collapse = ', ') 
  
  #mysplines = paste(mysplines, '+ ')
  
  if(length(stratalist) > 0) {
    mysplines = paste(splinelist, collapse = ' + ')
    mysplines = paste(mysplines, '+ ')
    mystrata = paste(' + strata(', stratas, ')')
    mycall = as.formula(paste('Surv(gest_weeks, preterm) ~', mysplines, myvars, mystrata))
    
  } else {
    
    mycall = as.formula(paste('Surv(gest_weeks, preterm) ~', splinelist))
  }
  
  print(mycall)
  mdl = coxph(mycall, data  = df)
  
  termp(mdl, fname)
  
  print(splinelist)
  print(pmspline)
  print(grep(pmspline, splinelist))
  if (length(splinelist) > 0) {
    #if (grepl(pmspline, splinelist) | pmspline == splinelist){
      print('plot')
      myplt = hazard_ratios(mdl, lbl, 'pm25mean', df)
    #}
  }
  write.csv(summary(mdl)$coefficients, file = output_name(paste('coeffs', lbl, sep = '')))
  return(list(mdl, myplt))
  }

## START #################

### read data files #####
births_raw = read.csv(births_fp)  # births data
births_add = read.csv(paste(fp_data, 'Jean_all_2020_file_102221.csv', sep = ''))
births_raw = births_combine(births_raw, births_add)

census_tracts = read.csv(census_tracts_fp)  # corrected census tracts
acs = read.csv(acs_fp)  # data from American Community Survey 2015-19
county_pop = read.csv(county_pop_fp)
county_area = read.csv(county_area_fp)
county_dict = read.csv(county_dict_fp)

# read exposure files #
#pocany = read.csv(pocany_fp)  # PM2.5 exposure data
#badair_days = read.csv(badair_days_fp)  # how many days was PM2.5 > 100ug/m3
#pm25_na = read.csv(pm25_na_fp)

### read data dictionaries #####
factor_defs = read.csv(factor_defs_fp)  # factor definitions (names, missing, default)
factor_comps = read.csv(factor_comps_fp)  # factor components
decodes = read.csv(decodes_fp)  # column renames
coldict = read.csv(coldict_fp)  # column dictionary

### process files ####
births_raw <- births_raw %>% dplyr::select(-PGB_CENSUS_BLOCK_193)  # remove truncated census tract variable
births_raw <- merge(births_raw, census_tracts)  # add correct census tract

# process the data, factorize where needed, determine sample sizes
keepcols <- decodes %>% filter(keep == 1) %>% dplyr::select(rawname) # drop unneeded columns
births <- births_raw %>% dplyr::select(keepcols$rawname)
keepcols <- decodes%>% filter(keep == 1)

# rename variables for readability
name_vec <- keepcols$rawname  # create a named vector to use with dplyr::rename
names(name_vec) <- keepcols$varname
births <- births %>% dplyr::rename(name_vec[name_vec %in% names(births)])

births <- births %>% mutate(female = case_when(male == 0 ~1, TRUE~0))  # add missing sex variable

births <- births %>% filter(gest_weeks > 0 & gest_weeks <42)

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


### calculate population density ###
county_area <- county_area[county_area$STATEFP == 6,]
county_density = merge(county_area, county_pop, by.x = 'NAMELSAD', by.y = 'County')
county_density = county_density %>% mutate(pop_density = Population/ALAND) %>% dplyr::select(COUNTYFP, pop_density)
births <- births %>% mutate(county_fips = as.numeric(substr(county, 3, 5)))

births <- merge(births, county_density, by.x = 'county_fips', by.y='COUNTYFP')

### determine pregnancy start date ##########
births <- births %>% mutate(baby_DOB = mdy(baby_DOB)) %>%
  mutate(start_date = ymd(baby_DOB) - weeks(gest_weeks))

births_tmp <- births %>% filter(start_date >= '2019-08-15' & start_date <= '2020-04-30')   # filter by date
births <- births_tmp

births <- births %>% mutate(
  tri1_start = start_date,
  tri2_start = start_date + make_difftime(weeks = 12), 
  tri3_start = tri2_start + make_difftime(weeks = 12),
  tri1_end = tri2_start - make_difftime(days = 1),
  tri2_end = tri3_start - make_difftime(days = 1),
  tri3_end = baby_DOB
  )

### write a coordinates file
exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB,
                                          tri1_start, tri1_end,
                                          tri2_start, tri2_end,
                                          tri3_start, tri3_end)

write.csv(exposure_calc, file = output_name('exposure_coordinates_final'), row.names = FALSE)


births <- births %>% mutate(covid_start_in = case_when(
  covid_start_date < start_date ~ 'tri_0',
  covid_start_date < tri2_start ~ 'tri_1',
  covid_start_date < tri3_start ~ 'tri_2',
  covid_start_date < baby_DOB ~ 'tri_3',
  TRUE ~ 'post_partum'))



### integrate exposure data #####
births_tmp = births

births_tmp = exposure_add(births_tmp)


births = births_tmp

# Start date of pregnancy
births_tmp <- births %>% mutate(year_mth_start = str_sub(start_date, 1, 7))
births_tmp <- births_tmp %>% mutate(year_mth_fact = as.factor(year_mth_start))

# to use the start year and month of the pregnancy, use a spline
datelist = sort(unique(births_tmp$year_mth_start))
births_tmp <- births_tmp %>% mutate(year_mth_ord = match(year_mth_start, datelist))
births = births_tmp

# did the pregnancy start after covid?
births <- births %>% mutate(after_covid = case_when(start_date > covid_start_date ~ TRUE,
                                                    TRUE ~ FALSE))

# if we choose to right-censor: calculate the time at event; if not preterm, use 37wks = 259 days
#births <- births %>% mutate(event_time = case_when(
 # gest_weeks > 37 ~ 37,
  #TRUE ~ gest_weeks
#))


### FILTER BY MONTH ####
#births <- births %>% filter(year_mth_ord <= 8) # 9 was chosen because it's the latest month with more than 2000 births

# add season variable
births <- births %>% mutate(season = ceil(year_mth_ord/3)) %>% mutate(season = as.factor(season))

births_tmp = births

### add in NICU distance ####
nicu = read.csv(paste(fp_data, "nicu_centroid_distances.csv", sep = ''))
nicu <- nicu %>% mutate(countystr = paste("0", STATEFP, str_pad(COUNTYFP, 3, 'left', '0'), sep = '')) %>% dplyr::select(countystr, HubDist) %>% 
  dplyr::group_by(countystr) %>% dplyr::summarize(min(HubDist))
colnames(nicu) = c('countystr', 'HubDist')
births <- merge(births, nicu, by.x = 'county', by.y = 'countystr')

### add in MSSA variables
#mssa = read.csv('~/covid/data/location_to_mssa.csv')
#mssa <- mssa %>% select(-c("longitude", 'latitude', 'start_date', 'baby_DOB'))
#births = merge(births, mssa, by = 'VS_unique')

# our term birth group is 39 - 41 weeks
births <- births %>% filter(gest_weeks < 37 | gest_weeks >= 39)

## END OF VARIABLE CREATION AND SAMPLE SELECTION ###


### write many nas file
many_nas = births %>% filter(na_days > 16) %>% dplyr::select(VS_unique, longitude, latitude, pm25mean, na_days)
write.csv(many_nas, file = '/Volumes/Padlock/covid/debug/nas/many_nas.csv', row.names = F)

  
### write a coordinates file
exposure_calc <- births %>% dplyr::select(VS_unique, longitude, latitude, start_date, baby_DOB,
                                          tri1_start, tri1_end,
                                          tri2_start, tri2_end,
                                          tri3_start, tri3_end)
write.csv(exposure_calc, file = output_name('exposure_coordinates'), row.names = FALSE)

# write a file of the processed columns
write.csv(colnames(births), "/Volumes/Padlock/covid/data/columns_processed.csv", row.names = FALSE)




### TABLE 1 ###############

mytab <- CreateTableOne(vars = c('race', 'pn_care', 'insurance', 'ipi', 
                                 'mom_age', 'inf_sex', 'edu', 
                                 'badair_binary', 'estimate_poverty'), data = births)
tbl1print = print(mytab)
write.csv(tbl1print, file = output_name('tbl1print'))





# frequency tables to check sample sizes

frq_pn <- births1 %>% tabyl(race, pn_care, insurance)

frq_bmi <- births %>% tabyl(race, bmi, insurance, mom_age)

frq_pn <- births512 %>% tabyl(race, pn_care, mom_age)


## ANALYSIS #######

# center and scale continuous variables
births <- births %>% mutate(poverty_scale = scale(estimate_poverty),
                            pm25_scale = pm25mean - mean(pm25mean),
                            pop_density_scale = scale(pop_density))
                            #hubdist_scale = scale(HubDist))
#births<- births %>% filter(year_mth_ord >2)

# split test and train
births_train = sample_frac(births, .75)
births_id <- as.numeric(rownames(births_train))
births_test <- births[-births_id]


# Basic Cox Proportional Hazards model, not adjusting for COVID status

### MODEL SELECTION #####

comparevars = c('ipi', 'mom_age', 'inf_sex', 'edu')
basevars = c('ipi', 'mom_age', 'inf_sex', 'edu', 'badair_binary', 'poverty_scale')
#basevars = c('ipi', 'mom_age', 'inf_sex', 'edu', 'poverty_scale')
basestrata = c('race', 'insurance', 'pn_care')
pmspline <<- c('ns(pm25mean, df = 3)')

seasonspline = c('pspline(year_mth_ord, df = 3)')
add_bad_air = c(TRUE, FALSE)
add_bad_air_count = c(TRUE, FALSE)


coxla = cox_cph("nola", birthsla, c(basevars), c(basestrata), c(pmspline))
coxla1 = cox_cph('la1', birthsla1, c(basevars), c(basestrata), c(pmspline))
coxnull = cox_cph("null", births_train, c(basevars, 'pm25mean'), c(basestrata), c())
coxempty = cox_cph("empty", births_train, c(''), c(''), c(pmspline))
births_train <- births_train %>% filter(pm25mean != 0)
# most basic model
cox0 = cox_cph("0", births_train, c(basevars), c(basestrata), c(pmspline))
cox05 = cox_cph("05", bl5, c(basevars), c(basestrata), c(pmspline))
cox01 = cox_cph("01", births_train, c(basevars), c(basestrata), c(pmspline))
# including a "season" spline
cox1 = cox_cph("1", births_train, c(basevars), c(basestrata), c(pmspline, seasonspline))

births$county = factor(births$county)
coxcounty = cox_cph("county", births_train, c(basevars, 'county'), c(basestrata), c(pmspline, seasonspline))

birthsba = births %>% filter(county %in% c('06001', '06013', '06041', '06075', '06081', '06085'))
coxba = cox_cph("ba", birthsba, c(basevars), c(basestrata), c(pmspline))
# comparing years
cox_comp17 = cox_cph("comp17", births_1718fire, c(comparevars), c(basestrata), c(pmspline))
cox_comp20 = cox_cph("comp20", births_train, c(comparevars), c(basestrata), c(pmspline))

cox_comp = cox_cph('cox_comp', births_train, comparevars, basestrata, c(pmspline, seasonspline))

# three separate terms for PM by trimester
cox2 = cox_cph("2", births_train, c(basevars, 'pm25meantri1', 'pm25meantri2', 'pm25meantri3'), c(basestrata), c(seasonspline))
write.csv(summary(cox2)$coefficients, file = output_name('coeffs_cox2'))

# separate terms for PM by trimester: three separate models
births1 = births %>% filter(season == 1)
births2 = births %>% filter(season == 2)
births3 = births %>% filter(season == 3)

cox2_tri1 = cox_cph("2_tri1", births1, c(basevars, 'pm25meantri1', 'pm25meantri2', 'pm25meantri3'), c(basestrata), c())
write.csv(summary(cox2)$coefficients, file = output_name('coeffs_cox2_tri1'))

cox2_tri2 = cox_cph("2_tri2", births2, c(basevars, 'pm25meantri1', 'pm25meantri2', 'pm25meantri3'), c(basestrata), c())
write.csv(summary(cox2)$coefficients, file = output_name('coeffs_cox2_tri2'))

cox2_tri3 = cox_cph("2_tri3", births3, c(basevars, 'pm25meantri1', 'pm25meantri2', 'pm25meantri3'), c(basestrata), c())
write.csv(summary(cox2)$coefficients, file = output_name('coeffs_cox2_tri3'))

# testing mssa variables
cox_pc = cox_cph('pc', births_train, c(basevars, 'PC_PHYS'), c(basestrata), c(pmspline, seasonspline))

cox_pcr = cox_cph('pcr', births_train, c(basevars, 'PC_PHYS_R'), c(basestrata), c(pmspline, seasonspline))

cox_pcsa = cox_cph('pcsa', births_train, c(basevars, 'PCSA_CIV'), c(basestrata), c(pmspline, seasonspline))

# separate by COVID status, model, plot
births_train_neg = births %>% filter(any_mat_covid == 0)
births_train_pos <- births %>% filter(any_mat_covid == 1)

cox1_neg = cox_cph("1_neg", births_train_neg, c(basevars), c(basestrata), c(pmspline, seasonspline))
cox1_pos = cox_cph("1_pos", births_train_pos, c(basevars), c(basestrata), c(pmspline, seasonspline))

prop_hr_plot = hazard_ratios_compare(cox1_neg[[1]], cox1_pos[[1]], 'tst', 'pm25mean', list('COVID-', 'COVID+')) 
prop_hr_plot

### RESULTS NUMBERS

# mean air
mean(births$pm25mean)

# range air
range(births$pm25mean)

# % of people with a bad air day
sum(births$badair_binary)/nrow(births)*100

# mean number of badair days among those with at least one
mean(births[births$badair_binary ==1, ]$badair_count)

# how many covid+ cases?
nrow(births[births$any_mat_covid == 1, ])

# overall ptb rate?
sum(births$preterm)/nrow(births)*100

# overall PTB rate by month
ptrate = births %>% dplyr::group_by(year_mth_fact) %>% dplyr::summarise(ptrate = 100*sum(preterm)/n())

#p1 <- ggMarginal(cox1_neg[[2]], type = 'boxplot')
#p1



## COVID ~ AIR
covidair = glm(any_mat_covid ~ pm25mean, data = births, family = 'binomial')
summary(covidair)

## LOGISTIC MODELS #######
logmdl = glm(preterm ~ pm25mean, data = birthsba, family = 'poisson')
logsum = tidy(logmdl)
logsum <- logsum %>% mutate(rr = exp(estimate), rrlow = exp(estimate - 1.96*std.error), rrhigh = exp(estimate + 1.96*std.error))
logsum


logreg = glm(preterm ~ pm25mean, data = birthsba, family = 'binomial')
logregs = tidy(logreg)
logregs <- logregs %>% mutate(rr = exp(estimate), rrlow = exp(estimate - 1.96*std.error), rrhigh = exp(estimate + 1.96*std.error))
logregs

### categorical pm25
births <- births %>% mutate(tile5 = ntile(pm25mean, 5))
table(births$tile5)
mean(births[births$tile5 == 1, ]$pm25mean)

tbl5 = epitable(births$tile5, births$preterm)
tbl5

tbl5rr = epitab(tbl5, method = 'riskratio')
tbl5rr

births_check = births
births_term = births[births$preterm == 0, ]
births_term <- births_term %>% mutate(tile10 = ntile(pm25mean, 10))

births_dec = births_term %>% dplyr::group_by(tile10) %>% 
  summarise(pct = n()/nrow(births_term),
            maxval = max(pm25mean))

deccut = births_dec$maxval

births_pt = births[births$preterm == 1, ]
births_pt <- births_pt %>% mutate(tile10 = 
                                    case_when(pm25mean < deccut[1] ~ 1,
                                              pm25mean < deccut[2] ~ 2,
                                              pm25mean < deccut[3] ~ 3,
                                              pm25mean < deccut[4] ~ 4,
                                              pm25mean < deccut[5] ~ 5,
                                              pm25mean < deccut[6] ~ 6,
                                              pm25mean < deccut[7] ~ 7,
                                              pm25mean < deccut[8] ~ 8,
                                              pm25mean < deccut[9] ~ 9,
                                              TRUE ~ 10))

births_pt_dec = births_pt %>% dplyr::group_by(tile10) %>%
  summarise(pct = n()/nrow(births_pt))

bt_plot = ggplot(births_pt_dec, aes(x = tile10, y = pct)) + 
  geom_col() + 
  geom_hline(yintercept = .1) + 
  ylab('Fraction of preterm births') + 
  xlab('PM deciles based on term births') + 
  theme_bw()
bt_plot

ggsave(bt_plot, filename = image_name('decile_plot'))



## PLOTS PT1 #############

# rural vs air quality
nicu_dist_plt = ggplot(data = births, aes(x = hubdist_scale, y = pm25mean)) + 
  geom_point() + 
  theme_bw()

nicu_dist_plt

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











cox_cph(births, c(basevars, 'season'), c(basestrata), c(pmspline, 'ns(estimate_poverty, df = 3)'))

births_18 = births %>% filter(mom_age %in% c("mom_age_18_34", "mom_age_gt34"))
cox_cph(births_18, basevars, c(basestrata, "season"), c(pmspline))

cox_cph(births_train, c(basevars, 'isrural'), basestrata, c(pmspline, seasonspline))

png(filename= image_name('termplot_full_noseason'))
cox_cph(births_train, basevars, basestrata, c(pmspline, seasonspline))
dev.off()


cox_cph(births_train, c(basevars, "hubdist_scale", 'season'), basestrata, c(pmspline))

cox_cph(births_train, c(basevars, "estimate_unemployment_rate"), basestrata, c(pmspline, seasonspline))

cox_cph(births_train, c(basevars, "HubDist"), basestrata, c(pmspline))

cox_cph(births_train, c(basevars, "HubDist", 'season'), basestrata, c(pmspline))

cox_cph(births_train, c(basevars, 'HubDist'), c(basestrata, 'season'), c(pmspline))

# build models
cox0 <- coxph(Surv(gest_weeks, preterm) ~  ns(pm25mean, df = 3), data = births_train)
summary(cox0)

tfn(cox0)

cox0a <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + pspline(year_mth_ord, df = 4) + 
                 ipi + mom_age  + inf_sex  + estimate_poverty + edu  +  
                 strata(race, insurance, pn_care), data = births_tmp)
tfn <- function(df) {
  termplot(df, term = 1, se = T, plot = T)
  
}

tfn(cox0a)

cox0b <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3)  + 
                 ipi + mom_age + estimate_poverty + edu + inf_sex + badair_count + season + 
                 strata(race, insurance, pn_care, after_covid), data = births_train)
tfn(cox0b)

coxseason <- coxph(Surv(gest_weeks, preterm) ~ pspline(year_mth_ord, df = 3), data = births_train )

tfn(coxseason)
summary(coxseason)

cox1 <- coxph(Surv(gest_weeks, preterm) ~  pm25mean + year_mth_ord +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox2 <- coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + year_mth_ord +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)

cox3 <- coxph(Surv(gest_weeks, preterm) ~  pm25mean + pspline(year_mth_ord, df = 4) +
                ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                strata(race, insurance, pn_care), data = births_train)
summary(cox3)

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

spline_termplot <- termplot(cox9, term = 1, se=T, plot = T)


# compare models
aic_cox = lapply(model_list, AIC)

min_aic_cox <- which.min(lapply(model_list, AIC))

# comparing PM25 to PM25 spline
exp((aic_cox[[4]] - aic_cox[[3]])/2)

# comparing linear + quadratic to spline
exp((aic_cox[[4]] - aic_cox[[5]])/2)
births_tmp = births %>% filter(edu_lt12 == TRUE)
log_like = lapply(model_list, logLik)
png(filename = image_name('cox_model_general_termplot'))
cox12_a <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + pspline(year_mth_ord, df = 4) + 
                 ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + pop_density_scale +
                 strata(race, insurance, pn_care, after_covid), data = births_tmp)
spline_termplot <- termplot(cox12_a, term = 1, se=T, plot = T)
dev.off()

cox12_b <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + pspline(year_mth_ord, df = 4) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   pop_density_scale + pop_density_scale^2 + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_b, term = 1, se=T, plot = T)

cox12_c <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + 
                   pspline(year_mth_ord, df = 4) + pspline(pop_density_scale, df = 4) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_c, term = 1, se=T, plot = T)

cox12_d <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + 
                   pspline(year_mth_ord, df = 4) + pspline(pop_density_scale, df = 3) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_d, term = 1, se=T, plot = T)

cox12_e <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + 
                   pspline(year_mth_ord, df = 4) + ntile(pop_density_scale, 10) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_e, term = 1, se=T, plot = T)

cox12_f <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + 
                   pspline(year_mth_ord, df = 4) + ns(pop_density_scale, 3) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_f, term = 1, se=T, plot = T)

cox12_g <- coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3) + 
                   pspline(year_mth_ord, df = 4) + ns(pop_density_scale, 4) + 
                   ipi + mom_age  + inf_sex  + estimate_poverty +  badair_binary + 
                   strata(race, insurance, pn_care, after_covid), data = births_train)
spline_termplot <- termplot(cox12_g, term = 1, se=T, plot = T)


# split data by covid status

cox_model <- function(df) {
  return(coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3)  + 
                     ipi + mom_age  + inf_sex  + estimate_poverty + pop_density_scale + edu + season + 
                     strata(race, insurance, pn_care), data = df))
  }

cox_model_nofire <- function(df) {
  return(coxph(Surv(gest_weeks, preterm) ~ ns(pm25mean, df = 3)  + 
                 ipi + mom_age  + inf_sex  + estimate_poverty  + pop_density_scale + edu + 
                 strata(race, insurance, pn_care), data = df))
}

fulldata_removed = cox_model(births)
termplot(fulldata_removed, term = 1, se = T, plot = T)

a = cox_model(births_fire)
termplot(a, term = 1, se = T, plot = T)

births_nofire = births %>% filter(baby_DOB<as.Date('2020-09-01'))
nofire_cox = cox_model_nofire(births_nofire)

termplot(nofire_cox, term = 1, se=T, plot = T)

births_train_neg = births_train %>% filter(any_mat_covid == 0)
cox_neg = cox_model(births_train_neg)

png(filename = image_name('neg_termplot'))
neg_termplot <- termplot(cox_neg, term = 1, se=T, plot = T)
dev.off()

births_train_pos <- births_train %>% filter(any_mat_covid == 1)

cox_pos = cox_model(births_train_pos)
png(filename = image_name('pos_termplot'))
pos_termplot <- termplot(cox_pos, term = 1, se=T, plot = T)
dev.off()

## subgroup analysis - before covid, separate from after
before_neg = births %>% filter(after_covid == FALSE, any_mat_covid == FALSE)
before_pos = births %>% filter(after_covid == FALSE, any_mat_covid == TRUE)

after_neg = births %>% filter(after_covid == TRUE, any_mat_covid == FALSE)
after_pos = births %>% filter(after_covid == TRUE, any_mat_covid == TRUE)

before_neg_cox = cox_model(before_neg)
before_pos_cox = cox_model(before_pos)

png(filename = image_name('pos_before_termplot'))
before_pos_termplot <- termplot(before_pos_cox, term = 1, se=T, plot = T)
dev.off()

png(filename = image_name('neg_before_termplot'))
before_neg_termplot <- termplot(before_neg_cox, term = 1, se=T, plot = T)
dev.off()

# after covid needs fewer strata probably, currently errors
after_neg_cox = cox_model(after_neg)
after_pos_cox = cox_model(after_pos)


pos_summ <- summary(cox_pos)
pos_coeff <- round(data.frame(pos_summ$coefficients), 2)
pos_ci <- round(data.frame(pos_summ$conf.int), 2)
pos_ci <- pos_ci %>% rownames_to_column()
colnames(pos_ci) <- c('varname', 'var_exp', 'var_nexp', 'lower', 'upper')

write.csv(pos_coeff, file = output_name('pos_model_coeff'), row.names = TRUE)
write.csv(pos_ci, file = output_name('pos_model_ci'), row.names  = TRUE)

hazard_ratios(cox_pos, 'pos', 'pm25mean')
hazard_forest(pos_ci, 'positive', 'positive')

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

### separate by covid status ############
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



### ARCHIVE #####


cox_coxph <- function(df) {
  coxmodel = coxph(Surv(gest_weeks, preterm) ~  pspline(pm25mean, df = 4) + pspline(year_mth_ord, df = 4) +
                     ipi + mom_age  + inf_sex  + estimate_poverty +  badair_count +
                     strata(race, insurance, pn_care), data = df)
  return(coxmodel)
}




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

covpos = births[births$any_mat_covid == 1,]
hist(covpos$pm25mean)

covpt = births %>% filter(any_mat_covid == 1 & preterm == 1)
hist(covpt$pm25mean)

covt = births %>% filter(any_mat_covid == 1 & preterm == 0)
hist(covt$pm25mean)
summary(covpt$pm25mean)
summary(covt$pm25mean)
