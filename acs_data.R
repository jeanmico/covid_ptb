library(tidycensus)
library(tidyverse)


setwd('/Users/student/covid/')

census_api_key('38df6ff2cf95698815a732683e38ca6ddb9a42cf')

# get list of acs variables
acs_vars <- read.csv('acs_vars.csv')
acs_vec <- as.vector(acs_vars$name)

# decode dataframes for renaming columns
decodes_est <- acs_vars %>% dplyr::select(name,rename)
decodes_est$codename <- paste('estimate',decodes_est$name,sep = '_')

# download acs variables
ca_acs <- get_acs(geography ="tract", variables = acs_vec, year=2019, state="CA", survey = 'acs5')

# add readable variable names
ca_acs <- merge(ca_acs, decodes_est, by.x = 'variable', by.y = 'name', all.x = TRUE)

# write output
write.csv(ca_acs, file = 'acs_data_5yr_2019.csv', row.names = FALSE)
