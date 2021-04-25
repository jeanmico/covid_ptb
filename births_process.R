library(tidyverse)

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
library(tableone)
mytab <- CreateTableOne(vars = c(fact_cols, 'wic'), data = births, factorVars = fact_cols)

mytab

# exploring sample and sample sizes


# determine start date of pregnancy



# plot pregnancy timeline with addl landmarks
#  when did the Bay Area shut down?
#  when did Los Angeles shut down?