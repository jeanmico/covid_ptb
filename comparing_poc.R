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
