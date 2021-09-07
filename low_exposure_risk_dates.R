
highlow = births %>% mutate(low = case_when(pm25mean < 11 ~ TRUE, TRUE ~FALSE))
table(highlow$low)

contingency_low <- epitable(highlow$low, highlow$preterm)
low_rr = epitab(contingency_low, method = 'riskratio')
low_rr


covid_start_date <<- ymd('2020-03-16')
# why do low patients have an increased risk of PTB?

# date of conception?
lowdate_box = ggplot(data = highlow, aes(x = low, y = start_date)) + 
  geom_boxplot() + 
  theme_bw()
lowdate_box

which.min(lowdf$start_date)
e_cnc_id = lowdf$VS_unique[which.min(lowdf$start_date)]
e_dob_id = lowdf$VS_unique[which.min(lowdf$baby_DOB)]
#t_cnc_id
#t_dob_id
#m_cnc_id = mean(lowdf$start_date)
#m_dob_id = mean(lowdf$baby_DOB)
#s_cnc_id
#s_dob_id
l_cnc_id = lowdf$VS_unique[which.max(lowdf$start_date)]
l_dob_id = lowdf$VS_unique[which.min(lowdf$baby_DOB)]

#vartype = c('earliest_cnc', 'earliest_dob', '25_start', '25_dob', 'mean_cnc','mean_dob', '75_cnc', '75_dob', 'latest_start', 'latest_dob')
vartype = c('second_quart', 'mean_cnc','third_quart')
start = c(quantile(lowdf$start_date, .25, type = 1)[[1]], mean(lowdf$start_date), quantile(lowdf$start_date, .75, type = 1)[[1]])
end = c(quantile(lowdf$baby_DOB, .25, type = 1)[[1]],
            mean(lowdf$baby_DOB),
            quantile(lowdf$baby_DOB, .75, type = 1)[[1]])
low = c("LOW", "LOW", "LOW")

low_timevals = data.frame(vartype, start, end, low)
low_timevals



vartype = c('second_quart', 'mean_cnc','third_quart')
start = c(quantile(highdf$start_date, .25, type = 1)[[1]], mean(highdf$start_date), quantile(highdf$start_date, .75, type = 1)[[1]])
end = c(quantile(highdf$baby_DOB, .25, type = 1)[[1]],
            mean(highdf$baby_DOB),
            quantile(highdf$baby_DOB, .75, type = 1)[[1]])
low = c('HIGH', 'HIGH', "HIGH")

high_timevals = data.frame(vartype, start, end, low)
high_timevals


timevals = rbind(low_timevals, high_timevals)
timevals$vartype = factor(timevals$vartype, levels = c('third_quart', 'mean_cnc', 'second_quart'))



lowhigh_timeplt = ggplot(data = timevals, aes(x = start, y = vartype)) + 
  geom_crossbar(aes(xmin = start, xmax = end), width = 0.2, fill = '#DDDDDD') + 
  geom_vline(xintercept = ymd('2020-03-16'), color = '#004488', size = 1) +  # lockdown in LA county
  #geom_vline(xintercept = ymd('2020-08-12'), color = '#BB5566', size = 1) +  # Lake fire in LA county
  #geom_vline(xintercept = ymd('2020-08-16'), color = '#DDAA33', size = 1) +  # CZU complex in Bay Area
  geom_rect(aes(xmin=ymd('2020-08-12'), xmax=ymd("2020-10-01")), ymin=0, ymax=4,
            color="#DDAA33", alpha=0.1, fill='#DDAA33') +
  theme_bw() + 
  facet_grid(rows = vars(low))

lowhigh_timeplt

ggsave(lowhigh_timeplt, file = '/Volumes/Padlock/covid/images/2021-09-06_lowhighPM_dates.png')


