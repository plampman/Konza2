#########################################################################
###
### Program name: Konza2 Purple Air Data
###
### Purpose: Analyze Purple Air data
###
### Author: Phinehas Lampman, plampman@uidaho.edu
###
### Last modified: 10/21/2022
###
#########################################################################

########### NOTES

## 


###########################################################################################
### PM2.5 EPA correction equation (EPA Tools and Resources Webinar 
### AirNow Fire and Smoke Map: Extension of the US-Wide Correction 
### for Purple PM2.5 Sensors)
# Low concentration <50 ug/m^3: PM2.5 = 0.52 * PM(ATM) - (0.086 * RH) + 5.75
# Mid concentration >50 ug/m^3 <229 ug/m^3: PM2.5 = 0.786 * PM(ATM) - (0.086 * RH) + 5.75
# High concentration >229 ug/m^3: PM2.5 = 0.69 * PM(ATM) + 8.84x10^-4 * PM(ATM)^2 + 2.97
###########################################################################################

library(tidyverse)
library(gt)
require(lubridate)
require(dplyr)
library(reshape2)


################### Sampling intervals

## C:\School_research\Konza2\PA

intervals <- read.table('C:/School_research/Konza2/PA/sampling_periods_konza2.txt' , header = TRUE, sep = ",", fill = TRUE) %>% ## add connection
  mutate(date_time_cdt = ymd_hms(date_time_cdt, tz = "CST6CDT"))

## BLUE

pm_blue <- read.csv('C:/School_research/Konza2/PA/blue_20220411.csv', header = TRUE) %>% ## 
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'), 
         local_time = with_tz(UTCDateTime, tzone = "CST6CDT"),
         ## Quality control of particulate matter using the two identical PMS5003 sensors
         difference = abs(pm2_5_atm - pm2_5_atm_b), 
         RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
         RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
         QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
         ## Calculating mixing ratio 
         temp_c = (5/9 * (current_temp_f - 32)),
         dp_c = (5/9 * (current_dewpoint_f - 32)),
         Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
         sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
         mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
         saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
         platform = "blue") %>%
  ## Averages of particulate matter between two sensors
  ## correcting PM 2.5 with EPA high concentration formula 
  rowwise() %>% 
  mutate(pm1_0_avg = mean(c(pm1_0_atm,pm1_0_atm_b), na.rm	= TRUE), 
         pm2_5_avg = mean(c(pm2_5_atm, pm2_5_atm_b), na.rm = TRUE), 
         pm10_0_avg = mean(c(pm10_0_atm, pm10_0_atm_b), na.rm = TRUE), 
         pm2_5_corr = as.integer(case_when(
           pm2_5_avg <= 50 ~ ((0.52 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 50 & pm2_5_avg <= 229 ~ ((0.786 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 229 ~ ((0.69 * pm2_5_avg) + ((8.84*10^-4) * (pm2_5_avg)^2) + 2.97)))) %>% 
  ungroup 


## Extracting sampling interval specific data
key1 <- c(
  '0'= "nothing", '1' = "Ambient 1",'2' = "nothing",'3' = "A2",'4' = "nothing",'5' = "Smoke 1",
  '6' = "nothing",'7' = 'S2','8'= "nothing",'9'= "S3", '10' = 'nothing','11' = "S4",
  '12' = "nothing",'13' = "S5",'14' = "nothing",'15' = "S6", '16' = "nothing", '17' = "S7",
  '18' = "nothing", '19' = "S8", '20' = "nothing")

key <- c(
  '0'= "nothing", '1' = "A",'2' = "nothing",'3' = "A",'4' = "nothing",'5' = "S",
  '6' = "nothing",'7' = 'S','8'= "nothing",'9'= "S", '10' = 'nothing','11' = "S",
  '12' = "nothing",'13' = "S",'14' = "nothing",'15' = "S", '16' = "nothing", '17' = "S",
  '18' = "nothing", '19' = "S", '20' = "nothing", '21' = "S", '22' = "nothing")

categories_blue <- pm_blue %>%
  mutate(int = findInterval(pm_blue$local_time, na.omit(intervals$date_time_cdt)),
         sample = key[as.character(int)]) %>%
  mutate(pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
         pm2_5_avg = ifelse(QC == "bad", NA_character_, pm2_5_avg),
         pm10_0_avg = ifelse(QC == "bad", NA_character_, pm10_0_avg),
         pm1_0_avg = ifelse(QC == "bad", NA_character_, pm1_0_avg))

## Filtering out rows that are outside of sampling intervals
samples_blue = categories_blue[categories_blue$sample != 'nothing',] 

samples_blueS = samples_blue[samples_blue$sample != "A",]


stats_blue <- samples_blue %>%
  group_by(sample) %>%
  summarise(
    "PM2.5 Corrected MEAN"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected SD" = sd(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MAX" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MIN" = min(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Avg. MEAN" =    mean(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. SD" = sd(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MAX" = max(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MIN" = min(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM1.0 Avg. MEAN" = mean(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. SD" = sd(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MAX" = max(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MIN" = min(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM10 Avg. MEAN" = mean(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. SD" = sd(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MAX" = max(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MIN" = min(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "TEMP (C) MEAN" = mean(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) SD" = sd(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MAX" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MIN" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "% RH MEAN" = mean(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH SD" = sd(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MIN" = min(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MAX" = max(as.numeric(unlist(current_humidity)), na.rm = T),
    "Mixing Ratio MEAN" = mean(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio SD" = sd(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MAX" = max(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MIN" = min(as.numeric(unlist(mixing_ratio)), na.rm = T),
  ) 

##write.csv(stats_blue,"C:/School_research/Inundation2/Blue_PA_Inundation2_Stats.csv")


ggplot(samples_blue, 
       aes(sample, as.numeric(unlist(pm2_5_corr)))) + geom_violin(scale = "width", fill = "lightskyblue1") + geom_jitter(height = 0, width = 0.1, size = 1) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7) +
  labs(x = "", y = 'PM 10 (ug/m3)') + scale_y_continuous(breaks = c(0, 1000, 2000, 3000), limits = c(0, 4000), labels = c(0, 1000, 2000, 3000))+
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))

samples_blueS = as.numeric(samples_blueS$pm2_5_corr)
samples_blueS = as.numeric(samples_blueS$pm10_0_avg)
samples_blueS = as.numeric(samples_blueS$current_temp_f)

blue_statsfig <- samples_blueS %>%
  pivot_longer(cols = c(pm2_5_corr, pm10_0_avg, current_temp_f), names_to = 'metric', values_to = 'long')


xorder_blue <- factor(blue_statsfig$metric, level = c("pm2_5_corr", "pm10_0_avg", "current_temp_f"))

ggplot(blue_statsfig, aes(x = xorder_blue, as.numeric(unlist(long)), fill = long)) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7, fill = "grey40") +
  labs(x = "", y = 'PM 2.5 (ug/m3)') + scale_x_discrete(labels = c("Blue", "Green", "Yellow")) +
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))


## YELLOW

pm_yellow <- read.csv('C:/School_research/Konza2/PA/yellow_20220411.csv', header = TRUE) %>% ## 
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'), 
         local_time = with_tz(UTCDateTime, tzone = "CST6CDT"),
         ## Quality control of particulate matter using the two identical PMS5003 sensors
         difference = abs(pm2_5_atm - pm2_5_atm_b), 
         RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
         RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
         QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
         ## Calculating mixing ratio 
         temp_c = (5/9 * (current_temp_f - 32)),
         dp_c = (5/9 * (current_dewpoint_f - 32)),
         Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
         sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
         mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
         saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
         platform = "Yellow") %>%
  ## Averages of particulate matter between two sensors
  ## correcting PM 2.5 with EPA high concentration formula 
  rowwise() %>% 
  mutate(pm1_0_avg = mean(c(pm1_0_atm,pm1_0_atm_b), na.rm	= TRUE), 
         pm2_5_avg = mean(c(pm2_5_atm, pm2_5_atm_b), na.rm = TRUE), 
         pm10_0_avg = mean(c(pm10_0_atm, pm10_0_atm_b), na.rm = TRUE), 
         pm2_5_corr = case_when(
           pm2_5_avg <= 50 ~ ((0.52 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 50 & pm2_5_avg <= 229 ~ ((0.786 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 229 ~ ((0.69 * pm2_5_avg) + ((8.84*10^-4) * (pm2_5_avg)^2) + 2.97))) %>% 
  ungroup 


categories_yellow <- pm_yellow %>%
  mutate(int = findInterval(pm_yellow$local_time, na.omit(intervals$date_time_cdt)),
         sample = key[as.character(int)]) %>%
  mutate(pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
         pm2_5_avg = ifelse(QC == "bad", NA_character_, pm2_5_avg),
         pm10_0_avg = ifelse(QC == "bad", NA_character_, pm10_0_avg),
         pm1_0_avg = ifelse(QC == "bad", NA_character_, pm1_0_avg))

## Filtering out rows that are outside of sampling intervals
samples_yellow = categories_yellow[categories_yellow$sample != 'nothing',] 

samples_yellowS = samples_yellow[samples_yellow$sample != "A",]


stats_yellow <- samples_yellow %>%
  group_by(sample) %>%
  summarise(
    "PM2.5 Corrected MEAN"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected SD" = sd(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MAX" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MIN" = min(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Avg. MEAN" =    mean(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. SD" = sd(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MAX" = max(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MIN" = min(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM1.0 Avg. MEAN" = mean(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. SD" = sd(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MAX" = max(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MIN" = min(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM10 Avg. MEAN" = mean(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. SD" = sd(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MAX" = max(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MIN" = min(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "TEMP (C) MEAN" = mean(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) SD" = sd(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MAX" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MIN" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "% RH MEAN" = mean(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH SD" = sd(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MIN" = min(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MAX" = max(as.numeric(unlist(current_humidity)), na.rm = T),
    "Mixing Ratio MEAN" = mean(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio SD" = sd(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MAX" = max(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MIN" = min(as.numeric(unlist(mixing_ratio)), na.rm = T),
  ) 

##write.csv(stats_yellow,"C:/School_research/Inundation2/Blue_PA_Inundation2_Stats.csv")


ggplot(samples_yellow, 
       aes(sample, as.numeric(unlist(pm2_5_corr)))) + geom_violin(scale = "width", fill = "yellow") + geom_jitter(height = 0, width = 0.1, size = 1) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7) +
  labs(x = 'Treatment', y = 'PM 2.5 (ug/m3)', title = 'Yellow Corrected PM2.5 Inundation 2') + 
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))

ggplot(samples_yellow, 
       aes(sample, as.numeric(unlist(pm2_5_corr)))) + geom_violin(scale = "width", fill = "yellow") + geom_jitter(height = 0, width = 0.1, size = 1) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7) +
  labs(x = "", y = 'PM 10 (ug/m3)') + scale_y_continuous(breaks = c(0, 1000, 2000, 3000), limits = c(0, 4000), labels = c(0, 1000, 2000, 3000))+
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))

## GREEN

pm_green <- read.csv('C:/School_research/Konza2/PA/green_20220411.csv', header = TRUE) %>% ## 
  mutate(UTCDateTime = ymd_hms(UTCDateTime, tz = 'UTC'), 
         local_time = with_tz(UTCDateTime, tzone = "CST6CDT"),
         ## Quality control of particulate matter using the two identical PMS5003 sensors
         difference = abs(pm2_5_atm - pm2_5_atm_b), 
         RPD = abs((pm2_5_atm - pm2_5_atm_b)/((pm2_5_atm + pm2_5_atm_b)/2)*100),
         RPD_diff = if_else(RPD > 70, 1, 0) + if_else(difference > 5, 1, 0),
         QC = case_when(RPD_diff == 2 ~ "bad", TRUE ~ "good"),
         ## Calculating mixing ratio 
         temp_c = (5/9 * (current_temp_f - 32)),
         dp_c = (5/9 * (current_dewpoint_f - 32)),
         Vapor_pressure = (6.11*10^((7.5*dp_c)/(237.7+dp_c))),
         sat_vapor = (6.11*10^((7.5*temp_c)/(237.7+temp_c))),
         mixing_ratio = (621.97*(Vapor_pressure/(pressure - Vapor_pressure))),
         saturated_mr = (621.97*(sat_vapor/(pressure - sat_vapor))),
         platform = "Green") %>%
  ## Averages of particulate matter between two sensors
  ## correcting PM 2.5 with EPA high concentration formula 
  rowwise() %>% 
  mutate(pm1_0_avg = mean(c(pm1_0_atm,pm1_0_atm_b), na.rm	= TRUE), 
         pm2_5_avg = mean(c(pm2_5_atm, pm2_5_atm_b), na.rm = TRUE), 
         pm10_0_avg = mean(c(pm10_0_atm, pm10_0_atm_b), na.rm = TRUE), 
         pm2_5_corr = case_when(
           pm2_5_avg <= 50 ~ ((0.52 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 50 & pm2_5_avg <= 229 ~ ((0.786 * pm2_5_avg) - (0.086 * current_humidity) + 5.75),
           pm2_5_avg > 229 ~ ((0.69 * pm2_5_avg) + ((8.84*10^-4) * (pm2_5_avg)^2) + 2.97))) %>% 
  ungroup 


categories_green <- pm_green %>%
  mutate(int = findInterval(pm_green$local_time, na.omit(intervals$date_time_cdt)),
         sample = key[as.character(int)]) %>%
  mutate(pm2_5_corr = ifelse(QC == "bad", NA_character_, pm2_5_corr),
         pm2_5_avg = ifelse(QC == "bad", NA_character_, pm2_5_avg),
         pm10_0_avg = ifelse(QC == "bad", NA_character_, pm10_0_avg),
         pm1_0_avg = ifelse(QC == "bad", NA_character_, pm1_0_avg))

## Filtering out rows that are outside of sampling intervals
samples_green = categories_green[categories_green$sample != 'nothing',] 

samples_greenS = samples_green[samples_green$sample != "A",]

stats_green <- samples_green %>%
  group_by(sample) %>%
  summarise(
    "PM2.5 Corrected MEAN"  = mean(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected SD" = sd(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MAX" = max(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Corrected MIN" = min(as.numeric(unlist(pm2_5_corr)), na.rm = T),
    "PM2.5 Avg. MEAN" =    mean(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. SD" = sd(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MAX" = max(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM2.5 Avg. MIN" = min(as.numeric(unlist(pm2_5_avg)), na.rm = T),
    "PM1.0 Avg. MEAN" = mean(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. SD" = sd(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MAX" = max(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM1.0 Avg. MIN" = min(as.numeric(unlist(pm1_0_avg)), na.rm = T),
    "PM10 Avg. MEAN" = mean(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. SD" = sd(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MAX" = max(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "PM10 Avg. MIN" = min(as.numeric(unlist(pm10_0_avg)), na.rm = T),
    "TEMP (C) MEAN" = mean(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) SD" = sd(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MAX" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "TEMP (C) MIN" = max(as.numeric(unlist(temp_c)), na.rm = T),
    "% RH MEAN" = mean(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH SD" = sd(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MIN" = min(as.numeric(unlist(current_humidity)), na.rm = T),
    "% RH MAX" = max(as.numeric(unlist(current_humidity)), na.rm = T),
    "Mixing Ratio MEAN" = mean(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio SD" = sd(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MAX" = max(as.numeric(unlist(mixing_ratio)), na.rm = T),
    "Mixing Ratio MIN" = min(as.numeric(unlist(mixing_ratio)), na.rm = T),
  ) 

## write.csv(stats_green,"C:/School_research/Inundation2/Blue_PA_Inundation2_Stats.csv")


ggplot(samples_green, 
       aes(sample, as.numeric(unlist(pm2_5_corr)))) + geom_violin(scale = "width", fill = "green") + geom_jitter(height = 0, width = 0.1, size = 1) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7) +
  labs(x = 'Treatment', y = 'PM 2.5 (ug/m3)', title = 'Green Corrected PM2.5 Inundation 2') + 
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))

ggplot(samples_green, 
       aes(sample, as.numeric(unlist(pm2_5_corr)))) + geom_violin(scale = "width", fill = "green") + geom_jitter(height = 0, width = 0.1, size = 1) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7) +
  labs(x = "", y = 'PM 2.5 (ug/m3)') + scale_y_continuous(breaks = c(0, 1000, 2000, 3000), limits = c(0, 4000), labels = c(0, 1000, 2000, 3000))+
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))




##JOINED

BYG = full_join(samples_blueS, samples_yellowS, by = "platform", keep = T) %>%
  full_join(samples_greenS, by = c('platform.y' = 'platform'), keep = T) %>%
  pivot_longer(cols = c(pm2_5_corr.x, pm2_5_corr.y, pm2_5_corr), names_to = 'Drone', values_to = 'pm2_5')


xorder <- factor(BYG$Drone, level = c("pm2_5_corr.x", "pm2_5_corr", "pm2_5_corr.y"))

ggplot(BYG, aes(x = xorder, as.numeric(unlist(pm2_5)), fill = Drone)) + 
  geom_boxplot(alpha = 0.7, fatten = 0.5, lwd=0.7, fill = "grey40") +
  labs(x = "", y = 'PM 2.5 (ug/m3)') + scale_x_discrete(labels = c("Blue", "Green", "Yellow")) +
  theme(legend.position = 'none') + theme(axis.text.x = element_text(angle = -55, hjust = 0)) + theme(axis.text = element_text(size = 12)) +
  theme(axis.title = element_text(size = 14, face = "bold")) + theme(title = element_text(size = 18, face = "bold")) + 
  theme(plot.background = element_rect(colour = "black", fill = NA, size = 2))



