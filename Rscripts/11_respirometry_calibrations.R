#### Respirometer calibrations

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)

#### step 1. read in the data

cal_20 <- read_csv("respirometry-data/respcalibrations/20C_calibration_may262016_Phase.csv")
cal_12 <- read_csv("respirometry-data/respcalibrations/12C_calibration_may132016_477_Phase.csv")
cal_16 <- read_csv("respirometry-data/respcalibrations/16C_calibration_may262016_Phase.csv")
cal_24 <- read_csv("respirometry-data/respcalibrations/24C_calibration_may262016_Phase.csv")
cal_12b <- read_csv("respirometry-data/12C_calibration_may302016_Phase.csv")




str(cal_16)

?read_csv
#### step 2. convert to long form

cal_20_long <- cal_20 %>% 
	gather(., "channel", "phase", 3:26) 

cal_20_long$oxygen_status <- ifelse(cal_20_long$channel %in% c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", "C6"), "anoxic", "full_oxy")

cal_12_long %>% 
	filter(Time_in_min > 140 & Time_in_min < 143) %>% 
	group_by(oxygen_status) %>% 
	mutate(mean_phase = mean(phase)) %>% 
	ggplot(data = ., aes(x = oxygen_status, y = mean_phase)) + geom_boxplot()

cal_12_longb %>% 
	filter(Time_in_min > 113) %>% 
	select(phase, oxygen_status) %>% 
group_by(oxygen_status) %>%
	filter(!is.na(phase)) %>% 
	summarise_each(funs(mean, std.error)) %>% View
	ggplot(data = ., aes(x = oxygen_status, y = mean)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + theme_bw() + ylab("phase") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))



cal_20_long %>% 
	filter(Time_in_min > 140 & Time_in_min < 143) %>% 
	filter(channel %in% c("B1", "B2", "B3", "B4", "B5", "B6", "D1", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(channel) %>% 
	mutate(mean_phase = mean(phase)) %>% 
	ggplot(data = ., aes(x = channel, y = mean_phase)) + geom_boxplot()


cal_12_long <- cal_12 %>% 
	gather(., "channel", "phase", 3:26) 
cal_12_long$oxygen_status <- ifelse(cal_12_long$channel %in% c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", "C6"), "anoxic", "full_oxy")

cal_12_longb <- cal_12b %>% 
	gather(., "channel", "phase", 3:26) 
cal_12_longb$oxygen_status <- ifelse(cal_12_longb$channel %in% c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", "C6"), "anoxic", "full_oxy")

cal_12_longb %>% View
cal_12_longb %>%
	filter(!is.na(phase)) %>% 
	# filter(Time_in_min >110) %>% 
	ggplot(data = ., aes(x = Time_in_min, y = phase)) + geom_point() + facet_wrap( ~ channel, scales = "free")

cal_16_long <- cal_16 %>% 
	gather(., "channel", "phase", 3:26) 
cal_16_long$oxygen_status <- ifelse(cal_16_long$channel %in% c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", "C6"), "anoxic", "full_oxy")

cal_16_long %>% 
	filter(Time_in_min > 140 & Time_in_min < 143) %>% 
	ggplot(data = ., aes(x = Time_in_min, y = phase)) + geom_point() + facet_wrap( ~ channel, scales = "free")

cal_24_long <- cal_24 %>% 
	gather(., "channel", "phase", 3:26) 
cal_24_long$oxygen_status <- ifelse(cal_24_long$channel %in% c("A1", "A2", "A3", "A4", "A5", "A6", "C1", "C3", "C4", "C5", "C6"), "anoxic", "full_oxy")

cal_24_long %>% 
	filter(Time_in_min > 140 & Time_in_min < 143) %>% 
	ggplot(data = ., aes(x = Time_in_min, y = phase)) + geom_point() + facet_wrap( ~ channel, scales = "free")


all_cals <- bind_rows(cal_24_long, cal_20_long, cal_16_long, cal_12_longb, .id = "temperature")

all_cals$temperature[all_cals$temperature == "1"] <- "24"
all_cals$temperature[all_cals$temperature == "2"] <- "20"
all_cals$temperature[all_cals$temperature == "3"] <- "16"
all_cals$temperature[all_cals$temperature == "4"] <- "12"
all_cals$temperature[all_cals$temperature == "4"] <- "12"

all_cals %>% 
	filter(channel != "C2") %>% 
	filter(!is.na(phase)) %>% 
	filter(Time_in_min > 110 & Time_in_min < 121) %>% 
	select(phase, oxygen_status, temperature) %>% 
	group_by(oxygen_status, temperature) %>%
	filter(!is.na(phase)) %>% 
	summarise_each(funs(mean, std.error)) %>%
	ggplot(data = ., aes(x = oxygen_status, y = mean)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error, group = temperature, color = temperature), width=.2) +
	geom_point(size = 6) + theme_bw() + ylab("body length, cm") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))


all_cals %>%
	# filter(temperature != "12") %>% 
	# filter(Time_in_min > 110 & Time_in_min < 115) %>% 
	filter(Time_in_min > 140 & Time_in_min < 143) %>% 
	select(phase, oxygen_status, temperature) %>% 
	group_by(oxygen_status, temperature) %>%
	filter(!is.na(phase)) %>% 
	summarise_each(funs(mean, std.error)) %>% View
	ggplot(data = ., aes(x = oxygen_status, y = mean, group = temperature, color = temperature)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.1) +
	geom_point(size = 3) + 
	geom_line() +
	theme_bw() + ylab("phase angle") + xlab("oxygen status of water") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))




	ggplot(data = ., aes(x = factor(oxygen_status), y = phase, fill = factor(temperature))) + geom_boxplot()