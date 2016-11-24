##### Phytoplankton growth rates
##### Joey Bernhardt
#### November 2016


# load libraries ----------------------------------------------------------

library(tidyverse)
library(minpack.lm)
library(broom)
library(gridExtra)
library(lubridate)
library(plotrix)

# load data ---------------------------------------------------------------

ptemp <- read_csv("data-processed/p_temp_processed.csv")
ptemp_algae <- read_csv("data-processed/p_temp_algae.csv") 



# data prep -----------------------------------------------------------

str(ptemp)

ptemp <- ptemp %>% 
	mutate(sample_date = mdy(sample_date))

ptemp$start.time <- ymd_hms("2016-03-30 14:15:43")
ptemp$time_since_innoc <- interval(ptemp$start.time, ptemp$sample_date)


ptemp2 <- ptemp %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))

### now with ptemp_algae

ptemp_algae$month_day <- NA
ptemp_algae$month_day[ptemp_algae$date == "MARCH29"] <- "2016-03-29"
ptemp_algae$month_day[ptemp_algae$date == "APRIL1"] <- "2016-04-01"
ptemp_algae$month_day[ptemp_algae$date == "APRIL5"] <- "2016-04-05"
ptemp_algae$month_day[ptemp_algae$date == "APRIL8"] <- "2016-04-08"
ptemp_algae$month_day[ptemp_algae$date == "APRIL12"] <- "2016-04-12"
ptemp_algae$month_day[ptemp_algae$date == "APRIL15"] <- "2016-04-15"
ptemp_algae$month_day[ptemp_algae$date == "APRIL19"] <- "2016-04-19"
ptemp_algae$month_day[ptemp_algae$date == "APRIL26"] <- "2016-04-26"
ptemp_algae$month_day[ptemp_algae$date == "MAY4"] <- "2016-05-04"

ptemp_algae <- ptemp_algae %>% 
	mutate(sample_date = ymd(month_day))

ptemp_algae$start.time <- ymd_hms("2016-03-28 14:15:43")
ptemp_algae$time_since_innoc <- interval(ptemp_algae$start.time, ptemp_algae$sample_date)


ptemp_algae2 <- ptemp_algae %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))



# inital plots ------------------------------------------------------------


ptemp2 %>% 
	ggplot(data = ., aes(x = days, y = algal_cell_concentration_cells_per_ml, group = unique_ID, color = factor(temperature))) + geom_line(aes(linetype = phosphorus_treatment), size = 2) +
	facet_wrap( ~ temperature) + scale_y_log10()


ptemp_algae %>% 
	filter(!grepl("C", replicate)) %>% 
	ggplot(data = ., aes(x = days, y = biovol, group = ID, color = factor(temp))) + geom_line(aes(linetype = P), size = 2) +
	facet_wrap( ~ temp, scales = "free") +
	scale_y_log10()

### find the max population size

ptemp_algae_summary <- ptemp_algae %>% 
	filter(!grepl("C", replicate)) %>% 
	group_by(temp, P, replicate) %>% 
	# ggplot(data = ., aes(x = factor(temp), y = biovol, fill = factor(P))) + geom_boxplot() +
	# scale_y_log10()
	summarise_each(funs(max, mean, std.error), biovol) %>%
	# ggplot(data = ., aes(x = temp, y = max, color = P)) + geom_point()
ggplot(data = ., aes(x = factor(temp), y = max, fill = factor(P))) + geom_boxplot() +
	scale_y_log10()

mod <- lm(log(max) ~ temp + P, data =ptemp_algae_summary)
summary(mod)

ptemp_algae_control_summary <- ptemp_algae %>% 
	filter(grepl("C", replicate)) %>% 
	group_by(temp, P, replicate) %>% 
	# ggplot(data = ., aes(x = factor(temp), y = biovol, fill = factor(P))) + geom_boxplot() +
	# scale_y_log10()
	summarise_each(funs(max, mean, std.error), biovol) %>% 
	ggplot(data = ., aes(x = factor(temp), y = max, fill = factor(P))) + geom_boxplot() +
	scale_y_log10()

ptemp_algae_control_summary$consumer <- "absent"
ptemp_algae_summary$consumer <- "present"

algae_summaries <- bind_rows(ptemp_algae_control_summary, ptemp_algae_summary)
write_csv(algae_summaries, "data-processed/algae_summaries.csv")



mod <- lm(log(max) ~ temp*P, data =ptemp_algae_control_summary)
summary(mod)


ptemp_algae %>% 
	filter(!grepl("C", replicate)) %>% 
	ggplot(data = ., aes(x = month_day, y = biovol, group = ID, color = factor(temp))) + geom_line(aes(linetype = P), size = 2) +
	facet_wrap( ~ temp) +
	scale_y_log10()

ptemp_algae %>% 
	filter(grepl("C", replicate)) %>% 
	ggplot(data = ., aes(x = month_day, y = biovol, group = ID, color = factor(temp))) + geom_point(size = 3) +
	geom_line(aes(linetype = P), size = 1) +
	facet_wrap( ~ temp) +
	scale_y_log10()


control_final_day <- ptemp_algae %>% 
	filter(grepl("C", replicate)) %>% 
	filter(month_day == "2016-05-04") %>% 
	group_by(temp, P) %>% 
	summarise_each(funs(mean, std.error), biovol) %>% 
	ggplot(data = ., aes(x = temp, y = mean, color = P)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error))
	

wdaph_final_day <- ptemp_algae %>% 
	filter(!grepl("C", replicate)) %>% 
	filter(month_day == "2016-05-04") %>% 
	group_by(temp, P) %>% 
	summarise_each(funs(mean, std.error), biovol) %>% 
	ggplot(data = ., aes(x = temp, y = mean, color = P)) + geom_point() +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error)) +
	scale_y_log10()


# logistic regressions ----------------------------------------------------

sum(is.na(ptemp2$algal_biovolume))
summary(ptemp2$algal_biovolume)

ptemp2 <- ptemp2 %>% 
	select(-time_since_innoc)

## get starting biovolumes

ptemp_algae2 %>% 
	filter(sample_date == "2016-03-29") %>% 
	summarise(mean_biovol = mean(biovol)) %>% View

## take out unique ID 37, 26, 27, 29, 45, 31, 33, 34, 35, 36

r20 <- ptemp2 %>%
	filter(phosphorus_treatment == "DEF") %>%
	filter(temperature == 24, replicate == 6) %>% 
	group_by(replicate) %>%
	do(tidy(nls(algal_cell_concentration_cells_per_ml ~ 200000 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000))))

ptemp_algae2 %>%
	# filter(phosphorus_treatment == "FULL") %>%
	# filter(!grepl('37|26|27|29|45|31|33|34|35|36', unique_ID)) %>% 
	group_by(ID, P, temp) %>%
	do(tidy(nls(biovol ~ 21435402 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=1000, minFactor=1/204800000)))) %>% 
	group_by(P, temp) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(data =., aes(x = temp, y = mean, color = P)) + geom_point() +
	geom_errorbar(aes(ymin = mean-std.error, ymax = mean + std.error))

growth_rates_r <- ptemp2 %>%
	group_by(unique_ID, phosphorus_treatment, temperature) %>%
	do(tidy(nls(algal_biovolume ~ 13830995 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=1000, minFactor=1/204800000)))) %>% 
	# filter(temperature != "24") %>% 
	mutate(inverse_temp = (-1/(.00008617*(temperature+273.15)))) %>%
	group_by(phosphorus_treatment) %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	summary
	
	
	growth_rates_r %>%
		filter(temperature != "24") %>%
		mutate(inverse_temp = (1/(.00008617*(temperature + 273.15)))) %>%
		ggplot(data = ., aes(x = inverse_temp, y = log(estimate), group = phosphorus_treatment, color = phosphorus_treatment)) + geom_point() +
		geom_smooth(method = "lm")	+ scale_x_reverse()
	
	
	