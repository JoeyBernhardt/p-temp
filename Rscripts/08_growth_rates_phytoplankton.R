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


# inital plots ------------------------------------------------------------


ptemp2 %>% 
	ggplot(data = ., aes(x = time_since_innoc_days, y = algal_cell_concentration_cells_per_ml, group = unique_ID, color = factor(temperature))) + geom_line(aes(linetype = phosphorus_treatment), size = 2) +
	facet_wrap( ~ temperature) + scale_y_log10()



# logistic regressions ----------------------------------------------------

sum(is.na(ptemp2$algal_biovolume))
summary(ptemp2$algal_biovolume)

ptemp2 <- ptemp2 %>% 
	select(-time_since_innoc)

## get starting biovolumes

ptemp2 %>% 
	filter(date == "01-Apr") %>% 
	summarise(mean_biovol = mean(algal_biovolume)) %>% View

## take out unique ID 37, 26, 27, 29, 45, 31, 33, 34, 35, 36

r20 <- ptemp2 %>%
	filter(phosphorus_treatment == "DEF") %>%
	filter(temperature == 24, replicate == 6) %>% 
	group_by(replicate) %>%
	do(tidy(nls(algal_cell_concentration_cells_per_ml ~ 200000 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000))))

r20 <- ptemp2 %>%
	# filter(phosphorus_treatment == "FULL") %>%
	filter(!grepl('37|26|27|29|45|31|33|34|35|36', unique_ID)) %>% 
	group_by(unique_ID, phosphorus_treatment, temperature) %>%
	do(tidy(nls(algal_biovolume ~ 13830995 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) %>% 
	group_by(phosphorus_treatment, temperature) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(data =., aes(x = temperature, y = mean, color = phosphorus_treatment)) + geom_point() +
	geom_errorbar(aes(ymin = mean-std.error, ymax = mean + std.error))


