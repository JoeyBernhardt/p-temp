library(tidyverse)
library(lubridate)
library(broom)
library(plotrix)
library(stringr)




ptemp <- read_csv("/Users/Joey/Documents/p-temp/data-processed/p_temp_processed.csv")
ptemp_algae <- read_csv("/Users/Joey/Documents/p-temp/data-processed/p_temp_algae.csv") 
algae_summaries <- read_csv("/Users/Joey/Documents/p-temp/data-processed/algae_summaries.csv")

### a bit of data prep
ptemp <- ptemp %>% 
	mutate(sample_date = mdy(sample_date))

ptemp$start.time <- ymd_hms("2016-03-30 14:15:43")
ptemp$time_since_innoc <- interval(ptemp$start.time, ptemp$sample_date)


ptemp <- ptemp %>% 
	mutate(time_since_innoc_days = time_since_innoc/ddays(1)) %>% 
	mutate(time_since_innoc_hours = time_since_innoc/dhours(1))


daphnia_growth_rates <- ptemp %>%
	# filter(phosphorus_treatment == "DEF") %>%
	group_by(temperature, phosphorus_treatment, replicate) %>%
	do(tidy(nls(daphnia_total ~ 1 * (1+a)^(time_since_innoc_days),
							data= .,  start=list(a=0.01),
							control = nls.control(maxiter=100, minFactor=1/204800000)))) 


### plot daphnia growth rates Ea
daphnia_growth_rates %>% 
	filter(temperature < 21) %>% 
	ungroup() %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "DEF", "low resource supply")) %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "FULL", "high resource supply")) %>% 
	mutate(estimate = ifelse(estimate <= 0, 0.0001, estimate)) %>%
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = estimate, color = phosphorus_treatment)) +
	geom_point(size = 4) + 
	geom_abline(slope = -0.9512561, intercept = 35.9033642) +
	scale_x_reverse() +
 theme_bw() + ylab("ln daphnia population growth rate") + theme_bw() +
	geom_smooth(method = "lm") +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 16)) +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16)) + xlab("temperature (1/kT)") 
	
ggsave("p-temp-figures_files/p-temp-daphnia-growth-rate-Ea.pdf")



daphnia_growth_rates %>% 
	filter(temperature < 21) %>% 
	ungroup() %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "DEF", "low resource supply")) %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "FULL", "high resource supply")) %>% 
	mutate(estimate = ifelse(estimate <= 0, 0.0001, estimate)) %>%
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	group_by(phosphorus_treatment) %>% 
	do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View



daphnia_growth_rates %>% 
	filter(temperature < 21) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	mutate(estimate = ifelse(estimate <= 0, 0.0001, estimate)) %>% 
	group_by(phosphorus_treatment) %>% 
do(tidy(lm(log(estimate) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	filter(term != "(Intercept)") %>% View
	ggplot(aes(x = term, y = estimate, color = phosphorus_treatment)) + geom_point() +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.1)


daphnia_growth_rates %>% 
	ungroup() %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "DEF", "low resource supply")) %>% 
	mutate(phosphorus_treatment = str_replace(phosphorus_treatment, "FULL", "high resource supply")) %>% 
	mutate(estimate = ifelse(estimate <= 0, 0.001, estimate)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	group_by(inverse_temp, phosphorus_treatment) %>% 
	mutate(estimate = log(estimate)) %>% 
	summarise_each(funs(mean, std.error), estimate) %>% 
	ggplot(data = ., aes(x = inverse_temp, y = mean, group = phosphorus_treatment, color = phosphorus_treatment)) + geom_point(size = 4) +
	geom_errorbar(aes(ymin = mean - std.error, ymax = mean + std.error), width = 0.1) + ylab("ln daphnia population growth rate") + theme_bw() + 
	scale_x_reverse() + 
	theme(legend.title=element_blank(),
legend.text = element_text(size = 18)) +
	theme(axis.text=element_text(size=16),
axis.title=element_text(size=16)) + xlab("temperature (1/kT)")

### Now onto abundances!

data <- read_csv("/Users/Joey/Documents/p-temp/p-temp-marcus/plotdata.csv")

##high resource prediction
data %>% 
	mutate(temperature = transformedtemp*-1) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	filter(inverse_temp < 37.12) %>% 
	gather(key = predicted_observed, value = abundance, 8:11) %>%
	mutate(type = ifelse(grepl("predicted", predicted_observed), "predicted", "observed")) %>% 
	mutate(trophic_level = ifelse(grepl("H$", predicted_observed), "consumer", "resource")) %>% 
	filter(trophic_level == "consumer") %>% 
	filter(Phosphorus == "FULL") %>% 
	filter(type == "observed") %>%
	mutate(log_abundance = log(abundance)) %>% 
	ggplot(aes(x = temperature, y = abundance)) + geom_point(size = 4) + 
	geom_abline(slope = 1.08, intercept = 46.183708, color = "#F8766D", size = 3) +
	geom_smooth(method = "lm", color = "black", size = 2) +
	# scale_shape(solid = FALSE) +
	ylab("ln consumer abundance at day 36") + theme_bw() + scale_x_reverse() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16)) + xlab("temperature (1/kT)") +ylim(0, 5)
ggsave("p-temp-figures_files/p-temp-daphnia-abundance-highP.pdf")



##lowresource prediction
data %>% 
	mutate(temperature = transformedtemp*-1) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	filter(inverse_temp < 37.12) %>% 
	gather(key = predicted_observed, value = abundance, 8:11) %>%
	mutate(type = ifelse(grepl("predicted", predicted_observed), "predicted", "observed")) %>% 
	mutate(trophic_level = ifelse(grepl("H$", predicted_observed), "consumer", "resource")) %>% 
	filter(trophic_level == "consumer") %>% 
	filter(Phosphorus == "DEF") %>% 
	filter(type == "observed") %>%
	mutate(log_abundance = log(abundance)) %>% 
	ggplot(aes(x = temperature, y = abundance)) + geom_point(size = 4) + 
	geom_abline(slope = 1.27572, intercept = 53.15588, color = "#00BFC4", size = 3) +
	geom_smooth(method = "lm", color = "black", size = 2) +
	# scale_shape(solid = FALSE) +
	ylab("ln consumer abundance at day 36") + theme_bw() + scale_x_reverse() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16)) + xlab("temperature (1/kT)") + ylim(0, 5)
ggsave("p-temp-figures_files/p-temp-daphnia-abundance-lowP.pdf")




data %>% 
	mutate(Phosphorus = str_replace(Phosphorus, "DEF", "low resource supply")) %>% 
	mutate(Phosphorus = str_replace(Phosphorus, "FULL", "high resource supply")) %>% 
	mutate(temperature = transformedtemp*-1) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	filter(inverse_temp < 37.12) %>% 
	gather(key = predicted_observed, value = abundance, 8:11) %>%
	mutate(type = ifelse(grepl("predicted", predicted_observed), "predicted", "observed")) %>% 
	mutate(trophic_level = ifelse(grepl("H$", predicted_observed), "consumer", "resource")) %>% 
	filter(trophic_level == "consumer") %>% 
	# filter(Phosphorus == "DEF") %>% 
	filter(type == "observed") %>%
	mutate(log_abundance = log(abundance)) %>% 
	ggplot(aes(x = temperature, y = abundance, color = Phosphorus)) + geom_point(size = 4) + 
	geom_smooth(method = "lm") +
	geom_abline(slope = 1.27572, intercept = 51.95588, color = "#00BFC4", size = 3) +
	geom_abline(slope = 1.08, intercept = 45.683708, color = "#F8766D", size = 3) +
	ylab("ln consumer abundance at day 36") + theme_bw() + scale_x_reverse() +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 18)) +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16)) + xlab("temperature (1/kT)")




data %>% 
	mutate(temperature = transformedtemp*-1) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	filter(inverse_temp < 37.12) %>% 
	gather(key = predicted_observed, value = abundance, 8:11) %>%
	mutate(type = ifelse(grepl("predicted", predicted_observed), "predicted", "observed")) %>% 
	mutate(trophic_level = ifelse(grepl("H$", predicted_observed), "consumer", "resource")) %>% 
	filter(trophic_level == "consumer") %>% 
	filter(Phosphorus == "DEF") %>% 
	filter(type == "predicted") %>%
	mutate(log_abundance = log(abundance)) %>% 
	do(tidy(lm(log_abundance ~ temperature, data = .), conf.int = TRUE)) %>% View


