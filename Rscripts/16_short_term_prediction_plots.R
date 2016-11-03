#### Plotting short term predictions and estimated activation energies


# load libraries ----------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(broom)


# read in data ------------------------------------------------------------

all_times <- read_csv("data-processed/CR_abundances_30days.csv")


# initial plots -----------------------------------------------------------


p_plot <-	all_times %>% 
	filter(temperature %in% c("12", "16", "20", "24")) %>% 
	ggplot(data = ., aes(x = time, y = P, color = factor(temperature))) +geom_point() +
	facet_wrap( ~ resource_level) +
	ggtitle("resource density") +
	scale_y_log10() +
	theme_minimal()

h_plot <-	all_times %>% 
	filter(temperature %in% c("12", "16", "20", "24")) %>% 
	ggplot(data = ., aes(x = time, y = H, color = factor(temperature))) +geom_point() +
	facet_wrap( ~ resource_level) +
	ggtitle("consumer density") +
	# scale_y_log10() +
	theme_minimal()
# Display both plots
grid.arrange(p_plot, h_plot, nrow = 2)
ggsave("p-temp-figures_files/CR-densities-over30days.png")


all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = H, group = resource_level, color = resource_level)) + geom_line(size = 3) +
	scale_x_reverse() + 
	scale_y_log10() +
	theme_minimal() + xlab("temperature (1/kT)") + ylab("log(consumer abundance)") + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 18)) +
	theme(legend.position="top")
ggsave("p-temp-figures_files/short_term_consumer_estimates_arrhenius.png", width = 8, height = 6)

all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = P, group = resource_level, color = resource_level)) + geom_line(size = 3) +
	scale_x_reverse() + 
	scale_y_log10() +
	theme_minimal() + xlab("temperature (1/kT)") + ylab("log(consumer abundance)") + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 18)) +
	theme(legend.position="top")
ggsave("p-temp-figures_files/short_term_producer_estimates_arrhenius.png", width = 8, height = 6)



### estimating activation energies

consumer_abundance_Ea <- all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	do(tidy(lm(log(H) ~ inverse_temp, data = .), conf.int = TRUE))

producer_abundance_Ea <- all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	do(tidy(lm(log(P) ~ inverse_temp, data = .), conf.int = TRUE))

