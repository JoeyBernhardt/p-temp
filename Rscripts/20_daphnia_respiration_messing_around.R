library(tidyverse)

metabolic_rates <- read_csv("data-processed/daphnia_metabolic_rates.csv")
mr <- read_csv("data-processed/metabolic_rates_all_unknown.csv")
plot_data <- read_csv("data-processed/daph_resp_extract.csv")


mr %>% 
	distinct(drymass, mean.slope, .keep_all = TRUE) %>% 
	mutate(respiration = microbe.corr.slope*-1) %>% 
	filter(temperature.x == 12) %>% 
	ggplot(aes(x = log(drymass), y = log(respiration))) + geom_point() +
	geom_smooth(method = "lm") 

mr %>% 
	distinct(drymass, mean.slope, .keep_all = TRUE) %>% 
	mutate(respiration = microbe.corr.slope*-1) %>% 
	# filter(temperature.x == 24) %>% 
	do(tidy(lm(log(respiration) ~ log(drymass), data = .), conf.int = TRUE)) %>% View


plot_data %>% 
	mutate(mass = exp(log_mass_mg)) %>% 
	mutate(respiration_rate = exp(log_respiration_mg_hour)) %>% 
	ggplot(aes(x = mass, y = respiration_rate)) + geom_point(size = 2) +
	geom_smooth(method = "lm", color = "black") +
	theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				axis.line = element_line(color="black"), 
				panel.border = element_rect(colour = "black", fill=NA, size=1))+
	theme(text = element_text(size=16, family = "Helvetica"))  + ylab('Log metabolic rate\n (mg oxygen/L/hour)') + xlab("Log mass (mg DW)") + 
	annotate("text", x = 0.0115, y= 0.08, label = "b = 1.44, CI (1.01, 1.89)", size = 6) + scale_y_log10(breaks = scales::pretty_breaks(n = 3))+ scale_x_log10(breaks = scales::pretty_breaks(n = 3))
ggsave("p-temp-figures_files/daphnia_metabolic_rate_per_mass.pdf", width = 5, height = 4)
ggsave("p-temp-figures_files/daphnia_metabolic_rate_per_mass.png", width = 5, height = 4)

?scale_y_continuous()

plot_data %>% 
	mutate(mass = exp(log_mass_mg)) %>% 
	mutate(respiration_rate = exp(log_respiration_mg_hour)) %>% 
	do(tidy(lm(log(respiration_rate) ~ log(mass), data = .), conf.int= TRUE)) %>% View