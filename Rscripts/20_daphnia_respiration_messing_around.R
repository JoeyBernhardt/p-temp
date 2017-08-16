

metabolic_rates <- read_csv("data-processed/daphnia_metabolic_rates.csv")
mr <- read_csv("data-processed/metabolic_rates_all_unknown.csv")


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