#### Daph Length weight

library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(plotrix)

lwr <- read_csv("data-raw/daph-length-weight-may31.csv")
str(lwr)

lwr %>% 
	filter(well_number != "c6" & well_number != "d4") %>% 
ggplot(data = ., aes(x = length, y = sample_weight)) + geom_point() +
	stat_smooth(method = "lm") +
	theme_bw()

mod <- lm(sample_weight ~ length, data = lwr)
summary(mod)

lwr %>% 
	filter(well_number != "c6" & well_number != "d4") %>% 
	lm(sample_weight ~ length, data = .) %>% 
	summary()
