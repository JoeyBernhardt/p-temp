### resp runs May 31 2016


library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
library(plotrix)

resp_16 <- read_csv("respirometry-data/may31resp/16_resp_may312016_Oxygen.csv")
resp_12 <- read_csv("respirometry-data/12C_daphnia_resp_June12016_Oxygen.csv")
resp_24 <- read_csv("respirometry-data/24C_resp_June12016_Oxygen.csv")
resp_20 <- read_csv("respirometry-data/may31resp/20C_resp_may312016_Oxygen.csv")


resp_16$temperature <- "16"
resp_12$temperature <- "12"
resp_20$temperature <- "20"
resp_24$temperature <- "24"

resp <- bind_rows(resp_24, resp_20, resp_16, resp_12)

weights <- read_csv("data-raw/daph-length-weight.csv")

weights <- weights %>% 
	mutate(length_mm = length/1000) %>% 
	mutate(drymass = 0.00402*((length_mm)^2.66)) %>% 
	unite("uniqueID", temperature, well_number, remove = FALSE) %>% View

?unite

resp<- resp %>% 
	select(-C2)

length(unique(weights$uniqueID))
length(unique(resp_long$uniqueID))



resp_long<- gather(resp, well_number, oxygen, 3:25) %>%
unite("uniqueID", temperature, well_number, remove = FALSE)


resp_weights <- left_join(weights, resp_long, by = "uniqueID", copy = TRUE)


resp_long %>%
	filter(temperature == "16") %>% 
	group_by(temperature) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	ggplot(data = ., aes(x = Time_in_min, y = oxygen, group = temperature, color = temperature)) + geom_point() + stat_summary(fun.y= "mean", geom = "point") +
	facet_wrap( ~ well_number) +
	geom_smooth(method = 'lm')

resp_weights %>%
	filter(temperature == "12") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2","C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D6")) %>%
	group_by(well_number) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	summarise(mean.slope2 = mean(mean.slope)) %>% 
	View


resp_long %>%
	filter(temperature == "24") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2","C3", "C5", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	summarise(mean.slope2 = mean(mean.slope)) %>% 
	View

resp_long %>%
	filter(temperature == "20") %>% 
	filter(well_number %in% c("B3", "B4", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	summarise(mean.slope2 = mean(mean.slope)) %>% 
	View

resp_long %>%
	filter(temperature == "16") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B4", "B5", "B6", "C3", "C4", "C6", "D1", "D2", "D3", "D4","D5", "D6")) %>%
	group_by(well_number) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	summarise(mean.slope2 = mean(mean.slope)) %>% 
	View

