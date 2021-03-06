
#### Load libraries
library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)

#### Import data
resp12 <- read_csv("respirometry-data/daph_resp_april29.csv")
resp20 <- read_csv("respirometry-data/daph_resp_april28.csv")
resp20_weights <- read_csv("respirometry-data/Daph_resp_weights_April28.csv")


#### Turn into long form
resp.long12 <- gather(resp12, channel, oxygen, 3:25)
resp.long20 <- gather(resp20, channel, oxygen, 3:25)

#### create new column with test temperature
resp.long20$temperature <- 20
resp.long12$temperature <- 12

#### fiddle with the time column so that it's reported in seconds, as a numeric variable
## this is going to be the dependent variable in my regressions
resp.long20$time <- as.numeric(resp.long20$`Time/Sec.`)
resp.long12$time <- as.numeric(resp.long12$`Time/Sec.`)


#### merge the two datasets
resp <- bind_rows(resp.long20, resp.long12)
resp$temperature <- as.factor(resp$temperature)

#### Plots
ggplot(data = resp.long12, aes(x = time, y = oxygen)) + geom_point() + facet_wrap( ~ channel)
ggplot(data = resp.long20, aes(x = time, y = oxygen)) + geom_point() +
	stat_summary(fun.y= "mean", geom = "point") +
	geom_smooth(method = 'lm') + 
	facet_wrap( ~ channel)


resp.long20 %>% 
ggplot(data = ., aes(x = `Time/Sec.`, y = oxygen)) + geom_point() + facet_wrap( ~ channel)


#### Density plot of de-oxygenated water, to compare the readings at 12 and 20C
resp %>% 
	ggplot(data = ., aes(x=oxygen)) + geom_density(aes(group=temperature, colour=temperature, fill=temperature), alpha=0.3) + facet_wrap( ~ channel, scales = "free")



#### Calculate slopes to get oxygen consumption (here for 20C run only)
control.slopes <- resp.long20 %>% 
	filter(channel %in% c("A3", "A4", "A5", "A6")) %>% ## I'm filtering out only the A row, b/c these were my plain COMBO wells
	group_by(channel) %>% 
	filter(time < 3600 & time > 300) %>% ## here I filter out the chunk of time between 5 minutes and one hour
	do(tidy(lm(oxygen ~ time, data = .), conf.int = TRUE)) %>% ## here I fit a linear model of oxygen concentration as a function of time, for each well
	filter(term != "(Intercept)") %>% ## get rid of the intercept estimate
	summarise(mean.slope = mean(estimate)) ## create a new column for the mean slope

mean.control.slope <- mean(control.slopes$mean.slope) ## calculate the mean control slope

#### Calculate and manipulate daphnia respiration slopes (want to eventually end up in units of mg O2/hr)
slopes <- resp.long20 %>% 
	filter(channel %in% c("D1", "D2", "D3", "D5", "B3", "B4", "B5", "B6", "C5")) %>% ## pull out the wells where it looks like the measurements worked
	group_by(channel) %>% 
	filter(time < 3600 & time > 300) %>% ## select time chunk
	do(tidy(lm(oxygen ~ time, data = .), conf.int = TRUE)) %>% ## fit linear models, grouped by channel
	filter(term != "(Intercept)") %>% ## get rid of intercept term
	## create new variable, called "microbe.corr.slope" in which I subtract the mean of the COMBO only slopes, correcting for microbial respiration
	## now units are mg O2 / L*s
	mutate(microbe.corr.slope = estimate - mean.control.slope) %>% 
	## now multiply the mg O2 /L*s by 0.0002L, which is the volume of the wells, and by 3600s, to convert my metric in seconds to hours
	## and finally, multiply by -1 to convert negative concentrations to a positive respiration value
	## units should now be mg O2/hr
	mutate(cons_per_hour = ((microbe.corr.slope * 0.0002 *3600) * -1)) 

#### Bring in the Daphnia weights, left join with the slope estimates
mass_slopes <- left_join(slopes, resp20_weights, by = "channel")

#### Get mass-specific respiration rates
## take the respiration rates (in mg O2/hr) and divide by daphnia weights
mass_slopes %>% 
	mutate(mass.corr.cons = (cons_per_hour/daph_minus_water)) %>% View


#### Visualize the slopes estimates in a coefficient plot
resp.long20 %>% 
	filter(channel %in% c("D1", "D2", "D3", "D5", "B3", "B4", "B5", "B6")) %>% 
	group_by(channel) %>% 
	filter(time < 3600 & time > 300) %>% 
	do(tidy(lm(oxygen ~ time, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)") %>% 
	ggplot(aes(y = estimate, x = factor(channel))) + geom_point() + 
	coord_flip() +
	geom_vline(xintercept = 0) +
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high)) 


mass_slopes %>% 
	filter(daph_minus_water < 0.3) %>% 
	ggplot(aes(x = daph_minus_water, y = cons_per_hour)) + geom_point()