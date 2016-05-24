### May 17 2016
### Final biomasses

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotrix)
library(MuMIn)


biomass <- read_csv("data-raw/daphnia_weights_may16.csv")
UniqueID <- read_csv("data-raw/P-TEMP-UniqueID-key.csv")

biomass <- biomass %>% 
	rename(., UniqueID = `Replicate number`)

biomass_ID <- left_join(biomass, UniqueID, by = "UniqueID") %>% 
	separate(TREATMENT_ID, c("treatment", "temperature", "replicate"))

biomass_ID$sample_weight[biomass_ID$sample_weight == "1.6343"]<-"3.0736"
biomass_ID$sample_weight[biomass_ID$sample_weight == "1.4393"]<-"0.00"
biomass_ID$sample_weight <- as.numeric(biomass_ID$sample_weight)
biomass_ID$temperature <- as.factor(biomass_ID$temperature)

biomass_ID %>% 
	# filter(UniqueID < 43) %>% 
	ggplot(data = ., aes(x = factor(temperature), y = sample_weight, fill = treatment)) +
	geom_boxplot() 

str(biomass_ID)

biomass_ID %>% 
	filter(UniqueID != "1") %>% 
	dplyr::select(sample_weight, treatment, temperature) %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean, median, sd, std.error)) %>%
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + ylab("biomass, mg DW") + xlab("temperature, C") +
	theme_bw() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))

biomass_ID$temperature <- as.numeric(biomass_ID$temperature)

mod1 <- lm(sample_weight ~ treatment + temperature, data = biomass_ID)
mod2 <- lm(sample_weight ~ treatment*temperature, data = biomass_ID)
summary(mod1)

model.sel(mod1, mod2)
