### consumer resource ratios
### May 25 2016

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotrix)
library(MuMIn)

# df <- read_csv("p_temp_Dominik/p_temp_daphnia_algae.csv")
df_raw <- read_csv("p_temp_algae.csv")

Unique_ID_key <- read_csv("data-raw/P-TEMP-UniqueID-key.csv")

Unique_ID_treat <- separate(Unique_ID_key, TREATMENT_ID, c("treatment", "temperature", "replicate")) %>% 
	mutate(UniqueID = as.character(UniqueID))

df <- df_raw %>% 
	rename(UniqueID = ID) %>%
	mutate(UniqueID = as.character(UniqueID)) %>% 
	filter(date == "MAY4")
	

df_all <- left_join(df, Unique_ID_treat, by = "UniqueID")



df_c <- df_all %>% 
	mutate(pp_carbon = biovol*0.103) %>% 
	mutate(UniqueID = as.integer(UniqueID)) %>% 
	filter(UniqueID < 49)

biomass <- read_csv("data-raw/daphnia_weights_may16.csv")

biomass <- biomass %>% 
	rename(., UniqueID = `Replicate number`) %>% 
	mutate(UniqueID = as.integer(as.character(UniqueID)))

df_c_daph <- left_join(df_c, biomass, by = "UniqueID") %>% 
	filter(UniqueID < 49)

write_csv(df_c_daph, "biomasses.csv")



df_c_daph <- df_c_daph %>% 
	mutate(daph_c = sample_weight*0.4) %>% 
	mutate(total_algal_c = pp_carbon*250) %>% 
	mutate(CRR_c = daph_c/total_algal_c) 

df_c_daph %>% 
	filter(temperature != "24") %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	lm(CRR_c ~ temperature, data = .) %>% 
	summary()

df_c_daph %>% 
	filter(UniqueID != "1") %>% 
	# filter(treatment.x == "FULL") %>% 
	filter(temperature <24) %>% 
	dplyr::select(CRR_c, temperature, treatment) %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean, median, sd, std.error)) %>% 
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + ylab("consumer:resource biomass, mgC") + xlab("temperature, C") +
	theme_bw() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/CRR.png")

df_c_daph %>% 
	filter(UniqueID != "1") %>% 
	# filter(treatment.x == "FULL") %>% 
	filter(temperature <24) %>% 
	dplyr::select(CRR_c, temperature) %>% 
	group_by(temperature) %>%
	summarise_each(funs(mean, median, sd, std.error)) %>% 
	ggplot(data = ., aes(temperature, y = mean)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + ylab("consumer:resource biomass, mgC") + xlab("temperature, C") +
	theme_bw() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/CRR-both-treatments.png")




df_c_daph %>% 
	filter(UniqueID != "1") %>% 
	# filter(treatment.x == "FULL") %>% 
	filter(temperature.x <24) %>% 
	dplyr::select(CRR_c, temperature.x) %>% 
	group_by(temperature.x) %>% 
	summarise_each(funs(mean, median, sd, std.error)) %>% 
	ggplot(data = ., aes(temperature.x, y = mean)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + ylab("CRR_c") + xlab("temperature, C") +
	theme_bw() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) + scale_y_log10()

df_c_daph$temperature.x <- as.numeric(df_c_daph$temperature.x)
mod <- lm(log(CRR_c) ~ temperature.x, data = df_c_daph)
summary(mod)

