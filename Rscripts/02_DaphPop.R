### DaphPop

library(readr)
library(tidyr)
library(dplyr)
library(ggplot2)


Daph_raw <- read_csv("data-processed/P-TEMP_DaphPop.csv")
UniqueID <- read_csv("P-TEMP-UniqueID-key.csv")

Daph_raw <- tbl_df(Daph_raw)

UniqueID_treat <- separate(UniqueID, TREATMENT_ID, c("treament", "temperature", "replicate"), remove = FALSE)

UniqueID_treat$UniqueID <- as.factor(UniqueID_treat$UniqueID)

Daph <- Daph_raw %>% 
	select(sample_date, uniqueID, daph_abundance)



Daph$UniqueID <- as.factor(as.character(Daph$uniqueID))

str(Daph_raw)

Daph <- left_join(Daph, UniqueID_treat, by = "UniqueID")

Daph %>% 
	filter(temperature.y != "20") %>% 
	ggplot(., aes(x = sample_date, y = daph_abundance, group = TREATMENT_ID, color = factor(treament))) + geom_line()

### April 25

library(devtools)
# devtools::install_github("jennybc/googlesheets")
library(googlesheets)
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)

DaphPop <- gs_title("P-TEMP_DaphPop")
Daph <- gs_read(DaphPop)

UniqueID <- read_csv("P-TEMP-UniqueID-key.csv")


UniqueID_treat <- separate(UniqueID, TREATMENT_ID, c("treatment", "temperature", "replicate"), remove = FALSE)

UniqueID_treat$UniqueID <- as.factor(UniqueID_treat$UniqueID)

Daph$UniqueID <- as.factor(Daph$uniqueID)
str(Daph)

Daph_pro <- left_join(Daph, UniqueID_treat, by = "UniqueID")

str(Daph_pro)
unique(Daph_pro$juveniles)

Daph_tot <- Daph_pro %>% 
	mutate(daph_tot = unknown + grown + juveniles) %>%
	mutate(sample_date = as.factor(sample_date))

Daph_tot$sample_date <- factor(Daph_tot$sample_date, levels = c("3/31/2016", "4/5/2016", "4/8/2016", "4/12/2016", "4/16/2016", "4/19/2016", "4/26/2016", "5/3/2016"))

hist(Daph_tot$daph_tot)

Daph_pro %>% 
	# filter(sample_date == "4/19/2016") %>% 
	ggplot(data = ., aes(x = temperature.y, y = juveniles)) + geom_boxplot() + facet_wrap( ~ sample_date)

Daph_tot %>% 
	# filter(sample_date %in% c("4/12/2016", "4/16/2016", "4/19/2016")) %>% 
	ggplot(data = ., aes(x = temperature.y, y = daph_tot)) + geom_boxplot() + facet_wrap( ~ sample_date)


View(Daph_tot)
Daph_tot %>% 
	# filter(sample_date %in% c("4/12/2016", "4/16/2016", "4/19/2016")) %>% 
	ggplot(., aes(x=sample_date, y=daph_tot, group=TREATMENT_ID, color=factor(temperature.y))) + 
	geom_point()

levels(factor(Daph_tot$treament))

daph_tot_stats <- Daph_tot %>% 
	group_by(sample_date, temperature.y) %>% 
	summarise_each(funs(mean, sd, se=sd(.)/sqrt(n())), daph_tot)
	
ggplot(data = Daph_tot, aes(x = factor(sample_date), y = daph_tot, fill = temperature.y)) + geom_boxplot() +theme_bw()


ggplot(data = daph_tot_stats, aes(x = sample_date, y = mean, group = temperature.y, color = temperature.y)) +
geom_line(size=1) +
	geom_point(size=3) + geom_pointrange()

Daph_pro %>% 
	filter(sample_date %in% c("4/12/2016", "4/16/2016", "4/19/2016", "5/3/2016")) %>% 
	filter(temperature.y == "12") %>% 
	ggplot(., aes(x=sample_date, y=juveniles, group=TREATMENT_ID, color=factor(temperature.y))) +
	geom_line(size=1) +
	geom_point(size=3)

april26 <- Daph_tot %>% 
	# mutate(daph_tot = as.numeric(daph_tot)) %>% 
	filter(temperature.y == "16") %>% 
	filter(sample_date == "4/19/2016") %>%
	select(treatment.y, daph_tot, temperature.y)


Daph_tot %>% 
	filter(sample_date == "5/3/2016") %>% 
	group_by(temperature.y) %>% 
	do(tidy(lm(daph_tot ~ treatment.y, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)") %>%  View

Daph_tot %>% 
	filter(sample_date == "5/3/2016" & temperature.y == 12) %>% 
	lm(data = ., daph_tot ~ treatment.y) %>% 
		summary()

Daph_tot <- plyr::revalue(Daph_tot$treatment.y, c(FULL = "High P", DEF = "Low P"))

Daph_tot$treatment.y <- as.factor(Daph_tot$treatment.y)
levels(Daph_tot$treatment.y)[levels(Daph_tot$treatment.y) == "DEF"] <- "Low P"
levels(Daph_tot$treatment.y)[levels(Daph_tot$treatment.y) == "FULL"] <- "High P"
str(Daph_tot)

Daph_tot %>% 
	dplyr::filter(sample_date == "5/3/2016") %>% 
	dplyr::select(daph_tot, treatment.y, temperature.y) %>% 
	group_by(temperature.y, treatment.y) %>%
	summarise_each(funs(mean,median, sd,std.error)) %>% 
	ggplot(data = ., aes(temperature.y, y = mean, group = treatment.y, color = treatment.y)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 2) + theme_bw() + ylab("population abundance") + xlab("temperature, C")



Daph_tot %>% 
	# filter(temperature.y != "20") %>% 
	ggplot(data = ., aes(x = factor(temperature.y), y = daph_tot, color = treatment.y)) + geom_boxplot() + facet_wrap( ~ sample_date, scales = "free")

april26$daph_tot <- as.numeric(april26$daph_tot)
	nut <- lm(daph_tot ~ treatment.y, data = april26)
	summary(nut)
	
	
	
	