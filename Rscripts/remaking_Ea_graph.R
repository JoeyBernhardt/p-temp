library(dplyr)
library(ggplot2)
library(lubridate)
library(readr)
library(broom)


algae <- read_csv("p_temp_algae.csv")
daphnia <- read_csv("p_temp_daphnia_algae.csv")
algae$P[algae$P == "DEF"] <- "Nutrient limited"
algae$P[algae$P == "FULL"] <- "Nutrient replete"
daphnia$P[daphnia$P == "DEF"] <- "Nutrient limited"
daphnia$P[daphnia$P == "FULL"] <- "Nutrient replete"
#
str(daphnia)
daphnia$sample_date <- as.character(daphnia$sample_date)
daphnia$sample_date <- mdy(daphnia$sample_date)
daphnia$uniqueID <- as.factor(daphnia$ID)


library(plyr)
ddply(df,"Category",transform,
			Growth=c(NA,exp(diff(log(Value)))-1))

## finagling witht the count data so I can plot it on log scale, so changing the 0s to 1s
daphnia$dapnia_count[daphnia$dapnia_count == 0] <- 1
daphnia$days[daphnia$days == 0] <- 3
daphnia$daydiff[daphnia$daydiff == 0] <- 3### fixing a point where original data had 0 days when it should be 3

daphnia <- daphnia %>% 
	group_by(uniqueID) %>%
	arrange(days) %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>%
	as.data.frame()

## finagling witht the count data so I can plot it on log scale, so changing the 0s to 1s
daphnia$dapnia_count[daphnia$dapnia_count == 0] <- 1
daphnia$days[daphnia$days == 0] <- 3
daphnia$daydiff[daphnia$daydiff == 0] <- 3### fixing a point where original data had 0 days when it should be 3


#### daphnia final population abundance
daphnia %>% 
	filter(date == "MAY4") %>% 
	group_by(uniqueID, P, inverse_temp) %>%
	# summarise(mean.daph = mean(dapnia_count)) %>%
	ggplot(data = ., aes(x = inverse_temp, y = log(dapnia_count), group = P, color = P)) + geom_point(size = 8) +
	geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	scale_colour_manual(values= c("darkolivegreen3", "cadetblue2")) + 
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	ylab("log daphnia population abundance") + xlab("temperature, 1/kt") +
	theme_bw() +
	theme(
		axis.ticks = element_line(color = "white"),
		panel.grid.minor = element_blank(), 
		panel.grid.major = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(color = "white"),
		plot.background = element_rect(fill = "transparent",colour = NA)) + 
	theme(legend.background = element_rect(fill = "transparent"), legend.margin = unit(1, "cm")) +
	theme(legend.text = element_text(size = 40, colour = "white")) +
	theme(axis.text.x=element_text(size=40, color = "white"),
				axis.text.y=element_text(size=40, color = "white"),
				axis.title=element_text(size=40,face="bold", color = "white"),
				legend.key = element_rect(fill = "transparent", colour = "transparent")) +
	theme(legend.position="top") + 
	theme(legend.title=element_blank()) +
	scale_x_reverse()

ggsave("daph_pop.png", bg = "transparent", width = 8, height = 8)


daphnia.max.abundance.Ea <- daphnia %>% 
	group_by(uniqueID, P, inverse_temp) %>% 
	summarise(max.daph.abundance = max(dapnia_count)) %>% 
	group_by(P) %>% 
	do(tidy(lm(log(max.daph.abundance) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)") 

daphnia.final.abundance.Ea <- daphnia %>% 
	filter(date == "MAY4") %>% 
	group_by(P) %>% 
	do(tidy(lm(log(dapnia_count) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)")

daphnia.max.abundance.Ea$rate <- "consumer max abundance"
daphnia.final.abundance.Ea$rate <- "consumer final abundance"

daph.growth <- daphnia %>% 
	select(dapnia_count, uniqueID) %>%
	as.data.frame() %>%
	group_by(uniqueID) %>% 
	mutate_each(funs(. / lag(.) - 1)) %>% 
	rename(growth.rate = dapnia_count) %>%
	as.data.frame() %>%
	select(-uniqueID) %>%
	bind_cols(daphnia, .) %>% 
	# mutate(growth_rate_days = growth.rate/daydiff) %>%
	select(growth.rate, uniqueID, P, inverse_temp) %>% 
	dplyr::group_by(uniqueID, P, inverse_temp) %>% 
	filter(!is.na(growth.rate)) %>% 
	summarise_each(., funs(max, mean)) %>%
	# 		rename(max_daph_population_growth = max,
	# 					 mean_daph_population_growth = mean)
	# 		group_by(P) %>% 
	# 			do(tidy(lm(log(mean) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	# 			filter(term != "(Intercept)") %>%
	# 		ggplot(data =., aes(x = P, y = estimate)) + geom_point() + 
	# 		geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width=.2)
	ggplot(data = ., aes(x = inverse_temp, y = log(mean), group = P, color = P)) + geom_point(size = 4) + 
	geom_smooth(method = "lm", se = FALSE) + theme_bw() + 
	scale_color_manual(values = c("royalblue1", "palegreen3")) +
	theme(axis.text=element_text(size=16), axis.title=element_text(size=16,face="bold")) +
	ylab("log max population growth rate") + xlab("temperature, 1/kt")
ggsave("p-temp-figures_files/figure-html/daph_max_pop_growth_arrhenius.png")

write_csv(daph.growth, "data-processed/daph.growth.rates.csv")		


daph.growth <- read_csv("data-processed/daph.growth.rates.csv")

daph.growth.Ea <- daph.growth %>% 
	group_by(P) %>% 
	do(tidy(lm(log(max_daph_population_growth) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)")


daph.growth.Ea$rate <- "consumer population growth rate"


### algae 

library(dplyr)
algae <- algae %>% 
	dplyr::rename(uniqueID = ID) %>%
	group_by(uniqueID) %>%
	arrange(days) %>% 
	mutate(inverse_temp = (1/(.00008617*(temp+273.15)))) %>%
	as.data.frame()


algae$treatment[algae$uniqueID < 49] <- "consumers present"
algae$treatment[algae$uniqueID > 48] <- "consumers absent"
algae <- tidyr::unite(algae, "Ptreatment", P, treatment, remove = FALSE)


write_csv(algae.growth, "data-processed/algae.growth.rates.csv")		



algae.growth.rates <- read_csv("data-processed/algae.growth.rates.csv")

algae.growth.rates$treatment[algae.growth.rates$uniqueID < 49] <- "grazed"
algae.growth.rates$treatment[algae.growth.rates$uniqueID > 48] <- "ungrazed"
algae.growth.rates$P[algae.growth.rates$P == "DEF"] <- "Low P"
algae.growth.rates$P[algae.growth.rates$P == "FULL"] <- "High P"


algae.growth.rates <- tidyr::unite(algae.growth.rates, "Ptreatment", P, treatment, remove = FALSE)

algae.growth.Ea <- algae.growth.rates %>%
	group_by(Ptreatment) %>% 
	do(tidy(lm(log(max_biovol_growth) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)")

algae.growth.Ea$rate <- "resource growth rate"	

algae.biol.Ea <- algae %>% 
	group_by(uniqueID, inverse_temp, P, treatment, Ptreatment) %>% 
	summarise(max.biovol = max(biovol)) %>% 
	group_by(Ptreatment) %>% 
	do(tidy(lm(log(max.biovol) ~ inverse_temp, data = .), conf.int = TRUE)) %>%
	filter(term != "(Intercept)")

algae.biol.Ea$rate <- "resource max abundance"


estimate <- -0.7114566
conf.low <- -0.9710273
conf.high <- -0.4518859
Ptreatment <- "High phosphorus"
term <- "inverse_temp"
rate <- "daphnia metabolic rate"

daph.mr <- data.frame(Ptreatment, term, estimate, conf.low, conf.high, rate)
str(daph.mr)
bind_rows(as.data.frame(daphnia.max.abundance.Ea), as.data.frame(daph.mr))



#### Activation energies ####

Ea <- dplyr::bind_rows(as.data.frame(algae.biol.Ea), as.data.frame(algae.growth.Ea), as.data.frame(daph.mr), as.data.frame(daphnia.max.abundance.Ea), as.data.frame(daph.growth.Ea), as.data.frame(daphnia.final.abundance.Ea))

write_csv(Ea, "Ea.csv")



Ea$group <- "community context"
Ea$group[Ea$Ptreatment == "High P_ungrazed"] <- "population context"
Ea$group[Ea$Ptreatment == "Low P_ungrazed"] <- "population context"


# Ea$Ptreatment[Ea$P == "High P"] <- "High P"
# Ea$Ptreatment[Ea$P == "Low P"] <- "Low P"
# Ea$Ptreatment[Ea$P == "FULL"] <- "High P"
# Ea$Ptreatment[Ea$P == "DEF"] <- "Low P"
# Ea$Ptreatment[Ea$P == "High P"] <- "High P"
# Ea$Ptreatment[Ea$P == "Low P"] <- "Low P"
# Ea$Ptreatment[Ea$P == "High phosphorus"] <- "High P"
# Ea$Ptreatment[Ea$P == "Low phosphorus"] <- "Low P"




Ea$Ptreatment[Ea$P == "High P"] <- "nutrient replete"
Ea$Ptreatment[Ea$P == "Low P"] <- "nutrient limited"
Ea$Ptreatment[Ea$P == "FULL"] <- "nutrient replete"
Ea$Ptreatment[Ea$P == "DEF"] <- "nutrient limited"
Ea$Ptreatment[Ea$P == "High P"] <- "nutrient replete"
Ea$Ptreatment[Ea$P == "Low P"] <- "nutrient limited"
Ea$Ptreatment[Ea$P == "High phosphorus"] <- "nutrient replete"
Ea$Ptreatment[Ea$P == "Low phosphorus"] <- "nutrient limited"


Ea$Ptreatment[Ea$Ptreatment == "High P_ungrazed"] <- "nutrient replete"
Ea$Ptreatment[Ea$Ptreatment == "High phosphorus"] <- "nutrient replete"
Ea$Ptreatment[Ea$Ptreatment == "Low P_ungrazed"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "High P_grazed"] <- "nutrient replete"
Ea$Ptreatment[Ea$Ptreatment == "Low P_grazed"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "Nutrient limited_consumers absent"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "Nutrient limited_consumers present"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "Nutrient replete_consumers absent"] <- "nutrient replete"
Ea$Ptreatment[Ea$Ptreatment == "Nutrient replete_consumers present"] <- "nutrient replete"
Ea$Ptreatment[Ea$Ptreatment == "Low P_grazed"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "Low phosphorus_consumer absent"] <- "nutrient limited"
Ea$Ptreatment[Ea$Ptreatment == "Low phosphorus_consumer present"] <- "nutrient limited"
Ea$Ptreatment[Ea$P == "Nutrient limited"] <- "nutrient limited"
Ea$Ptreatment[Ea$P == "Nutrient replete"] <- "nutrient replete"

Ea$group <- as.factor(Ea$group)


pd <- position_dodge(width = 0.4)

Ea$rate


Ea$rate <- factor(Ea$rate, c("daphnia metabolic rate", "consumer population growth rate", "consumer max abundance", "resource growth rate", "resource max abundance", "consumer max abundance", "consumer final abundance"))


Ea$rate[is.na(Ea$rate)] <- "resource growth rate"

levels(Ea$rate)

Ea$rate <- as.character(Ea$rate)



## activation energies plot
Ea %>% 
	filter(rate != "consumer final abundance") %>% 
	ggplot(data =., aes(x = factor(rate), y = estimate, group = Ptreatment, color = Ptreatment)) + geom_point(aes(shape = group), size = 10, position = pd) + 
	geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width=.2, position = pd) +
	theme_bw() + 
	scale_colour_manual(values= c("darkolivegreen3", "cadetblue2")) +
	theme(axis.text.x=element_text(size=16, angle = 25, hjust = 1), axis.title=element_text(size=16,face="bold")) +
	ylab("activation energy") + xlab("rate") +
	theme(
		# panel.border = element_blank(),
		# legend.key = element_blank(),
		# axis.ticks = element_line(color = "white"),
		# axis.text.y = element_blank(),
		# axis.text.x = element_blank(),
		# panel.grid = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.grid.major = element_blank(),
		panel.background = element_blank(),
		# panel.border = element_rect(color = "white"),
		plot.background = element_rect(fill = "transparent",colour = NA)) + 
	theme(legend.background = element_rect(fill = "transparent"), legend.margin = unit(1, "cm")) +
	theme(legend.text = element_text(size = 20)) +
	theme(axis.text.x=element_text(size=30),
				axis.text.y=element_text(size=30),
				axis.title=element_text(size=30,face="bold"),
				legend.key = element_rect(fill = "transparent", colour = "transparent")) 

ggsave("activation_energies.png", bg = "transparent", width = 20, height = 8)



### second try!!!

Ea_mod <- read_csv("Ea_mod.csv")

View(Ea_mod)
str(Ea_mod)	
	

Ea_mod %>%
	filter(rate != "4. consumer final abundance") %>% 
	ggplot(data =., aes(x = factor(rate), y = estimate, group = Ptreatment, color = Ptreatment)) +
	geom_point(aes(shape = Consumer_treatment), size = 10, position = pd) + 
	geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width=.2, position = pd) +
	theme_bw() + 
	scale_colour_manual(values= c("chartreuse4", "darkblue")) +
	theme(axis.text.x=element_text(size=16, angle = 25, hjust = 1), axis.title=element_text(size=16,face="bold")) +
	ylab("activation energy (eV)") + xlab("rate") +
	theme(
		# panel.border = element_blank(),
		# legend.key = element_blank(),
		# axis.ticks = element_line(color = "white"),
		# axis.text.y = element_blank(),
		# axis.text.x = element_blank(),
		# panel.grid = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.grid.major = element_blank(),
		panel.background = element_blank(),
		# panel.border = element_rect(color = "white"),
		plot.background = element_rect(fill = "transparent",colour = NA)) + 
	theme(legend.background = element_rect(fill = "transparent"), legend.margin = unit(1, "cm")) +
	theme(legend.text = element_text(size = 20)) +
	theme(legend.position = "none") + 
	theme(axis.text.x=element_text(size=30),
				axis.text.y=element_text(size=30),
				axis.title=element_text(size=30,face="bold"),
				legend.key = element_rect(fill = "transparent", colour = "transparent")) 

ggsave("activation_energies2.png", bg = "transparent", width = 20, height = 8)
	
	
	