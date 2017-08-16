### resp runs May 31 2016


library(readr)
library(ggplot2)
library(tidyr)
library(dplyr)
library(broom)
library(plotrix)
library(stringr)

resp_16 <- read_csv("respirometry-data/daphrespinmg/16C_daph_resp_mg_Oxygen.csv")
resp_12 <- read_csv("respirometry-data/daphrespinmg/12C_daph_resp_mg_Oxygen.csv")
resp_24 <- read_csv("respirometry-data/daphrespinmg/24C_daph_resp_may31_mg_Oxygen.csv")
resp_20 <- read_csv("respirometry-data/daphrespinmg/20C_daph_resp_mg_Oxygen.csv")


resp_16$temperature <- "16"
resp_12$temperature <- "12"
resp_20$temperature <- "20"
resp_24$temperature <- "24"

resp <- bind_rows(resp_24, resp_20, resp_16, resp_12)

weights <- read_csv("data-raw/daph-length-weight.csv")

weights <- weights %>% 
	mutate(length_mm = length/1000) %>% 
	mutate(drymass = 0.00402*((length_mm)^2.66)) %>% 
	unite("uniqueID", temperature, well_number, remove = FALSE) %>%
	select(-temperature) %>% 
	select(-well_number) %>% 
	filter(row_number() < 77) ## since I'm just using the may 31 24C run, not the June 1, I'm getting rid of the June 1 24C rows


resp <- resp %>% 
	select(-C2)

resp_long <- gather(resp, well_number, oxygen, 3:25) %>%
unite("uniqueID", temperature, well_number, remove = FALSE)

### change the weights df to get rid of the C after temp, which is preventing merging

weights <- weights %>% 
	mutate(uniqueID = str_replace(uniqueID, "^16C", "16")) %>% 
	mutate(uniqueID = str_replace(uniqueID, "^12C", "12")) %>% 
	mutate(uniqueID = str_replace(uniqueID, "^20C", "20")) %>% 
	mutate(uniqueID = str_replace(uniqueID, "^24C", "24"))


resp_weights <- left_join(resp_long, weights, by = "uniqueID")

resp_weights <- resp_weights %>% 
	filter(uniqueID != "20_A3")
	

resp_weights %>%
	filter(temperature == "24") %>% 
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



resp_weights %>%
	filter(temperature == "12") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2","C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D6")) %>%
	group_by(well_number, drymass) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	mutate(met.rate = mean.slope*-1) %>%
# lm(log(met.rate) ~ log(drymass), data = .) %>% 
# summary() %>% 
	ggplot(data =., aes(x = drymass, y = met.rate)) + geom_point() +
	geom_smooth(method = "lm")





resp_weights %>%
	filter(temperature == "24") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C2","C3", "C5", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% View
	summarise(mean.slope2 = mean(mean.slope)) %>% 
	View

#### checking for relationship between mass and metabolic rate

	resp_weights %>%
		filter(temperature == "24") %>% 
		filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C2","C3", "C5", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
		group_by(well_number, drymass) %>% 
		filter(Time_in_min > 10 & Time_in_min < 100) %>% 
		do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
		filter(term != "(Intercept)") %>% 
		summarise(mean.slope = mean(estimate)) %>% 
		mutate(met.rate = mean.slope*-1) %>%
# 		lm(log(met.rate) ~ log(drymass), data = .) %>% 
# 		tidy(., conf.int = TRUE) %>% View
		ggplot(data =., aes(x = log(drymass), y = log(met.rate))) + geom_point() +
		geom_smooth(method = "lm")

	
	
	
	
	
	
	
	
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

### checking mass dependence of mr
resp_weights %>%
	filter(temperature == "20") %>% 
	filter(well_number %in% c("B3", "B4", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number, drymass) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>%
	mutate(met.rate = mean.slope*-1) %>%
	lm(log(met.rate) ~ log(drymass), data = .) %>% 
	summary()


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

### checking for mass dependence of metabolic rate
resp_weights %>%
	filter(temperature == "16") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B4", "B5", "B6", "C3", "C4", "C6", "D1", "D2", "D3", "D4","D5", "D6")) %>%
	group_by(well_number, drymass) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>% 
	mutate(met.rate = mean.slope*-1) %>%
	ggplot(data =., aes(x = drymass, y = met.rate)) + geom_point()
	lm(log(met.rate) ~ log(drymass), data = .) %>% 
	summary()




#### control slopes

str(resp_weights)
280*2

control.slopes <- resp_weights %>% 
	dplyr::filter(well_number %in% c("A2", "A3", "A4")) %>% ## I'm filtering out only the A row, b/c these were my plain COMBO wells
	group_by(well_number, temperature) %>% 
	filter(Time_in_min > 5) %>% ## here I filter out the chunk of time between 5 minutes and one hour
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% ## here I fit a linear model of oxygen concentration as a function of time, for each well
	filter(term != "(Intercept)") %>% ## get rid of the intercept estimate
	summarise(mean.slope = mean(estimate)) %>% ## create a new column for the mean slope
	group_by(temperature) %>% 
	summarise(mean.control.slope = mean(mean.slope))

control.slopes.mat <- as.matrix(control.slopes)

control.slope.16 <- control.slopes.mat[2,2]
control.slope.16 <- as.numeric(control.slope.16)

control.slope.20 <- control.slopes.mat[3,2]
control.slope.20 <- as.numeric(control.slope.20)

control.slope.24 <- control.slopes.mat[4,2]
control.slope.24 <- as.numeric(control.slope.24)

control.slope.12 <- control.slopes.mat[1,2]
control.slope.12 <- as.numeric(control.slope.12)


#### 16C slopes ####
slopes.16 <- resp_weights %>%
	filter(temperature == "16") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B4", "B5", "B6", "C3", "C4", "C6", "D1", "D2", "D3", "D4","D5", "D6")) %>%
	group_by(well_number, temperature, uniqueID) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>%
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>%
	mutate(microbe.corr.slope = mean.slope - control.slope.16) 


slopes.16.weights <- left_join(slopes.16, resp_weights, by = "uniqueID")

slopes.16.weights <- as.data.frame(slopes.16.weights)

slopes.16.weights %>% 
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) * 0.0002 *60 * -1) %>% View
	
#### 24C slopes

slopes.24 <- resp_weights %>%
	filter(temperature == "24") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6","C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number, temperature, uniqueID) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>%
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>%
	mutate(microbe.corr.slope = mean.slope - control.slope.16) 


slopes.24.weights <- left_join(slopes.24, resp_weights, by = "uniqueID")

slopes.24.weights <- as.data.frame(slopes.24.weights)

slopes.24.weights %>% 
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) * 0.0002 *60 * -1) %>%
	ggplot(data =., aes(cons_per_hour)) + geom_density()

#### 12C slopes ####

slopes.12 <- resp_weights %>%
	filter(temperature == "12") %>% 
	filter(well_number %in% c("A5", "A6", "B1", "B2", "B3", "B4", "B5", "B6", "C1", "C2","C3", "C4", "C5", "C6", "D1", "D2", "D3", "D4", "D6")) %>%
	group_by(well_number, temperature, uniqueID) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>%
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>%
	mutate(microbe.corr.slope = mean.slope - control.slope.16) 


slopes.12.weights <- left_join(slopes.12, resp_weights, by = "uniqueID")

slopes.12.weights <- as.data.frame(slopes.12.weights)

slopes.12.weights %>% 
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) * 0.0002 *60 * -1) %>%
	ggplot(data =., aes(cons_per_hour)) + geom_density()


#### 20C slopes ####

slopes.20 <- resp_weights %>%
	filter(temperature == "20") %>% 
	filter(well_number %in% c("B3", "B4", "C6", "D2", "D3", "D4", "D5", "D6")) %>%
	group_by(well_number, temperature, uniqueID) %>% 
	filter(Time_in_min > 10 & Time_in_min < 100) %>%
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	summarise(mean.slope = mean(estimate)) %>%
	mutate(microbe.corr.slope = mean.slope - control.slope.16) 


slopes.20.weights <- left_join(slopes.20, resp_weights, by = "uniqueID")

slopes.20.weights <- as.data.frame(slopes.20.weights)

slopes.20.weights %>% 
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) * 0.0002 *60 * -1) %>%
	ggplot(data =., aes(cons_per_hour)) + geom_density()

all.but.12 <- bind_rows(slopes.20.weights, slopes.24.weights, slopes.16.weights)
all <- bind_rows(slopes.20.weights, slopes.24.weights, slopes.16.weights, slopes.12.weights)
write_csv(all, "data-processed/metabolic_rates_all_unknown.csv")

all.but.12 %>% 
	filter(temperature.x == "24") %>% View

resp.mass <- all %>%
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	filter(!is.na(mass)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) *60 * 0.0002* -1) %>% 
	select(temperature.x, cons_per_hour, well_number.x, uniqueID) %>% 
	group_by(temperature.x, well_number.x, uniqueID) %>% 
	summarise_each(funs(mean,median, sd,std.error)) %>%
	ungroup() %>% 
	mutate(temperature = as.numeric(temperature.x)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+ 273.15)))) %>% 
ggplot(data = ., aes(y = mean, x = inverse_temp)) +
	# geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 4, color = "blue", alpha = 0.5) +
	geom_smooth(method = "lm", color = "blue", size = 2) +
	scale_x_reverse() + 
	theme_bw() + ylab("mass-normalized oxygen flux") + xlab("temperature (1/kT)") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) + 
	theme(axis.text.y   = element_text(size=20),
				axis.text.x   = element_text(size=20),
				axis.title.y  = element_text(size=20),
				axis.title.x  = element_text(size=20),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.ticks = element_line(size = 1)) +
	theme(panel.border = element_blank(), axis.line = element_line(colour="black", size=1, lineend="square"))


ggsave("daphnia_metabolic_rate.png")

### convert the 12C values which are in air sat to torr

resp.mass$mean[resp.mass$uniqueID == "12_A5"] <- 0.12
resp.mass$mean[resp.mass$uniqueID == "12_A6"] <- 0.19
resp.mass$mean[resp.mass$uniqueID == "12_B1"] <- 0.15
resp.mass$mean[resp.mass$uniqueID == "12_B2"] <- 0.11
resp.mass$mean[resp.mass$uniqueID == "12_B3"] <- 0.11
resp.mass$mean[resp.mass$uniqueID == "12_B4"] <- 0.12
resp.mass$mean[resp.mass$uniqueID == "12_B5"] <- 0.13
resp.mass$mean[resp.mass$uniqueID == "12_B6"] <- 0.16
resp.mass$mean[resp.mass$uniqueID == "12_C1"] <- 0.12
resp.mass$mean[resp.mass$uniqueID == "12_C3"] <- 0.21
resp.mass$mean[resp.mass$uniqueID == "12_C4"] <- 0.12
resp.mass$mean[resp.mass$uniqueID == "12_C5"] <- 0.09
resp.mass$mean[resp.mass$uniqueID == "12_C6"] <- 0.07
resp.mass$mean[resp.mass$uniqueID == "12_D1"] <- 0.11
resp.mass$mean[resp.mass$uniqueID == "12_D2"] <- 0.18
resp.mass$mean[resp.mass$uniqueID == "12_D3"] <- 0.09
resp.mass$mean[resp.mass$uniqueID == "12_D4"] <- 0.14
resp.mass$mean[resp.mass$uniqueID == "12_D6"] <- 0.08


str(resp.mass)

resp.mass$temperature.x <- as.numeric(as.character(resp.mass$temperature.x))
resp.mass$temp.k <- resp.mass$temperature.x + 273.15

resp.mass <- resp.mass %>% 
	mutate(temp.inv = (1/(.00008617*(temperature.x+ 273.15))))


write_csv(resp.mass, "data-processed/daphnia_metabolic_rates.csv")

tidy(lm(log(mean) ~ inverse_temp, data = resp.mass), conf.int = TRUE) %>% 
	View

ggplot(aes(x = temp.inv, y = log(mean)), data = resp.mass) + geom_point()



estimate.m(resp.mass$temperature.x, resp.mass$mean)

resp.mass <- as.data.frame(resp.mass)

resp.mass$temperature.x <- as.numeric(resp.mass$temperature.x)

str(resp.mass)

resp.mass %>% 
	mutate(inv_temp = (1/(.00008617*(temperature.x + 273.15)))) %>% 
ggplot(data = ., aes(x = inv_temp, y = log(mean))) +
			 	geom_point(size = 4) + geom_smooth(method = "lm") + theme_bw() +
	xlab("temperature, 1/kt") + ylab("log mass-normalized oxygen flux, torr") +
# 	theme(axis.text=element_text(size=16, color = "white"),
# 				axis.title=element_text(size=16,face="bold", color = "white")) +
theme_bw() +
	theme(
		# panel.border = element_blank(),
		# legend.key = element_blank(),
		axis.ticks = element_line(color = "white"),
		# axis.text.y = element_blank(),
		# axis.text.x = element_blank(),
		# panel.grid = element_blank(),
		panel.grid.minor = element_blank(), 
		panel.grid.major = element_blank(),
		panel.background = element_blank(),
		panel.border = element_rect(color = "white"),
		plot.background = element_rect(fill = "transparent",colour = NA)) + 
	theme(axis.text=element_text(size=24, color = "white"),
				axis.title=element_text(size=20,face="bold", color = "white")) + 
	scale_x_reverse()

ggsave("test.png", bg = "transparent")





resp.mass <- resp.mass %>% 
	mutate(temp.kelvin = (1/(.00008617*(temperature.x+ 273.15))))


### convert torr to mg/l

library(purrr)
resp.k <- resp.mass %>% 
	dplyr::select(temperature.x, mean, uniqueID, temp.kelvin) %>% 
	by_row(~ convert_torr_to_mg_per_litre(.x$temperature.x, .x$mean), .to = "oxygen_in_mg") %>%
	mutate(oxygen_in_mg = as.numeric(as.character(oxygen_in_mg))) %>%
	ggplot(data = ., aes(x = temp.kelvin, y = log(oxygen_in_mg))) + geom_point() + geom_smooth(method = "lm")

mod <- lm(log(oxygen_in_mg) ~ temp.kelvin, data = resp.k)
summary(mod)

tidy(mod, conf.int = TRUE) %>% 
	View



resp.mass %>% 
	select(temp.kelvin, mean, temperature.x) %>% 
	group_by(temperature.x) %>% 
summarise_each(funs(mean,median, sd,std.error)) %>% 
	ggplot(data = ., aes(x = temp.kelvin_mean, y = mean_mean)) +
	geom_errorbar(aes(ymin=mean_mean-mean_std.error, ymax=mean_mean+mean_std.error), width=.2) +
	geom_point(size = 2) + theme_bw() + ylab("oxygen flux, torr") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) + scale_y_log10()



mod <- lm(data = resp.mass, log(mean) ~ temp.kelvin)
summary(mod)
tidy(mod, conf.int = TRUE) %>% 
	View



all$temperature.x <- as.numeric(all$temperature.x)

 all <- all %>% 
	mutate(mass = 0.00402*(length_mm^2.66)) %>% 
	filter(!is.na(mass)) %>% 
	mutate(cons_per_hour = (microbe.corr.slope/mass) *60 * 0.0002* -1) %>% 
	select(temperature.x, cons_per_hour, well_number.x) %>% 
	group_by(temperature.x, well_number.x) %>% 
	ggplot(data = ., aes(y = cons_per_hour, x = temperature.x)) +
	# geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 2) + theme_bw() + ylab("oxygen flux, torr") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	geom_smooth(method = "lm")



oxygen_in_mg <- c("0.009", "0.020", "0.024", "0.029")
std.err.mg <- c("0.000", "0.000", "0.000", "0.000")
temperature <- c("12", "16", "20", "24")

df <- data.frame(temperature, oxygen_in_mg, std.err.mg)

str(df)
df$std.err.mg <- as.numeric(as.character(df$std.err.mg))

df %>% 
	ggplot(data = ., aes(x = temperature, y = oxygen_in_mg)) +
	geom_errorbar(aes(ymin=oxygen_in_mg-std.err.mg, ymax=oxygen_in_mg+std.err.mg), width=.2) +
	geom_point(size = 6) + theme_bw() + 
	geom_smooth(method = "lm") +
	ylab("oxygen flux, mg/l*mgDW") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))


estimate.m<-function(tvalues,rvalues)
{    
	
	t.length<-length(tvalues)
	r.length<-length(rvalues)
	
	#Boltzmann's constant in eV
	k=.00008617
	
	#convert celsius values to Kelvin
	tvalues<-tvalues+273.15
	
	#convert kelvin values to 1/kT
	xval<-1/(k*tvalues)
	
	#convert response values to log(response) values
	yval<-log(rvalues)
	
	num.sum<-0
	denom.sum<-0
	x.mean<-mean(xval)
	y.mean<-mean(yval)
	
	#calculates activation energy
	
	for(i in 1:t.length)
	{
		num<-(xval[i]-x.mean)*(yval[i]-y.mean)
		denom<-(xval[i]-x.mean)*(xval[i]-x.mean)
		num.sum<-num.sum+num
		denom.sum<-denom.sum+denom
	}
	
	m<-num.sum/denom.sum  
	
	#calculates reciprocal temperature range (rtr) i.e., (1/kTmin)-(1/kTmax)
	kTmax.inv<-min(xval)
	kTmin.inv<-max(xval)
	rtr<-kTmin.inv-kTmax.inv
	
	#calculates reciprocal temperature range-location (rtr.l) i.e., mean(1/kT)
	rtr.l<-mean(xval)
	
	#calculates reciprocal temperature range-spread (rtr.s) i.e., sum(xi-xbar)^2
	
	rtr.s<-0
	
	for(j in 1:t.length)
	{
		rtr.s.part<-(xval[j]-x.mean)*(xval[j]-x.mean)
		rtr.s<-rtr.s+rtr.s.part
	}
	
	#calculates standard error
	resid<-residuals(lm(yval~xval))
	se.num<-sum(resid*resid)
	se.denom<-rtr.s*(t.length-2)
	se<-sqrt(se.num/se.denom)
	
	#return calculated values to user
	print(c("Estimated activation energy (m)",m))
	print(c("Reciprocal temperature range",rtr))
	print(c("Reciprocal temperature range-location",rtr.l))
	print(c("Reciprocal temperature range-spread",rtr.s))
	print(c("Standard error (i.e., precision)",se))
}


resp.mass %>% 
	group_by(well_number.x, temperature.x) %>% 
	summarise(mean.cons = mean(cons_per_hour)) %>% View

resp.mass$temperature.x <- as.numeric(resp.mass$temperature.x)
estimate.m(resp.mass$temperature.x, resp.mass$cons_per_hour)

str(df)
df$temperature <- as.numeric(as.character(df$temperature))
df$oxygen_in_mg <- as.numeric(as.character(df$oxygen_in_mg))
estimate.m(df$temperature, df$oxygen_in_mg)

mod <- lm(oxygen_in_mg ~ temperature, data = df)
summary(mod)


#### Microbial respiration ####

mresp_16 <- read_csv("respirometry-data/16C_microbe_resp_June22016_Oxygen.csv")
mresp_12 <- read_csv("respirometry-data/12C_microbe_resp_June22016_Oxygen.csv")
mresp_24 <- read_csv("respirometry-data/24C_microbe_resp_June222016_Oxygen.csv")
mresp_20 <- read_csv("respirometry-data/20C_microbe_resp_June22016_Oxygen.csv")


mresp_16$temperature <- "16"
mresp_12$temperature <- "12"
mresp_20$temperature <- "20"
mresp_24$temperature <- "24"

mresp <- bind_rows(mresp_24, mresp_20, mresp_16, mresp_12) %>% 
	select(-C2)
mresp_long <- gather(mresp, well_number, oxygen, 3:25) %>%
	unite("uniqueID", temperature, well_number, remove = FALSE)

mresp_long %>% 
	filter(temperature == "24") %>% 
	filter(Time_in_min > 10) %>% 
	ggplot(data =., aes(x = Time_in_min, y = oxygen)) + geom_point() +
	geom_smooth(method = "lm") +
	facet_wrap( ~ well_number)

well_IDs <- read_csv("data-raw/microbe-resp-well-IDs-June-2-2016.csv")

m_resp <- left_join(mresp_long, well_IDs, by = "well_number") %>% 
	rename(uniquewell = uniqueID) %>% 
	rename(uniqueID = UniqueID) %>% 
	mutate(uniqueID = as.character(uniqueID))

sep_June2 <- read_csv("data-processed/ptemp_summaries_June2.csv")

sep_June2 <- sep_June2 %>% 
	mutate(uniqueID = as.character(uniqueID))

m_resp_biovol <- left_join(m_resp, sep_June2, by = "uniqueID")

m_resp_biovol %>% 
	mutate(treatment = str_replace(treatment, "DEF", "low phosporus")) %>% 
	mutate(treatment = str_replace(treatment, "FULL", "high phosporus")) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', uniquewell)) %>% 
	group_by(uniqueID, temperature, treatment) %>% 
	filter(Time_in_min > 20) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(temperature = as.numeric(temperature)) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
	# dplyr::arrange(conf.high) %>% 
	# group_by(temperature, treatment) %>% View
	# summarise(mean.slope = mean(estimate)) %>%
	ggplot(data = ., aes(x = inverse_temp, y = log(estimate*-1), group = treatment, color = treatment)) +
	geom_point(size = 4) +
	geom_smooth(method = "lm") + scale_x_reverse() + ylab("log microbial oxygen flux") + xlab("inverse temperature (1/kT)") +
	theme(axis.text.y   = element_text(size=20),
				axis.text.x   = element_text(size=20),
				axis.title.y  = element_text(size=20),
				axis.title.x  = element_text(size=20),
				panel.background = element_blank(),
				panel.grid.major = element_blank(), 
				panel.grid.minor = element_blank(),
				axis.line = element_line(colour = "black"),
				axis.ticks = element_line(size = 1)) +
	theme(panel.border = element_blank(), axis.line = element_line(colour="black", size=1, lineend="square"))

ggsave("microbial_respiration.png")


	m_resp_biovol %>% 
		filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
		filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', uniquewell)) %>% 
		group_by(uniqueID, temperature, treatment) %>% 
		filter(Time_in_min > 20) %>% 
		do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
		filter(term != "(Intercept)") %>% 
		as.data.frame() %>% 
		group_by(treatment) %>% 
		mutate(temperature = as.numeric(temperature)) %>% 
		mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>% 
		do(tidy(lm(log(estimate*-1) ~ inverse_temp, data = .), conf.int = TRUE)) %>% View
	
	
	
	
	
	
m_resp_biov <- m_resp_biovol %>% 
	mutate(total_biovol_pwell = (cell_count*biovolume)/5) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', uniquewell)) %>% 
	group_by(uniqueID, temperature, treatment, total_biovol_pwell) %>% 
	filter(Time_in_min > 20) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(temperature = as.numeric(as.character(temperature))) %>% 
	mutate(mass_corr_slope = (estimate/total_biovol_pwell)*-1) %>% 
	mutate(inverse.temp = (1/(.00008617*(temperature+273.15)))) %>%
	# group_by(inverse.temp, treatment) %>% 
	# summarise(mean.slope = mean(mass_corr_slope), std.err.slope = std.error(mass_corr_slope)) %>% 
	ggplot(data = ., aes(x = inverse.temp, y = log(mass_corr_slope), group = treatment, color = treatment)) +
	geom_point() +
	geom_smooth(method = "lm")
geom_errorbar(aes(ymin=mean.slope - std.err.slope, ymax = mean.slope + std.err.slope), width=.2)


m_resp_biov %>% 
	group_by(treatment) %>% 
	summarise(mean = mean(total_biovol_pwell), std.err = std.error(total_biovol_pwell)) %>% 
	ggplot(data = ., aes(x = treatment, y = mean)) + geom_point() +
	geom_errorbar(aes(ymin=mean- std.err, ymax = mean + std.err), width=.2)



library(purrr)
str(m_resp)
as.data.frame(m_resp)

m_resp %>% 
	by_row(.d = ., ..f = convert_torr_to_mg_per_litre, temperature, oxygen_in_torr = oxygen, .to = "oxygen_in_mg") %>% View

mo <- m_resp %>% 
	select(temperature, oxygen)

mo$oxygen_in_mg <- mo %>% 
	pmap(convert_torr_to_mg_per_litre)

str(as.data.frame(mo))

m_resp2 <- bind_cols(m_resp, mo)




View(mo)
#### This didnt work, come back to it later
m_resp %>% 
	dplyr::select(temperature, oxygen, uniqueID) %>% 
	by_row(~ convert_torr_to_mg_per_litre(.x$temperature, .x$oxygen), .to = "oxygen in mg") %>%
	filter(temperature == 12) %>% View

?by_row
