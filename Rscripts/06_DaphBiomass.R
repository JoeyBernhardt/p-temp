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
# biomass_ID$temperature <- as.factor(biomass_ID$temperature)

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

#### Activation energies ####
biomassEa <- biomass_ID %>% 
	rename(., rvalues = sample_weight, tvalues = temperature) %>% 
	filter(UniqueID != 1 & UniqueID < 43) %>% 
	filter(tvalues < 24) %>% 
	select(rvalues, tvalues) 

biomassEa$tvalues <- as.numeric(biomassEa$tvalues)
biomassEa$rvalues <- as.numeric(biomassEa$rvalues)


#assumptions: the Boltzmann-Arrhenius model is a good descriptor of the temperature response
#curve.  residuals of the BA plot (ln rate vs. 1/kT) are uncorrelated and homoscedastic

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

estimate.m(biomassEa$tvalues, biomassEa$rvalues)


####

#Boltzmann's constant in eV
k=.00008617

#convert celsius values to Kelvin
biomassEa$tvalues<-biomassEa$tvalues+273.15

#convert kelvin values to 1/kT
biomassEa$xval<-1/(k*biomassEa$tvalues)

#convert response values to log(response) values
biomassEa$yval<-log(biomassEa$rvalues)

plot(biomassEa$xval, biomassEa$yval)
mod <- lm(biomassEa$yval ~ biomassEa$xval)
summary(mod)
plot(mod)

