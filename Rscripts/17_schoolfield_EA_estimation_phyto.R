#### Fitting ptemp data to schoolfield model to estimate activation energies
#### Adapted from code by SSB
#### November 7 2016

# load libraries ----------------------------------------------------------

library(tidyverse)
library(minpack.lm)
library(broom)
library(gridExtra)
library(lubridate)

# load data ---------------------------------------------------------------

ptemp <- read_csv("data-processed/p_temp_processed.csv")
ptemp_algae <- read_csv("data-processed/p_temp_algae.csv") 



# initial plots -----------------------------------------------------------

str(ptemp)

ptemp <- ptemp %>% 
	mutate(sample_date = mdy(sample_date))

ptemp %>% 
	# group_by(unique_ID) %>% 
ggplot(data = ., aes(x = sample_date, y = algal_biovolume, group = unique_ID, color = factor(temperature))) + geom_line(aes(linetype = phosphorus_treatment), size = 2) +
 facet_wrap( ~ temperature)



# get data in order -------------------------------------------------------

current_dataset <- ptemp %>% 
	filter(date == "04-May", phosphorus_treatment == "FULL") %>%
	select(algal_biovolume, temperature) %>% 
	mutate(K = temperature + 273.15) %>% 
	rename(OriginalTraitValue = algal_biovolume) %>% 
	select(-temperature)


current_dataset_def <- ptemp %>% 
	filter(date == "04-May", phosphorus_treatment == "DEF") %>%
	select(algal_biovolume, temperature) %>% 
	mutate(K = temperature + 273.15) %>% 
	rename(OriginalTraitValue = algal_biovolume) %>% 
	select(-temperature)



current_dataset <- ptemp_algae %>% 
	filter(date == "MAY4", P == "FULL") %>%
	filter(grepl("C", replicate)) %>% 
	select(biovol, temp) %>% 
	mutate(K = temp + 273.15) %>% 
	rename(OriginalTraitValue = biovol) %>% 
	select(-temp)
	
current_dataset_def <- ptemp_algae %>% 
	filter(date == "MAY4", P == "DEF") %>%
	filter(grepl("C", replicate)) %>% 
	select(biovol, temp) %>% 
	mutate(K = temp + 273.15) %>% 
	rename(OriginalTraitValue = biovol) %>% 
	select(-temp)


current_dataset_def$OriginalTraitValue[current_dataset_def$OriginalTraitValue == 0] <- 1



# inital plots ------------------------------------------------------------

full <- ggplot(data = current_dataset, aes(x = K, y = OriginalTraitValue)) + geom_point() + scale_y_log10()
def <- ggplot(data = current_dataset_def, aes(x = K, y = OriginalTraitValue)) + geom_point() + ylim(0, 3 * 10^9)

grid.arrange(full, def, ncol = 2)


#### assign Tref as GlobalEnv
# T_ref is the standardization temperature (in K). 
# This needs to be any value below the peak of the curve.
assign("Tref", 285.15, envir = .GlobalEnv) 


# get starting values -----------------------------------------------------
GetE <- function(tmp, rate, T.p, k=8.62e-5)
{
	# Estimate starting value for E, taking linear regression using the rise part
	# of the curve only.
	# ~~~ Parameters ~~~
	# tmp  : temperature data (in K).
	# rate : rate data corresponding to temperature above.
	# T.p  : temperature at which rate peaks, used as a cutoff point.
	# k    : Boltzmann constant.
	
	tmp.w <- which(tmp <= T.p)
	if (length(tmp.w) > 1)
	{
		m <- lm(log(rate[tmp.w]) ~ I(1 / (k * (tmp[tmp.w]))))
		return(abs(summary(m)$coefficients[2, 1]))
	} else
	{
		return(0.6)
	}
}

GetB0 <- function(tmp, rate)
{
	# Estimate starting value for the normalising constant.
	# ~~~ Parameters ~~~
	# tmp   : temperature data (in K).
	# rate  : rate data corresponding to temperature above.
	# T.ref : estimate normalising constant at this temperature (in K).
	
	if (min(tmp,na.rm=TRUE) > Tref)
	{
		return(log(min(rate[1],na.rm=TRUE)))
	} else
	{
		return(log(max(rate[which(tmp <= Tref)],na.rm=TRUE)))
	}
}


GetTpk <- function(tmp, rate)
{
	# Temperature at which the rate is maximised (estimate of T.peak).
	# ~~~ Parameters ~~~
	# tmp  : Temperature data (in K).
	# rate : Rate data corresponding to temperature above.
	
	return(max(tmp[which.max(rate)]))
}


# define schoolfield function ---------------------------------------------

Schoolfield <- function(B0, E, E_D, T_h, temp) {
	
	# Boltzmann's constant. Units imply that E and E_D are in eV.
	k <- 8.62e-5
	
	# B0 is the normalization constant.    
	# E is the activation energy.
	# E_D is the de-activation energy.    
	# T_h is the temperature at which the rate-limiting enzyme 
	# is 50% active and 50% denatured due to high temperature.
	
	# return(B0 - E/k * (1/temp - 1/Tref) - log(1 + exp((E_D/k) * (1/T_h - 1/temp)))) #Schoolfied model in original form (not with T_pk as an explicit parameter)
	
	return(B0 + log(exp((-E/k) * ((1/temp) - (1/Tref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_h - 1/temp))))) ## T_pk as an explicit parameter. FITS BETTER
	
}

# get starting values -----------------------------------------------------

T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)


# get schoolfield fit ---------------------------------------------------------


schoolfield_nls_full <- nlsLM(
	log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
	start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
	lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
	upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
	data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024))

schoolfield_nls_def <- nlsLM(
	log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
	start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
	lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
	upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
	data=current_dataset_def, control=list(minFactor=1 / 2^16, maxiter=1024))

full_est <- tidy(schoolfield_nls_full) %>% 
	mutate(phosphorus = "replete")

def_est <- tidy(schoolfield_nls_def) %>% 
	mutate(phosphorus = "deficient")

all_estimates <- bind_rows(full_est, def_est)

write_csv(all_estimates, "data-processed/schoolfield_Ea_estimates_phyto.csv")

## plot the activation energy estimates for the P replete and deficient treatments
all_estimates %>% 	
filter(term == "E") %>%
ggplot(aes(x = phosphorus, y = estimate, group = phosphorus)) + geom_point() +
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error))


### bootstrapping
bootstrapped_estimates_100 <- current_dataset%>% 
	bootstrap(10) %>%
do(tidy(nlsLM(
	log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
	start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
	lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
	upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
	data=., control=list(minFactor=1 / 2^16, maxiter=1024)))) 

write_csv(bootstrapped_estimates_100, "data-processed/schoolfield_boot_estimates.csv")

## plot these babies!

tmp_temps <- seq(min(
	floor(current_dataset$K)), 
	ceiling(max(current_dataset$K)
	), length = 200)

	plot_sf <- function(df) {
		exp(Schoolfield(
		df$estimate[df$term == "E"],
		df$estimate[df$term == "B0"],
		df$estimate[df$term == "E_D"],
		df$estimate[df$term == "T_h"],
		tmp_temps))
	}

	tmp_model_boot <- bootstrapped_estimates_100 %>% 
		split(.$replicate) %>% 
		map_df(.f = plot_sf) %>% 
		gather(key = replicate, value = TraitValue, 1:100) %>% 
		mutate(Temperature = rep(tmp_temps - 273.15, 100))


ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																						 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
											alpha = 0.7, pch = 21) + 
	geom_line(data = tmp_model_boot, 
						aes(x = Temperature, y = TraitValue, group = replicate), colour = "#1b9e77",
						lwd = 1.3, alpha=.2) + 
# geom_line(data = ModelToPlotS, 
# 					aes(x = Temperature, y = TraitValue), colour = "black", 
# 					lwd = 1.3) +
	ggtitle("Phosphorus replete") +
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("daphnia population abundance @ day 30") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16))

## estimating confidence intervals using the quantile method
alpha = .05
bootnls %>% group_by(term) %>%
	summarize(low=quantile(estimate, alpha / 2),
																				 high=quantile(estimate, 1 - alpha / 2))

ggplot(bootnls, aes(estimate)) + geom_histogram() + facet_wrap(~ term, scales="free")



# Generate predictions from the model fit (non bootstrapped) ---------------------------------

tmp_temps <- seq(min(
	floor(current_dataset$K)), 
	ceiling(max(current_dataset$K)
	), length = 200)

tmp_model <- exp(Schoolfield(
	coef(schoolfield_nls_full)["B0"],
	coef(schoolfield_nls_full)["E"],
	coef(schoolfield_nls_full)["E_D"],
	coef(schoolfield_nls_full)["T_h"],
	tmp_temps
))
tmp_model

ModelToPlotS <- data.frame(
	Temperature = tmp_temps - 273.15, 
	TraitValue = tmp_model
)

DataToPlot <- data.frame(
	Temperature = current_dataset$K - 273.15, 
	TraitValue = current_dataset$OriginalTraitValue
)
DataToPlot <- na.omit(DataToPlot)



##### DEF 

tmp_model_def <- exp(Schoolfield(
	coef(schoolfield_nls_def)["B0"],
	coef(schoolfield_nls_def)["E"],
	coef(schoolfield_nls_def)["E_D"],
	coef(schoolfield_nls_def)["T_h"],
	tmp_temps
))
tmp_model

ModelToPlotS_def <- data.frame(
	Temperature = tmp_temps - 273.15, 
	TraitValue = tmp_model_def
)

DataToPlot_def <- data.frame(
	Temperature = current_dataset_def$K - 273.15, 
	TraitValue = current_dataset_def$OriginalTraitValue
)
DataToPlot_def <- na.omit(DataToPlot_def)



# plot it! ----------------------------------------------------------------

full_plot <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																						 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
											alpha = 0.7, pch = 21) + 
	geom_line(data = ModelToPlotS, 
						aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
						lwd = 1.3) +                           
	ggtitle("Phosphorus replete") +
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("daphnia population abundance @ day 30") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16))

#### DEF
def_plot <- ggplot() + geom_point(data = DataToPlot_def, aes(x = Temperature, 
																						 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
											alpha = 0.7, pch = 21) + 
	geom_line(data = ModelToPlotS_def, 
						aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
						lwd = 1.3) +                           
	ggtitle("Phosphorus deficient") +
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("daphnia population abundance @ day 30") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16)) + scale_y_log10()

grid.arrange(full_plot, def_plot, ncol = 2)

