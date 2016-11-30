## P-TEMP
## messing around with activation energy estimation
## Nov 7 2016
## JB


# load libraries ----------------------------------------------------------

library(tidyverse)
library(minpack.lm)

# read in data ------------------------------------------------------------

ptemp <- read_csv("data-processed/p_temp_processed.csv")

data_full <- data %>% 
	filter(Phosphorus == "FULL")

ggplot(data = data_full, aes(x = temp, y = r)) + geom_point()

ptemp_final <- ptemp %>% 
	filter(date == "04-May", phosphorus_treatment == "FULL") %>% 
	select(unique_ID, daphnia_total, temperature) %>% 
	mutate(tmp = temperature + 273.15)


Data <- ptemp_final %>% 
	mutate(rate = daphnia_total) %>% 
	select(-temperature) %>% 
	select(-unique_ID) %>% 
	select(-daphnia_total)


# get data in order -------------------------------------------------------

current_dataset <- ptemp %>% 
	filter(date == "04-May", phosphorus_treatment == "DEF") %>% 
	select(daphnia_total, temperature) %>% 
	mutate(K = temperature + 273.15) %>% 
	rename(OriginalTraitValue = daphnia_total) %>% 
	select(-temperature)

current_dataset$OriginalTraitValue[current_dataset$OriginalTraitValue == 0] <- 1

## If there are negative values, substract the minimum value
MinVal <- NA
if (min(current_dataset$OriginalTraitValue)<=0){
	MinVal <- min(current_dataset$OriginalTraitValue)
	current_dataset$OriginalTraitValue <-current_dataset$OriginalTraitValue - MinVal
	current_dataset <-current_dataset[-which(current_dataset$OriginalTraitValue==0),]}



#### assign Tref as GlobalEnv
# T_ref is the standardization temperature (in K). 
# This needs to be any value below the peak of the curve.
assign("Tref", 285.15, envir = .GlobalEnv) 

 

# Estimate STARTING VALUES for the nls ------------------------------------

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




# Schoolfield fitting -----------------------------------------------------

Schoolfield <- function(B0, E, E_D, T_h, temp) {
	
	# Boltzmann's constant. Units imply that E and E_D are in eV.
	k <- 8.62e-5
	
	# B0 is the normalization constant.    
	# E is the activation energy.
	# E_D is the de-activation energy.    
	# T_h is the temperature at which the rate-limiting enzyme 
	# is 50% active and 50% denatured due to high temperature.
	
	#     return(B0 - E/k * (1/temp - 1/Tref) - log(1 + exp((E_D/k) * (1/T_h - 1/temp)))) #Schoolfied model in original form (not with T_pk as an explicit parameter)
	
	return(B0 + log(exp((-E/k) * ((1/temp) - (1/Tref)))/(1 + (E/(E_D - E)) * exp(E_D/k * (1/T_h - 1/temp))))) ## T_pk as an explicit parameter. FITS BETTER
	
}


B0_sch <- c()
E_sch <- c()
E_D_sch <- c()	
T_h_sch <- c()
T_pk_sch <- c()
P_pk_sch <- c()
AIC_sch <- c()
r_sq_sch <- c()

T.h.st  <- GetTpk(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)
E.st    <- GetE(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue, T.p=T.h.st)
B.st <- GetB0(tmp=current_dataset$K, rate=current_dataset$OriginalTraitValue)



# schoolfield -------------------------------------------------------------

schoolfield_nls <- NA
	schoolfield_nls <- nlsLM(
		log(OriginalTraitValue) ~ Schoolfield(B0, E, E_D, T_h, temp = K), 
		start=c(B0 = B.st, E = E.st, E_D = 4*E.st, T_h=T.h.st),
		lower=c(B0=-Inf,   E=0,    E.D=0, Th=0),
		upper=c(B0=Inf,    E=Inf,  E.D=Inf, Th=273.15+150),
		data=current_dataset, control=list(minFactor=1 / 2^16, maxiter=1024))


if(!is.na(schoolfield_nls[1])) 
{ 
	
	# Collect the parameter estimates..
	# if (!is.na(MinVal)){ ## Add MinVal if it was substracted
	# 	B0Bug <- c(B0Bug,TRUE)
	# 	B0_sch <- c(B0_sch, (coef(schoolfield_nls)["B0"]+MinVal))
	# }else {
	# 	B0Bug <- c(B0Bug,FALSE)
	# 	B0_sch <- c(B0_sch, coef(schoolfield_nls)["B0"])
	# }
	E_sch <- c(E_sch, coef(schoolfield_nls)["E"])
	E_D_sch <- c(E_D_sch, coef(schoolfield_nls)["E_D"])
	T_h_sch <- c(T_h_sch, coef(schoolfield_nls)["T_h"])
	AIC_sch<- c(AIC_sch, AIC(schoolfield_nls))
	
	# Calculate the R squared value as: 1 - (rss/tss)
	rss <- sum((exp(predict(schoolfield_nls)) - 
								current_dataset$OriginalTraitValue)^2, 
						 na.rm = TRUE)
	tss <- sum(
		(current_dataset$OriginalTraitValue - 
		 	mean(current_dataset$OriginalTraitValue, na.rm = TRUE))^2, 
		na.rm = TRUE)
	
	if ( tss != 0 )
	{
		r_sq_sch <- c(r_sq_sch, 1 - (rss/tss))
	} else
	{
		r_sq_sch <- c(r_sq_sch, 1)
	}
	
	# Calculate the peak of the curve and its 
	# corresponding temperature value.
	curr_prediction <- predict(schoolfield_nls)
	for (j in 1:length(curr_prediction))
	{
		# If we found the maximum performance, exit the loop.
		if (curr_prediction[j] == max(curr_prediction))
		{
			break
		}
	}
	
	T_pk_sch <- c(T_pk_sch, current_dataset$K[j])
	# if (!is.na(MinVal)){ ## Add MinVal if it was substracted
	# 	P_pkBug <- c(P_pkBug,TRUE)
	# P_pk_sch <- c(P_pk_sch, (curr_prediction[j]+MinVal))
	# }else {
	# 	P_pkBug <- c(P_pkBug,FALSE)
	# 	P_pk_sch <- c(P_pk_sch, curr_prediction[j])
	# }
}	
	
	
	##############################
	# Plotting Schoolfield's fit #
	##############################
	
	# Create a name for the output file using:
	#	- the original id number
	#   - the species name
	#   - the model
	# output_name <- paste(
	# 	current_dataset$FinalID[1], 
	# 	current_dataset$Consumer[1], 
	# 	'Schoolfield',
	# 	sep = "_"
	# )
	# 
	# 
	# # Remove any characters that won't look good in a file name,
	# # using a regular expression.
	# output_name <- gsub("[^\\w|\\s](|)", "", output_name, perl=TRUE)
	# 
	# # Convert spaces to underscores.
	# output_name <- gsub("\\s+", "_", output_name, perl=TRUE)
	# 
	# # CHANGE THIS to set an alternative output directory.
	# outdir <- "./"
	
	# Generate predictions from the model fit...
	tmp_temps <- seq(min(
		floor(current_dataset$K)), 
		ceiling(max(current_dataset$K)
		), length = 200)
	
	tmp_model <- exp(Schoolfield(
		coef(schoolfield_nls)["B0"],
		coef(schoolfield_nls)["E"],
		coef(schoolfield_nls)["E_D"],
		coef(schoolfield_nls)["T_h"],
		tmp_temps
	))
	
	ModelToPlotS <- data.frame(
		Temperature = tmp_temps - 273.15, 
		TraitValue = tmp_model
	)
	
	# Prepare the data points of the original values.
	DataToPlot <- data.frame(
		Temperature = current_dataset$K - 273.15, 
		TraitValue = current_dataset$OriginalTraitValue
	)
	DataToPlot <- na.omit(DataToPlot)
	

# plot it! ----------------------------------------------------------------

	p_full <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																										y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
														 alpha = 0.7, pch = 21) + 
		geom_line(data = ModelToPlotS, 
							aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
							lwd = 1.3) +                           
		ggtitle("Phosphorus replete") +
		xlab(expression(paste("Temperature (", degree, C, ")"))) + 
		ylab("daphnia population abundance @ day 30") +
		theme_bw() + theme(plot.title = element_text(size = 16), 
											 axis.title = element_text(size = 16)) +
		ylim(0,250)
		annotate("text", size = 3, label=             
						 	paste("R^2","sch=", sprintf("%.2f", r_sq_sch),"\nE sch=", format(coef(schoolfield_nls)["E"], digits = 3),"\nAIC sch=",format(AIC(schoolfield_nls),digits=3)), 
						 x = min(DataToPlot[, "Temperature"]),
						 y = mean(DataToPlot[, "TraitValue"]),
						 hjust=0,
						 fontface = 3)
	
	ggsave("p-temp-figures_files/schoolfield_def_p.png")

	p_def <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																												 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
																	alpha = 0.7, pch = 21) + 
		geom_line(data = ModelToPlotS, 
							aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
							lwd = 1.3) +                           
		ggtitle("Phosphorus deficient") +
		xlab(expression(paste("Temperature (", degree, C, ")"))) + 
		ylab("daphnia population abundance @ day 30") +
		theme_bw() + theme(plot.title = element_text(size = 16), 
											 axis.title = element_text(size = 16)) +
		ylim(0,250)
	
	
	
library(gridExtra)	
grid.arrange(p_full, p_def, ncol = 2)
ggsave("p-temp-figures_files/schoolfield_plots.png")
	
# looking at model results ------------------------------------------------
library(broom)
coef(schoolfield_nls)
tidy(schoolfield_nls, conf.int = TRUE)

tidy()#Deficient phosphorus
# B0           E         E_D         T_h 
# 0.2420816   4.7971742   5.3311550 294.2349310 

# > coef(schoolfield_nls) # Full phospohorus
# B0          E        E_D        T_h 
# 1.552015   3.416349   6.019757 294.212115 
	