### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)

### Data frame ###

# Read the data from the consumer free controls, and store as object
controldata <- read.csv("C:\\Users\\Matt\\My Documents\\repo\\p-temp\\data-processed\\p_temp_algae.csv",
	stringsAsFactors = FALSE,
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data manipulation ###

# Rename column "P" in data frame to "Phosphorus". This is necessary to avoid
# confusion in the next step.
controldata <- rename(controldata, Phosphorus = P)

# Rename columns in data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
controldata <- rename(controldata, P = biovol)

# Isolate only the data involving controls, labelled in the original data with "CX"; X is some number
controldata <- filter(controldata, grepl("C", replicate))

# Reorder rows in data frame by resource treatment, temperature, and ID
controldata <- arrange(controldata, Phosphorus, temp, ID)

# Split data frame into multiple indexed data frames based on their ID
controldata <- split(controldata, f = controldata$ID)

### Model Construction ###

## Declare Arrhenius Function ##

# Create an Arrhenius function to transform metabolic rates based on
# temperature. Here, "T" is the temperature. "E" is the activation energy
# constant. Its structure is identical to the one used in O'Connor et al.
# 2011

# Declare constants for Arrhenius function
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
BasalTemperature <- 12 + 273.15 # the "base temperature" that determines the basal metabolic rate; 12 was the lowest temp used during the experiment.

# Declare Arrhenius function
arrhenius <- function(T,E){
	output <- exp(E * ((T + 273.15) - BasalTemperature) / (Boltz * (T + 273.15) * BasalTemperature))
return(output)
}

## Declare parameters ##

# Declare the parameters to be used in the dynamical models
LowResourceParameters <- c(r = 1, K = 5, Er = 0.32, EK = -0.32)

# MAKE TEMP A PARAMETER!!!

CRmodel <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
			dp <- arrhenius(temp, Er) * r * P * (1 - (P / (arrhenius(temp, EK) * 10^K))) # K is scaled exponentially to assist the PORT algorithm
		list(c(dp))
		})
	},
	parms = LowResourceParameters,
	times = c(from = 0, to = 35, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 100000),
	solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
		)

fittedparms <- c("r", "K")
temp <- 10
controlfit <- function(data){
		
		CRmodel <- CRmodel
		init(CRmodel) <- c(P = data$P[1])
		obstime <- data$days
		yobs <- select(data, P)
		temp <- data$temp[1]

		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = c(r = 0.1, K = 3),
  		control = list(trace = T)
		)
		
		Phosphorus <- data$Phosphorus[1]
		r <- coef(fittedCRmodel)[1]
		K <- coef(fittedCRmodel)[2]
		ID <- data$ID[1]
		output <- data.frame(ID, Phosphorus, temp, r, K)
		return(output)
}

# If you would like to fit parameters for all of the control replicates, and
# output all of the results together as a dataframe, use:
# map_df(controldata, controlfit)

# If you would like to plot the model fit for a single replicate X, please use:
# plotsinglefit(controldata[['X']]

plotsinglefit <- function(data){
		
		CRmodel <- CRmodel
		init(CRmodel) <- c(P = data$P[1])
		obstime <- data$days
		yobs <- select(data, P)
		temp <- data$temp[1]

		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = c(r = 0.1, K = 3),
  		control = list(trace = T)
		)

	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)

	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=40, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))

	observeddata <- data.frame(obstime, yobs)
	simulateddata <- ysim

	# Plot the results of our model fitting.
	biol_plot <- ggplot() +
		geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) +
		geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) +
		labs(x = "Time (days)", y = "Algal Biovolume")

	output <- biol_plot
	return(output)
}
