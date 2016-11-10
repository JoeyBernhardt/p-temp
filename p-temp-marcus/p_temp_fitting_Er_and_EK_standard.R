### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)
# Load the gridExtra package for conveniently arranging ggplot objects
library(gridExtra)

### Data frame ###

# Read the data from the consumer free controls, and store as object
controldata <- read.csv(file = file.path("data-processed", "p_temp_algae.csv"),
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

# Remove "weird" replicates

controldata <- filter(controldata, ID != 49 & ID != 51 & ID != 54 & ID != 57)

# Divide the data by their different resource treatments.
DefPdata <- filter(controldata, Phosphorus == "DEF")
FullPdata <- filter(controldata, Phosphorus == "FULL")

# Split the above two dataframes into two indexed lists of data frames based on their ID
DefPdata <- split(DefPdata, f = DefPdata$ID)
FullPdata <- split(FullPdata, f = FullPdata$ID)

# Split entire control dataset into multiple indexed data frames based on their ID
controldata <- split(controldata, f = controldata$ID)

### Model Construction ###

## Arrhenius Function ##

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

## Consumer Resource Model ##

# Create a new odeModel object. This represents the "base" Lotka-Volterra
# consumer resource dynamics model that we would like to eventually fit, with
# added metabolic effects due to temperature.

# dp is the differential equation for the phytoplankton population dynamics. P
# refers to the population density. Both the intrinsic growth rate r and the
# carrying capacity K are subject to metabolic scaling by the Arrhenius
# function defined above

# Declare the parameters to be used in the dynamical models

# Activation energies used in O'Connor et al.
Parameters <- c(r = 1, K = 5, Er = 0.32, EK = -0.32)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.15, K = 10 ^ 6)
UpperBound <- c(r = 2, K = 10 ^ 11) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

CRmodel <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
			dp <-  r * P * (1 - (P /  K))
		list(c(dp))
		})
	},
	parms = Parameters,
	times = c(from = 0, to = 35, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 100000),
	solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
		)

# These vectors simply contain strings, which are used to facilitate parameter
# assignment in the below code. While seemingly clumsy, it appears to be a
# necessary step.
fittedparms <- c("r", "K") # for assigning fitted parameter values to fittedCRmodel

## Model Fitting Function ##

# The following function is intended to be used with map_df() on the nested
# dataframe called "controldata". It takes a single dataframe of various
# observations for a control replicate, and outputs a dataframe consisting of
# the replicate ID, the Phosphorus treatment, the temperature, and the
# parameter estimates for r and K. It can also be used to output parameter
# values for a single replicate. To do this call controlfit(controldata[['X']],
# where "X" is the replicate's ID number.

rKfit <- function(data){
		
		init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
		parms(CRmodel) <- Parameters
		temp <- data$temp[1]

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.

		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = T)
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' ID, Phosphorus level, temperature, and the
		# fitted parameters. "truer" and "trueK" are the fitted
		# parameters, but scaled using the appropriate arrhenius
		# transform.	
		Phosphorus <- data$Phosphorus[1]
		r <- coef(fittedCRmodel)[1]
		K <- coef(fittedCRmodel)[2]
		ID <- data$ID[1]
		transformedtemp <- -1/(Boltz * (temp + 273.15))
		output <- data.frame(ID, Phosphorus, temp, transformedtemp, r, K)
		return(output)
}

# Here we fit the values for r and K, using the 12C replicates, and then
# prepare these data to be used to fit the activation energies further below in
# the code.

# Fit r and K
rKdefdata <- map_df(DefPdata, rKfit)
rKfulldata <- map_df(FullPdata, rKfit)
rKdata <- rbind(rKdefdata, rKfulldata)

## Plotting Functions ##

plotsinglefit <- function(data){

		init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
		
		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.
		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = fittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = T)
		)

	# To display the fitted results we need to create a new OdeModel object. Here
	# we duplicate CRmodel and then alter it to use our new fitted parameters.
	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[fittedparms] <- coef(fittedCRmodel)

	# set model parameters to fitted values and simulate again
	times(plotfittedCRmodel) <- c(from=0, to=40, by=1)
	ysim <- out(sim(plotfittedCRmodel, rtol = 1e-9, atol = 1e-9))

	# Form observed data into a dataframe; the simulated data are already in a dataframe
	observeddata <- data.frame(obstime, yobs)
	simulateddata <- ysim

	# Plot the results of our model fitting.
	biol_plot <- ggplot() +
		geom_point(data = observeddata, aes(x = obstime, y = yobs, color = "observed")) + # Observed data are points
		geom_line(data = simulateddata, aes(x = time, y = P, color = "simulated")) + # Simulated data are in a continuous line
		labs(x = "Time (days)", y = "Algal Biovolume")
	
	# Output the results as a ggplot2 object
	output <- biol_plot
	return(output)
}


r_plot <- ggplot(data = rKdata, aes(x = transformedtemp, y = log(r), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(r) Values") +
		labs(x = "-1/kT", y = "log(r)")

K_plot <- ggplot(data = rKdata, aes(x = transformedtemp, y = log(K), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(K) Values") +
		labs(x = "-1/kT", y = "log(K)")

grid.arrange(r_plot, K_plot, nrow = 2)

r_model <- lm(log(r) ~ transformedtemp, data = rKdata)
summary(r_model)

K_model <- lm(log(K) ~ transformedtemp, data = rKdata)
summary(K_model)

## Strange replicates

# 49: either bizarre measurement error or severe demographic noise
# 51: appears to be measurement error at t=20 days. Impossible to fit logistic growth to resulting pattern.
# 54: dynamics appear to undergo strange bifurcation. Chaos?
# 57: dynamics appear periodic with high amplitude; possibly chaotic
# 63,64,65,66,67,68,70,71,73, and 74: None of these appear to reach equilibrium density.

###!! HOW TO USE THIS SCRIPT: !!###

# 1. Simply run all of the code, and it will produce four dataframes of interest. 

# Dataframes of the fitted r's and K's, grouped by phosphorus treatment:

# rKdefdata
# rKfulldata

# 2. To visually investigate the fit of a single replicate with ID = "X", please use:
# plotErEKfit(controldata[['X']])
