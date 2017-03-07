### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)

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

# Remove strange replicates
controldata <- filter(controldata, ID != 49 & ID != 51 & ID != 54 & ID != 57)

# Isolate only the data involving the 12C controls, dividing into different resource treatments.
TwelveDefPdata <- filter(controldata, temp == 12, Phosphorus == "DEF")
TwelveFullPdata <- filter(controldata, temp == 12, Phosphorus == "FULL")

# Split the 12C control dataframes into indexed data frames based on their ID
TwelveDefPdata <- split(TwelveDefPdata, f = TwelveDefPdata$ID)
TwelveFullPdata <- split(TwelveFullPdata, f = TwelveFullPdata$ID)

# Isolate and divide the non-12C data into different dataframes for each resource treatment
DefPdata <- filter(controldata, temp != 12, Phosphorus == "DEF")
FullPdata <- filter(controldata, temp != 12, Phosphorus == "FULL")

# Split non-12C dataframes into multiple indexed data frames based on their ID
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
Parameters <- c(r = 1, K = 5, Er = 0.32, EK = -0.32, temp = 12)

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.15, K = 10 ^ 6)
UpperBound <- c(r = 2, K = 10 ^ 11) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

CRmodel <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
			dp <- arrhenius(temp, Er) * r * P * (1 - (P / (arrhenius(temp, EK) * K)))
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
TempName <- c("temp") # for assigning temperature values to CRmodel
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
		parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.

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
		temp <- tail(parms(CRmodel), n=1)
		output <- data.frame(ID, Phosphorus, temp, r, K)
		return(output)
}

# Here we fit the values for r and K, using the 12C replicates, and then
# prepare these data to be used to fit the activation energies further below in
# the code.

# Fit r and K for 12C replicates
rKdefdata <- map_df(TwelveDefPdata, rKfit)
rKfulldata <- map_df(TwelveFullPdata, rKfit)

# Extract fitted values for r and K, and take their arithmetic means
defr <- mean(rKdefdata$r)
defK <- mean(rKdefdata$K)

fullr <- mean(rKfulldata$r)
fullK <- mean(rKfulldata$K)

## Declare function for fitting Er and EK ##

# Declare parameters to be used for fitting Er and EK. r and K taken from above steps. We investigate the Phosphorus "FULL" and
# "DEF" cases separately.
defParameters <- c(r = defr, K = defK, Er = 0.2, EK = 0.2, temp = 12)
fullParameters <- c(r = fullr, K = fullK, Er = 0.1, EK = 0.1, temp = 12)

# Declare the parameters to be used as the bounds for the fitting algorithm
ErEKLowerBound <- c(Er = -1, EK = -1)
ErEKUpperBound <- c(Er = 1, EK = 1)

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ErEKParamScaling <- 1 / ErEKUpperBound

# Declare fitted parameters for ErEKfit
ErEKfittedparms <- c("Er", "EK")

ErEKfit <- function(data){
				
		init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
		parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.

		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = ErEKfittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = ErEKLowerBound, upper = ErEKUpperBound, scale.par = ErEKParamScaling,
  		control = list(trace = T)
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' ID, Phosphorus level, temperature, and the
		# fitted parameters. "truer" and "trueK" are the fitted
		# activation energies, but scaled using the appropriate arrhenius
		# transform.
		Phosphorus <- data$Phosphorus[1]
		Er <- coef(fittedCRmodel)[1]
		EK <- coef(fittedCRmodel)[2]
		ID <- data$ID[1]
		temp <- tail(parms(CRmodel), n=1)
		output <- data.frame(ID, Phosphorus, temp, Er, EK)
		return(output)
}

# Generating dataframes with fitted Er and EK #

# Phosphorus deficient case #

parms(CRmodel) <- defParameters

fittedErEKdefdata <- map_df(DefPdata, ErEKfit)

# Phosphorus rich case #
parms(CRmodel) <- fullParameters

fittedErEKfulldata <- map_df(FullPdata, ErEKfit)

## Plotting Functions ##

# Here, intErEKfit is just an intermediate "helper function" - please do not
# call it. It is used by plotErEKfit for visualizing the model fit of a single
# replicate. To plot something with ID "X", please call:
# plotErEKfit(controldata[['X']])

intErEKfit <- function(data){

		init(CRmodel) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to
		parms(CRmodel)[TempName] <- data$temp[1] # Set the temperature parameter in CRmodel to whatever our control replicate used.
		
		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.
		
		fittedCRmodel <- fitOdeModel(CRmodel, whichpar = ErEKfittedparms, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = ErEKLowerBound, upper = ErEKUpperBound, scale.par = ErEKParamScaling,
  		control = list(trace = T)
		)

	# To display the fitted results we need to create a new OdeModel object. Here
	# we duplicate CRmodel and then alter it to use our new fitted parameters.
	plotfittedCRmodel <- CRmodel
	parms(plotfittedCRmodel)[ErEKfittedparms] <- coef(fittedCRmodel)

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

plotErEKfit <- function(data){
	if(data$Phosphorus[1] == "DEF"){
		parms(CRmodel) <- defParameters
		output <- intErEKfit(data)
		return(output)
		} else {
		parms(CRmodel) <- fullParameters
		output <- intErEKfit(data)
		return(output)
		}
}

###!! HOW TO USE THIS SCRIPT: !!###

# 1. Simply run all of the code, and it will produce four dataframes of interest. 

# Dataframes of the fitted r's and K's, grouped by phosphorus treatment:

# rKdefdata
# rKfulldata

# dataframes of the fitted Er's and EK's, grouped by phosphorus treatment:

# fittedErEKdefdata
# fittedErEKfulldata

# 2. To visually investigate the fit of a single replicate with ID = "X", please use:
plotErEKfit(controldata[["74"]])
str(controldata)
# Please keep in mind that step (2.) above is only designed to fit THE ACTIVATION
# ENERGIES of single replicates, and then plot the fitting result. It won't
# make any sense to call it on data from any of the 12C replicates.


### Joey playing around
mean(fittedErEKfulldata$Er)
mean(fittedErEKdefdata$Er)
