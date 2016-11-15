### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)
# Load the gridExtra package for conveniently arranging ggplot objects
library(gridExtra)

### Data Frame ###

# Read the data from the consumer free controls, and store as object
pdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data Manipulation ###

# Rename columns in data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
pdata <- rename(pdata, P = algal_biovolume)
pdata <- rename(pdata, H = daphnia_total)

# Reorder rows in data frame by ID, and then days
pdata <- arrange(pdata, unique_ID, days)

# Split entire control dataset into multiple indexed data frames based on their ID
pdata <- split(pdata, f = pdata$unique_ID)

### Model Construction ###

# Here we construct the consumer resource model. First the Arrhenius function
# is declared, and it is then used in the creation of the dynamical model. The
# resulting system of equations, called "CRModel", can be used to produce
# theoretical predictions.

# Create an Arrhenius function to transform metabolic rates based on
# temperature. Here, "T" is the temperature. "E" is the activation energy
# constant. Its structure is identical to the one used in O'Connor et al.
# 2011

## Arrhenius Function ##

# First we declare constants used in the Arrhenius function

Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
BasalTemperature <- 12 + 273.15 # the "base temperature" that determines the basal metabolic rate; 12 was the lowest temp used during the experiment.

# Declare Arrhenius function
arrhenius <- function(T,E){
	output <- exp(E * ((T + 273.15) - BasalTemperature) / (Boltz * (T + 273.15) * BasalTemperature))
return(output)
}

## Consumer Resource Model ##

# Declare the parameters to be used in the dynamical models #
Parameters <- c(r = 1, K = 10 ^ 8, a = 3, b = 10 ^ 5, eps = 0.1, m = 0.2, Er = 0.32, EK = -0.32, Ea = 0.65, Em = 0.65, temp = 12)

# This vector simply contains strings; they are used to tell the function
# "fitOdeModel" which parameters it is supposed to fit
FittedParameters <- c("r", "K", "a", "b", "eps", "m")

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.15, K = 10 ^ 6, a = 0, b = 0, eps = 0, m = 0)
UpperBound <- c(r = 2, K = 10 ^ 11, a = 1000, b = 10 ^ 7, eps = 10, m = 10) 

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

# Create a new odeModel object. This represents the "base" Lotka-Volterra
# consumer resource dynamics model that we would like to eventually fit, with
# added metabolic effects due to temperature.

# dp and dh are the differential equations for producers and heterotrophs,
# respectively. Likewise, P and H refer to the population densities. For
# producers, both the intrinsic growth rate r and the carrying capacity K are
# subject to metabolic scaling. For heterotrophs, metabolic scaling applies to
# the attack rate a, and the intrinsic mortality rate m. "temp" refers to the
# temperature of the modelled system, used for determining the arrhenius
# scaling factor.

# If you would like to view the model output, you can use the following
# commands afer hp has been evaluated:
# hp <- sim(hp)
# plot(hp)

CRmodel <- new("odeModel",
	main = function (time, init, parms) {
			with(as.list(c(init, parms)), {
		dp <- arrhenius(temp, Er) * r * P * (1 - (P / (arrhenius(temp, EK) * K))) - arrhenius(temp, Ea) * a * H * (P / (P + b))
		dh <- arrhenius(temp, Ea) * a * eps * H * (P / (P + b)) - arrhenius(temp, Em) * m * H
		list(c(dp, dh))
		})
	},
	parms = Parameters,
	times = c(from = 0, to = 35, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 100000),
	solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
		)

create_model <- function(temp){
		model <- CRmodel
		parms(model)["temp"] <- temp
		output <- model 
return(output)
}

### FUNCTIONS ###

# This section contains two functions, one for fitting, and one for plotting
# the fit of a single replicate.

## 1. Model Fitting Function ##

# The following function is intended to be used with map_df() on nested
# dataframes like "rKdefdata". It takes a single dataframe of various
# observations for a control replicate, and outputs a dataframe consisting of
# the replicate ID, the Phosphorus treatment, the temperature, and the
# parameter estimates for r and K. It can also be used to output parameter
# values for a single replicate. To do this call rKfit(controldata[['X']],
# where "X" is the replicate's ID number.

pfit <- function(data){

		temp <- data$temperature[1]
		model <- create_model(temp)
		init(model) <- c(P = data$P[1]) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P, H) # The Y values of the observed data points we are fitting our model to

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.

		fittedCRmodel <- fitOdeModel(model, whichpar = FittedParameters, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = T)
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' ID, Phosphorus level, temperature, and the
		# fitted parameters. "transformedtemp" is simply -1/kT, where
		# k is the Boltzmann constant and T is temperature. This measure
		# will be used for plotting purposes.

		ID <- data$unique_ID[1]
		Phosphorus <- data$phosphorus_treatment[1]
		r <- coef(fittedCRmodel)["r"]
		K <- coef(fittedCRmodel)["K"]
		a <- coef(fittedCRmodel)["a"]
		b <- coef(fittedCRmodel)["eps"]
		m <- coef(fittedCRmodel)["m"]
		
		output <- data.frame(ID, Phosphorus, temp, r, K, a, b, m)
		return(output)
}

## 2. Plotting Function ##

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

### Output Data ###

# Dataframes of the fitted r's and K's, grouped by phosphorus treatment:
rKdata <- map_df(controldata, rKfit)

### Output activation energies ###

# We regress log(fitted parameter) onto -1/kT; here the activation energies are the slopes of each OLS model.

r_model <- lm(log(r) ~ transformedtemp, data = rKdata)
K_model <- lm(log(K) ~ transformedtemp, data = rKdata)

# Create data frame containing the fitted activation energies
r <- last(coef(r_model))
K <- last(coef(K_model))

Eadata <- data.frame(r, K)

### Plot results of entire dataset ###

# Plot for fitted r values for all temperature treatments
r_plot <- ggplot(data = rKdata, aes(x = transformedtemp, y = log(r), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(r) Values") +
		labs(x = "-1/kT", y = "log(r)")

# Plot for fitted K values for all temperature treatments
K_plot <- ggplot(data = rKdata, aes(x = transformedtemp, y = log(K), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(K) Values") +
		labs(x = "-1/kT", y = "log(K)")

# Display both plots together
grid.arrange(r_plot, K_plot, nrow = 2)
