### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)
# Load the gridExtra package for conveniently arranging ggplot objects
library(gridExtra)
# Load the knitr package for producing tables and graphics with markdown
library(knitr)

### Data Frames ###

# Read the data from the consumer free controls, and store as object
pdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

# Read the day zero data from the consumer free controls, and store as object
initdata <- read.csv(file = file.path("data-processed", "march29_cell_data.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data Manipulation ###

# Remove unneeded columns from "initdata" dataframe
initdata <- initdata[-c(1, 5, 7)] # column 1 is "filename"; column 5 is "date"; column 7 is "start time"

# Rename columns in "initdata" to correspond to those in pdata
initdata <- rename(initdata, phosphorus_treatment = nutrient_level, 
				     volume_cell = cell_volume,
				     algal_cell_concentration_cells_per_ml = cell_density
			)
# Calculate the initial algal biovolume, and also add in some new columns corresponding to initial values for days and daphnia.
initdata <- mutate(initdata, algal_biovolume = volume_cell * algal_cell_concentration_cells_per_ml,
				     days = 0,
				     daphnia_total = 10)

# Vertically merge "pdata" and "initdata" data frames
pdata <- bind_rows(pdata, initdata)

# Rename columns in the "pdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
pdata <- rename(pdata, P = algal_biovolume)
pdata <- rename(pdata, H = daphnia_total)
pdata <- rename(pdata, Phosphorus = phosphorus_treatment)

scalefactor <- 250
pdata <- mutate(pdata, P = P * scalefactor)

# Reorder rows in data frame by treatment ID, and then by days
pdata <- arrange(pdata, Phosphorus, temperature, replicate, days)

# Create a column for boltzmann-transformed temperatures
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
pdata <- mutate(pdata, transformedtemp = -1/(Boltz * (temperature + 273.15)))

# Split entire dataset into multiple indexed data frames based on their ID

full12data <- filter(pdata, Phosphorus == "FULL" & temperature == 12)

full16data <- filter(pdata, Phosphorus == "FULL" & temperature == 16)

full20data <- filter(pdata, Phosphorus == "FULL" & temperature == 20)

full24data <- filter(pdata, Phosphorus == "FULL" & temperature == 24)

def12data <- filter(pdata, Phosphorus == "DEF" & temperature == 12)

def16data <- filter(pdata, Phosphorus == "DEF" & temperature == 16)

def20data <- filter(pdata, Phosphorus == "DEF" & temperature == 20)

def24data <- filter(pdata, Phosphorus == "DEF" & temperature == 24)

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

# The Boltzmann constant "Boltz" has already been declared earlier in the code
BasalTemperature <- 12 + 273.15 # the "base temperature" that determines the basal metabolic rate; 12 was the lowest temp used during the experiment.

# Declare Arrhenius function
arrhenius <- function(T,E){
	output <- exp(E * ((T + 273.15) - BasalTemperature) / (Boltz * (T + 273.15) * BasalTemperature))
return(output)
}

## Consumer Resource Model ##

# Declare the parameters to be used in the dynamical models

Parameters <- c(r = 2.67, K = 92647516291, a = 0.27, eps = 5.834479e-11, m = 0.11)

# WORKED FOR DEF16!!!
# Parameters <- c(r = 2.5, K = 1e13, a = 0.15, eps = 0.00000000008, m = 0.1)


# This vector simply contains strings; they are used to tell the function
# "fitOdeModel" which parameters it is supposed to fit
FittedParameters <- c("r", "K", "a", "eps", "m")

# Declare the parameters to be used as the bounds for the fitting algorithm
LowerBound <- c(r = 0.5, K = 1e9, a = 0.1, eps = 0.000000000005, m = 0.001)
UpperBound <- c(r = 6, K = 1e12, a = 10, eps = 0.0000001, m = 1)

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
# command afer CRmodel has been evaluated:
# plot(sim(CRmodel))

CRmodel <- new("odeModel",
	main = function (time, init, parms) {
			with(as.list(c(init, parms)), {
		dp <-  r * P * (1 - (P / K)) - a * H * P
		dh <-  a * eps * H * P - m * H
		list(c(dp, dh))
		})
	},
	parms = Parameters,
	times = c(from = 0, to = 36, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 1e6, H = 10),
	solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
		)

### FUNCTIONS ###

# This section contains two functions, one for fitting, and one for plotting
# the fit of a single replicate.

## 1. Model Fitting Function ##

# The following function is intended to be used with map_df() on nested
# dataframes like "pdata". It takes a single dataframe of various observations
# for a replicate, and outputs a dataframe consisting of the replicate ID, the
# Phosphorus treatment, the temperature, and the parameter estimates for our
# differential equation. It can also be used to output parameter values for a
# single replicate. To do this call rKfit(controldata[['X']], where "X" is the
# replicate's ID number.


pfit <- function(data) {
	
		dayzerodata <- filter(data, days == 0)
		model <- CRmodel
		init(model) <- c(P = mean(dayzerodata$P[1:6]), H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P, H) # The Y values of the observed data points we are fitting our model to

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "lower" is a vector containing the lower bound constraints
		# for the parameter values. This may need tweaking.

		fittedmodel <- fitOdeModel(model, whichpar = FittedParameters, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = TRUE),
			    rtol = 1e-10,
			    atol = 1e-10,
			    maxsteps = 5000
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' ID, Phosphorus level, temperature, and the
		# fitted parameters.
	
		coef(fittedmodel)

		simmodel <- model
		# set model parameters to fitted values and simulate
		parms(simmodel)[FittedParameters] <- coef(fittedmodel)
		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		simdata <- mutate(simdata, r = coef(fittedmodel)["r"],
						   K = coef(fittedmodel)["K"],
						   a = coef(fittedmodel)["a"],
						   eps = coef(fittedmodel)["eps"],
						   m = coef(fittedmodel)["m"],
						   temp = data$temperature[1],
						   transformedtemp = data$transformedtemp[1]
					)
		return(simdata)
}

targetdata <- full20data
fitteddata <- pfit(targetdata)

### Calibration function ###

# Declare the parameters to be used in the dynamical models
SimParameters <- c(r = 4.05, K = 113051101738, a = 0.24029643, eps = 5.812694e-06, m = 0.009410279)


simfit <- function(data){

		dayzerodata <- filter(data, days == 0)
		simmodel <- CRmodel
		init(simmodel) <- c(P = mean(dayzerodata$P[1:6]), H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
		parms(simmodel) <- SimParameters

		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		
		return(simdata)
}

simulateddata <- simfit(targetdata)

### PLOTS ###

prod_plot <- ggplot() +
		geom_point(data = targetdata, aes(x = days, y = P)) +
		geom_line(data = fitteddata, aes(x = time, y = P), color = "red")


het_plot <- ggplot() +
		geom_point(data = targetdata, aes(x = days, y = H)) +
		geom_line(data = fitteddata, aes(x = time, y = H), color = "red")

grid.arrange(prod_plot, het_plot, ncol=2)

# initial parameter seedings
# for def16: 2.676643 92647516291 0.1523916 5.834479e-11 0.1133629

# even better for def16!:
#  4 178420773594 0.2292639 1.174038e-11 0.01743577

# calculate ssq

obstime <- targetdata$days # The X values of the observed data points we are fitting our model to
yobs <- select(targetdata, P, H) # The Y values of the observed data points we are fitting our model to

findinit <- function(num){

rand_r <- runif(1, min = 0, max = 5)
rand_K <- runif(1, min = 1e9, max = 1e12)
rand_a <- runif(1, min = 0, max = 10)
rand_eps <- runif(1, min = 0, max = 1e-5)
rand_m <- runif(1, min = 0, max = 0.5)

SearchParameters <- c(r = rand_r, K = rand_K, a = rand_a, eps = rand_eps, m = rand_m)
ssq <- ssqOdeModel(SearchParameters, CRmodel, obstime, yobs)

output <- data.frame(ssq, rand_r, rand_K, rand_a, rand_eps, rand_m)
return(output)
}

seq <- seq(1,1000,1)
findinitdata <- map_df(seq, findinit)
findinitdata
