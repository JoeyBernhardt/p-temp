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
controldata <- read.csv(file = file.path("data-processed", "p_temp_algae.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

# Read the day zero data from the consumer free controls, and store as object
dayzerodata <- read.csv(file = file.path("data-processed", "march29_cell_data.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

# Read the data from the consumer free controls, and store as object
ptempdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data Manipulation ###

## controldata ##

# Rename columns in data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
controldata <- rename(controldata, phosphorus = P,
					     P = biovol,
					     temperature = temp
			)

# Isolate only the data involving controls, labelled in the original data with "CX"; X is some number
controldata <- filter(controldata, grepl("C", replicate))

# Create a column for boltzmann-transformed temperatures; this is useful for plotting purposes later
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
controldata <- mutate(controldata, transformedtemperature = -1/(Boltz * (temperature + 273.15)))

# Reorder rows in data frame by resource treatment, temperature, and ID
controldata <- arrange(controldata, phosphorus, temperature, ID)

# Create a column that lists the treatment type of each replicate, which depends on its phosphorus and temperature combination
controldata <- mutate(controldata, treatment = paste(controldata$phosphorus, controldata$temperature, sep=""))

# Split data frame into multiple indexed data frames based on their ID
controldata <- split(controldata, f = controldata$treatment)

## dayzerodata ##

# Remove unneeded columns from "dayzerodata" dataframe
dayzerodata <- dayzerodata[-c(1, 5, 7)] # column 1 is "filename"; column 5 is "date"; column 7 is "start time"

# Rename columns in "dayzerodata" to correspond to those in pdata
dayzerodata <- rename(dayzerodata, phosphorus_treatment = nutrient_level, 
				     volume_cell = cell_volume,
				     algal_cell_concentration_cells_per_ml = cell_density
			)

# Calculate the initial algal biovolume, and also add in some new columns corresponding to day zero values for days and daphnia.
dayzerodata <- mutate(dayzerodata, algal_biovolume = volume_cell * algal_cell_concentration_cells_per_ml,
				     days = 0,
				     daphnia_total = 10)

## ptempdata ##

# Vertically merge "pdata" and "dayzerodata" data frames
ptempdata <- bind_rows(ptempdata, dayzerodata)

# Rename columns in the "pdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
ptempdata <- rename(ptempdata, P = algal_biovolume)
ptempdata <- rename(ptempdata, H = daphnia_total)
ptempdata <- rename(ptempdata, Phosphorus = phosphorus_treatment)

# Here we scale the phytoplankton density by 250, in order to get the total biovolume present in the beaker.
# We want to model the total algal biovolume because we have modelled the total Daphnia density.
scalefactor <- 250
ptempdata <- mutate(ptempdata, P = P * scalefactor)

# Create a column for boltzmann-transformed temperatures; this is useful for plotting purposes later. "Boltz" is defined
# earlier in the code, for a similar transformation of controldata
ptempdata <- mutate(ptempdata, transformedtemperature = -1/(Boltz * (temperature + 273.15)))

# Reorder rows in data frame by treatment ID, and then by days
ptempdata <- arrange(ptempdata, Phosphorus, temperature, replicate, days)

# Create a column that lists the treatment type of each replicate, which depends on its phosphorus and temperature combination
ptempdata <- mutate(ptempdata, treatment = paste(ptempdata$Phosphorus, ptempdata$temperature, sep=""))

# Split entire dataset into multiple indexed data frames based on their treatment
ptempdata <- split(ptempdata, f = ptempdata$treatment)

### FITTING CONSUMER FREE CONTROLS ###

ControlParameters <- c(r = 2, K = 1000000)

ControlFittedParameters <- c("r", "K")

r_min <- 0.01
r_max <- 0.5

K_min <- 1e12
K_max <- 1e17

ControlLowerBound <- c(r = r_min, K = K_min)

ControlUpperBound <- c(r = r_max, K = K_max)

ControlParamScaling <- 1 / ControlUpperBound

CONTROLmodel <- new("odeModel",
	main = function (time, init, parms) {
			with(as.list(c(init, parms)), {
		dp <-  r * P * (1 - (P / K))
		list(c(dp))
		})
	},
	parms = ControlParameters,
	times = c(from = 0, to = 36, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 1e6),
	solver = "lsoda" #lsoda will be called with tolerances of 1e-9, as seen directly below. Default tolerances are both 1e-6. Lower is more accurate.
		)

controlfit <- function(data, num_replicates) {

repfit <- function(data) {

PORTfit <- function(num) {

		# First we initialize the simecol model object; we are using a CONTROLmodel here, which will be modified further
		model <- CONTROLmodel

		# Here we randomly generate our initial parameter settings by sampling from uniform distributions, with minima and maxima 
		# corresponding to the bounds declared earlier in the code. These are then assigned to our simecol model object.
		parms(model)[ControlFittedParameters] <- c(r = runif(1, min = r_min, max = r_max),
								K = runif(1, min = K_min, max = K_max)
								)

		# Here we isolate all of the day zero data, and then use that to set our initial model conditions. The initial phytoplankton density 
		# is set to the mean biovolume taken from the first measurement day.
		dayzerodata <- filter(data, days == 0)
		init(model) <- c(P = mean(dayzerodata$P))

		# Here we make sure that the times are correct
		times(model) <- c(from = 0, to = max(data$days), by = 0.1)
		
		# Here we declare our observed data in a form that fitOdeModel (used below) can understand.
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P) # The Y values of the observed data points we are fitting our model to

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences (SSQ) between the experimental data and our modelled data.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "LowerBound" and "UpperBound" are vectors containing the lower and upper bound constraints
		# for the parameter values. "ParamScaling" assists the function with fitting when the parameters take different orders of magnitude.

		fittedmodel <- fitOdeModel(model, whichpar = ControlFittedParameters, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = ControlLowerBound, upper = ControlUpperBound, scale.par = ControlParamScaling,
  		control = list(trace = TRUE),
			    rtol = 1e-9,
			    atol = 1e-9,
			    maxsteps = 5000
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' SSQ, phosphorus level, temperature, transformed temperature, and the values of 
		# the fitted parameters.

		ssq <- ssqOdeModel(coef(fittedmodel), model, obstime, yobs) # simulates another CONTROLmodel, but now using our fitted parameters, and outputs the SSQ.
		treatment <- data$treatment[1]
		phosphorus <- data$phosphorus[1]
		temperature <- data$temperature[1]
		transformedtemperature <- data$transformedtemperature[1]
		r <- coef(fittedmodel)["r"]
		K <- coef(fittedmodel)["K"]

		# Create one-row dataframe of the above vectors
		simulation_data <- data.frame(ssq, treatment, phosphorus, temperature, transformedtemperature, r, K)
		
		# Return above dataframe as the function's output
		return(simulation_data)
}
replicates <- seq(1, num_replicates, 1) # number of fitting attempts

outputdf <- map_df(replicates, PORTfit) # use map_df here; DO NOT USE A FOR-LOOP!!!

return(outputdf)
}

outputdf <- map_df(data, repfit)

return(outputdf)
}

rawfittedcontroldata <- controlfit(controldata, 10)

plotbest3controls <- function(phosphorus, temp) {

x <- paste(phosphorus, temp, sep = "")
obsdata <- controldata[[x]]

fitdata <- filter(rawfittedcontroldata, treatment == x)

fitdata <- filter(fitdata, ssq != 0)
fitdata <- arrange(fitdata, ssq)
fitdata <- fitdata[1:3,] # use best 3 for now

fitdata$repnumber <- rownames(fitdata)
fitdata <- split(fitdata, f = fitdata$repnumber)

# objective is to write a function that operates on a single data frame of parameters and other info, and then produce the density estimates
# using the sim function in simecol. Then use map_df on this guy.

innerfunction <- function(xdata) {

fittedr <- xdata$r
fittedK <- xdata$K

SimParameters <- c(r = fittedr, K = fittedK)

		dayzerodata <- filter(obsdata, days == 0)
		simmodel <- CONTROLmodel
		init(simmodel) <- c(P = mean(dayzerodata$P)) # Set initial model conditions to the biovolume taken from the first measurement day
		times(simmodel) <- c(from = 0, to = max(controldata[[x]]$days), by = 0.1)
		parms(simmodel) <- SimParameters

		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		simdata <- mutate(simdata, repnumber = xdata$repnumber)

return(simdata)
}

best3simdata <- map_df(fitdata, innerfunction)

output_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = P)) +
		geom_line(data = best3simdata, aes(x = time, y = P, color = repnumber))

return(output_plot)

}

plotbest3controls("DEF", 24)

### FITTING POPULATIONS WITH CONSUMERS ###

# Here, we construct the consumer resource model. The
# resulting system of equations, called "CRModel", can be used to produce
# theoretical predictions.

# Declare the parameter values which will be used in the dynamical models. This is just a dummy set 
# of parameters so that the CRmodel object below will properly initialize. The model parameters will 
# be set in a non-arbitary fashion later on in the code; details on this can be found further down.

# For reference:
# r = intrinsic growth rate of phytoplankton
# K = carrying capacity of phytoplankton
# a = attack rate of Daphnia
# eps = transfer efficiency for Daphnia
# m = intrinsic mortality rate for Daphnia

Parameters <- c(r = 2, K = 1000000, a = 0.5, eps = 1e-10, m = 0.1)

# This vector simply contains strings; they are used to tell the function
# "fitOdeModel" which parameters it is supposed to fit
FittedParameters <- c("r", "K", "a", "eps", "m")

## Declare the parameters to be used as the bounds for the fitting algorithm. ##

# These bounds will be used in two different ways:

# (1) They determine the minimum and maximum ranges for the parameter values which the fitting algorithm is allowed to fit.

# (2) They determine the entire range of the sample space from which our initial parameter values will be drawn from.

# N.B. Note that because the parameter scaling vector, "ParamScaling" (located below), depends on "UpperBound", changing the upper bound 
# of any of the parameter values can have unintended effects on the quality of fits. Thus, it recommended not to set the upper bounds on parameters 
# too high; don't increase them unless you can think of a biological/mathematical reason for doing so!

r_min <- 0.1
r_max <- 10

K_min <- 1e9
K_max <- 1e13

a_min <- 0
a_max <- 1

eps_min <- 1e-13
eps_max <- 1e-11

m_min <- 0.0001
m_max <- 0.1

LowerBound <- c(r = r_min, K = K_min, a = a_min, eps = eps_min, m = m_min)
UpperBound <- c(r = r_max, K = K_max, a = a_max, eps = eps_max, m = m_max)

# Declare the "step size" for the PORT algorithm. 1 / UpperBound is recommended
# by the simecol documentation.
ParamScaling <- 1 / UpperBound

# Create a new odeModel object. This represents the "base" Lotka-Volterra
# consumer resource dynamics model that we would like to eventually fit.

# dp and dh are the differential equations for producers and heterotrophs,
# respectively. Likewise, P and H refer to the population densities.

# If you would like to view the model output (density vs. time plots) for troubleshooting purposes, you can use the following
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

### Fitting Functions ###

## 1. Model Fitting Function ##

# The following function is intended to be used with map_df() on nested
# dataframes like "pdata". It takes a single dataframe of various observations
# for a replicate, and outputs a dataframe consisting of the replicate ID, the
# Phosphorus treatment, the temperature, and the parameter estimates for our
# differential equation. It can also be used to output parameter values for a
# single replicate. To do this call rKfit(controldata[['X']], where "X" is the
# replicate's ID number.

fitall <- function(data, num_replicates) {

repfit <- function(data) {

PORTfit <- function(num) {

		# First we initialize the simecol model object
		model <- CRmodel

		# Here we randomly generate our initial parameter settings by sampling from uniform distributions, with minima and maxima 
		# corresponding to the bounds declared earlier in the code. These are then assigned to our simecol model object.
		parms(model)[FittedParameters] <- c(r = runif(1, min = r_min, max = r_max),
								K = runif(1, min = K_min, max = K_max),
								a = runif(1, min = a_min, max = a_max),
								eps = runif(1, min = eps_min, max = eps_max),
								m = runif(1, min = m_min, max = m_max)
								)

		# Here we isolate all of the day zero data, and then use that to set our initial model conditions. The initial Daphnia density is always 10,
		# regardless of treatment group. The initial phytoplankton density is set to the mean biovolume taken from the first measurement day.
		dayzerodata <- filter(data, days == 0)
		init(model) <- c(P = mean(dayzerodata$P), H = 10)
		
		# Here we declare our observed data in a form that fitOdeModel (used below) can understand.
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P, H) # The Y values of the observed data points we are fitting our model to

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences (SSQ) between the experimental data and our modelled data.
		
		# The PORT algorithm is employed to perform the model fitting, analogous to O'Connor et al.
		# "LowerBound" and "UpperBound" are vectors containing the lower and upper bound constraints
		# for the parameter values. "ParamScaling" assists the function with fitting when the parameters take different orders of magnitude.

		fittedmodel <- fitOdeModel(model, whichpar = FittedParameters, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = TRUE),
			    rtol = 1e-9,
			    atol = 1e-9,
			    maxsteps = 5000
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' SSQ, phosphorus level, temperature, transformed temperature, and the values of 
		# the fitted parameters.

		ssq <- ssqOdeModel(coef(fittedmodel), model, obstime, yobs) # simulates another CRmodel, but now using our fitted parameters, and outputs the SSQ.
		treatment <- data$treatment[1]
		phosphorus <- data$Phosphorus[1]
		temperature <- data$temperature[1]
		transformedtemperature <- data$transformedtemperature[1]
		r <- coef(fittedmodel)["r"]
		K <- coef(fittedmodel)["K"]
		a <- coef(fittedmodel)["a"]
		eps <- coef(fittedmodel)["eps"]
		m <- coef(fittedmodel)["m"]
		
		# Create one-row dataframe of the above vectors
		simulation_data <- data.frame(ssq, treatment, phosphorus, temperature, transformedtemperature, r, K, a, eps, m)
		
		# Return above dataframe as the function's output
		return(simulation_data)
}
replicates <- seq(1, num_replicates, 1) # number of fitting attempts

outputdf <- map_df(replicates, PORTfit) # use map_df here; DO NOT USE A FOR-LOOP!!!

return(outputdf)
}

outputdf <- map_df(data, repfit)

return(outputdf)
}


### PRODUCING DATA ###

# rawfitteddata <- fitall(ptempdata, 2)

fitteddata <- rawfitteddata %>%
filter(ssq !=0) %>%
arrange(treatment, ssq) %>%
group_by(treatment) %>%
slice(1) %>%
ungroup

### OUTPUTTING DATA ###

# write.csv(fitteddata, "fitteddata06_2017_24_01.csv")
# write.csv(rawfitteddata, "rawfitteddata06_2017_24_01.csv")

### INPUTTING DATA ###

rawfitteddata <- read.csv(file = file.path("p-temp-marcus", "outputs", "rawfitteddata05_2017_22_01.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

rawfitteddata <- mutate(rawfitteddata, treatment = paste(rawfitteddata$phosphorus, rawfitteddata$temperature, sep=""))

fitteddata <- read.csv(file = file.path("p-temp-marcus", "outputs", "fitteddata05_2017_22_01.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

# use this function on 
plotbest3 <- function(phosphorus, temp) {

x <- paste(phosphorus, temp, sep = "")
obsdata <- ptempdata[[x]]

fitdata <- filter(rawfitteddata, treatment == x)

fitdata <- filter(fitdata, ssq != 0)
fitdata <- arrange(fitdata, ssq)
fitdata <- fitdata[1:3,] # use best 3 for now

fitdata$repnumber <- rownames(fitdata)
fitdata <- split(fitdata, f = fitdata$repnumber)

# objective is to write a function that operates on a single data frame of parameters and other info, and then produce the density estimates
# using the sim function in simecol. Then use map_df on this guy.

innerfunction <- function(xdata) {

fittedr <- xdata$r
fittedK <- xdata$K
fitteda <- xdata$a
fittedeps <- xdata$eps
fittedm <- xdata$m

SimParameters <- c(r = fittedr, K = fittedK, a = fitteda, eps = fittedeps, m = fittedm)

		dayzerodata <- filter(obsdata, days == 0)
		simmodel <- CRmodel
		init(simmodel) <- c(P = mean(dayzerodata$P), H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
		parms(simmodel) <- SimParameters

		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		simdata <- mutate(simdata, repnumber = xdata$repnumber)

return(simdata)
}

best3simdata <- map_df(fitdata, innerfunction)

prod_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = P)) +
		geom_line(data = best3simdata, aes(x = time, y = P, color = repnumber))


het_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = H)) +
		geom_line(data = best3simdata, aes(x = time, y = H, color = repnumber))

output_plot <- arrangeGrob(prod_plot, het_plot, ncol=2)

return(output_plot)

}

def12plot <- plotbest3("DEF", 12)
# ggsave("def12plot.png", plot = def12plot)

def16plot <- plotbest3("DEF", 16)
# ggsave("def16plot.png", plot = def16plot)

def20plot <- plotbest3("DEF", 20)
# ggsave("def20plot.png", plot = def20plot)

def24plot <- plotbest3("DEF", 24)
# ggsave("def24plot.png", plot = def24plot)

full12plot <- plotbest3("FULL", 12)
# ggsave("full12plot.png", plot = full12plot)

full16plot <- plotbest3("FULL", 16)
# ggsave("full16plot.png", plot = full16plot)

full20plot <- plotbest3("FULL", 20)
# ggsave("full20plot.png", plot = full20plot)

full24plot <- plotbest3("FULL", 24)
# ggsave("full24plot.png", plot = full24plot)

### PLOTS ###

fittedr_plot <- ggplot(data = rawfitteddata, aes(x = transformedtemperature, y = log(r), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(r) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log intrinsic growth rate (r)")
fittedr_plot
# ggsave("fittedr_2017_24_01.png", plot = last_plot())

r_full_model <- lm(data = filter(rawfitteddata, phosphorus == "FULL"), log(r) ~ transformedtemperature)
confint(r_full_model)

r_def_model <- lm(data = filter(rawfitteddata, phosphorus == "DEF"), log(r) ~ transformedtemperature)
confint(r_def_model)

fittedK_plot <- ggplot(data = rawfitteddata, aes(x = transformedtemperature, y = log(K), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(K) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log carrying capacity (K)")
fittedK_plot
# ggsave("fittedK_2017_24_01.png", plot = last_plot())

fitteda_plot <- ggplot(data = rawfitteddata, aes(x = transformedtemperature, y = log(a), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(a) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log attack rate (a)")
fitteda_plot
# ggsave("fitteda_2017_24_01.png", plot = last_plot())