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
dayzerodata <- read.csv(file = file.path("data-processed", "march29_cell_data.csv"), #file.path() is used for cross-platform compatibility
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data Manipulation ###

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

# Vertically merge "pdata" and "dayzerodata" data frames
pdata <- bind_rows(pdata, dayzerodata)

# Rename columns in the "pdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
pdata <- rename(pdata, P = algal_biovolume)
pdata <- rename(pdata, H = daphnia_total)
pdata <- rename(pdata, Phosphorus = phosphorus_treatment)

# Here we scale the phytoplankton density by 250, in order to get the total biovolume present in the beaker.
scalefactor <- 250
pdata <- mutate(pdata, P = P * scalefactor)

# Reorder rows in data frame by treatment ID, and then by days
pdata <- arrange(pdata, Phosphorus, temperature, replicate, days)

# Create a column for boltzmann-transformed temperatures
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
pdata <- mutate(pdata, transformedtemperature = -1/(Boltz * (temperature + 273.15)))

# Split entire dataset into multiple data frames based on their phosphorus and temperature combinations

full12data <- filter(pdata, Phosphorus == "FULL" & temperature == 12)
full16data <- filter(pdata, Phosphorus == "FULL" & temperature == 16)
full20data <- filter(pdata, Phosphorus == "FULL" & temperature == 20)
full24data <- filter(pdata, Phosphorus == "FULL" & temperature == 24)

def12data <- filter(pdata, Phosphorus == "DEF" & temperature == 12)
def16data <- filter(pdata, Phosphorus == "DEF" & temperature == 16)
def20data <- filter(pdata, Phosphorus == "DEF" & temperature == 20)
def24data <- filter(pdata, Phosphorus == "DEF" & temperature == 24)

### Model Construction ###

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

K_min <- 1e11
K_max <- 1e14

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

portfit <- function(data) {

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
		# "lower" and upper are vectors containing the lower and upper bound constraints
		# for the parameter values.

		fittedmodel <- fitOdeModel(model, whichpar = FittedParameters, obstime, yobs,
 		debuglevel = 0, fn = ssqOdeModel,
   		method = "PORT", lower = LowerBound, upper = UpperBound, scale.par = ParamScaling,
  		control = list(trace = TRUE),
			    rtol = 1e-10,
			    atol = 1e-10,
			    maxsteps = 5000
		)
		
		# Here we create vectors to be used to output a dataframe of
		# the replicates' SSQ, phosphorus level, temperature, transformed temperature, and the values of 
		# the fitted parameters.

		ssq <- ssqOdeModel(coef(fittedmodel), model, obstime, yobs) # simulates another CRmodel, but now using our fitted parameters, and outputs the SSQ.
		phosphorus <- data$Phosphorus[1]
		temperature <- data$temperature[1]
		transformedtemperature <- data$transformedtemperature[1]
		r <- coef(fittedmodel)["r"]
		K <- coef(fittedmodel)["K"]
		a <- coef(fittedmodel)["a"]
		eps <- coef(fittedmodel)["eps"]
		m <- coef(fittedmodel)["m"]
		
		# Create one-row dataframe of the above vectors
		simulation_data <- data.frame(ssq, phosphorus, temperature, transformedtemperature, r, K, a, eps, m)
		
		# Return above dataframe as the function's output
		return(simulation_data)
}

repfit <- function(data, num){

replicates <- seq(1, num, 1)
simulation_data <- data.frame(ssq = double(),
					phosphorus = integer(),
					temperature = double(),
					transformedtemperature = double(),
					r = double(),
			  		K = double(), 
                	   		a = double(), 
			   		eps = double(),
					m = double(),
                	   		stringsAsFactors = FALSE)

for (i in replicates) {
simulation_output <- portfit(data)
simulation_data <- rbind(simulation_data, simulation_output)
}

return(simulation_data)
}

trimdata <- function(data){

data <- filter(data, ssq != 0)
data <- arrange(data, ssq)
data <- data[1,]

return(data)
}

num <- 10
rawfitteddef12data <- repfit(def12data, num)
rawfitteddef16data <- repfit(def16data, num)
rawfitteddef20data <- repfit(def20data, num)
rawfitteddef24data <- repfit(def24data, num)

rawfittedfull12data <- repfit(full12data, num)
rawfittedfull16data <- repfit(full16data, num)
rawfittedfull20data <- repfit(full20data, num)
rawfittedfull24data <- repfit(full24data, num)

fitteddef12data <- trimdata(rawfitteddef12data)
fitteddef16data <- trimdata(rawfitteddef16data)
fitteddef20data <- trimdata(rawfitteddef20data)
fitteddef24data <- trimdata(rawfitteddef24data)

fittedfull12data <- trimdata(rawfittedfull12data)
fittedfull16data <- trimdata(rawfittedfull16data)
fittedfull20data <- trimdata(rawfittedfull20data)
fittedfull24data <- trimdata(rawfittedfull24data)

rawfitteddata <- bind_rows(rawfitteddef12data,
				   rawfitteddef16data,
				   rawfitteddef20data,
				   rawfitteddef24data,
				   rawfittedfull12data,
				   rawfittedfull16data,
				   rawfittedfull20data,
				   rawfittedfull24data)

fitteddata <- bind_rows(fitteddef12data,
				fitteddef16data,
				fitteddef20data,
				fitteddef24data,
				fittedfull12data,
				fittedfull16data,
				fittedfull20data,
				fittedfull24data)

# Output data as csv

# write.csv(fitteddata, "fitteddata05_2017_22_01.csv")
# write.csv(rawfitteddata, "rawfitteddata05_2017_22_01.csv")

# use this function on the "rawfittedx" subsets
plotbest3 <- function(obsdata, fitdata) {

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

best10simdata <- map_df(fitdata, innerfunction)

prod_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = P)) +
		geom_line(data = best10simdata, aes(x = time, y = P, color = repnumber))


het_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = H)) +
		geom_line(data = best10simdata, aes(x = time, y = H, color = repnumber))

output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)

return(output_plot)

}

plotbest3(def12data, rawfitteddef12data)

### PLOTS ###

fittedr_plot <- ggplot(data = fitteddata, aes(x = transformedtemperature, y = log(r), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(r) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log intrinsic growth rate (r)")
fittedr_plot

r_model <- lm(data = fitteddata, log(r) ~ transformedtemperature)
confint(r_model)

fittedK_plot <- ggplot(data = fitteddata, aes(x = transformedtemperature, y = log(K), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(K) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log carrying capacity (K)")
fittedK_plot

fitteda_plot <- ggplot(data = fitteddata, aes(x = transformedtemperature, y = log(a), color = phosphorus)) +
        geom_point() +
        geom_smooth(method = lm) +
        ggtitle("Fitted log(a) Values") +
        labs(x = "inverse temperature (-1/kT)", y = "log attack rate (a)")
fitteda_plot
