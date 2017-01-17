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

# This vector simply contains strings; they are used to tell the function
# "fitOdeModel" which parameters it is supposed to fit
FittedParameters <- c("r", "K", "a", "eps", "m")

# Declare the parameters to be used as the bounds for the fitting algorithm

r_min <- 0.5
r_max <- 8

K_min <- 1e10
K_max <- 1e12

a_min <- 0.1
a_max <- 1

eps_min <- 1e-12
eps_max <- 1e-10

m_min <- 0.01
m_max <- 0.1

LowerBound <- c(r = r_min, K = K_min, a = a_min, eps = eps_min, m = m_min)
UpperBound <- c(r = r_max, K = K_max, a = a_max, eps = eps_max, m = m_max)

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

portfit <- function(data) {

		dayzerodata <- filter(data, days == 0)
		model <- CRmodel # initial the simecol model object
		parms(model)[FittedParameters] <- c(r = runif(1, min = r_min, max = r_max),
								K = runif(1, min = K_min, max = K_max),
								a = runif(1, min = a_min, max = a_max),
								eps = runif(1, min = eps_min, max = eps_max),
								m = runif(1, min = m_min, max = m_max)
								)
		init(model) <- c(P = mean(dayzerodata$P), H = 10) # Set initial model conditions to the mean biovolume taken from the first measurement day
		obstime <- data$days # The X values of the observed data points we are fitting our model to
		yobs <- select(data, P, H) # The Y values of the observed data points we are fitting our model to

		# Below we fit a CRmodel to the replicate's data. The optimization criterion used here is the minimization of the sum of
		# squared differences between the experimental data and our modelled data. This
		# is fairly standard, although alternatives do exist.
		
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
		# the replicates' ID, Phosphorus level, temperature, and the
		# fitted parameters.

		# set model parameters to fitted values and simulate

		ssq <- ssqOdeModel(coef(fittedmodel), model, obstime, yobs)
		phosphorus <- data$Phosphorus[1]
		temperature <- data$temperature[1]
		transformedtemperature <- data$transformedtemperature[1]
		r <- coef(fittedmodel)["r"]
		K <- coef(fittedmodel)["K"]
		a <- coef(fittedmodel)["a"]
		eps <- coef(fittedmodel)["eps"]
		m <- coef(fittedmodel)["m"]

		simulation_data <- data.frame(ssq, phosphorus, temperature, transformedtemperature, r, K, a, eps, m)

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

num <- 50
fitteddef12data <- trimdata(repfit(def12data, num))
fitteddef16data <- trimdata(repfit(def16data, num))
fitteddef20data <- trimdata(repfit(def20data, num))
fitteddef24data <- trimdata(repfit(def24data, num))

fittedfull12data <- trimdata(repfit(full12data, num))
fittedfull16data <- trimdata(repfit(full16data, num))
fittedfull20data <- trimdata(repfit(full20data, num))
fittedfull24data <- trimdata(repfit(full24data, num))

fitteddata <- bind_rows(fitteddef12data,
				fitteddef16data,
				fitteddef20data,
				fitteddef24data,
				fittedfull12data,
				fittedfull16data,
				fittedfull20data,
				fittedfull24data)

write.csv(fitteddata, "fitteddata.csv")

targetdata <- def16data

#def24data
# SimParameters <- c(r = 8.000000, K = 231170612575, a = 0.2096601, eps = 4.818179e-11, m = 0.02274694)
#def16data
SimParameters <- c(r = 6.0269286, K = 2.040687e+11, a = 0.3919607, eps = 5.000000e-12, m = 0.01026190)
simfit <- function(data){

		dayzerodata <- filter(data, days == 0)
		simmodel <- CRmodel
		init(simmodel) <- c(P = mean(dayzerodata$P), H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
		parms(simmodel) <- SimParameters

		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		
		return(simdata)
}

simulateddata <- simfit(targetdata)

### PLOTS ###

prod_plot <- ggplot() +
		geom_point(data = targetdata, aes(x = days, y = P)) +
		geom_line(data = simulateddata, aes(x = time, y = P), color = "red")


het_plot <- ggplot() +
		geom_point(data = targetdata, aes(x = days, y = H)) +
		geom_line(data = simulateddata, aes(x = time, y = H), color = "red")

grid.arrange(prod_plot, het_plot, ncol=2)

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
