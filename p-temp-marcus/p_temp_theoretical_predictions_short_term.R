### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)
# Load the gridExtra package for conveniently arranging ggplot objects
library(gridExtra)

### Objectives ###

# Here we use a for-loop in order to iteratively model the short-term
# population dynamics at different temperatures. This script can be
# conceptually broken down into three sections:

# Pre-loop
# 1. Declare Arrhenius function
# 2. Generate an empty "master" dataframe
# 3. Declare parameters used in simulated models
# 4. Declare temperature range that the for-loop will iterate over

# For-loop
# At each step i of the for-loop, the following happens:
# 1. Simulate CR dynamics using the temperature i.
# 2. Output the simulation results as a dataframe
# 3. Subset the above dataframe to obtain a smaller one corresponding to the population densities on "day" 30.
# 4. Merge this dataframe into the master data frame.

# Post-loop
# Generate plot using ggplot2.

### Pre-Loop Setup ###

# Create an Arrhenius function to transform metabolic rates based on
# temperature. Here, "T" is the temperature. "E" is the activation energy
# constant. Its structure is identical to the one used in O'Connor et al.
# 2011

# Declare constants for Arrhenius function

Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
BasalTemperature <- 12 + 273.15 # the "base temperature" that determines the basal metabolic rate

# Declare Arrhenius function
arrhenius <- function(T,E){
	output <- exp(E * ((T + 273.15) - BasalTemperature) / (Boltz * (T + 273.15) * BasalTemperature))
return(output)
}

# Here we generate an empty dataframe which we will subsequently fill with data
# on the low-resource short-term dynamics, using the for-loop.
LowResourceDF <- data.frame(temp = double(),
			   P = double(), 
                	   H = double(), 
                	   stringsAsFactors = FALSE)

# Here we do the same as above, only for the high-resource dynamics.
HighResourceDF <- data.frame(temp = double(),
			   P = double(), 
                	   H = double(), 
                	   stringsAsFactors = FALSE)

# Declare the parameters to be used in the dynamical models 
LowResourceParameters <- c(r = 1, K = 10000000, a = 3, b = 500000, eps = 0.2, m = 0.2, Er = 0.32, EK = -0.32, Ea = 0.65, Em = 0.65)

HighResourceParameters <- replace(LowResourceParameters, "K", LowResourceParameters["K"] * 10)

# Declare the temperatures over which we will be simulating our consumer-resource model.
TempRange <- seq(12, 24, 0.5) # We are principally interested in these temperatures: 12, 16, 20, 24

### For-Loops ###

# Here we declare both for-loops; one each for the low and high resource dynamics

# Low Consumer Resources
for (i in TempRange) {
temp <- i
CRmodel <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
	dp <- arrhenius(temp, Er) * r * P * (1 - (P / (arrhenius(temp, EK) * K))) - arrhenius(temp, Ea) * a * H * (P / (P + b))
	dh <- arrhenius(temp, Ea) * a * eps * H * (P / (P + b)) - arrhenius(temp, Em) * m * H
	list(c(dp, dh))
	})
	},
	parms = LowResourceParameters,
	times = c(from = 0, to = 35, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 500000, H = 10),
	solver = "lsoda" # We use lsoda here because it was used in the O'Connor paper. Are there better methods?
)
CRmodeloutput <- out(sim(CRmodel)) # Output the simulated dynamics of the above model as a dataframe
CRmodeloutput <- filter(CRmodeloutput, time == 30) # Select only the row corresponding to "day" 30
CRmodeloutput <- select(CRmodeloutput, P:H) # Remove the "time" column
CRmodeloutput <- mutate(CRmodeloutput, temp = i) # Add a column corresponding to the temperature i
LowResourceDF <- rbind(LowResourceDF, CRmodeloutput) # Merge dataframe into "LowResourceDF" dataframe
}

# High Consumer Resources
for (i in TempRange) {
temp <- i
CRmodel <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
	dp <- arrhenius(temp, Er) * r * P * (1 - (P / (arrhenius(temp, EK) * K))) - arrhenius(temp, Ea) * a * H * (P / (P + b))
	dh <- arrhenius(temp, Ea) * a * eps * H * (P / (P + b)) - arrhenius(temp, Em) * m * H
	list(c(dp, dh))
	})
	},
	parms = HighResourceParameters,
	times = c(from = 0, to = 35, by = 0.1), # the time interval over which the model will be simulated.
	init = c(P = 500000, H = 10),
	solver = "lsoda" # We use lsoda here because it was used in the O'Connor paper. Are there better methods?
)
CRmodeloutput <- out(sim(CRmodel)) # Output the simulated dynamics of the above model as a dataframe
CRmodeloutput <- filter(CRmodeloutput, time == 30) # Select only the row corresponding to "day" 30
CRmodeloutput <- select(CRmodeloutput, P:H) # Remove the "time" column
CRmodeloutput <- mutate(CRmodeloutput, temp = i) # Add a column corresponding to the temperature i
HighResourceDF <- rbind(HighResourceDF, CRmodeloutput) # Merge dataframe into "HighResourceDF" dataframe
}

### Plotting ###

# Generate plot of the simulated data under low-resource conditions
low_resource_short_term_plot <- ggplot(data = LowResourceDF) + # declare dataframe
	geom_point(aes(x = temp, y = P, colour = "Producer")) +
	geom_point(aes(x = temp, y = H, colour = "Heterotroph")) +
	scale_y_log10() +
	ggtitle("Theoretical Consumer Resource Density at 30 Days; Low Resources") +
	labs(x = "Temperature", y = "Density")

# Generate plot of the simulated data under high-resource conditions
high_resource_short_term_plot <- ggplot(data = HighResourceDF) + # declare dataframe
	geom_point(aes(x = temp, y = P, colour = "Producer")) +
	geom_point(aes(x = temp, y = H, colour = "Heterotroph")) +
	ggtitle("Theoretical Consumer Resource Density at 30 Days; High Resources") +
	labs(x = "Temperature", y = "Density")

# Display both plots
grid.arrange(low_resource_short_term_plot, high_resource_short_term_plot, ncol = 2)


#### Joey's plots (just playing around here!)

HighResourceDF %>% 
	filter(temp < 15) %>% 
ggplot(data = .,  aes(x = temp, y = P)) + geom_point() +
	scale_y_log10()
ggplot(data = HighResourceDF, aes(x = temp, y = H)) + geom_point() +
	scale_y_log10()
