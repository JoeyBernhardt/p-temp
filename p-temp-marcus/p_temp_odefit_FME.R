### R Packages ###

# Load the FME package for fitting ordinary differential equation models
library(FME)
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

# Here we scale the phytoplankton density by 250, in order to get the total biovolume present in the beaker.
# We want to model the total algal biovolume because we have modelled the total Daphnia density.
scalefactor <- 1000000
controldata <- mutate(controldata, P = P / scalefactor)

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
ptempdata <- mutate(ptempdata, P = P / scalefactor)

# Create a column for boltzmann-transformed temperatures; this is useful for plotting purposes later. "Boltz" is defined
# earlier in the code, for a similar transformation of controldata
ptempdata <- mutate(ptempdata, transformedtemperature = -1/(Boltz * (temperature + 273.15)))

# Reorder rows in data frame by treatment ID, and then by days
ptempdata <- arrange(ptempdata, Phosphorus, temperature, replicate, days)

# Create a column that lists the treatment type of each replicate, which depends on its phosphorus and temperature combination
ptempdata <- mutate(ptempdata, treatment = paste(ptempdata$Phosphorus, ptempdata$temperature, sep=""))

# Split entire dataset into multiple indexed data frames based on their treatment
ptempdata <- split(ptempdata, f = ptempdata$treatment)

pdata <- ptempdata[["DEF16"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, 36, 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P), H = 10)

# works well for def16: c(r = 4.8, K = 585, a = 0.25, eps = 0.002, m = 0.01)
# works well for full16: c(r = 3.14, 882, a = 0.18, eps = 0.004, m = 0.01)
# works well for def20: c(r = 4.9, K = 997, a = 0.23, eps = 0.002, m = 0.03)
# works SOMEWHAT well for full20: c(r = 2.1, K = 928, a = 0.12, eps = 0.005, m = 0.09)
# model_parameters <- c(r = 2, K = 500, a = 0.5, eps = 0.001, m = 0.01)

model_parameters <- c(r = 4.8, K = 585, a = 0.25, eps = 0.002, m = 0.01)
lower_parameters <- c(0, 100, 0.1, 0.001, 0)
upper_parameters <- c(6, 1000, 1, 0.01, 0.1)

CRmodel <- function (times, initial_state, model_parameters) {
		with(as.list(c(initial_state, model_parameters)), {
	dP <-  r * P * (1 - (P /  K)) - a * H * P 
	dH <-  a * eps * H * P - m * H
	list(c(dP, dH))
	})
	}


Lsoda <- ode(initial_state, times, CRmodel, model_parameters) # lsoda is default method
dfsoda <- data.frame(Lsoda)

ModelCost2 <- function(parms) {
odemodel <- ode(initial_state, times, CRmodel, parms)
cost <- modCost(odemodel, fittingdata, scaleVar = TRUE) # Note that scaleVar is true here...
return(cost) # object of class modCost
}

Fit <- modFit(f = ModelCost2,
		  p = model_parameters,
		  lower = lower_parameters,
		  upper = upper_parameters,
		  method = c("Marq")
		  )

fitsoda <- ode(initial_state, times, CRmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

prod_plot <- ggplot() + # declare ggplot object
	geom_line(data = fitsodadf, aes(x = times, y = P, colour = "red")) +
	geom_point(data = fittingdata, aes(x = time, y = P))
	
het_plot <- ggplot() + 
	geom_line(data = fitsodadf, aes(x = times, y = H, colour = "red")) +
	geom_point(data = fittingdata, aes(x = time, y = H))

output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)

var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
MCMC <- modMCMC(f = ModelCost2, p = Fit$par, 
					  lower = lower_parameters,
					  upper = upper_parameters,
					  niter = 5000, jump = cov0,
                                var0 = var0, wvar0 = 0.1,
					  updatecov = 50)

Sfun <- sensFun(ModelCost2, model_parameters)
ident <- collin(Sfun)

