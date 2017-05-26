#### R Packages ####

# Load the FME package for fitting ordinary differential equation models
library(FME)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)
# Load the gridExtra package for conveniently arranging ggplot objects
library(gridExtra)
# Load the knitr package for producing tables and graphics with markdown
library(knitr)

#### Data Frames ####

# Read the day zero data, and store as an object
initialdata <- read.csv(file = file.path("data-processed", "march29_cell_data.csv"), # file.path() is used for cross-platform compatibility
												strip.white = TRUE, # remove leading and trailing spaces from character string entries
												na.strings = c("NA","") # treat empty fields as missing
												)

# Read the data, and store as an object
rawdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"), # file.path() is used for cross-platform compatibility
											strip.white = TRUE, # remove leading and trailing spaces from character string entries
											na.strings = c("NA","") # treat empty fields as missing
											)

#### Data Manipulation ####

# Remove unneeded columns from "initialdata" dataframe
initialdata <- initialdata[-c(1, 5, 7)] # column 1 is "filename"; column 5 is "date"; column 7 is "start time"

# Rename columns in "initialdata" to correspond to those in "rawdata"
initialdata <- rename(initialdata,
											phosphorus_treatment = nutrient_level,
											volume_cell = cell_volume,
				     					algal_cell_concentration_cells_per_ml = cell_density
											)

# Calculate the initial algal biovolume, and also add in some new columns corresponding to the day zero values for days and daphnia.
# Here, the algal biovolume is calculated by multiplying the mean cell volume by the algal concentration.
initialdata <- mutate(initialdata,
											algal_biovolume = volume_cell * algal_cell_concentration_cells_per_ml,
				     					days = 0,
				     					daphnia_total = 10
											)

# Vertically merge "rawdata" and "initialdata" data frames
rawdata <- bind_rows(rawdata, initialdata)

# Rename columns in the "ptempdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary before we can initialize the fitting process.
ptempdata <- rename(rawdata,
										P = algal_biovolume,
										H = daphnia_total,
										phosphorus = phosphorus_treatment
									  )

# Here we scale the phytoplankton density downwards by six orders of magnitude.
scalefactor <- 10 ^ 6
ptempdata <- mutate(ptempdata,
										P = P / scalefactor
										)

# Create a column for boltzmann-transformed temperatures; this is useful for plotting purposes later.
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
ptempdata <- mutate(ptempdata,
										transformedtemperature = -1/(Boltz * (temperature + 273.15))
										)

# Reorder rows in "ptempdata" by treatment (phosphorus x temperature), replicate ID, and then by days
ptempdata <- arrange(ptempdata,
										 phosphorus,
										 temperature,
										 replicate,
										 days
										 )

# Create a column that lists the treatment type of each replicate, which depends on its phosphorus and temperature combination
ptempdata <- mutate(ptempdata,
										treatment = paste(ptempdata$phosphorus, ptempdata$temperature, sep="")
										)

# Split entire dataset into multiple indexed data frames based on their treatment
ptempdata <- split(ptempdata, f = ptempdata$unique_ID)

#### Objects ####

rkfit <- function(pdata){
# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, 36, 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P), H = 10)

model_parameters <- c(r = 4, K = 505, a = 0.11, eps = 0.008, m = 0.03)
lower_parameters <- c(0, 50, 0.1, 0.001, 0)
upper_parameters <- c(5, 1000, 1, 0.5, 0.1)

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
		  method = c("Nelder-Mead") # options are "Marq" or "Nelder-Mead"
		  )

fitsoda <- ode(initial_state, times, CRmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

return(Fit$par)
}

results <- map_df(ptempdata, rkfit)


