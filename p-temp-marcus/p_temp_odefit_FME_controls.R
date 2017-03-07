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

# Read the data from the consumer free controls, and store as an object
rawdata <- read.csv(file = file.path("data-processed", "p_temp_algae.csv"), # file.path() is used for cross-platform compatibility
												strip.white = TRUE, # remove leading and trailing spaces from character string entries
												na.strings = c("NA","") # treat empty fields as missing
												)

#### Data Manipulation ####

# Rename columns in the "ptempdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary before we can initialize the fitting process.
controldata <- rename(rawdata,
										phosphorus = P,
										P = cellcon, # used biovol before
										temperature = temp
									  )

# Isolate only the data involving controls, labelled in the original data with "CX"; X is some number
controldata <- filter(controldata, grepl("C", replicate))

# Here we scale the phytoplankton density downwards by six orders of magnitude.
scalefactor <- 10 ^ 3 # use 10 ^ 6 for biovol
controldata <- mutate(controldata,
										P = P / scalefactor
										)

# Create a column for boltzmann-transformed temperatures; this is useful for plotting purposes later.
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
controldata <- mutate(controldata,
										transformedtemperature = -1/(Boltz * (temperature + 273.15))
										)

# Reorder rows in "ptempdata" by treatment (phosphorus x temperature), and then replicate ID
controldata <- arrange(controldata,
										 phosphorus,
										 temperature,
										 ID
										 )

# Create a column that lists the treatment type of each replicate, which depends on its phosphorus and temperature combination
controldata <- mutate(controldata,
										treatment = paste(controldata$phosphorus, controldata$temperature, sep="")
										)

# Split entire dataset into multiple indexed data frames based on their treatment
controldata <- split(controldata, f = controldata$treatment)

#### Objects ####

model_times <- seq(0, 36, 0.1)
initial_state <- c(P = 100)
model_parameters <- c(r = 1, K = 1000)

Pmodel <- function(model_times, initial_state, model_parameters) {
	with(as.list(c(initial_state, model_parameters)), {
		dP <-  r * P * (1 - (P /  K))
		list(c(dP))
	})
}

#### Functions ####

pdata <- controldata[["DEF16"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, max(pdata$days), 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P))

model_parameters <- c(r = 0.5, K = 1000)
lower_parameters <- c(0, 10)
upper_parameters <- c(1, 4000)

Pmodel <- function (times, initial_state, model_parameters) {
		with(as.list(c(initial_state, model_parameters)), {
	dP <-  r * P * (1 - (P /  K))
	list(c(dP))
	})
	}


ModelCost2 <- function(parms) {
odemodel <- ode(initial_state, times, Pmodel, parms)
cost <- modCost(odemodel, fittingdata, scaleVar = TRUE) # Note that scaleVar is true here...
return(cost) # object of class modCost
}

Fit <- modFit(f = ModelCost2,
		  p = model_parameters,
		  lower = lower_parameters,
		  upper = upper_parameters,
		  method = c("Marq") # options are "Marq" or "Nelder-Mead"
		  )

fitsoda <- ode(initial_state, times, Pmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

prod_plot <- ggplot() + # declare ggplot object
	geom_line(data = fitsodadf, aes(x = times, y = P, colour = "red")) +
	geom_point(data = fittingdata, aes(x = time, y = P)) +
  ggtitle("Simulated Algal Biovolume") +
	labs(x = "Days", y = "Algal Biovolume") +
	theme(legend.position = "none")
prod_plot

set.seed(7) # set pseudorandom seed; for reproducibility
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
MCMC <- modMCMC(f = ModelCost2, p = Fit$par, 
					  lower = lower_parameters,
					  upper = upper_parameters,
					  niter = 50000, jump = cov0,
            var0 = var0, wvar0 = 0.1,
					  updatecov = 50)
summary(as.mcmc(MCMC$pars))

sink("full12control_MCMC_parameter_summary.txt")
summary(as.mcmc(MCMC$pars))
sink()

size_plot <- ggplot(data = rawdata, aes(x = temp, y = vol_cell, color = P)) + # declare ggplot object
	geom_point() +
	geom_smooth(method = lm) +
	ggtitle("Cell Volume vs. temperature") +
	labs(x = "Temperature", y = "Cell volume") +
	theme(legend.position = "none")
size_plot

vol.lm <- lm(data = rawdata, vol_cell ~ temp)
summary(vol.lm)

size_box <- boxplot(data = rawdata, vol_cell ~ temp)