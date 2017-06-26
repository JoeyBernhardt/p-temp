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
# The units of the phytoplankton density in the raw data are cells/mL.
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

# Split "ptempdata" dataframe into multiple indexed data frames based on their treatment
ptempdata.treatment <- split(ptempdata, f = ptempdata$treatment)

# Split "ptempdata" dataframe into multiple indexed data frames based on their replicate ID
ptempdata.replicate <- split(ptempdata, f = ptempdata$replicate)


pdata <- ptempdata.treatment[["FULL12"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

# works well for def12: 4, 337, 3.901555e-01 4.001129e-03 3.700928e-02 USE MARQ FOR DEF12
# works well for full12: 5.654957e+00 2.025641e+02 3.683304e-01 1.261719e-03 1.009373e-03 USE NM

times <- seq(0, 36, 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P), H = 10)

temp = pdata$temperature[1]

model_times <- seq(0, 36, 0.1)
initial_state <- c(P = mean(dayzerodata$P), H = 10)
model_parameters <- c(r = 5.65, K = 200, a = 0.36, eps = 0.0012, m = 0.03)
lower_parameters <- c(0, 50, 0.1, 0.0001, 0)
upper_parameters <- c(6, 1000, 1, 0.1, 0.1)

CRmodel <- function(model_times, initial_state, model_parameters) {
	with(as.list(c(initial_state, model_parameters)), {
		dP <-  r * P * (1 - (P / K)) - a * H * P 
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

fitfull12<- Fit$par

#### Functions ####

pdata <- ptempdata.treatment[["FULL16"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, 36, 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P), H = 10)

temp = pdata$temperature[1]

## Arrhenius Function ##

# Create an Arrhenius function to transform metabolic rates based on
# temperature. Here, "T" is the temperature. "E" is the activation energy
# constant. Its structure is identical to the one used in O'Connor et al.
# 2011

# Declare constants for Arrhenius function
Boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
BasalTemperature <- 12 + 273.15 # the "base temperature" that determines the basal metabolic rate; 12 was the lowest temp used during the experiment.

arrhenius <- function(E){
	output <- exp(E * ((temp + 273.15) - BasalTemperature) / (Boltz * (temp + 273.15) * BasalTemperature))
	return(output)
}

model_times <- seq(0, 36, 0.1)
initial_state <- c(P = mean(dayzerodata$P), H = 10)
model_parameters <- c(Er = 0.1, EK = 0.5, Ea = 0.5, Eeps = 0.1, Em = 0.5)
lower_parameters <- c(-1, -1, -1, -1, -1)
upper_parameters <- c(1, 1, 1, 1, 1)

basevec <- fitfull12
r0 <- basevec[1]
K0 <- basevec[2]
a0 <- basevec[3]
eps0 <- basevec[4]
m0 <- basevec[5]

arrvec <- sapply(model_parameters, function(x) arrhenius(x))
vecdf <- data.frame(arrvec,basevec)
effvec <- transmute(vecdf, effvec = arrvec*basevec)

CRmodel <- function(model_times, initial_state, model_parameters) {
	with(as.list(c(initial_state, model_parameters)), {
		dP <-  r0*arrhenius(Er) * P * (1 - (P / (K0*arrhenius(EK)))) - a0*arrhenius(Ea) * H * P 
		dH <-  a0*arrhenius(Ea) * eps0*arrhenius(Eeps) * H * P - m0*arrhenius(Em) * H
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
		  method = c("Marq") # options are "Marq" or "Nelder-Mead"
		  )

fitsoda <- ode(initial_state, times, CRmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

fitfull16 <- Fit$par

# Gather results and plot
	
def16 <- c(Er = -0.1452267, EK = 0.8430171)
def20 <- c(Er = 0.3043211, EK = 0.5400571)
def24 <- c(Er = 0.59267223, EK = 0.59975538)

Er <- c(-0.1452267, 0.3043211, 0.59267223)
EK <- c(0.8430171, 0.5400571, 0.59975538)
phos <- rep("DEF", 3)
def <- data.frame(Er,EK, phos)

full16 <- c(Er = 0.1507282, EK = 0.4899299)
full20 <- c(Er = 0.59997975, EK = 0.07442368)
full24 <- c(Er = 0.4999998, EK = 0.6301322)

Er <- c(0.1507282, 0.59997975, 0.4999998)
EK <- c(0.4899299, 0.07442368, 0.63011322)
phos <- rep("FULL", 3)

full <- data.frame(Er,EK,phos)

res.df <- rbind(def,full)

write.csv(res.df, "activation_energies_reworked_2017_26_06.csv")

Er_plot <- ggplot(res.df, aes(phos, Er))+
	geom_boxplot(outlier.shape=1) +    
	ggtitle("Effect of Nutrient Availability on the Activation Energy \nof Population Intrinsic Growth Rate") + # set title
	labs(x = "Phosphorus Treatment", y = expression("Activation Energy of Intrinsic Growth rate r")) 
Er_plot

EK_plot <- ggplot(res.df, aes(phos, EK))+
	geom_boxplot(outlier.shape=1) +    
	ggtitle("Effect of Nutrient Availability on the Activation Energy \nof Population Carrying Capacity") + # set title
	labs(x = "Phosphorus Treatment", y = expression("Activation Energy of Carrying Capacity K"))
EK_plot

#### Plotting Functions ####

prod_plot <- ggplot() + # declare ggplot object
	geom_line(data = fitsodadf, aes(x = times, y = P, colour = "red")) +
	geom_point(data = fittingdata, aes(x = time, y = P)) +
  ggtitle("Simulated Algal Biovolume") +
	labs(x = "Days", y = "Algal Biovolume") +
	theme(legend.position = "none")
	
het_plot <- ggplot() + 
	geom_line(data = fitsodadf, aes(x = times, y = H, colour = "red")) +
	geom_point(data = fittingdata, aes(x = time, y = H)) +
  ggtitle("Simulated Daphnia Density") +
	labs(x = "Days", y = "Total Daphnia Density") +
  theme(legend.position = "none")

output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)


#### MCMC ####

set.seed(7) # set pseudorandom seed; for reproducibility
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
MCMC <- modMCMC(f = ModelCost2, p = Fit$par, 
								lower = lower_parameters,
								upper = upper_parameters,
								niter = 10000, jump = cov0,
								var0 = var0, wvar0 = 0.1,
								updatecov = 50)
summary(as.mcmc(MCMC$pars))

sink("full16_MCMC_parameter_summary.txt")
summary(as.mcmc(MCMC$pars))
sink()
