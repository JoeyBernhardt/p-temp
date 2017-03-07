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
ptempdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"), # file.path() is used for cross-platform compatibility
											strip.white = TRUE, # remove leading and trailing spaces from character string entries
											na.strings = c("NA","") # treat empty fields as missing
											)

#### Data Manipulation ####

# Remove unneeded columns from "initialdata" dataframe
initialdata <- initialdata[-c(1, 5, 7)] # column 1 is "filename"; column 5 is "date"; column 7 is "start time"

# Rename columns in "initialdata" to correspond to those in "ptempdata"
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

# Vertically merge "ptempdata" and "initialdata" data frames
ptempdata <- bind_rows(ptempdata, initialdata)

# Rename columns in the "ptempdata" data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary before we can initialize the fitting process.
ptempdata <- rename(ptempdata,
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
ptempdata <- split(ptempdata, f = ptempdata$treatment)

#### Objects ####

model_times <- seq(0, 36, 0.1)
initial_state <- c(P = 100, H = 10)
model_parameters <- c(r = 1, K = 1000, a = 1, eps = 0.1, m = 0.1)

PHmodel <- function(model_times, initial_state, model_parameters) {
	with(as.list(c(initial_state, model_parameters)), {
		dP <-  r * P * (1 - (P /  K)) - a * H * P 
		dH <-  a * eps * H * P - m * H
		list(c(dP, dH))
	})
}

#### Functions ####

SANNfit <- function(data){
	    
			model_times <- seq(0, max(data$days), 0.1)
			
			# Here we isolate all of the day zero data, which is then used to set our initial model conditions. The initial Daphnia density is always 10,
			# regardless of treatment group. The initial phytoplankton density is set to the mean biovolume taken from the first measurement day.
			dayzerodata <- filter(data, days == 0)
			initial_state <- c(P = mean(dayzerodata$P), H = 10)
			
			SANN_parameters <- c(r = 5, K = 500, a = 5, eps = 0.002, m = 0.01)
			SANN_lower_bound <- c(0, 50, 0, 0.001, 0)
			SANN_upper_bound <- c(10, 1500, 10, 0.01, 0.5)

			obsdata <- data %>%
								 select(days, P, H) %>%
								 rename(time = days)
			
			SANNcost <- function(parms) {
				odemodel <- ode(initial_state, model_times, PHmodel, parms)
				cost <- modCost(odemodel, obsdata, scaleVar = TRUE) # Note that scaleVar is true here...
				return(cost) # returns object of class modCost
																	}

			SANNfit <- modFit(f = SANNcost,
								        p = SANN_parameters,
								        lower = SANN_lower_bound,
												upper = SANN_upper_bound,
												method = c("SANN"),
									    	control = list(maxit = 5,
																	     temp = 50,
																	     tmax = 100)
												)
			
			SANNsim <- ode(initial_state, model_times, PHmodel, SANNfit$par)
			SANNfitdata <- data.frame(SANNsim)
			
			return(SANNfitdata)
}

sfitdf <- SANNfit(ptempdata[["FULL20"]])

obsdata <- ptempdata[["FULL20"]] %>%
	select(days, P, H) %>%
	rename(time = days)

			prod_plot <- ggplot() + # declare ggplot object
				geom_line(data = sfitdf, aes(x = model_times, y = P, colour = "red")) +
				geom_point(data = obsdata, aes(x = time, y = P)) +
				ggtitle("Simulated Algal Biovolume") +
				labs(x = "Days", y = "Algal Biovolume") +
				theme(legend.position = "none")
			
			het_plot <- ggplot() + 
				geom_line(data = sfitdf, aes(x = model_times, y = H, colour = "red")) +
				geom_point(data = obsdata, aes(x = time, y = H)) +
				ggtitle("Simulated Daphnia Density") +
				labs(x = "Days", y = "Total Daphnia Density") +
				theme(legend.position = "none")
			
			output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)

pdata <- ptempdata[["FULL24"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, 36, 0.1)
dayzerodata <- filter(pdata, days == 0)
initial_state <- c(P = mean(dayzerodata$P), H = 10)

# works well for def12: 5.524790e+00 3.378748e+02 3.901555e-01 4.001129e-03 3.700928e-02 USE MARQ FOR DEF12
# works well for full12: 5.654957e+00 2.025641e+02 3.683304e-01 1.261719e-03 1.009373e-03 USE NM
# works well for def16: c(r = 4.8, K = 585, a = 0.25, eps = 0.002, m = 0.01)
# works well for full16: c(r = 3.14, 882, a = 0.18, eps = 0.004, m = 0.01)
# works well for def20: c(r = 4.9, K = 997, a = 0.23, eps = 0.002, m = 0.03)
# works well for def24: c(r = 4.1, K = 700, a = 0.3, eps = 0.008, m = 0.02) WEIRD STUFF HAPPENING HERE. USE 10^5 for MCMC
# works SOMEWHAT well for full20: c(r = 2.1, K = 928, a = 0.12, eps = 0.005, m = 0.09)
# model_parameters <- c(r = 2, K = 500, a = 0.5, eps = 0.001, m = 0.01)

model_parameters <- c(r = 4.1, K = 700, a = 0.3, eps = 0.008, m = 0.02)
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
		  method = c("Nelder-Mead") # options are "Marq" or "Nelder-Mead"
		  )

fitsoda <- ode(initial_state, times, CRmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

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

set.seed(7) # set pseudorandom seed; for reproducibility
var0 <- Fit$var_ms_unweighted
cov0 <- summary(Fit)$cov.scaled * 2.4^2/5
MCMC <- modMCMC(f = ModelCost2, p = Fit$par, 
					  lower = lower_parameters,
					  upper = upper_parameters,
					  niter = 50000, jump = cov0,
            var0 = var0, wvar0 = 0.1,
					  updatecov = 50)
summary(MCMC)

sink("full24_MCMC_parameter_summary.txt")
summary(as.mcmc(MCMC$pars))
sink()

Sfun <- sensFun(ModelCost2, model_parameters)
ident <- collin(Sfun)


temp <- c(12, 16, 20, 24)
transformedtemp <- -1/(Boltz * (temp + 273.15))

defr <- c(5.428, 5.034, 5.130, 4.409)
defK <- c(556.0, 688.5, 629.1, 659.0)
defa <- c(0.4625, 0.2999, 0.1269, 0.3450)

fullr <- c(3.891, 4.432, 5.045, 5.408)
fullK <- c(423.3, 262.9, 495.5, 597.9)
fulla <- c(0.5348, 0.2062, 0.1343, 0.3168)

slope.df <- data.frame(temp,transformedtemp,defr,defK,defa,fullr,fullK,fulla)


defrlm <- lm(data = slope.df, log(defr) ~ transformedtemp)
summary(defrlm)

defr_plot <- ggplot(data = slope.df, aes(x = transformedtemp, y = log(defr))) + 
	geom_point() +
	geom_smooth(method = lm, col = "red") +
	ggtitle("Estimated values of r - Phosphorus Deficient") +
	labs(x = "Temperature (-1/kT)", y = "log(r)") +
	theme(legend.position = "none")
defr_plot

fullrlm <- lm(data = slope.df, log(fullr) ~ transformedtemp)
summary(fullrlm)

fullr_plot <- ggplot(data = slope.df, aes(x = transformedtemp, y = log(fullr))) + 
	geom_point() +
	geom_smooth(method = lm, col = "red") +
	ggtitle("Estimated values of r - Phosphorus Rich") +
	labs(x = "Temperature (-1/kT)", y = "log(r)") +
	theme(legend.position = "none")
fullr_plot

defKlm <- lm(data = slope.df, log(defK) ~ transformedtemp)
summary(defKlm)

defK_plot <- ggplot(data = slope.df, aes(x = transformedtemp, y = log(defK))) + 
	geom_point() +
	geom_smooth(method = lm, col = "red") +
	ggtitle("Estimated values of K - Phosphorus Deficient") +
	labs(x = "Temperature (-1/kT)", y = "log(K)") +
	theme(legend.position = "none")
defK_plot

fullKlm <- lm(data = slope.df, log(fullK) ~ transformedtemp)
summary(fullKlm)

fullK_plot <- ggplot(data = slope.df, aes(x = transformedtemp, y = log(fullK))) + 
	geom_point() +
	geom_smooth(method = lm, col = "red") +
	ggtitle("Estimated values of K - Phosphorus Rich") +
	labs(x = "Temperature (-1/kT)", y = "log(K)") +
	theme(legend.position = "none")
fullK_plot

# density plots

pdata <- ptempdata[["DEF24"]]

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

prod_plot <- ggplot() + # declare ggplot object
	geom_point(data = fittingdata, aes(x = time, y = P)) +
	ggtitle("Simulated Algal Biovolume") +
	labs(x = "Days", y = "Algal Biovolume") +
	theme(legend.position = "none")

het_plot <- ggplot() + 
	geom_point(data = fittingdata, aes(x = time, y = H)) +
	ggtitle("Simulated Daphnia Density") +
	labs(x = "Days", y = "Total Daphnia Density") +
	theme(legend.position = "none")

output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)