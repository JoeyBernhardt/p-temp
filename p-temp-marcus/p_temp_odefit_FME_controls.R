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
										P = biovol, # used biovol before
										temperature = temp
									  )

# Isolate only the data involving controls, labelled in the original data with "CX"; X is some number
controldata <- filter(controldata, grepl("C", replicate))

# Here we scale the phytoplankton density downwards by six orders of magnitude.
scalefactor <- 10 ^ 6 # use 10 ^ 6 for biovol
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

pdata <- controldata[["FULL16"]]

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

sink("full16control_MCMC_parameter_summary.txt")
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

size_box <- boxplot(data = rawdata, vol_cell ~ temp)

## k0 is the ref K
## m is body mass
## ss is the temperature size slope (i.e. change in size due to change in temp)
## t is the temperature
## k is boltzman's constant
## EP is the activation energy of photosynthesis


KMT2 <- function(k0, m, EM, EP, k, t) k0*(m^(-3/4))*exp((3*EM - 4*EP)/(4*k*t))

## draw K curve with TSR (in black)
curve(KMT2(k0 = 10000, m = 100, EP = -0.30, EM = -0.4, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

## now add curve for K without TSR, shown in blue
curve(KMT2(k0 = 10000, m = 100, EP = -0.32, EM = 0, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", add = TRUE, col = "blue")


KMT3 <- function(k0, m, EM, EP, k, t) k0*((m*exp((-EM/(k*t))))^(-3/4))*exp((-EP/(k*t)))

## draw K curve with TSR (in black)
curve(KMT3(k0 = 10000, m = 100, EP = -0.32, EM = -0.32, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

## now add curve for K without TSR, shown in blue
curve(KMT3(k0 = 10000, m = 100, EP = -0.32, EM = 0, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", add = TRUE, col = "blue")


Barr <- function(EP,k,t, m) (m^(-3/4))*exp(-EP/(k*t))

curve(Barr(m = 500, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15, to=273.15+10, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

JMT <- function(k0, m, ss, EP, k, t) k0*((m + ss*(t-273.15))^(-3/4))*exp(-EP/(k*t))

## draw K curve with TSR (in black)
curve(JMT(k0 = 10000, m = 500, ss = -15, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

## now add curve for K without TSR, shown in blue
curve(JMT(k0 = 10000, m = 500, ss = 0, EP = -0.32, k = 8.62 * 10^(-5), x), from=273.15, to=273.15+20, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", add = TRUE, col = "blue")


DMT <- function(k0, m, EM, k, t)  k0*m*(exp(-EM/(k*t)))^(-3/4)

#curve(mmt(k0 = 10000, m = 100, EP = -0.32, EM = 0, k = 8.62 * 10^(-5), x), from=273.15-20, to=273.15+50, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y")

KMT <- function(k0, m, ss, EP, k, t) k0*((m + ss*(t-273.15) + 30)^(-3/4))*exp(-EP/(k*t))

## draw K curve with TSR (in green)
curve(KMT(k0 = 100, m = 10, ss = (-3/100)*10, EP = -0.32, k = 8.62 * 10^(-5), x), from=270, to=310, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "green")

## now add curve for K without TSR, shown in red
curve(KMT(k0 = 100, m = 10, ss = 0, EP = -0.32, k = 8.62 * 10^(-5), x), from=270, to=310, xlab="Temperature (Kelvins)", ylab="log(K)", log = "y", col = "red", add = TRUE)

require(ggplot2)
KMT <- function(x) 6.5*((1000 + ((-2/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT2 <- function(x) 6.5*((1000 + ((0/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))
KMT3 <- function(x) 6.5*((1000 + ((-2.5/100)*1000)*(x-278.15))^(-3/4))*exp(0.32/(8.62 * 10^(-5)*x))


p <- ggplot(data = data.frame(x = 0), mapping = aes(x = x))

p + 
	stat_function(fun = KMT, color = "black", size = 2) +  stat_function(fun = KMT2, color = "cadetblue", size = 2) +
	stat_function(fun = KMT3, color = "green", size = 2) +
	xlim(273.15 + 5, 273.15 + 32) +
	scale_y_continuous(trans = "log", breaks = 5) + theme_bw() +
	xlab("temperature (kelvin)") + ylab("ln carrying capacity (K)") +
	theme(text = element_text(size=20))