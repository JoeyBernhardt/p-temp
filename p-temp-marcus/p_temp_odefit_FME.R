### R Packages ###

# Load the FME package for fitting differential equation models
library(FME)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)

# Read the data from the consumer-resource treatments, and store them as an object
pdata <- read.csv(file = file.path("data-processed", "p_temp_processed.csv"),
	stringsAsFactors = FALSE,
	strip.white = TRUE,
	na.strings = c("NA","") )

### Data manipulation ###

# In this section we are extracting only the data for the first replicate
# mesocosm, and formatting the data for fitting a dynamical model. The
# fitOdeModel function in the simecol package requires very specific object
# names.

# If you would like to fit the model using a different replicate, simply change
# the value of unique_ID

uniqueid <- 1
pdata <- filter(pdata, unique_ID == uniqueid)

# Arrange data by sampling day, in ascending order
pdata <- arrange(pdata, days)

# Rename columns in data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
pdata <- rename(pdata, P = algal_biovolume)
pdata <- rename(pdata, H = daphnia_total)

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

fittingdata <- data.frame(obstime, yobs)
fittingdata <- rename(fittingdata, time = obstime)

times <- seq(0, 40, 0.1)
initial_state <- c(P = yobs$P[1], H = 10)
model_parameters <- c(r = 1, K = 10 ^ 9, a = 3, b = 5000000, eps = 0.1, m = 0.2)
lower_parameters <- c(0, 10 ^ 5, 0, 0, 0, 0)
upper_parameters <- c( 5, 10 ^ 20, 5, 10 ^ 10, 5, 5)

CRmodel <- function (times, initial_state, model_parameters) {
		with(as.list(c(initial_state, model_parameters)), {
	dP <-  r * P * (1 - (P /  K)) - a * H * (P / (P + b))
	dH <-  a * eps * H * (P / (P + b)) -  m * H
	list(c(dP, dH))
	})
	}


Lsoda <- ode(initial_state, times, CRmodel, model_parameters) # lsoda is default method
dfsoda <- data.frame(Lsoda)

ModelCost2 <- function(parms) {
odemodel <- ode(initial_state, times, CRmodel, parms)
cost <- modCost(odemodel, fittingdata, scaleVar = TRUE)
return(cost) # object of class modCost
}

Fit <- modFit(f = ModelCost2, p = model_parameters, lower = lower_parameters, upper = upper_parameters)

fitsoda <- ode(initial_state, times, CRmodel, Fit$par)
fitsodadf <- data.frame(fitsoda)

producer_plot <- ggplot() + # declare ggplot object
	geom_line(data = fitsodadf, aes(x = times, y = P, colour = "Producer")) +
	geom_point(data = fittingdata, aes(x = time, y = P, colour = "Producer"))
	
hetero_plot <- ggplot() + 
	geom_line(data = fitsodadf, aes(x = times, y = H, colour = "Heterotroph")) +
	geom_point(data = fittingdata, aes(x = time, y = H, colour = "Heterotroph"))

