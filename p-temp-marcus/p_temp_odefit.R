### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)

### Data frame

# Read the data and store them as an object
pdata <- read.csv("C:\\Users\\Matt\\Google Drive\\Work\\joey_project\\raw_data\\p_temp_processed.csv",
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

# Insert a new row into our data frame, representing the initial conditions for each mesocosm.
day_three_cell_volume <- pdata$volume_cell[1]
day_one_cell_concentration <- 10 ^ 5
initial_algal_biovolume <- day_three_cell_volume * day_one_cell_concentration
add_row(pdata, unique_ID = uniqueid, daphnia_total = 10, algal_biovolume = initial_algal_biovolume, days = 1, .before = 1)

# Rename columns in data frame to correspond to the names of the stocks in our
# dynamical model. This is necessary to invoke the fitOdeModel function.
pdata <- rename(pdata, P = algal_biovolume)
pdata <- rename(pdata, H = daphnia_total)

# Extract from the above subset only what we require to fit our model
obstime <- pdata$days # this is the time interval over which we are fitting our model
yobs <- select(pdata, P, H)

### Model Construction ###

# Here we construct the consumer resource model. First the Arrhenius function
# is declared, and it is then used in the creation of the dynamical model. The
# resulting system of equations, called "hp", can be used to produce
# theoretical predictions. It has been seeded with some parameters in an effort
# to resemble the initial state of your experimental setup, but these need
# tweaking.

# Create an Arrhenius function to transform metabolic rates based on
# temperature. Its structure is identical to the one used in Mary's paper.

boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
temperature <- 25 + 273.15 # the temperature of the modelled system
btemperature <- 21 + 273.15 # the "base temperature" that determines the basal metabolic rate

arrhenius <- function(E){
	output <- exp(E * (temperature - btemperature) / (boltz * temperature * btemperature))
return(output)
}

# Create new odeModel object. This represents the "base" Lotka-Volterra
# consumer resource dynamics model that we would like to eventually fit, with
# added metabolic effects due to temperature.

# dp and dh are the differential equations for producers and heterotrophs,
# respectively. Likewise, P and H refer to the population densities. For
# producers, both the intrinsic growth rate r and the carrying capacity K are
# subject to metabolic scaling. For heterotrophs, metabolic scaling applies to
# the attack rate a, and the intrinsic mortality rate m.

# If you would like to view the model output, you can use the following
# commands afer hp has been evaluated:
# hp <- sim(hp)
# plot(hp)

hp <- new("odeModel",
	main = function (time, init, parms) {
		with(as.list(c(init, parms)), {
	dp <- arrhenius(Er) * r * P * (1 - (P / (arrhenius(EK) * K))) - arrhenius(Ea) * a * H * (P / (P + b))
	dh <- arrhenius(Ea) * a * eps * H * (P / (P + b)) - arrhenius(Em) * m * H
	list(c(dp, dh))
	})
	},
	parms = c(r = 1, K = 10000000, a = 3, b = 500000, eps = 0.1, m = 0.2, Er = 0.32, EK = -0.32, Ea = 0.65, Em = 0.65),
	times = c(from = 0, to = 400, by = 0.1),
	init = c(P = initial_algal_biovolume, H = 10),
	solver = "rk4"
)

### Model fitting ###

# Here we attempt to fit the consumer resource model we constructed above,
# using experimental data from the first replicate. Instructions for fitting
# the model to other replicates can be found above in the "Data Manipulation"
# section.

# The optimization criterion used here is the minimization of the sum of
# squared differences between the experimental data and our modelled data. This
# is fairly standard, although alternatives do exist. The specific algorithm we
# are using to do this is called PORT. This is a poor choice for fitting this
# particular model, but it evaluates quickly and is thus useful for debugging
# the code.

# The authors of the 2011 O'Connor et al. Am Nat paper used the BOBYQA
# algorithm. Currently, this algorithm doesn't function correctly!!! I have
# contacted Ben Gilbert to discuss this issue, but I have a number of ideas
# why.

# Declare a vector containing the parameters we would like to fit.
fittedparms <- c("r", "K", "a", "b", "eps", "m")

# Call fitOdeModel using the PORT algorithm.
fittedhp <- fitOdeModel(hp, whichpar = fittedparms, obstime, yobs,
  debuglevel = 0, fn = ssqOdeModel,
  method = "PORT", lower = 0,
  control = list()
)

# If you would like to output the fitted parameters, you can use:
# coef(fittedhp)

### Display model fitting results ###

# To display the fitted results we need to create a new OdeModel object. Here
# we duplicate hp and then alter it to use our new fitted parameters.
displayfittedhp <- hp
parms(displayfittedhp)[fittedparms] <- coef(fittedhp)

# set model parameters to fitted values and simulate again
times(displayfittedhp) <- c(from=0, to=40, by=1)
ysim <- out(sim(displayfittedhp))

# Plot the results of our model fitting. If you're reading this, sorry for the
# shitty plots. I'll use ggplot2 next time!
par(mfrow=c(2,1))
plot(obstime, yobs$P,
	ylim = range(yobs$P, ysim$P),
	xlab = "Days",
	ylab = "Algal Biovolume"
)
lines(ysim$time, ysim$P, col="red")

plot(obstime, yobs$H,
	ylim = range(yobs$H, ysim$H),
	xlab = "Days",
	ylab = "Daphnia Density"
)
lines(ysim$time, ysim$H, col="red")
