### R Packages ###

# Load the simecol package for fitting ordinary differential equation models
library(simecol)
# Load the tidyverse package for improved data manipulation and plotting functions
library(tidyverse)

### log-scaled Coexistence Equilibria functions ###

# The overall goal of this script is to produce plots of equilibrium dynamics
# for producers and heterotrophs, analogous to those found in fig 1 of O'Connor
# et al. (2011).

### Declare parameters ###

# In O'Connor et al. (2011) they used the parameters fit using the lowest
# temperature. In our case that corresponds to 12 C or 285.15 K.

# !!! Unfortunately, because the model fitting is not currently working, I
# will, in the interim, be using the parameters used to set the initial
# conditions for the differential equations themselves. !!!

Ep <- -0.32
Eh <- 0.65
boltz <- 8.62 * 10 ^ (-5) # Boltzmann constant
r <- 0.1
a <- 2.1
m <- 2
eps <- 1
b <- 50000000
K <- 10 ^ 7

# Function for determining equilibrium density for producers.

Pfunc <- function(T){
	output <- log( (m * b) / (eps * a - m))
return(output)
}

# Function for determining equilibrium density for heterotrophs.

Hfunc <- function(T){
	output <- ((Eh - Ep) / ( boltz * T)) + log(
		(
		(b * eps * r) *
		K * (eps * a - m) - b * m * exp( (-Ep / (boltz * T))) /
		(K * (eps * a - m) ^ 2)
		)
	)
return(output)
}


boltzscale <- function(T){
	output <- 1 / (boltz * T)
return(output)
}

temp <- seq(285.15, 297.15, 0.5)
scaledtemp <- sapply(temp, boltzscale)
logPhat <- rep(Pfunc(290.15), length(temp))
logHhat <- sapply(temp, Hfunc)

predictdata <- data.frame(temp, scaledtemp, logPhat, logHhat)

predict_plot <- ggplot(data = predictdata) + # declare data
	geom_line(aes(x = temp, y = logPhat, colour = "Producer")) + 
	geom_line(aes(x = temp, y = logHhat, colour = "Heterotroph")) +
	ggtitle("Theoretical Consumer-Resource Density") +
	labs(x="Temperature",y="Density") 
predict_plot
