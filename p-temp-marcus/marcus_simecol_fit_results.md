# Fitting results

### Data

Here we read in the data used for all of the below plots:
```r
#file.path() is used for cross-platform compatibility
plotdata <- read.csv(file = file.path("p-temp-marcus", "outputs", "rawfitteddata05_2017_22_01.csv"),
	strip.white = TRUE,
	na.strings = c("NA","") )
```

### Plotting Simulated Density Curves vs. Observed Data

Here we show plots of the "best" three simulation results vs. the observed experimental data, for each treatment combination of phosphorus and temperature.

The goal here is to assess the quality of our model fits.

##### 12 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/def12plot.png" width="1000">

##### 12 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/full12plot.png" width="1000">

##### 16 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/def16plot.png" width="1000">

##### 16 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/full16plot.png" width="1000">

##### 20 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/def20plot.png" width="1000">
##### 20 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/full20plot.png" width="1000">

##### 24 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/def24plot.png" width="1000">
##### 24 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/full24plot.png" width="1000">

```r
plotbest3 <- function(phosphorus, temp) {

x <- paste(phosphorus, temp, sep = "")
obsdata <- ptempdata[[x]]

fitdata <- filter(rawfitteddata, treatment == x)

fitdata <- filter(fitdata, ssq != 0)
fitdata <- arrange(fitdata, ssq)
fitdata <- fitdata[1:3,] # use best 3 for now

fitdata$repnumber <- rownames(fitdata)
fitdata <- split(fitdata, f = fitdata$repnumber)

# objective is to write a function that operates on a single data frame of parameters and other info, and then produce the density estimates
# using the sim function in simecol. Then use map_df on this guy.

innerfunction <- function(xdata) {

fittedr <- xdata$r
fittedK <- xdata$K
fitteda <- xdata$a
fittedeps <- xdata$eps
fittedm <- xdata$m

SimParameters <- c(r = fittedr, K = fittedK, a = fitteda, eps = fittedeps, m = fittedm)

		dayzerodata <- filter(obsdata, days == 0)
		simmodel <- CRmodel
		init(simmodel) <- c(P = mean(dayzerodata$P), H = 10) # Set initial model conditions to the biovolume taken from the first measurement day
		parms(simmodel) <- SimParameters

		simdata <- out(sim(simmodel, rtol = 1e-10, atol = 1e-10))
		simdata <- mutate(simdata, repnumber = xdata$repnumber)

return(simdata)
}

best3simdata <- map_df(fitdata, innerfunction)

prod_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = P)) +
		geom_line(data = best3simdata, aes(x = time, y = P, color = repnumber))


het_plot <- ggplot() +
		geom_point(data = obsdata, aes(x = days, y = H)) +
		geom_line(data = best3simdata, aes(x = time, y = H, color = repnumber))

output_plot <- grid.arrange(prod_plot, het_plot, ncol=2)

return(obsdata)

}
```


### Estimating Activation Energies for Fitted Parameter Values

Here we estimate the activation energies for 3 different fitted parameters: **r**, **K**, and **a**. In general, the current fitting implementation likely underestimated **r** by a small amount, and also underestimated **K**, possibly by as much as an order of magnitude in some cases. This is due to the influence of the half-saturation constant (**b**) on the fitting process. You can think of **b** as being a kind of "pseudo-carrying capacity", because a higher **b** depresses the _effective_ attack rate, which then increases the _effective_ carrying capacity.

From the fitting algorithm's perspective, fitting a higher **b** is very similar to fitting a higher **K**, which is why I think it screwed up here. I have attempted to address this problem by coding in specific constraints as to how high **b** is allowed to be fit compared to **K**. The allowed gap between them was simply too narrow on the previous try.

Immediately below we can see the plots for the fitted **r** values. Note the large amount of alternatively high and low values for the two middle temperatures. This is because, for some of these replicates, the observed phytoplankton population dynamics appear to be periodic. Periodic dynamics in the producer are typically associated with higher **r** values, but the fitting algorithm was not able to successfully "figure out" that the dynamics were periodic in some of these cases, and subsequently fit very low values for **r**. We expect the fit to improve on our next few tries, by reducing the influence of **b**, as stated above.

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedr_2017_24_01.png" width="600">

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedK_2017_24_01.png" width="600">

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fitteda_2017_24_01.png" width="600">
