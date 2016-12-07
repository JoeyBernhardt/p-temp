# Fitting results

### Data

Here we read in the data used for all of the below plots:
```r
#file.path() is used for cross-platform compatibility
plotdata <- read.csv(file = file.path("p-temp-marcus", "plotdata.csv"),
	strip.white = TRUE,
	na.strings = c("NA","") )
```

### Plotting Day 36 Densities for Predicted vs. Observed Data

Here we show plots of the simulation results vs. the experimental data taken on day 36 (the final observation).
The major point here is to assess the quality of our model fits. I have outputted some (possibly) relevant summary statistics for the linear regressions below each plot. Overall the fits look surprisingly good, especially for the Daphnia!

For the phytoplankton, there is a noticeable outlier, which I believe is an artefact of the fitting process. I'm going to take another look at this individual replicate. In general, some of the fits are off due to the fitting algorithm underestimating **K**. This should improve in the next fitting attempt (more details in the fitting section), so we should expect the new simulation results to do a better job of predicting the experimental data (and we'll see a higher R^2).

```r
# Plot the results of our model fitting.
	producer_plot <- ggplot(data = plotdata, aes(x = logpredictedfinalP, y = logobservedfinalP, color = Phosphorus)) +
		geom_point() + # predicted data
		geom_smooth(method = lm, col = "red") +
		labs(x = "log(Predicted phytoplankton density)", y = "log(observed Phytoplankton density)") +
		ggtitle("Phytoplankton Densities: Observed vs. Predicted")
	producer_plot
```

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/Phyto_predicted_vs_observedplot.png" width="600">

```r
P_model <- lm(logobservedfinalP ~ logpredictedfinalP, data = plotdata)
summary(P_model)
# adjusted R^2 of 0.4257, p-value of 3.04 * 10^-7
```
For the Daphnia, we can see some obvious outliers where the predicted values are much higher than the observed values; these are concentrated in the bottom left corner of the below plot.

```r
# Plot the results of our model fitting for daphnia
	hetero_plot <- ggplot(data = plotdata, aes(x = logpredictedfinalH, y = logobservedfinalH, color = Phosphorus)) +
		geom_point() + # predicted data
		geom_smooth(method = lm, col = "red") +
		labs(x = "log(Predicted Daphnia density)", y = "log(observed Daphnia density)") +
		ggtitle("Daphnia Densities: Observed vs. Predicted")
	hetero_plot
```

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/Daphnia_predicted_vs_observedplot.png" width="600">

```r
H_model <- lm(logobservedfinalH ~ logpredictedfinalH, data = plotdata)
summary(H_model)
# adjusted R^2 of 0.7706, p-value of 2 * 10^-16
```
### Estimating Activation Energies for Fitted Parameter Values

Here we estimate the activation energies for 3 different fitted parameters: **r**, **K**, and **a**. In general, the current fitting implementation likely underestimated **r** by a small amount, and also underestimated **K**, possibly by as much as an order of magnitude in some cases. This is due to the influence of the half-saturation constant (**b**) on the fitting process. You can think of **b** as being a kind of "pseudo-carrying capacity", because a higher **b** depresses the _effective_ attack rate, which then increases the _effective_ carrying capacity.

From the fitting algorithm's perspective, fitting a higher **b** is very similar to fitting a higher **K**, which is why I think it screwed up here. I have attempted to address this problem by coding in specific constraints as to how high **b** is allowed to be fit compared to **K**. The allowed gap between them was simply too narrow on the previous try.

Immediately below we can see the plots for the fitted **r** values. Note the large amount of alternatively high and low values for the two middle temperatures. This is because, for some of these replicates, the observed phytoplankton population dynamics appear to be periodic. Periodic dynamics in the producer are typically associated with higher **r** values, but the fitting algorithm was not able to successfully "figure out" that the dynamics were periodic in some of these cases, and subsequently fit very low values for **r**. We expect the fit to improve on our next few tries, by reducing the influence of **b**, as stated above.

```r
fittedr_plot <- ggplot(data = plotdata, aes(x = transformedtemp, y = log(r), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(r) Values") +
		labs(x = "-1/kT", y = "log(r)")
fittedr_plot
ggsave("fittedr_plot.png", plot = last_plot())

r_model <- lm(log(r) ~ transformedtemp, data = plotdata)
summary(r_model) #slope of 0.2375, 95% CI: (-0.17, 0.65), p-value = 0.253
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedr_plot.png" width="600">


Here we can see that **K** appears to be underestimated, at least for the replicates in the highest temperature treatment.
```r
fittedK_plot <- ggplot(data = plotdata, aes(x = transformedtemp, y = log(K), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(K) Values") +
		labs(x = "-1/kT", y = "log(K)")
fittedK_plot
ggsave("fittedK_plot.png", plot = last_plot())

K_model <- lm(log(K) ~ transformedtemp, data = plotdata)
summary(K_model) # slope of -0.1918, 95% CI: (-0.76, 0.38), p-value = 0.50
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedK_plot.png" width="600">


```r
fitteda_plot <- ggplot(data = plotdata, aes(x = transformedtemp, y = log(a), color = Phosphorus)) +
		geom_point() +
		geom_smooth(method = lm, col = "red") +
		ggtitle("Fitted log(a) Values") +
		labs(x = "-1/kT", y = "log(a)")
fitteda_plot
ggsave("fitteda_plot.png", plot = last_plot())

a_model <- lm(log(a) ~ transformedtemp, data = plotdata)
summary(a_model) # slope of 0.2092, 95% CI: (-0.93, 1.35), p-value = 0.71
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fitteda_plot.png" width="600">
