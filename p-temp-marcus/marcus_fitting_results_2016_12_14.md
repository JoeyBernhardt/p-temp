# Fitting results

### Data

Here we read in the data used for all of the below plots:
```r
#file.path() is used for cross-platform compatibility
fittedpdata <- read.csv(file = file.path("p-temp-marcus", "fittedpdata6.csv"),
	strip.white = TRUE,
	na.strings = c("NA","") )
```
### Estimating Activation Energies for Fitted Parameter Values

Here we estimate the activation energies for 3 different fitted parameters: **r**, **K**, and **a**.

```r
fittedr_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(r), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(r) Values") +
        labs(x = "-1/kT", y = "log(r)")
fittedr_plot

r_model <- lm(log(r) ~ transformedtemp, data = fittedpdata)
summary(r_model) 
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedr_plot2.png" width="600">


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

Some of the fitted **a's** are almost certainly off; this is due to the influence of the transfer efficiency **e**. For some replicates a very, very high **a** was fit (one order of magnitude larger than in the other replicates), but the fitted **e** for the same replicate was an order of magnitude smaller than what was produced by the other fittings. This is just a matter of tightening up the parameter constraints for the fitting. Improvements have already been made to the code to address this issue for the next step.


<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fitteda_plot.png" width="600">
