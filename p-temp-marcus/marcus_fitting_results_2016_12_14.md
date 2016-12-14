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
Fitted activation energy for **r**: 0.1474; 95% confidence intervals:

```r
confint(r_model)
                      2.5 %     97.5 %
(Intercept)     -20.8155358 30.5050648
transformedtemp  -0.4963514  0.7911925
```

```r
fittedK_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(K), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(K) Values") +
        labs(x = "-1/kT", y = "log(K)")
fittedK_plot
ggsave("fittedK_plot2.png", plot = last_plot())

K_model <- lm(log(K) ~ transformedtemp, data = fittedpdata)
summary(K_model)
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedK_plot2.png" width="600">

```r
fitteda_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(a), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(a) Values") +
        labs(x = "-1/kT", y = "log(a)")
fitteda_plot
ggsave("fitteda_plot2.png", plot = last_plot())

a_model <- lm(log(a) ~ transformedtemp, data = fittedpdata)
summary(a_model)
```

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fitteda_plot2.png" width="600">
