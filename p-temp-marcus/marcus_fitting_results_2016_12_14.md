# Fitting results

## Data

Here we read in the data used for all of the below plots:
```r
#file.path() is used for cross-platform compatibility
fittedpdata <- read.csv(file = file.path("p-temp-marcus", "fittedpdata6.csv"),
	strip.white = TRUE,
	na.strings = c("NA","") )
```
## Estimating Activation Energies for Fitted Parameter Values

Here we estimate the activation energies for 3 different fitted parameters: **r**, **K**, and **a**.

### Estimating r
```r
fittedr_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(r), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(r) Values") +
        labs(x = "-1/kT", y = "log(r)")
fittedr_plot
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedr_plot3.png" width="600">

#### Activation energy for r: Phosphorus Rich
```r
# Finding Ea for phosphorus rich treatment
r_fullP_model <- lm(log(r) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "FULL"))
summary(r_fullP_model)
confint(r_fullP_model)
```
**Ea**: 0.5350; 95% confidence intervals:

```r
confint(r_fullP_model)
                     2.5 %    97.5 %
(Intercept)     -7.4565840 48.288173
transformedtemp -0.1643051  1.234233
```

#### Activation energy for r: Phosphorus Poor
```r
# Finding Ea for phosphorus poor treatment
r_defP_model <- lm(log(r) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "DEF"))
summary(r_defP_model)
confint(r_defP_model)
```
**Ea**: 0.1434; 95% confidence intervals:

```r
confint(r_defP_model)
                      2.5 %    97.5 %
(Intercept)     -29.3847414 40.460383
transformedtemp  -0.7227137  1.029578

```
### Estimating K
```r
fittedK_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(K), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(K) Values") +
        labs(x = "-1/kT", y = "log(K)")
fittedK_plot
ggsave("fittedK_plot2.png", plot = last_plot())
```
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fittedK_plot3.png" width="600">

#### Activation energy for K: Phosphorus Rich
```r
# Finding Ea for phosphorus rich treatment
K_fullP_model <- lm(log(K) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "FULL"))
summary(K_fullP_model)
confint(K_fullP_model)
```

**Ea**: -0.2115; 95% confidence intervals:

```r
confint(K_fullP_model)
                    2.5 %     97.5 %
(Intercept)     -24.27721 42.4069370
transformedtemp  -1.04801  0.6249789
```
#### Activation energy for K: Phosphorus Poor
```r
# Finding Ea for phosphorus poor treatment
K_defP_model <- lm(log(K) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "DEF"))
summary(K_defP_model)
confint(K_defP_model)
```

**Ea**: -0.3097; 95% confidence intervals:

```r
confint(K_defP_model)
                    2.5 %     97.5 %
(Intercept)     -32.60623 42.2947284
transformedtemp  -1.24927  0.6298635
```

### Estimating a
```r
fitteda_plot <- ggplot(data = fittedpdata, aes(x = transformedtemp, y = log(a), color = Phosphorus)) +
        geom_point() +
        geom_smooth(method = lm, col = "red") +
        ggtitle("Fitted log(a) Values") +
        labs(x = "-1/kT", y = "log(a)")
fitteda_plot
ggsave("fitteda_plot2.png", plot = last_plot())
```

<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/fitteda_plot3.png" width="600">

#### Activation energy for a: Phosphorus Rich
**Ea**: -0.2606; 95% confidence intervals:

```r
confint(a_model)
                    2.5 %     97.5 %
(Intercept)     8.8368465 23.8860429
transformedtemp 0.1669927  0.5445506
```
#### Activation energy for a: Phosphorus Poor
**Ea**: -0.2606; 95% confidence intervals:
