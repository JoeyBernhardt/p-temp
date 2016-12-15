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
### Interpreting the estimates for r

Currently, the fitting function is still struggling to fit **r** properly. Visually, we can think of **r** as primarily responsible for the _initial_ slope of the growth curve for phytoplankton. A higher **r** will lead to faster population growth, and thus a steeper, more positive slope. From looking at the shape of the plots of experimental data, and then comparing them to the analytical model using various parameter settings, own guess is that the mean _true_ value of **r** falls somewhere between 0.2 and 0.8 for our experimental data. I hope to see fits that fall within this range in a few days, if all goes well.

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
### Interpreting the estimates for K

For both the low and high phosphorus treatments we see what appears to be a non-linear trend, with the fitted carrying capacity **K** steadily increasing with temperature until around 20 degrees C, and then dropping. When fitting a linear regression, this pattern manifests itself as a negative slope (and thus a negative activation energy), but given the apparent non-linear trend, linear regression is not appropriate.

I also recall seeing a similar pattern in the experimental data itself, when we plot the maximum observed abundances for phytoplankton. Given the perceived shape of the curve, attempting non-linear regression with a function such as the schoolfield model seems to be the best course of action.

The mean K's do appear to be larger for the full-phosphorus treatments, and this can also be seen by calling `t.test(data = fittedpdata, K~Phosphorus)`, however the difference is **not significant** in this case (for now). _As an aside, it is worth noting that we can also call a t-test on the log-transformed data; but this is in some sense a very different (but still valid) kind of test. In this case the t-test would be comparing the **geometric** means, which do have meaning when dealing with dynamical systems that can exhibit periodic cycling, such as ours._

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
```r
a_fullP_model <- lm(log(a) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "FULL"))
summary(a_fullP_model)
confint(a_fullP_model)
```

**Ea**: 0.2862; 95% confidence intervals:

```r
confint(a_fullP_model)
                      2.5 %     97.5 %
(Intercept)      1.78058145 25.4574317
transformedtemp -0.01080157  0.5832091
```
#### Activation energy for a: Phosphorus Poor
```r
# Finding Ea for phosphorus poor treatment
a_defP_model <- lm(log(a) ~ transformedtemp, data = filter(fittedpdata, Phosphorus == "DEF"))
summary(a_defP_model)
confint(a_defP_model)
```

**Ea**: 0.4253; 95% confidence intervals:

```r
confint(a_defP_model)
                    2.5 %     97.5 %
(Intercept)     8.7756074 29.4321582
transformedtemp 0.1662212  0.6844579
```

### Interpreting the estimates for a
