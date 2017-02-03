# Fitting Results

Here we used fitted **r** values which were estimated from the consumer free controls. Below are the results of the best 3 fits ("best" is determined by lowest sum of squared deviations).

We did not pass on fitted **K** values from the consumer free controls here, as it was determined that **K** could not be properly estimated from the consumer-free controls. Because the observed data for the consumer-free controls appear to show that the phytoplankton are experiencing exponential growth during the entire course of the experiment, it is very difficult to estimate K with high precision. Attempting to do so resulted in obtaining fits with almost identical SSQs, but where the estimated value of **K** varied by up to 7 orders of magnitude.

### Code

The code used to generate the fits in seen below, as well as all related plotting functions, can be found at:
https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/p_temp_simecol_fit_with_controls.R

### Plotting Simulated Density Curves vs. Observed Data

##### 20 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/densityplot_def20_2017_FEB_03.png" width="900">
##### 20 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/densityplot_def24_2017_FEB_03.png" width="900">

##### 24 Degrees, Phosphorus Deficient
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/densityplot_full20_2017_FEB_03.png" width="900">
##### 24 Degrees, Phosphorus Rich
<img src="https://github.com/JoeyBernhardt/p-temp/blob/master/p-temp-marcus/plots/densityplot_full24_2017_FEB_03.png" width="900">
