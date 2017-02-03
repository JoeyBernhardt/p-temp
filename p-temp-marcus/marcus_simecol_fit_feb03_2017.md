# Fitting Results

Here we used fitted **r** values which were estimated from the consumer free controls, in order to fit our data from the consumer-resource treatments. Below are the results of the best 3 fits ("best" is determined by lowest sum of squared deviations). We show both the 20C and 24C treatments here, as we had lower confidence in the **r** estimates from the 16C and 12C consumer-free controls. Basically, if this approach does not work well for 20C and 24C, I have a hard time imagining it working for the other treatments.

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
