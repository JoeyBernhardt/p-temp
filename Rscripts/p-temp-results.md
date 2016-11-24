# P-TEMP results




Load libraries and data

```r
library(tidyverse)
library(minpack.lm)
library(broom)
library(gridExtra)
library(lubridate)
library(plotrix)
library(stringr)

# load data ---------------------------------------------------------------

ptemp <- read_csv("/Users/Joey/Documents/p-temp/data-processed/p_temp_processed.csv")
ptemp_algae <- read_csv("/Users/Joey/Documents/p-temp/data-processed/p_temp_algae.csv") 
algae_summaries <- read_csv("/Users/Joey/Documents/p-temp/data-processed/algae_summaries.csv")
```


Phytoplankton populations over time



Phytoplankton populations over time (daphnia present)

```r
ptemp_algae %>% 
	filter(!grepl("C", replicate)) %>% 
	ggplot(data = ., aes(x = month_day, y = biovol, group = ID, color = factor(temp))) + geom_line(aes(linetype = P), size = 2) +
	facet_wrap( ~ temp) +
	scale_y_log10() +
	ylab("phytoplankton biovolume") +
	xlab("date") +
	ggtitle("Phytoplankton dynamics, with daphnia")
```

![](p-temp-results_files/figure-html/unnamed-chunk-3-1.png)<!-- -->

Phytoplankton populations over time (daphnia absent)

```r
ptemp_algae %>% 
	filter(grepl("C", replicate)) %>% 
	ggplot(data = ., aes(x = month_day, y = biovol, group = ID, color = factor(temp))) + geom_line(aes(linetype = P), size = 2) +
	facet_wrap( ~ temp) +
	scale_y_log10() +
	ylab("phytoplankton biovolume") +
	xlab("date") +
	ggtitle("Phytoplankton dynamics, without daphnia")
```

![](p-temp-results_files/figure-html/unnamed-chunk-4-1.png)<!-- -->

Phytoplankton maximum biovolume, over the entire experiment

```r
algae_summaries %>% 
	mutate(consumer = str_replace(consumer, "present", "daphnia present")) %>% 
		mutate(consumer = str_replace(consumer, "absent", "daphnia absent")) %>% 
	ggplot(data = ., aes(x = factor(temp), y = max, fill = factor(P))) + geom_boxplot() +
	# scale_y_log10() +
	facet_wrap( ~ consumer) +
	ylab("phytoplankton max biovolume")
```

![](p-temp-results_files/figure-html/unnamed-chunk-5-1.png)<!-- -->

Activation energies of phytoplankton densitiies, a la Schoolfield


```r
current_dataset <- algae_summaries %>% 
	filter(P == "FULL", consumer == "absent") %>%
	select(max, temp) %>% 
	mutate(K = temp + 273.15) %>% 
	rename(OriginalTraitValue = max) %>% 
	select(-temp)

current_dataset_def <- algae_summaries %>% 
	filter(P == "DEF", consumer == "absent") %>%
	select(max, temp) %>% 
	mutate(K = temp + 273.15) %>% 
	rename(OriginalTraitValue = max) %>% 
	select(-temp)
```


term      estimate      std.error    statistic     p.value  phosphorus 
-----  -----------  -------------  -----------  ----------  -----------
B0       18.700827   2.250786e-01   83.0857648   0.0000000  replete    
E         2.329546   6.108948e-01    3.8133334   0.0041319  replete    
E_D      31.523532   1.012739e+07    0.0000031   0.9999976  replete    
T_h     296.354656   1.715220e+05    0.0017278   0.9986591  replete    
B0       18.858266   4.731711e-01   39.8550650   0.0000000  deficient  
E         1.233073   1.284253e+00    0.9601481   0.3620638  deficient  
E_D      30.348379   2.945930e+07    0.0000010   0.9999992  deficient  
T_h     296.420100   4.542636e+05    0.0006525   0.9994936  deficient  

Plot the activation energies, daphnia absent


```r
all_estimates %>% 	
filter(term == "E") %>%
ggplot(aes(x = phosphorus, y = estimate, group = phosphorus)) + geom_point() +
	geom_errorbar(aes(ymin = estimate - std.error, ymax = estimate + std.error))
```

![](p-temp-results_files/figure-html/unnamed-chunk-8-1.png)<!-- -->

Generate predictions from the model fit (non bootstrapped) ---------------------------------

```r
 tmp_temps <- seq(min(
	floor(current_dataset$K)), 
	ceiling(max(current_dataset$K)
	), length = 200)

tmp_model <- exp(Schoolfield(
	coef(schoolfield_nls_full)["B0"],
	coef(schoolfield_nls_full)["E"],
	coef(schoolfield_nls_full)["E_D"],
	coef(schoolfield_nls_full)["T_h"],
	tmp_temps
))


ModelToPlotS <- data.frame(
	Temperature = tmp_temps - 273.15, 
	TraitValue = tmp_model
)

DataToPlot <- data.frame(
	Temperature = current_dataset$K - 273.15, 
	TraitValue = current_dataset$OriginalTraitValue
)
DataToPlot <- na.omit(DataToPlot)



##### DEF 

tmp_model_def <- exp(Schoolfield(
	coef(schoolfield_nls_def)["B0"],
	coef(schoolfield_nls_def)["E"],
	coef(schoolfield_nls_def)["E_D"],
	coef(schoolfield_nls_def)["T_h"],
	tmp_temps
))


ModelToPlotS_def <- data.frame(
	Temperature = tmp_temps - 273.15, 
	TraitValue = tmp_model_def
)

DataToPlot_def <- data.frame(
	Temperature = current_dataset_def$K - 273.15, 
	TraitValue = current_dataset_def$OriginalTraitValue
)
DataToPlot_def <- na.omit(DataToPlot_def)
```


plot them!

```r
full_plot <- ggplot() + geom_point(data = DataToPlot, aes(x = Temperature, 
																						 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
											alpha = 0.7, pch = 21) + 
	geom_line(data = ModelToPlotS, 
						aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
						lwd = 1.3) +                           
	ggtitle("Phosphorus replete") +
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("max phytoplankton abundance") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16))

#### DEF
def_plot <- ggplot() + geom_point(data = DataToPlot_def, aes(x = Temperature, 
																						 y = TraitValue), size = 3, col = "black", bg = "lightcyan2", 
											alpha = 0.7, pch = 21) + 
	geom_line(data = ModelToPlotS_def, 
						aes(x = Temperature, y = TraitValue), colour = "#1b9e77", 
						lwd = 1.3) +                           
	ggtitle("Phosphorus deficient") +
	xlab(expression(paste("Temperature (", degree, C, ")"))) + 
	ylab("max phytoplankton abundance") +
	theme_bw() + theme(plot.title = element_text(size = 16), 
										 axis.title = element_text(size = 16)) + scale_y_log10()

grid.arrange(full_plot, def_plot, ncol = 2)
```

![](p-temp-results_files/figure-html/unnamed-chunk-10-1.png)<!-- -->

