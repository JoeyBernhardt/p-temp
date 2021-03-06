# plots_fun



```r
# load libraries ----------------------------------------------------------
library(tidyverse)
library(gridExtra)
library(broom)
```



```r
# read in data ------------------------------------------------------------

all_times <- read_csv("CR_abundances_30days.csv")
```


```r
all_times %>% 
	filter(temperature %in% c("12", "16", "20", "24")) %>% 
	ggplot(data = ., aes(x = time, y = P, color = factor(temperature))) +geom_point() +
	facet_wrap( ~ resource_level) +
	ggtitle("resource density") +
	scale_y_log10() +
	theme_minimal()
```

![](ea_plots_files/figure-html/unnamed-chunk-3-1.png)<!-- -->



```r
all_times %>% 
	filter(temperature %in% c("12", "16", "20", "24")) %>% 
	ggplot(data = ., aes(x = time, y = H, color = factor(temperature))) +geom_point() +
	facet_wrap( ~ resource_level) +
	ggtitle("consumer density") +
	# scale_y_log10() +
	theme_minimal()
```

![](ea_plots_files/figure-html/unnamed-chunk-4-1.png)<!-- -->


```r
all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = H, group = resource_level, color = resource_level)) + geom_line(size = 3) +
	scale_x_reverse() + 
	scale_y_log10() +
	theme_minimal() + xlab("temperature (1/kT)") + ylab("log(consumer abundance)") + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 18)) +
	theme(legend.position="top")
```

![](ea_plots_files/figure-html/unnamed-chunk-5-1.png)<!-- -->


```r
all_times %>% 
	filter(time == 30) %>% 
	mutate(inverse_temp = (1/(.00008617*(temperature+273.15)))) %>%
	ggplot(aes(x = inverse_temp, y = P, group = resource_level, color = resource_level)) + geom_line(size = 3) +
	scale_x_reverse() + 
	scale_y_log10() +
	theme_minimal() + xlab("temperature (1/kT)") + ylab("log(producer abundance)") + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	theme(legend.title=element_blank(),
				legend.text = element_text(size = 18)) +
	theme(legend.position="top")
```

![](ea_plots_files/figure-html/unnamed-chunk-6-1.png)<!-- -->
