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


