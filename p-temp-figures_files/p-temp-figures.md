# P-TEMP_figures

```r
library(readr)
library(tidyr)
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
## 
## The following objects are masked from 'package:stats':
## 
##     filter, lag
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
library(purrr)
```

```
## 
## Attaching package: 'purrr'
## 
## The following object is masked from 'package:dplyr':
## 
##     order_by
```

```r
library(ggplot2)
```

```
## Warning: package 'ggplot2' was built under R version 3.2.4
```

```r
#### April 5 ####

ptemp_sep_29 <- read_csv("ptemp_summaries_March29.csv")
fnams <- list.files("/Users/Joey/Documents/p-temp/run-summaries-April5", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

#### Step 3: create a df with the dataset ID and the cell count ####
ptemp_summaries_april5 <- fnams %>%  
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	select(-V1)


ptemp_sep_april5 <- separate(ptemp_summaries_april5, dataset, c("UniqueID", "date", "try"), extra = "drop")


sep_april5 <- separate(ptemp_sep_april5, UniqueID, c("x", "Unique_ID"), sep = 1) %>% 
	select(-1)

Unique_ID_key <- read_csv("P-TEMP-UniqueID-key.csv")

Unique_ID_key <- separate(Unique_ID_key, TREATMENT_ID, c("treatment", "temperature", "replicate"))
sep_april5 <- sep_april5 %>% 
	rename(UniqueID = Unique_ID)
sep_april5$UniqueID <- as.integer(sep_april5$UniqueID)
	
april_5_cellcount <- left_join(sep_april5, Unique_ID_key, by = "UniqueID")

ggplot(april_5_cellcount, aes(x = temperature, y = cell_count, group = treatment, color = factor(treatment))) + geom_point()
```

![](p-temp-figures_files/figure-html/unnamed-chunk-1-1.png) 

```r
april_5_cellcount %>% 
	filter(replicate %in% c("1", "2", "3", "4", "5", "6")) %>% 
	ggplot(., aes(x = as.factor(temperature), y = cell_count, fill = factor(treatment), geom = "boxplot")) +
	geom_boxplot()
```

![](p-temp-figures_files/figure-html/unnamed-chunk-1-2.png) 

```r
ggsave("cellcount_t2.png")
```

```
## Saving 7 x 5 in image
```

```r
ptemp_sep_29_u <- unite(ptemp_sep_29, treatment_ID, treatment:replicate, remove = FALSE)

april_5_cellcount_u <- unite(april_5_cellcount, treatment_ID, treatment:replicate, remove = FALSE)

sep_29_u <- ptemp_sep_29_u %>% 
	select(treatment_ID, treatment, temperature, date, cell_count)

april_5_u <- april_5_cellcount_u %>% 
	select(treatment_ID, treatment, temperature, date, cell_count)

sep_29_u$temperature <- as.integer(sep_29_u$temperature)
april_5_u$temperature <- as.integer(april_5_u$temperature)


cell_counts <- bind_rows(sep_29_u, april_5_u)

cell_counts$date <- factor(cell_counts$date, c("MARCH24","MARCH29", "APRIL5"))
	
cell_counts$date[cell_counts$date == "MARCH24"] <- "MARCH29"


cell_counts %>% 
	# filter(temperature != 20) %>% 
ggplot(., aes(x=date, y=cell_count, group=treatment_ID, color=factor(temperature))) +
	geom_line(size=1) +
	geom_point(size=3)
```

![](p-temp-figures_files/figure-html/unnamed-chunk-1-3.png) 

```r
ggsave("cell_count_t1t2w-20.png")
```

```
## Saving 7 x 5 in image
```

```r
cell_counts %>% 
	filter(temperature %in% c("12", "16")) %>% 
	ggplot(., aes(x=date, y=cell_count, group=treatment_ID, color=factor(temperature))) +
	geom_line(size=1) +
	geom_point(size=3)
```

![](p-temp-figures_files/figure-html/unnamed-chunk-1-4.png) 

```r
cell_counts %>% 
	filter(temperature %in% c("12", "16")) %>% 
ggplot(., aes(x = date, y = cell_count, fill = factor(temperature), geom = "boxplot")) +
	geom_boxplot()
```

![](p-temp-figures_files/figure-html/unnamed-chunk-1-5.png) 

```r
ggsave("cellcounts_12.16.png")
```

```
## Saving 7 x 5 in image
```
