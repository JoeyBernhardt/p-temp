library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)

## cp **/*summary.csv /Users/Joey/Desktop/run-summaries

fnams <- list.files("/Users/Joey/Documents/p-temp/run-summaries", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

ptemp_summaries <- fnams %>%  
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	select(-V1)

View(ptemp_summaries)

write_csv(ptemp_summaries, "ptemp_summaries.csv")

ptemp_sep <- separate(ptemp_summaries, dataset, c("treatment", "temperature", "replicate", "date"))

View(ptemp_sep)

ggplot(ptemp_sep, aes(y = cell_count, x = temperature, color = treatment)) + geom_point()

