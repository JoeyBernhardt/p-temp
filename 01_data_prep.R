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

### March 30

fnams <- list.files("/Users/Joey/Documents/p-temp/run-summaries-March30", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

ptemp_summaries_March30 <- fnams %>%  
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	select(-V1)

View(ptemp_summaries_March30)

ptemp_sep_30 <- separate(ptemp_summaries_March30, dataset, c("treatment", "temperature", "replicate", "date", "try", "type"))

write_csv(ptemp_sep_30, "ptemp_sep_30.csv")
## edited out R1s where there was an R2

p_temp_march30 <- read_csv("ptemp_sep_30_edits.csv")

march30 <- p_temp_march30 %>% 
	select(cell_count, temperature, replicate, date, treatment)

march29 <- ptemp_sep %>% 
	select(cell_count, temperature, replicate, date, treatment)
	
all <- rbind(march29, march30)

all <- unite(all, uniqueID, treatment, temperature, replicate, sep = "_", remove = FALSE)

all$date[all$date == "MARCH24"] <- "MARCH29" 
?unite


ggplot(all, aes(x = date, y = cell_count, group = uniqueID, color = factor(all$treatment))) + geom_line() +ylim(0,20000)
