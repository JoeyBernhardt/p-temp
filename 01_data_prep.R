library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)


#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Desktop/run-summaries

#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/Joey/Documents/p-temp/run-summaries-March29", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

#### Step 3: create a df with the dataset ID and the cell count ####
ptemp_summaries_march29 <- fnams %>%  
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

ptemp_sep_29 <- separate(ptemp_summaries_march29, dataset, c("treatment", "temperature", "replicate", "date"), extra = "drop")


#### Step 4: write out the df to a csv ####
write_csv(ptemp_sep_29, "ptemp_summaries_March29.csv")



ggplot(ptemp_sep_29, aes(y = cell_count, x = temperature, color = treatment)) + geom_point()

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



ggplot(all, aes(x = date, y = cell_count, group = uniqueID, color = factor(all$treatment))) + geom_line() +ylim(0,20000)


### April 2

fnams <- list.files("/Users/Joey/Documents/p-temp/run-summaries-April2", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

ptemp_summaries_April2 <- fnams %>%  
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	select(-V1)

ptemp_sep_April2 <- separate(ptemp_summaries_April2, dataset, c("UniqueID", "treatment", "temperature", "replicate", "date", "try", "type"))

ggplot(ptemp_sep_April2, aes(x = temperature, y = cell_count, group = treatment, color = factor(treatment))) + geom_point()


#### April 5 ####


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

View(ptemp_summaries_april5)

ptemp_sep_april5 <- separate(ptemp_summaries_april5, dataset, c("UniqueID", "date", "try"), extra = "drop")

?separate

sep_april5 <- separate(ptemp_sep_april5, UniqueID, c("x", "Unique_ID"), sep = 1) %>% 
	select(-1)

Unique_ID_key <- read_csv("P-TEMP-UniqueID-key.csv")

Unique_ID_key <- separate(Unique_ID_key, TREATMENT_ID, c("treatment", "temperature", "replicate"))
sep_april5 <- sep_april5 %>% 
	rename(UniqueID = Unique_ID)
sep_april5$UniqueID <- as.integer(sep_april5$UniqueID)
	
april_5_cellcount <- left_join(sep_april5, Unique_ID_key, by = "UniqueID")

ggplot(april_5_cellcount, aes(x = temperature, y = cell_count, group = treatment, color = factor(treatment))) + geom_point()

april_5_cellcount %>% 
	filter(replicate %in% c("1", "2", "3", "4", "5", "6")) %>% 
	ggplot(., aes(x = as.factor(temperature), y = cell_count, fill = factor(treatment), geom = "boxplot")) +
	geom_boxplot()
ggsave("cellcount_t2.png")

ptemp_sep_29_u <- unite(ptemp_sep_29, treatment_ID, treatment:replicate, remove = FALSE)

april_5_cellcount_u <- unite(april_5_cellcount, treatment_ID, treatment:replicate, remove = FALSE)

sep_29_u <- ptemp_sep_29_u %>% 
	select(treatment_ID, treatment, temperature, date, cell_count)

april_5_u <- april_5_cellcount_u %>% 
	select(treatment_ID, treatment, temperature, date, cell_count)

cell_counts <- bind_rows(sep_29_u, april_5_u)


cell_counts %>% 
	group_by(treatment_ID) %>% 
	
	
	cell_counts$date <- factor(cell_counts$date, c("MARCH24","MARCH29", "APRIL5"))
	
cell_counts$date[cell_counts$date == "MARCH24"] <- "MARCH29"


cell_counts %>% 
	# filter(temperature != 20) %>% 
ggplot(., aes(x=date, y=cell_count, group=treatment_ID, color=temperature)) +
	geom_line(size=1) +
	geom_point(size=3)
ggsave("cell_count_t1t2w-20.png")


cell_counts %>% 
	filter(temperature %in% c("12", "16")) %>% 
	ggplot(., aes(x=date, y=cell_count, group=treatment_ID, color=temperature)) +
	geom_line(size=1) +
	geom_point(size=3)


cell_counts %>% 
	filter(temperature %in% c("12", "16")) %>% 
ggplot(., aes(x = date, y = cell_count, fill = factor(temperature), geom = "boxplot")) +
	geom_boxplot()
ggsave("cellcounts_12.16.png")
	
	
	
	


