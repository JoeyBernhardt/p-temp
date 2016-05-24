library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)


#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Desktop/run-summaries

#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/Joey/Documents/P-TEMP-FLOWCAM/P-TEMP-MAY4/csvtotal", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

#### Step 3: create a df with the dataset ID and the cell count ####
ptemp_summaries_may4 <- fnams %>%  
	lapply(FUN = function(p) read.csv(p)) %>% 
	as.data.frame(.) %>% 
	filter(List.File == "Particles / ml") %>% 
	select(- starts_with("List")) %>%
	t(.) %>%
	as.data.frame() %>%
	mutate(dataset = rownames(.)) %>%
	mutate(cell_count = as.numeric(as.character(V1))) %>%
	select(-V1)

ptemp_summaries_may4 <- ptemp_summaries_may4 %>% 
	filter(dataset != "X1.MAY4.R1_copy_copy.lst") %>% 
	filter(dataset != "X39.MAY4.R1.lst")


	

ptemp_sep_may4 <- separate(ptemp_summaries_may4, dataset, c("replicate", "date", "copy"), extra = "drop")
ptemp_sep_may4 <- separate(ptemp_sep_may4, replicate, c("x", "Unique_ID"), sep = 1)
ptemp_sep_may4 <- ptemp_sep_may4 %>% 
	select(-1)



#### Step 4: write out the df to a csv ####
write_csv(ptemp_sep_may4, "data-processed/ptemp_summaries_may4.csv")


Unique_ID_key <- read_csv("data-raw/P-TEMP-UniqueID-key.csv")

Unique_ID_key <- separate(Unique_ID_key, TREATMENT_ID, c("treatment", "temperature", "replicate"))
Unique_ID_key <- rename(Unique_ID_key, Unique_ID = UniqueID)
Unique_ID_key$Unique_ID <- as.character(Unique_ID_key$Unique_ID)

str(Unique_ID_key)
may4_cells <- left_join(ptemp_sep_may4, Unique_ID_key, by = "Unique_ID")

may4_cells %>% 
	filter(Unique_ID < 49) %>% 
	dplyr::select(cell_count, treatment, temperature) %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean,median, sd,std.error)) %>% 
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	ylim(0,250000) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) + 
	geom_point(size = 2) + theme_bw() + ylab("cell_count") + xlab("temperature, C")
