#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/Joey/Desktop/run-summaries

#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/Joey/Desktop/PK-TEMP-JUNE2016/JULY6-MORNING_processed/summaries", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

library(stringr)
?str_extract
copy_summaries <- str_extract(pattern = "copy", string = fnams)
copy_summaries


#### Step 3: create a df with the dataset ID and the cell count ####
July6_copies<- copy_summaries %>% 
	lapply(FUN = function(p) read.csv(p)) %>% View
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

