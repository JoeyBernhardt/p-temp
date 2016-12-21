

### bringing in march 29 flow cam cell data
### last updated by JB Dec 21 2016

cell_files <- list.files("/Users/Joey/Documents/p-temp/data-processed/run-summaries-March29", full.names = TRUE)  ## find out the names of all the files in data-summary, use full.names to get the relative path for each file


names(cell_files) <- cell_files %>% 
	gsub(pattern = ".csv$", replacement = "")


#### Step 3: read in all the files!

march29_cells <- map_df(cell_files, read_csv, col_names = FALSE,
										.id = "file_name")
march29_cell_data <- march29_cells %>% 
	rename(obs_type = X1,
				 value = X2) %>% 
	filter(obs_type %in% c("List File", "Start Time", "Particles / ml", "Volume (ABD)")) %>%
	spread(obs_type, value) %>% 
	separate(`List File`, into = c("nutrient_level", "temperature", "replicate", "date"), sep = "-") %>%
	rename(start_time = `Start Time`,
				 cell_density = `Particles / ml`,
				 cell_volume = `Volume (ABD)`)

write_csv(march29_cell_data, "data-processed/march29_cell_data.csv")
