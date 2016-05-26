library(readxl)
library(dplyr)
library(purrr)
library(readr)
library(tidyr)
library(ggplot2)
library(plotrix)


x <- read_csv("daph_size_csv/268-4-size-April19.csv")

fnams <- list.files("daph_size_data/May3_daph_size", full.names = TRUE) ## find out the names of all the files in data-raw, use full.names to get the relative path for each file

## Step 2: Name the the vector with filenames (i.e. fnams) with the basename (i.e file name without the directory) of the data files.
## then get rid of the ".csv" part
names(fnams) <- fnams %>% 
	basename %>% # basename removes all of the path up to and including the last path separator
	tools::file_path_sans_ext(.) ## remove .csv part of the file name

## Step 3. Read in every file, using the relative path
## keep the filename in a nice column called "dataset".
daph_size <- fnams %>% 
	map_df(read.csv, .id = "dataset") ## use the map function to read in the csvs listed in fnams
## .id argument allows us to create a variable giving the name of the dataset

daph_size_sep <- daph_size %>% 
	separate(dataset, c("photo", "UniqueID", "date", "type")) %>% 
	select(-X)

Unique_ID_key <- read_csv("data-raw/P-TEMP-UniqueID-key.csv")
Unique_ID_key <- separate(Unique_ID_key, TREATMENT_ID, c("treatment", "temperature", "replicate"))
Unique_ID_key$UniqueID <- as.factor(Unique_ID_key$UniqueID)
daph_size_sep$UniqueID <- as.factor(daph_size_sep$UniqueID)

Daph <- left_join(daph_size_sep, Unique_ID_key, by = "UniqueID")

Daph.size <- Daph %>% 
	group_by(UniqueID, temperature, date, treatment) %>% 
	summarise(mean.length = mean(Length)) 

ggplot(data = Daph.size, aes(x = factor(temperature), y = mean.length)) + geom_boxplot()


Daph %>% 
	dplyr::select(Length, treatment, temperature) %>%
	filter(temperature != "12") %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean,median, sd,std.error)) %>% 
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + theme_bw() + ylab("body length, cm") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/body_size.png")


Daph %>% 
	filter(temperature != "12") %>% 
	mutate(temperature = as.numeric(temperature)) %>% 
	# filter(treatment == "DEF") %>% 
	lm(Length ~ temperature + treatment, data = .) %>% 
	summary()



### top quartile

Daph25 <- Daph %>% 
	group_by(UniqueID, temperature, date, treatment) %>% 
	summarise(top = quantile(Length, probs = 0.75),
						avg = mean(Length),
						n = n())
Daph25 <- as.data.frame(Daph25)

Daph25 %>% 
	dplyr::select(top, treatment, temperature) %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean,median, sd,std.error)) %>%
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + theme_bw() + ylab("ln(body length)") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) +
	scale_y_continuous(trans = "log")


Daph25 %>% 
	ggplot(data = ., aes(x = factor(temperature), y = top)) + geom_boxplot()

mod <- lm(top ~ temperature + treatment, data = Daph25)
summary(mod)


### top 10
Daphtop10 <- Daph %>% 
	group_by(UniqueID, temperature, date, treatment) %>% 
	# filter(temperature != "20") %>% 
	arrange(desc(Length)) %>% 
	top_n(., 10, wt = Length) %>% 
	summarise(top10 = mean(Length),
						n = n()) %>% 
	ggplot(data = ., aes(x = factor(temperature), y = top10)) + geom_boxplot()


mod <- lm(top2 ~ temperature, data = Daphtop10)
summary(mod)


Daph.size %>% 
	# filter(UniqueID != "28") %>% 
	filter(date == "April12") %>% 
ggplot(data = ., aes(x = factor(temperature), y = mean.length)) + geom_boxplot()
ggsave("april_19_daphsize.png")
	
Daph.size %>% 
	# filter(UniqueID != "28") %>% 
	filter(date == "April12") %>% 
	ggplot(data = ., aes(x = factor(temperature), y = mean.length, group = treatment, color = treatment)) + geom_point()

Daph.max <- Daph %>% 
	group_by(UniqueID, temperature) %>% 
	summarise(max.length = max(Length)) 

Daph.max %>% 
	# filter(UniqueID != "28") %>% 
	ggplot(data = ., aes(x = factor(temperature), y = max.length)) + geom_boxplot()

mod <- lm(max.length ~ temperature, data = Daph.max)
summary(mod)


