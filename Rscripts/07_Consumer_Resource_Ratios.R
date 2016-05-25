### consumer resource ratios
### May 25 2016

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(plotrix)
library(MuMIn)

df <- read_csv("p_temp_Dominik/p_temp_daphnia_algae.csv")

df <- df %>% 
	mutate(CRR = dapnia_count/biovol) %>% 
	rename(UniqueID = ID) %>%
	mutate(UniqueID = as.factor(UniqueID))
	

df <- left_join(df, UniqueID_treat, by = "UniqueID")

df %>% 
	filter(date == "MAY4") %>% 
	filter(CRR < 0.003) %>% 
ggplot(., aes(x = temperature, y = CRR)) + geom_point() +
	ylim(0,0.0001)

df %>% 
	filter(date == "MAY4") %>% 
	filter(CRR < 0.003) %>% 
	dplyr::select(CRR, treatment, temperature) %>% 
	group_by(temperature, treatment) %>%
	summarise_each(funs(mean, median, sd, std.error)) %>%
	ggplot(data = ., aes(temperature, y = mean, group = treatment, color = treatment)) +
	geom_errorbar(aes(ymin=mean-std.error, ymax=mean+std.error), width=.2) +
	geom_point(size = 6) + ylab("CRR") + xlab("temperature, C") +
	theme_bw() +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
