library(dplyr)
library(ggplot2)

p_temp <- read.csv("p_temp_algae.csv")
p_temp_Daphnia <- read.csv("p_temp_daphnia_algae.csv")

#calculate growth rates
datevec<-c('MARCH29','APRIL1','APRIL5','APRIL8',
					 'APRIL12','APRIL15','APRIL19',
					 'APRIL26','MAY4')

realdatevec<-c('2016-03-29','2016-04-01','2016-04-05','2016-04-08',
					 '2016-04-12','2016-04-15','2016-04-19',
					 '2016-04-26','2016-05-04')

p_temp$date_edit <- NA
p_temp_Daphnia$date_edit <- NA
for (i in 1:length(datevec)){
	p_temp$date_edit[p_temp$date==datevec[i]] <- as.Date(realdatevec[i])
	p_temp_Daphnia$date_edit[p_temp_Daphnia$date==datevec[i]] <- as.Date(realdatevec[i])
}

#sorting data
p_temp<-p_temp %>%
	group_by(ID) %>%
	arrange(date) 

p_temp_Daphnia<-p_temp_Daphnia %>%
	group_by(ID) %>%
	arrange(date) 

#calculate growth rate (difference between recent and previous timepoint/number of days between sampling)
#Biovol alge
p_temp$voldiff <- ave(p_temp$biovol, p_temp$ID, FUN=function(x) c(0, diff(x)))
p_temp$daydiff <- ave(p_temp$date_edit, p_temp$ID, FUN=function(x) c(0, diff(x)))
p_temp$growthrate <- p_temp$voldiff/p_temp$daydiff

#abundance Daphnia
p_temp_Daphnia$abu_diff <- ave(p_temp_Daphnia$daphnia_ab, p_temp_Daphnia$ID, FUN=function(x) c(0, diff(x)))
p_temp_Daphnia$daydiff <- ave(p_temp_Daphnia$date_edit, p_temp_Daphnia$ID, FUN=function(x) c(0, diff(x)))
p_temp_Daphnia$growthrate <- p_temp_Daphnia$abu_diff/p_temp_Daphnia$daydiff

p_temp %>% 
	#filter(temp %in% c('20') | ID <49) %>% 
	filter(ID <49) %>% 
	ggplot(., aes(x = date, y = growthrate, fill = factor(temp), geom = "boxplot")) +
	geom_boxplot()

p_temp_Daphnia %>% 
	#filter(temp %in% c('20') | ID <49) %>% 
	filter(date != 'APRIL1') %>% 
	ggplot(., aes(x = date, y = growthrate, fill = factor(temp), geom = "boxplot")) +
	geom_boxplot()

