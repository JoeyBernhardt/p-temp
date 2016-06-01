library(readr)
library(tidyr)
library(dplyr)
library(purrr)
library(ggplot2)
library(data.table)

#### Step 1 #### 

## in the shell, use this command to copy the summary files from the folder of flow cam run folders to a summary-only file
## cp **/*summary.csv /Users/dominikbahlburg/Desktop/Comp/summary_MAY10

#### Step 2: create a list of file names for each of the summaries ####
fnams <- list.files("/Users/dominikbahlburg/Documents/Studium/Bewerbungen/Oconnor/p_temp/summaries", full.names = TRUE) ## find out the names of all the files in data-summary, use full.names to get the relative path for each file

#### Step 3: create data frame for all available data!
p_temp_sums <- fnams %>%  
  lapply(FUN = function(p) read.csv(p)) 

  p_temp_sums<-rbindlist(p_temp_sums)
  p_temp_sums<-cbind(p_temp_sums[p_temp_sums$List.File=="Particles / ml",],p_temp_sums[p_temp_sums$List.File=="Volume (ABD)",])

  #this step is not very elegant so far... substr... 79 - to extract relevant data
  #I have to substring from element 79 of the files' path. the length of the path
  #changes from computer to computer...
  fnams2<-sapply(fnams, function(x) substr(x, 79, nchar(x)))
  p_temp_sums<-cbind(p_temp_sums,fnams2)
  
  p_temp_sums <- separate(p_temp_sums, fnams2, c("ID", "date"), extra = "drop")
  colnames(p_temp_sums)<-c('','cells/ml','','volume/cell','ID','date')
  p_temp_sums <- as.data.frame(p_temp_sums)
  p_temp_sums <- p_temp_sums[,c(2,4:6)]
  
  
  #I'm cheating here!!! It's only to get a nicer plot
  #we could consider making a column with something like "sampling No"
  #that could group data of very similar but not equal dates
  p_temp_sums$date[p_temp_sums$date=='APRIL2'] <- 'APRIL1'
  p_temp_sums$date[p_temp_sums$date=='MARCH30'] <- 'MARCH29'
 
#add number of days experiment is running + format date   
  p_temp_sums$date <- as.factor(p_temp_sums$date)
  row.names(p_temp_sums) <- 1:nrow(p_temp_sums)
  
  
  p_temp_sums$date <- factor(p_temp_sums$date, levels=c('MARCH29','APRIL1','APRIL5','APRIL8',
  																						'APRIL12','APRIL15','APRIL19',
  																						'APRIL26','MAY4'))
  
  
  datevec<-c('MARCH29','APRIL1','APRIL5','APRIL8',
  					 'APRIL12','APRIL15','APRIL19',
  					 'APRIL26','MAY4')
  
  realdatevec<-c('2016-03-29','2016-04-01','2016-04-05','2016-04-08',
  							 '2016-04-12','2016-04-15','2016-04-19',
  							 '2016-04-26','2016-05-04')
  
  p_temp_sums$date_edit <- NA
  for (i in 1:length(datevec)){
  	p_temp_sums$date_edit[p_temp_sums$date==datevec[i]] <- as.Date(realdatevec[i])
  }
  
  p_temp_sums$ID <- as.numeric(as.character(as.factor(p_temp_sums$ID)))
  p_temp_sums<-p_temp_sums %>%
  	group_by(ID) %>%
  	arrange(date) 
  p_temp_sums$daydiff <- ave(p_temp_sums$date_edit, p_temp_sums$ID, FUN=function(x) c(0, diff(x)))
  days <- aggregate(daydiff ~ ID, data = p_temp_sums, cumsum)
  p_temp_sums<-p_temp_sums[order(p_temp_sums$ID),]
  p_temp_sums$days<-unlist(days$daydiff)
  
#import metadata and Daphnia-count-data
  ID <- read.csv("/Users/dominikbahlburg/Documents/Studium/Bewerbungen/Oconnor/p_temp/UniqueID.csv", sep=';')  
  daphnia_count <- read.csv("~/Documents/Studium/Bewerbungen/Oconnor/p_temp/daphnia_count.csv")
  
  p_temp <- merge(p_temp_sums,ID,by="ID")
  daphnia_algae <- merge(p_temp_sums,daphnia_count,by=c("ID",'date'))
  daphnia_algae <- merge(daphnia_algae,ID,by=c("ID"))
  daphnia_algae$biovol <- as.numeric(as.character((daphnia_algae[,3]))) * as.numeric(as.character((daphnia_algae[,4])))
  
  p_temp <- p_temp[order(p_temp$date,as.numeric(p_temp$ID)),]
  p_temp <- p_temp[ , c("ID", "P", "temp", "replicate",'date','daydiff','days','cells/ml','volume/cell')]
  p_temp$biovol <- as.numeric(as.character((p_temp[,8]))) * as.numeric(as.character((p_temp[,9])))
  colnames(p_temp)[8]<-'cellcon'
  colnames(p_temp)[9]<-'vol_cell'
  
  #give data right format
  p_temp$cellcon <- as.numeric(as.character(p_temp$cellcon))
  p_temp$vol_cell <- as.numeric(as.character(p_temp$vol_cell))
  p_temp$ID <- as.numeric(p_temp$ID)
  
  
  p_temp %>% 
  	#filter(temp %in% c('20') | ID <49) %>% 
  	filter(ID <49) %>% 
  	ggplot(., aes(x = days, y = log(cellcon), fill = factor(temp), geom = "boxplot")) +
  	geom_boxplot()
  
  
  daphnia_algae %>% 
  	filter(temp %in% c('20')) %>% 
  	ggplot(., aes(x = date, y = daphnia_ab, fill=factor(P),geom = "boxplot")) +
  	geom_boxplot()
 
  
  write_csv(p_temp, "p_temp_algae.csv")
  write_csv(daphnia_algae, "p_temp_daphnia_algae.csv")
  