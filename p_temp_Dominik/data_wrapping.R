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
  
#import metadata and Daphnia-count-data
  ID <- read.csv("/Users/dominikbahlburg/Documents/Studium/Bewerbungen/Oconnor/p_temp/UniqueID.csv", sep=';')  
  daphnia_count <- read.csv("~/Documents/Studium/Bewerbungen/Oconnor/p_temp/daphnia_count.csv")
  
  p_temp <- merge(p_temp_sums,ID,by="ID")
  daphnia_algae <- merge(p_temp_sums,daphnia_count,by=c("ID",'date'))
  daphnia_algae <- merge(daphnia_algae,ID,by=c("ID"))
  daphnia_algae$biovol <- as.numeric(as.character((daphnia_algae[,3]))) * as.numeric(as.character((daphnia_algae[,4])))
  
  p_temp <- p_temp[order(p_temp$date,as.numeric(p_temp$ID)),]
  p_temp <- p_temp[ , c("ID", "P", "temp", "replicate",'date','cells/ml','volume/cell')]
  p_temp$biovol <- as.numeric(as.character((p_temp[,6]))) * as.numeric(as.character((p_temp[,7])))
  colnames(p_temp)[6]<-'cellcon'
  colnames(p_temp)[7]<-'vol_cell'
  
  #I'm cheating here!!! It's only to get a nicer plot
  #we could consider making a column with something like "sampling No"
  #that could group data of very similar but not equal dates
  p_temp$date[p_temp$date=='APRIL2'] <- 'APRIL1'
  p_temp$date[p_temp$date=='MARCH30'] <- 'MARCH29'
  
  p_temp$date <- as.factor(p_temp$date)
  row.names(p_temp) <- 1:nrow(p_temp)
  
  
  p_temp$date <- factor(p_temp$date, levels=c('MARCH29','APRIL1','APRIL5','APRIL8',
  																						'APRIL12','APRIL15','APRIL19',
  																						'APRIL26','MAY4'))
  
  colnames(daphnia_algae)[6]<-'daphnia_ab' 
  daphnia_algae$date <- factor(daphnia_algae$date, levels=c('APRIL1','APRIL5','APRIL8',
  																													'APRIL12','APRIL15','APRIL19',
  																													'APRIL26','MAY4'))
  
  
  #give data right format
  p_temp$cellcon <- as.numeric(as.character(p_temp$cellcon))
  p_temp$vol_cell <- as.numeric(as.character(p_temp$vol_cell))
  p_temp$ID <- as.numeric(p_temp$ID)
  
  
  p_temp %>% 
  	#filter(temp %in% c('20') | ID <49) %>% 
  	filter(ID <49) %>% 
  	ggplot(., aes(x = date, y = log(cellcon), fill = factor(temp), geom = "boxplot")) +
  	geom_boxplot()
  
  
  daphnia_algae %>% 
  	filter(temp %in% c('20')) %>% 
  	ggplot(., aes(x = date, y = daphnia_ab, fill=factor(P),geom = "boxplot")) +
  	geom_boxplot()
 
  
  write_csv(p_temp, "p_temp_algae.csv")
  write_csv(daphnia_algae, "p_temp_daphnia_algae.csv")
  