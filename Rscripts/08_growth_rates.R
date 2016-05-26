#-----------------------------------------------------------------------------------#
#multiplot function for better presentation of results
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
#-----------------------------------------------------------------------------------#

library(dplyr)
library(ggplot2)


p_temp <- read.csv("/Users/dominikbahlburg/p_temp_res/p-temp/p_temp_algae.csv")
p_temp_Daphnia <- read.csv("/Users/dominikbahlburg/p_temp_res/p-temp/p_temp_daphnia_algae.csv")

#
p_temp$date <- factor(p_temp$date, levels=c('MARCH29','APRIL1','APRIL5','APRIL8',
                                            'APRIL12','APRIL15','APRIL19',
                                            'APRIL26','MAY4'))
p_temp_Daphnia$date <- factor(p_temp_Daphnia$date, levels=c('APRIL1','APRIL5','APRIL8',
                                            'APRIL12','APRIL15','APRIL19',
                                            'APRIL26','MAY4'))


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


#Plot preparation
#Algae without controls
p_temp$temp <- as.factor(p_temp$temp)
sums_algae <- p_temp[,c(1,2,3,5,12)]
sums_algae <- sums_algae[sums_algae$ID<49,]

sums_algae <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = sums_algae[,c(2:5)],
                                        FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_algae <- na.omit(sums_algae)

#Algae controls
sums_algae_controls <- p_temp[,c(1,2,3,5,12)]
sums_algae_controls <- sums_algae_controls[sums_algae_controls$ID>48,]

sums_algae_controls <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = sums_algae[,c(2:5)],
																							FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_algae_controls <- na.omit(sums_algae_controls)

#Daphnia
p_temp_Daphnia$temp <- as.factor(p_temp_Daphnia$temp)
sums_daphnia <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = p_temp_Daphnia[,c(2,10,11,17)],
                                        FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_daphnia <- na.omit(sums_daphnia)

#plot algae without controls
plots_algae<-list()
for (k in 1:length(levels(sums_algae$temp))){
  p<-ggplot(data =sums_algae[sums_algae$temp==levels(sums_algae$temp)[k],], aes(date, y = growthrate.mn, group = P, color = P)) +
  geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
  geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('') + ggtitle(paste('Temperature:',levels(sums_algae$temp)[k],'ºC'))+
  geom_hline(aes(yintercept=0))+
  theme(title = element_text(size=20),
          legend.text = element_text(size = 22),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 20))
  plots_algae[[k]]<-p
}
png(filename = "growth_rate_algae.png",
    width = 1400, height = 1000, units = "px",
    pointsize = )
multiplot(plotlist=plots_algae,cols=2)
dev.off()

#plot algae controls
plots_algae_controls<-list()
for (k in 1:length(levels(sums_algae_controls$temp))){
	p<-ggplot(data =sums_algae_controls[sums_algae_controls$temp==levels(sums_algae_controls$temp)[k],], aes(date, y = growthrate.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('') + ggtitle(paste('Temperature:',levels(sums_algae$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_algae_controls[[k]]<-p
}
png(filename = "growth_rate_algae_ctrls.png",
		width = 1400, height = 1000, units = "px",
		pointsize = )
multiplot(plotlist=plots_algae_controls,cols=2)
dev.off()


#plot daphnia
plots_daphnia<-list()
for (k in 1:length(levels(sums_daphnia$temp))){
  p<-ggplot(data =sums_daphnia[sums_daphnia$temp==levels(sums_daphnia$temp)[k],], aes(date, y = growthrate.mn, group = P, color = P)) +
    geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
    geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('') + ggtitle(paste('Temperature:',levels(sums_daphnia$temp)[k],'ºC'))+
    geom_hline(aes(yintercept=0))+
    theme(title = element_text(size=20),
          legend.text = element_text(size = 22),
          axis.text.y = element_text(size = 16),
          axis.text.x = element_text(size = 16),
          axis.title.y = element_text(size = 20))
  plots_daphnia[[k]]<-p
}
png(filename = "growth_rate_daphnia.png",
    width = 1400, height = 1000, units = "px",
    pointsize = )
multiplot(plotlist=plots_daphnia,cols=2)
dev.off()
