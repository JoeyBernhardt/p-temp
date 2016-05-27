#plots of mean algae biovolume over time

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

p_temp_Daphnia <- read.csv("/Users/dominikbahlburg/p_temp_res/p-temp/p_temp_daphnia_algae.csv")

p_temp <- read.csv("/Users/dominikbahlburg/p_temp_res/p-temp/p_temp_algae.csv")
p_temp$date <- factor(p_temp$date, levels=c('MARCH29','APRIL1','APRIL5','APRIL8',
																						'APRIL12','APRIL15','APRIL19',
																						'APRIL26','MAY4'))
p_temp$temp <- as.factor(p_temp$temp)
biovol_mean <- p_temp[p_temp$ID<49,c(1,2,3,5,8)]
biovol_ctrl_mean <- p_temp[p_temp$ID>48,c(1,2,3,5,8)]
daphnia_abu <- p_temp_Daphnia[,c(2,10,11,6:9)]
daphnia_abu$temp <- as.factor(daphnia_abu$temp)
daphnia_abu$date <- factor(daphnia_abu$date, levels=c('APRIL1','APRIL5','APRIL8',
																						'APRIL12','APRIL15','APRIL19',
																						'APRIL26','MAY4'))
#for samples with Daphnia -> n=11 and n=12 for 12ºC&FULL/DEF on MARCH29 because those were flow cammed at March29 and March30 (which I merged as one time point)
biovol_mean <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = biovol_mean[,c(2:5)],
																							FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))

#for controls -> n=11 and n=12 for 12ºC&FULL/DEF on MARCH29 because those were flow cammed at March29 and March30 (which I merged as one time point)
biovol_ctrl_mean <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = biovol_ctrl_mean[,c(2:5)],
																							 FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
biovol_ctrl_mean<-na.omit(biovol_ctrl_mean)

#for Daphnia:
daphnia_mean <- as.data.frame(as.list(aggregate(. ~ date+temp+P,data = daphnia_abu[,c(1:7)],
																										FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))


#plot biovol with daphnia
plots_biovol<-list()
for (k in 1:length(levels(biovol_mean$temp))){
	p<-ggplot(data =biovol_mean[biovol_mean$temp==levels(biovol_mean$temp)[k],], aes(date, y = biovol.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=biovol.mn-biovol.se, ymax=biovol.mn+biovol.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("biovolume") + xlab('') + ggtitle(paste('Temperature:',levels(biovol_mean$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		ylim(0,max(biovol_mean$biovol.mn+biovol_mean$biovol.se)+200)+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_biovol[[k]]<-p
}
png(filename = "biovol_grazed_algae.png",
		width = 1800, height = 1200, units = "px",
		pointsize = )
multiplot(plotlist=plots_biovol,cols=2)
dev.off()

#plot controls
plots_ctrl_biovol<-list()
for (k in 1:length(levels(biovol_ctrl_mean$temp))){
	p<-ggplot(data =biovol_ctrl_mean[biovol_ctrl_mean$temp==levels(biovol_ctrl_mean$temp)[k],], aes(date, y = biovol.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=biovol.mn-biovol.se, ymax=biovol.mn+biovol.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("biovolume") + xlab('') + ggtitle(paste('Temperature:',levels(biovol_ctrl_mean$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		ylim(0,max(biovol_ctrl_mean$biovol.mn+biovol_ctrl_mean$biovol.se)+200)+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_ctrl_biovol[[k]]<-p
}
png(filename = "biovol_controls_algae.png",
		width = 1600, height = 1100, units = "px",
		pointsize = )
multiplot(plotlist=plots_ctrl_biovol,cols=2)
dev.off()




###

plots_daphnia<-list()
for (k in 1:length(levels(daphnia_mean$temp))){
	p<-ggplot(data =daphnia_mean[daphnia_mean$temp==levels(daphnia_mean$temp)[k],], aes(date, y = daphnia_ab.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=daphnia_ab.mn-daphnia_ab.se, ymax=daphnia_ab.mn+daphnia_ab.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("abundance") + xlab('') + ggtitle(paste('Temperature:',levels(biovol_ctrl_mean$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		ylim(0,max(daphnia_mean$daphnia_ab.mn+daphnia_mean$daphnia_ab.se))+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_daphnia[[k]]<-p
}
png(filename = "daphnia_abund.png",
		width = 1600, height = 1100, units = "px",
		pointsize = )
multiplot(plotlist=plots_daphnia,cols=2)
dev.off()
