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


p_temp <- read.csv("p_temp_algae.csv")
p_temp_Daphnia <- read.csv("p_temp_daphnia_algae.csv")

#
p_temp$date <- factor(p_temp$date, levels=c('MARCH29','APRIL1','APRIL5','APRIL8',
																						'APRIL12','APRIL15','APRIL19',
																						'APRIL26','MAY4'))
p_temp_Daphnia$date <- factor(p_temp_Daphnia$date, levels=c('APRIL1','APRIL5','APRIL8',
																														'APRIL12','APRIL15','APRIL19',
																														'APRIL26','MAY4'))


#calculate absolute growth rates (difference between recent and previous timepoint/number of days between sampling)
#Biovol alge
p_temp <- p_temp[order(p_temp$ID),]
p_temp$voldiff <- ave(p_temp$biovol, p_temp$ID, FUN=function(x) c(0, diff(x)))
p_temp$growthrate <- p_temp$voldiff/p_temp$daydiff
p_temp <- na.omit(p_temp)

#abundance Daphnia
p_temp_Daphnia<-p_temp_Daphnia %>%
	group_by(ID) %>%
	arrange(date) 
p_temp_Daphnia$abu_diff <- ave(p_temp_Daphnia$dapnia_count, p_temp_Daphnia$ID, FUN=function(x) c(0, diff(x)))
p_temp_Daphnia$growthrate <- p_temp_Daphnia$abu_diff/p_temp_Daphnia$daydiff
p_temp_Daphnia <- na.omit(p_temp_Daphnia)


library(plotrix)
library(broom)
p_temp_Daphnia %>%
	dplyr::select(rel_growth, P, temp) %>%
	filter(rel_growth > 0) %>% 
	filter(temp != "24") %>% 
	select(-ID) %>%
	group_by(P, temp) %>%
	mutate(temp.num = as.numeric(as.character(temp))) %>% 
	mutate(inverse.temp = 1/(.00008617*(temp.num+273.15))) %>%
	select(-(starts_with("ID"))) %>%
	summarise_each(funs(mean,max, median, sd,std.error)) %>%
	ggplot(data = ., aes(inverse.temp_max, y = rel_growth_mean, group = P, color = P)) +
	geom_errorbar(aes(ymin=rel_growth_mean-rel_growth_std.error, ymax=rel_growth_mean+rel_growth_std.error), width=.2) +
	geom_point(size = 6) + theme_bw() + ylab("daphnia relative growth/day") + xlab("temperature, C") +
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold")) + scale_y_log10()


ptemp_daph_inverse <- p_temp_Daphnia %>%
	dplyr::select(rel_growth, P, temp) %>%
	filter(rel_growth > 0) %>% 
	filter(temp != "24") %>% 
	select(-ID) %>%
	group_by(P, temp) %>%
	mutate(temp.num = as.numeric(as.character(temp))) %>% 
	mutate(inverse.temp = 1/(.00008617*(temp.num+273.15))) %>%
	select(-(starts_with("ID"))) %>%
	mutate(log.growth = log(rel_growth)) %>% 
	ggplot(data =., aes(x = inverse.temp, y = log.growth, group = P, color = P)) + geom_point() +
	geom_smooth(method = "lm")
ggsave("daph.growth.arrhenius.png")

mod <- lm(log(rel_growth) ~ inverse.temp, data = ptemp_daph_inverse)
summary(mod)

ptemp_daph_inverse %>% 
	group_by(P) %>%
	do(tidy(lm(log(rel_growth) ~ inverse.temp, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>%
	ggplot(data = ., aes(x = P, y = estimate)) + geom_point() + 
	geom_errorbar(aes(ymin = conf.low, ymax = conf.high, width=.2))
ggsave("activation_energy_daph_population_growth.png")
#relative growth rates
#Algae
p_temp$rel_growth <- p_temp$voldiff/p_temp$biovol

#relative growth Daphnia
p_temp_Daphnia$rel_growth <- p_temp_Daphnia$abu_diff/p_temp_Daphnia$dapnia_count

#Plot preparation
#Algae without controls
p_temp$temp <- as.factor(p_temp$temp)
sums_algae <- p_temp[,c(1,2,3,7,12,13)]
sums_algae <- sums_algae[sums_algae$ID<49,]
sums_algae <- as.data.frame(as.list(aggregate(. ~ days+temp+P,data = sums_algae[,c(2:6)],
																							FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_algae <- na.omit(sums_algae)

#Algae controls
sums_algae_controls <- p_temp[,c(1,2,3,7,12,13)]
sums_algae_controls <- sums_algae_controls[sums_algae_controls$ID>48,]

sums_algae_controls <- as.data.frame(as.list(aggregate(. ~ days+temp+P,data = sums_algae_controls[,c(2:6)],
																											 FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_algae_controls <- na.omit(sums_algae_controls)

#Daphnia
p_temp_Daphnia$temp <- as.factor(p_temp_Daphnia$temp)
p_temp_Daphnia<-do.call(data.frame,lapply(p_temp_Daphnia, function(x) replace(x, is.infinite(x),NA)))
sums_daphnia <- as.data.frame(as.list(aggregate(. ~ days+temp+P,data = p_temp_Daphnia[,c(7,13,14,18,19)],
																								FUN=function(x) c(mn =mean(x), n=length(x), se=sd(x)/length(x)))))
sums_daphnia <- na.omit(sums_daphnia)

#plot absolute growth rates
#plot algae without controls
plots_algae<-list()
for (k in 1:length(levels(sums_algae$temp))){
	p<-ggplot(data =sums_algae[sums_algae$temp==levels(sums_algae$temp)[k],], aes(days, y = growthrate.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('days') + ggtitle(paste('Temperature:',levels(sums_algae$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		ylim(min(sums_algae$growthrate.mn)-max(sums_algae$growthrate.se),max(sums_algae$growthrate.mn+sums_algae$growthrate.se))+
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
	p<-ggplot(data =sums_algae_controls[sums_algae_controls$temp==levels(sums_algae_controls$temp)[k],], aes(days, y = growthrate.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('days') + ggtitle(paste('Temperature:',levels(sums_algae$temp)[k],'ºC'))+
		geom_hline(aes(yintercept=0))+
		ylim(min(sums_algae_controls$growthrate.mn)-max(sums_algae_controls$growthrate.se),max(sums_algae_controls$growthrate.mn+sums_algae_controls$growthrate.se))+
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
	p<-ggplot(data =sums_daphnia[sums_daphnia$temp==levels(sums_daphnia$temp)[k],], aes(days, y = growthrate.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=growthrate.mn-growthrate.se, ymax=growthrate.mn+growthrate.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('days') + ggtitle(paste('Temperature:',levels(sums_daphnia$temp)[k],'ºC'))+
		ylim(-16,16)+
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

#------------------------------------------------------------------------------------------------------------------#
#plot relative growth rates
#Daphnia + Algae
plots_rel_growth<-list()
for (k in 1:length(levels(sums_algae$temp))){
	p<-ggplot(data =sums_algae[sums_algae$temp==levels(sums_algae$temp)[k],], aes(days, y = rel_growth.mn, group = P, colour=P)) +
		geom_line(size=1.5,linetype="dotted") + 
		geom_line(data =sums_daphnia[sums_daphnia$temp==levels(sums_daphnia$temp)[k],], 
							aes(date, y = rel_growth.mn, group = P, colour=P),size=1.5)+
		#geom_hline(aes(yintercept=0))+
		theme_bw() + ylab("relative growth rate") + xlab('') + ggtitle(paste('Temperature:',levels(sums_algae$temp)[k],'ºC'))+	
		ylim(-85,10)+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_rel_growth[[k]]<-p
}
png(filename = "relative_growth.png",
		width = 1400, height = 1000, units = "px",
		pointsize = )
multiplot(plotlist=plots_rel_growth,cols=2)
dev.off()

#Daphnia only
plots_rel_growth_daphnia<-list()
for (k in 1:length(levels(sums_daphnia$temp))){
	p<-ggplot(data =sums_daphnia[sums_daphnia$temp==levels(sums_daphnia$temp)[k],], aes(days, y = rel_growth.mn, group = P, color = P)) +
		geom_errorbar(aes(ymin=rel_growth.mn-rel_growth.se, ymax=rel_growth.mn+rel_growth.se), width=.2) +
		geom_point(size = 6) + theme_bw() + ylab("growth rate") + xlab('') + ggtitle(paste('Temperature:',levels(sums_daphnia$temp)[k],'ºC'))+
		ylim(min(sums_daphnia$rel_growth.mn)-max(sums_daphnia$rel_growth.se),max(sums_daphnia$rel_growth.se+sums_daphnia$rel_growth.mn))+
		geom_hline(aes(yintercept=0))+
		theme(title = element_text(size=20),
					legend.text = element_text(size = 22),
					axis.text.y = element_text(size = 16),
					axis.text.x = element_text(size = 16),
					axis.title.y = element_text(size = 20))
	plots_rel_growth_daphnia[[k]]<-p
}
png(filename = "relative_growth_Daphnia.png",
		width = 1400, height = 1000, units = "px",
		pointsize = )
multiplot(plotlist=plots_rel_growth_daphnia,cols=2)
dev.off()



#########
#get maximum growth rates for different temperatures
#algae
max_algae <- as.data.frame(as.list(aggregate(. ~ sums_algae$temp+sums_algae$P,data = sums_algae[,c(4,7)],
																						 FUN=function(x) c(max_growth =max(x)))))
growth_days_vec<-"NA"
rel_growth_days_vec<-"NA"
for (i in 1:length(max_algae$growthrate.mn)){
	growth_days_vec[i]<-as.character(sums_algae$days[sums_algae$growthrate.mn==max_algae$growthrate.mn[i]])
	rel_growth_days_vec[i]<-as.character(sums_algae$days[sums_algae$rel_growth.mn==max_algae$rel_growth.mn[i]])
}
png(filename = "max_growth_rates.png",
		width = 1000, height = 750, units = "px")
ggplot(data = max_algae, aes(sums_algae.temp, y = log(growthrate.mn),label=paste(growth_days_vec,'days'), group = sums_algae.P, color = sums_algae.P)) +
	geom_point(size = 15) + theme_bw() + ylab("log(growth rate)") + xlab('temperature') +
	geom_text(size=8, nudge_y = 0.25)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()

png(filename = "max_relative_growth_rates.png",
		width = 1000, height = 750, units = "px")
ggplot(data = max_algae, aes(sums_algae.temp, y = rel_growth.mn, group = sums_algae.P,label=paste(growth_days_vec,'days'), color = sums_algae.P)) +
	geom_point(size = 15) + theme_bw() + ylab("relative growth rate") + xlab('temperature') +
	geom_text(size=10, nudge_y = 0.02)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()

#####################################################################################################################
#####################################################################################################################
#controls
max_controls <- as.data.frame(as.list(aggregate(. ~ sums_algae_controls$temp+sums_algae_controls$P,data = sums_algae_controls[,c(4,7)],
																								FUN=function(x) c(max_growth =max(x)))))
growth_days_ctrl<-"NA"
rel_growth_days_ctrl<-"NA"
for (i in 1:length(max_controls$growthrate.mn)){
	growth_days_ctrl[i]<-as.character(sums_algae_controls$days[sums_algae_controls$growthrate.mn==max_controls$growthrate.mn[i]])
	rel_growth_days_ctrl[i]<-as.character(sums_algae_controls$days[sums_algae_controls$rel_growth.mn==max_controls$rel_growth.mn[i]])
}
png(filename = "max_growth_rates_controls.png",
		width = 1300, height = 750, units = "px")
ggplot(data = max_controls, aes(sums_algae_controls.temp, y = log(growthrate.mn),label=paste(growth_days_ctrl,'days'), group = sums_algae_controls.P, color = sums_algae_controls.P)) +
	geom_point(size = 15) + theme_bw() + ylab("log(growth rate)") + xlab('temperature') +
	geom_text(size=8, nudge_y = 0.25)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()

png(filename = "max_relative_growth_rates_controls.png",
		width = 1300, height = 750, units = "px")
ggplot(data = max_controls, aes(sums_algae_controls.temp, y = rel_growth.mn,label=paste(rel_growth_days_ctrl,'days'), group = sums_algae_controls.P, color = sums_algae_controls.P)) +
	geom_point(size = 15) + theme_bw() + ylab("relative growth rate") + xlab('temperature') +
	geom_text(size=8, nudge_y = 0.02)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()

#####################################################################################################################
#####################################################################################################################
#Daphnia
max_daphnia <- as.data.frame(as.list(aggregate(. ~ sums_daphnia$temp+sums_daphnia$P,data = sums_daphnia[,c(4,7)],
																							 FUN=function(x) c(max_growth =max(x)))))
growth_days_daphnia<-"NA"
rel_growth_days_daphnia<-"NA"
for (i in 1:length(max_daphnia$growthrate.mn)){
	growth_days_daphnia[i]<-as.character(sums_daphnia$days[sums_daphnia$growthrate.mn==max_daphnia$growthrate.mn[i]])
	rel_growth_days_daphnia[i]<-as.character(sums_daphnia$days[sums_daphnia$rel_growth.mn==max_daphnia$rel_growth.mn[i]])
}
png(filename = "max_growth_rates_daphnia.png",
		width = 1300, height = 750, units = "px")
ggplot(data = max_daphnia, aes(sums_daphnia.temp, y = growthrate.mn,label=paste(growth_days_daphnia,'days'), group = sums_daphnia.P, color = sums_daphnia.P)) +
	geom_point(size = 15) + theme_bw() + ylab("growth rate") + xlab('temperature') +
	geom_text(size=8, nudge_y = .75)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()

png(filename = "max_relative_growth_rates_daphnia.png",
		width = 1300, height = 750, units = "px")
ggplot(data = max_daphnia, aes(sums_daphnia.temp, y = rel_growth.mn,label=paste(rel_growth_days_daphnia,'days'), group = sums_daphnia.P, color = sums_daphnia.P)) +
	geom_point(size = 15) + theme_bw() + ylab("relative growth rate") + xlab('temperature') +
	geom_text(size=8, nudge_y = .05)+
	theme(axis.title=element_text(size=30),
				axis.text=element_text(size=25),
				legend.title=element_text(size=25),
				legend.text=element_text(size=25))
dev.off()


#### Arrhenius plots of max growth rates, Daphnia

#Boltzmann's constant in eV
k=.00008617

#convert celsius values to Kelvin
max_daphnia$tvalues<-as.numeric(as.character(max_daphnia$sums_daphnia.temp))
max_daphnia$tvalues<-max_daphnia$tvalues+273.15

#convert kelvin values to 1/kT
max_daphnia$xval<-1/(k*max_daphnia$tvalues)

#convert response values to log(response) values
max_daphnia$yval<-log(max_daphnia$rel_growth.mn)



max_daphnia %>% 
group_by(sums_daphnia.P) %>%
	 ggplot(data = ., aes(x = xval, y = yval, group = sums_daphnia.P, color = sums_daphnia.P)) +
	geom_point(size = 3) + theme_bw() + ylab("ln(daphnia max population growth rate)") + xlab("temperature, 1/kt") +
	stat_summary(fun.y= "mean", geom = "point") +
	geom_smooth(method = 'lm') + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/daphnia_max_growth_arrhenius.png")

mod1 <- lm(yval ~ xval*sums_daphnia.P, data = max_daphnia)
summary(mod1)

#### Arrhenius plots of max growth rates, algae controls

#Boltzmann's constant in eV
k=.00008617

#convert celsius values to Kelvin
max_controls$tvalues<-as.numeric(as.character(max_controls$sums_algae_controls.temp))
max_controls$tvalues<-max_controls$tvalues+273.15

#convert kelvin values to 1/kT
max_controls$xval<-1/(k*max_controls$tvalues)

#convert response values to log(response) values
max_controls$yval<-log(max_controls$rel_growth.mn)



max_controls %>% 
	group_by(sums_algae_controls.P) %>%
	ggplot(data = ., aes(x = xval, y = yval, group = sums_algae_controls.P, color = sums_algae_controls.P)) +
	geom_point(size = 3) + theme_bw() + ylab("ln(algae max population growth rate)") + xlab("temperature, 1/kt") +
	stat_summary(fun.y= "mean", geom = "point") +
	geom_smooth(method = 'lm') + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/algae_controls_max_growth_arrhenius.png")

mod1 <- lm(yval ~ xval, data = max_controls)
summary(mod1)

