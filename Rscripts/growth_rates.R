#calculate growth rates
datevec<-c('MARCH29','APRIL1','APRIL5','APRIL8',
					 'APRIL12','APRIL15','APRIL19',
					 'APRIL26','MAY4')

realdatevec<-c('2016-03-29','2016-04-01','2016-04-05','2016-04-08',
					 '2016-04-12','2016-04-15','2016-04-19',
					 '2016-04-26','2016-05-04')

p_temp$datediff <- NA
for (i in 1:length(datevec)){
	p_temp$datediff[p_temp$date==datevec[i]] <- as.Date(realdatevec[i])
}

#sorting data
p_temp<-p_temp %>%
	group_by(ID) %>%
	arrange(date) 

#calculate numerator (biovol difference between recent and previous timepoint)
p_temp$voldiff <- ave(p_temp$biovol, p_temp$ID, FUN=function(x) c(0, diff(x)))
p_temp$daydiff <- ave(p_temp$datediff, p_temp$ID, FUN=function(x) c(0, diff(x)))
p_temp$growthrate <- p_temp$voldiff/p_temp$daydiff




