#### arrhenius plots


daph <- read_csv("p_temp_daphnia_algae.csv")
body_sizes <- read_csv("Daphnia_body_sizes_May4.csv")

daph_mass <- body_sizes %>% 
	group_by(UniqueID) %>% 
	summarise(mean.length = mean(Length)) %>%
	mutate(mean.size = (0.0042*(mean.length^2.66))) %>% 
	rename(ID = UniqueID)

daph_final <- daph %>% 
	filter(date == "MAY4") %>% 
	rename(daphnia_final = daphnia_ab)

daph_final_mass <- left_join(daph_final, daph_mass, by = "ID")

daph_initial <- daph %>% 
	filter(date == "APRIL1") %>% 
	rename(daphnia_initial = daphnia_ab)


daph_r <- left_join(daph_initial, daph_final, by = "ID")
daph_r$daphnia_final[daph_r$daphnia_final == "0"] <- "1"
daph_r$daphnia_final <- as.numeric(daph_r$daphnia_final)

daph_r_mass <- left_join(daph_r, daph_mass, by = "ID")

	
mean(daph_mass$mean.size)
daph_r_mass$mean.size[is.na(daph_r_mass$mean.size)] <- "0.00002389045"
daph_r_mass$mean.size <- as.numeric(daph_r_mass$mean.size)
summary(daph_r_mass$mean.size)


daph_r_Ea <- daph_r_mass %>% 
	mutate(growth_rate = daphnia_final/daphnia_initial) %>%
	mutate(mass_corr_r = growth_rate*(mean.size^0.25)) %>% 
	filter(temp.y < 24) %>% 
	rename(., rvalues = mass_corr_r, tvalues = temp.y) %>%
	select(rvalues, tvalues, P.y) %>% 
	filter(!is.na(rvalues))

# biomass_ID$sample_weight[biomass_ID$sample_weight == "1.6343"]<-"3.0736"
daph_r_Ea$rvalues[daph_r_Ea$rvalues == "0"] <- "1"

daph_r_Ea$tvalues <- as.numeric(daph_r_Ea$tvalues)
daph_r_Ea$rvalues <- as.numeric(daph_r_Ea$rvalues)
View(daph_r_Ea)

estimate.m(daph_r_Ea$tvalues, daph_r_Ea$rvalues)

#Boltzmann's constant in eV
k=.00008617

#convert celsius values to Kelvin
daph_r_Ea$tvalues<-daph_r_Ea$tvalues+273.15

#convert kelvin values to 1/kT
daph_r_Ea$xval<-1/(k*daph_r_Ea$tvalues)

#convert response values to log(response) values
daph_r_Ea$yval<-log(daph_r_Ea$rvalues)


hist(daph_r_Ea$yval)

daph_r_Ea %>% 
	# filter(P.y == "FULL") %>% 
	group_by(P.y) %>% 
	ggplot(data = ., aes(x = xval, y = yval, group = P.y, color = P.y)) +
	geom_point(size = 3) + theme_bw() + ylab("population growth rate") + xlab("temperature, 1/kt") +
	stat_summary(fun.y= "mean", geom = "point") +
	geom_smooth(method = 'lm') + 
	theme(axis.text=element_text(size=16),
				axis.title=element_text(size=16,face="bold"))
ggsave("p-temp-figures_files/figure-html/growthrate_arrhenius.png")


daph_r_Ea %>% 
	group_by(P.y) %>%
	do(tidy(lm(yval ~ xval, data = .), conf.int = TRUE)) %>% View