#### Microbial respiration ####

mresp_16 <- read_csv("respirometry-data/16C_microbe_resp_June22016_Oxygen.csv")
mresp_12 <- read_csv("respirometry-data/12C_microbe_resp_June22016_Oxygen.csv")
mresp_24 <- read_csv("respirometry-data/24C_microbe_resp_June222016_Oxygen.csv")
mresp_20 <- read_csv("respirometry-data/20C_microbe_resp_June22016_Oxygen.csv")


mresp_16$temperature <- "16"
mresp_12$temperature <- "12"
mresp_20$temperature <- "20"
mresp_24$temperature <- "24"

mresp <- bind_rows(mresp_24, mresp_20, mresp_16, mresp_12) %>% 
	select(-C2)
mresp_long <- gather(mresp, well_number, oxygen, 3:25) %>%
	unite("unique_well_ID", temperature, well_number, remove = FALSE)



well_IDs <- read_csv("data-raw/microbe-resp-well-IDs-June-2-2016.csv")

m_resp <- left_join(mresp_long, well_IDs, by = "well_number") %>% 
	rename(uniqueID = UniqueID) %>% 
	mutate(uniqueID = as.character(uniqueID))

sep_June2 <- sep_June2 %>% 
	mutate(uniqueID = as.character(uniqueID))

m_resp_biovol <- left_join(m_resp, sep_June2, by = "uniqueID") %>% 
	filter(well_number != "A1") %>%
	mutate(biovol_per_well = ((cell_count*biovolume)/5))


# 12C ---------------------------------------------------------------------


## Check out the slopes
m_resp_biovol %>% 
	filter(temperature == "12") %>% 
	filter(Time_in_min > 10) %>% 
	ggplot(data =., aes(x = Time_in_min, y = oxygen)) + geom_point() +
	geom_smooth(method = "lm") +
	facet_wrap( ~ well_number)


control.slopes.12 <- m_resp_biovol %>% 
	filter(temperature == 12) %>% 
	filter(grepl('A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	summarise(mean.control.slope.12 = mean(estimate))

control.slope.12 <- control.slopes.12$mean.control.slope.12

slopes.12 <- m_resp_biovol %>% 
	filter(temperature == 12) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number, treatment, unique_well_ID, biovol_per_well) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(microbe.corr.slope = estimate - control.slope.12) %>%
	filter(microbe.corr.slope < 0)



# 16C ---------------------------------------------------------------------

## Check out the slopes
m_resp_biovol %>% 
	filter(temperature == "16") %>% 
	filter(Time_in_min > 10) %>% 
	ggplot(data =., aes(x = Time_in_min, y = oxygen)) + geom_point() +
	geom_smooth(method = "lm") +
	facet_wrap( ~ well_number)


control.slopes.16 <- m_resp_biovol %>% 
	filter(temperature == 16) %>% 
	filter(grepl('A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	summarise(mean.control.slope.16 = mean(estimate))

control.slope.16 <- control.slopes.16$mean.control.slope.16

slopes.16 <- m_resp_biovol %>% 
	filter(temperature == 16) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number, treatment, unique_well_ID, biovol_per_well) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(microbe.corr.slope = estimate - control.slope.12) %>%
	filter(microbe.corr.slope < 0)


# 20C ---------------------------------------------------------------------


## Check out the slopes
m_resp_biovol %>% 
	filter(temperature == "20") %>% 
	filter(Time_in_min > 10) %>% 
	ggplot(data =., aes(x = Time_in_min, y = oxygen)) + geom_point() +
	geom_smooth(method = "lm") +
	facet_wrap( ~ well_number)


control.slopes.20 <- m_resp_biovol %>% 
	filter(temperature == 20) %>% 
	filter(grepl('A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	summarise(mean.control.slope.20 = mean(estimate))

control.slope.20 <- control.slopes.20$mean.control.slope.20

slopes.20 <- m_resp_biovol %>% 
	filter(temperature == 20) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number, treatment, unique_well_ID, biovol_per_well) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(microbe.corr.slope = estimate - control.slope.20) %>%
	filter(microbe.corr.slope < 0)


# 24C ---------------------------------------------------------------------

## Check out the slopes
m_resp_biovol %>% 
	filter(temperature == "24") %>% 
	filter(Time_in_min > 10) %>% 
	ggplot(data =., aes(x = Time_in_min, y = oxygen)) + geom_point() +
	geom_smooth(method = "lm") +
	facet_wrap( ~ well_number)


control.slopes.24 <- m_resp_biovol %>% 
	filter(temperature == 24) %>% 
	filter(grepl('A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	summarise(mean.control.slope.24 = mean(estimate))

control.slope.24 <- control.slopes.24$mean.control.slope.24

slopes.24 <- m_resp_biovol %>% 
	filter(temperature == 24) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(Time_in_min > 10) %>% 
	group_by(well_number, treatment, unique_well_ID, biovol_per_well) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(microbe.corr.slope = estimate - control.slope.24) %>%
	filter(microbe.corr.slope < 0)


slopes.12$temperature <- 12
slopes.16$temperature <- 16
slopes.20$temperature <- 20
slopes.24$temperature <- 24

all.slopes <- bind_rows(slopes.12, slopes.16, slopes.20, slopes.24) %>% 
	mutate(microbe.corr.slope = -1*microbe.corr.slope) %>% 
	mutate(mass.corr.slope = microbe.corr.slope/biovol_per_well)




all.slopes %>% 
filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', unique_well_ID)) %>% 
ggplot(data = ., aes(x = factor(temperature), y = log(mass.corr.slope), fill = treatment)) + geom_boxplot()

all.slopes %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', unique_well_ID)) %>% 
	ggplot(data = ., aes(x = temperature, y = log(estimate*-1), group = treatment, color = treatment)) + geom_point() +
	geom_smooth(method = "lm")

all.slopes %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', unique_well_ID)) %>% 
	ggplot(data = ., aes(x = factor(temperature), y = log(biovol_per_well), fill = treatment)) + geom_boxplot() +
	geom_smooth(method = "lm")



all.slopes.mg <- all.slopes %>% 
	# dplyr::select(temperature, mass.corr.slope, uniqueID) %>% 
	purrr::by_row(~ convert_torr_to_mg_per_litre(.x$temperature, .x$mass.corr.slope), .to = "mass.corr.slope.in.mg") %>% 
	tidyr::unnest(mass.corr.slope.in.mg) %>% 
	as.data.frame()



all.slopes.mg %>% 
	mutate(inverse.temp = (1/(.00008617*(temperature+273.15)))) %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', unique_well_ID)) %>%
group_by(treatment) %>% 
		do(tidy(lm(log(mass.corr.slope) ~ inverse.temp, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	# as.data.frame() %>%
	ggplot(data = ., aes(x = treatment, y = estimate)) + geom_point() +
	geom_errorbar(aes(ymin=conf.low, ymax = conf.high), width=.2)

	


m_resp_biovol %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', uniquewell)) %>% 
	group_by(uniqueID, temperature, treatment) %>% 
	filter(Time_in_min > 20) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	# dplyr::arrange(conf.high) %>% 
	# group_by(temperature, treatment) %>% View
	# summarise(mean.slope = mean(estimate)) %>%
	ggplot(data = ., aes(x = temperature, y = estimate, group = treatment, color = treatment)) +
	geom_point() +
	geom_smooth(method = "lm")
geom_errorbar(aes(ymin=estimate-std.error, ymax=estimate + std.error), width=.2)



m_resp_biov <- m_resp_biovol %>% 
	mutate(total_biovol_pwell = (cell_count*biovolume)/5) %>% 
	filter(!grepl('A1|A2|A3|A4|A5', well_number)) %>% 
	filter(!grepl('24_A6|24_B5|24_D4|20_A6|20_B3|20_B4|20_D2|20_D3|16_A6|16_B1|16_B5|16_C1|16_C3|16_C4', uniquewell)) %>% 
	group_by(uniqueID, temperature, treatment, total_biovol_pwell) %>% 
	filter(Time_in_min > 20) %>% 
	do(tidy(lm(oxygen ~ Time_in_min, data = .), conf.int = TRUE)) %>% 
	filter(term != "(Intercept)") %>% 
	as.data.frame() %>%
	mutate(temperature = as.numeric(as.character(temperature))) %>% 
	mutate(mass_corr_slope = (estimate/total_biovol_pwell)*-1) %>% 
	mutate(inverse.temp = (1/(.00008617*(temperature+273.15)))) %>%
	# group_by(inverse.temp, treatment) %>% 
	# summarise(mean.slope = mean(mass_corr_slope), std.err.slope = std.error(mass_corr_slope)) %>% 
	ggplot(data = ., aes(x = inverse.temp, y = log(mass_corr_slope), group = treatment, color = treatment)) +
	geom_point() +
	geom_smooth(method = "lm")
geom_errorbar(aes(ymin=mean.slope - std.err.slope, ymax = mean.slope + std.err.slope), width=.2)


m_resp_biov %>% 
	group_by(treatment) %>% 
	summarise(mean = mean(total_biovol_pwell), std.err = std.error(total_biovol_pwell)) %>% 
	ggplot(data = ., aes(x = treatment, y = mean)) + geom_point() +
	geom_errorbar(aes(ymin=mean- std.err, ymax = mean + std.err), width=.2)



library(purrr)
str(m_resp)
as.data.frame(m_resp)

m_resp %>% 
	by_row(.d = ., ..f = convert_torr_to_mg_per_litre, temperature, oxygen_in_torr = oxygen, .to = "oxygen_in_mg") %>% View

mo <- m_resp %>% 
	select(temperature, oxygen)

mo$oxygen_in_mg <- mo %>% 
	pmap(convert_torr_to_mg_per_litre)

str(as.data.frame(mo))

m_resp2 <- bind_cols(m_resp, mo)




View(mo)
#### This didnt work, come back to it later
m_resp %>% 
	dplyr::select(temperature, oxygen, uniqueID) %>% 
	by_row(~ convert_torr_to_mg_per_litre(.x$temperature, .x$oxygen), .to = "oxygen in mg") %>%
	filter(temperature == 12) %>% View

?by_row
