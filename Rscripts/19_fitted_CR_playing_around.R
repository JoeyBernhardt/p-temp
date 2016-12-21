### messing around with the fitted parameter estimates Dec 13 2016


library(tidyverse)


data <- read_csv("fittedpdata5.csv")

summary(data$K)
data %>% 
	filter(K < 10000000000) %>% 
	# filter(r > 0.1, r < 5) %>% 
	group_by(Phosphorus) %>% 
ggplot(data = ., aes(x = temp, y = log(a), color = Phosphorus)) + geom_point() + 
	geom_smooth(method = "lm")
