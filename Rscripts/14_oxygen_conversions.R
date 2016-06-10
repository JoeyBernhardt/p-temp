##### script for converting torr to mg/l of oxygen

## at 12C, oxygen solubility is 2.1415 umol/l/torr
## at 16C, oxygen solubility is 1.9705 umol/l/torr
## at 20C, oxygen solubility is 1.8230 umol/l/torr
## at 24C, oxygen solubility is 1.6993 umol/l/torr

convert_torr_to_mg_per_litre <- function(temperature, oxygen_in_torr) {
	if (temperature == 16) {
		return(oxygen_in_torr*0.0000019705*32*1000)
	} else if (temperature == 20) { 
		return(oxygen_in_torr*0.0000018239*32*1000)
	} else if (temperature == 24) {
		return(oxygen_in_torr*0.0000016993*32*1000)
	} else if (temperature == 12) {
		return(oxygen_in_torr*0.0000021415*32*1000)
	}
	}

#### if mo is a data frame with 2 columns, one with temperature, the other with oxygen in torr
#### using the pmap function in purrr will create another column with the oxygen values converted to mg/l
mo$oxygen_in_mg <- mo %>% 
	pmap(convert_torr_to_mg_per_litre)



