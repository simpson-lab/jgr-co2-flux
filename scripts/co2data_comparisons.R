## received complete time series of data from kerri (.../fromkerri/): pco2foremmasept2015
##    and pco2foremmatidy - the latter has been saved as csv
## will now see if my data reproduce her results

## get necessary data from previous processing
if (file.exists("../data/private/gasFlux.rds")) {
  co2emma <- readRDS("../data/private/gasFlux.rds")
} else {
  print("run pressuremanipulations.R")
}

## read kerri's data and comment on parameters
co2kerri <- read.csv("../data/private/kerri_co2flux_results_complete.csv")
names(co2kerri) <- c("Lake", "Date", "pH", "Alkalinity", "Temperature", "Conductivity", "pCO2atm", "Pressure",
                     "Altitude", "Ionstrength", "pk1", "pk2", "ao", "a1", "a2", "pkh", "DIC", "CO2uM", "CO2uatm",
                     "CO2eq", "CO2log", "CO2sat", "alpha", "Flux", "FluxEnh")
co2kerri <- transform(co2kerri, Date = as.POSIXct(as.character(Date), format = "%d-%m-%Y"))
co2kerri <- transform(co2kerri, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")))

co2kerrisub <- subset(co2kerri, select = c('Lake', 'Date', 'Year', 'Month', 'Flux', 'FluxEnh'))
# Other params such as Salinity and Wind were originally on a different sheet and therefore not on this
#   file
# alk is calculated based on DIC
# pCO2atm is 370
# Pressure is 101.325 (sea level pressure)
# Altitude is 509.3
# Ionstrength gives values where cond missing so there are a few instances where the IF clause is invoked
#     (i.e. mean taken), but these should disappear once I date match
# DIC entered, not calculated
# CO2uM is the gasExchange.R one (co2 <- dic*ao)
# CO2uatm is ibid (pco2 <- co2/(10^-pkh)
# CO2eq is calculated with the IF clause that defaults to alt.. I used Pressure
# Wind: is the annual average value of each site (though in fact is constant 4.06 m/s for this practice run)
# Salinity: 0 entered for all rows. Not calculated from cond relationship either

## Grab my flux calculations and check whether dates can be matched
original <- readRDS("../data/private/gasFlux.rds") # (pco2atm = 370) (see pressuremanipulations.R)
joined <- merge(co2kerrisub, original, by = c("Lake", "Date"), all.x = TRUE)

## play with some basic diagnostics
with(joined, plot(FluxEnh ~ co2Flux, xlab = "my flux", ylab = "kerri's flux")) # y is from kerri, x from me
abline(0,1)
with(joined, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)"))
# so here, the difference is salinity, pressure (rather than altitude), and wind
with(joined, max(FluxEnh - co2Flux, na.rm = TRUE)) # [1] 1029.797
with(joined, which(FluxEnh - co2Flux == max(FluxEnh - co2Flux, na.rm = TRUE)))
# [1] 303; in the warm July day of 2009 when pH recorded = 12.2 !!!???
with(joined, plot(FluxEnh - co2Flux ~ pH, xlab = "pH", ylab = "difference(kerri's - mine)"))
# what??? pH > 13???
with(joined, which(pH >= 13))
# joined[850,] on a September day in 2011.
## 12.2 and 13.5 also in lab notebooks; sensor error.

## need to source function that will run through the calculations
source("../functions/gasExchange.R")

## start messing around with parameters
## ====================================
##    Part 1: create complete similarity
##    created new gasExchangealt.R to use alt rather than kpa and suppress salt calculation
source("../functions/gasExchangealt.R")
identical <- transform(joined,
                    co2Flux = gasExchangealt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC,
                                          pco2atm = rep(370, times = nrow(joined)),
                                          alt = rep(509.3, times = nrow(joined)),
                                          wind = rep(4.06, times = nrow(joined)),
                                          salt = rep(0, times = nrow(joined))))
with(identical, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(identical, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)"))
# assuming R rounding errors here

##    Part 2: introduce salt (all else same)
##    created new gasexchangesalt.R to phase back the salt calculation and use salt
source("../functions/gasExchangesalt.R")
salty <- transform(joined,
                       co2Flux = gasExchangesalt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC,
                                                pco2atm = rep(370, times = nrow(joined)),
                                                alt = rep(509.3, times = nrow(joined)),
                                                wind = rep(4.06, times = nrow(joined)),
                                                salt = Salinity))
with(salty, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(salty, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)", col = Lake))
with(salty, which(FluxEnh - co2Flux == max(FluxEnh - co2Flux, na.rm = TRUE)))
# one "outlier" here is a normal day in 2001, the other one of the high pH days;
#     salt is used in equations r1 and r2
take <- with(salty, which(abs(FluxEnh - co2Flux) >= 10))
salty[take,]

##    Part 3: introduce wind (all else same)
windy <- transform(joined,
                       co2Flux = gasExchangealt(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC,
                                                pco2atm = rep(370, times = nrow(joined)),
                                                alt = rep(509.3, times = nrow(joined)),
                                                wind = Wind,
                                                salt = rep(0, times = nrow(joined))))
with(windy, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(windy, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)", col = Lake))
with(windy, which(FluxEnh - co2Flux == max(FluxEnh - co2Flux, na.rm = TRUE)))
# again max difference at the pH anomaly. funnily #303 seems to produce most differences, #850 not so

##    Part 4: replace alt with Pressure (all else same)
pressed <- transform(joined,
                       co2Flux = gasExchange(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC,
                                                pco2atm = rep(370, times = nrow(joined)),
                                                kpa = Pressure,
                                                wind = rep(4.06, times = nrow(joined)),
                                                salt = rep(0, times = nrow(joined))))
with(pressed, plot(FluxEnh ~ co2Flux))
abline(0,1)
with(pressed, plot(FluxEnh - co2Flux ~ Date, xlab = "Date", ylab = "difference(kerri's - mine)"))
# using pressure rather than altitude incurs the least changes; less than 10 units of flux.

## Save some output for later
## saveRDS(, "data/private/.rds")

