## adapted scripts from bottomoftheheap/...canadian climate data, from 24.05.2016,
##    which is nowadays flexible for both hourly and daily data (see git archives)

## source necessary functions:
source("../functions/getData.R")
source("../functions/genURLS.R")

## my station data frame

## nb all stations can cover my year range apart from outlook which has hourly
## only from 1996... there are randomly monthly data going way back but the
## monthly interface does not show pressure... if we want to request data it
## costs.... $100!!!! http://climate.weather.gc.ca/newPriceAnnounce_e.html

#maximum example now
stations <- data.frame(StationID = c(3002, 51441, 3062, 44203, 48977,
                       3318, 2926, 2925),
                       start = c(1994, 2013, 1994, 2005, 2011, 1996, 1994, 1996),
                       end = c(2013, 2015, 2005, 2011, 2015, 2015, 1996, 2015))
## see weatherstations.csv in /RScripts; these will cover lakes b,c,d,k,l,p,ww.

met <- getData(stations, folder = "../data/weatherdata", timeframe = 'hourly')
saveRDS(met, "../data/met.rds") # save for next script
