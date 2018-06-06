## outputting co2 exchange with all data combined as one table
## get necessary data from database queries

if (file.exists("../data/pressuredata.rds")) {
  pressuredata <- readRDS("../data/pressuredata.rds")
} else {
  source("../scripts/weathermanipulations.R")
  pressuredata <- readRDS("../data/pressuredata.rds")
  }

if (file.exists("../data/windsdata.rds")) {
  windsdata <- readRDS("../data/windsdata.rds")
} else {
  source("../scripts/weathermanipulations.R")
  windsdata <- readRDS("../data/windsdata.rds")
}
if (file.exists("../data/relhumdata.rds")) {
  relhumdata <- readRDS("../data/relhumdata.rds")
} else {
  source("../scripts/weathermanipulations.R")
  relhumdata <- readRDS("../data/relhumdata.rds")
}

if (file.exists("../data/private/fluxquery1.csv")) { # DIC, pH, wind
  surfd <- read.csv("../data/private/fluxquery1.csv")
} else {
  stop("get fluxquery1.csv from dropbox")
}

if (file.exists("../data/private/fluxquery2.csv")) { # elevation and lake
  # abbreviation data
  elevd <- read.csv("../data/private/fluxquery2.csv")
} else {
  stop("get fluxquery2.csv from dropbox")
}

if (file.exists("../data/private/fluxquery3.csv")) { # T, cond, salinity
  tcsd <- read.csv("../data/private/fluxquery3.csv")
} else {
  stop("get fluxquery3.csv from dropbox")
}

if (file.exists("../data/superstationz.csv")) { # links superstations with lakes
  superstations <- read.csv("../data/superstationz.csv")
} else {
  stop("get superstationz.csv from emma")
}

## are we really removing all NA data? THERE'S LOTS! complete.cases, na.omit
## NO - Leave NA's alone so we can see where they are at this stage
remNA <- FALSE
##
surfset <- surfd
#
tcsset <- tcsd
## Deal with NAs?
if (remNA) {
  tcsset <- droplevels(na.omit(tcsset))
  surfset <- droplevels(na.omit(surfset))
}

## create date objects for later
tcsset <- transform(tcsset, Date = as.Date(as.character(Date), format = "%d-%m-%Y", tz="Canada/Saskatchewan"))
surfset <- transform(surfset, Date = as.Date(as.character(Date), format = "%d-%m-%Y", tz="Canada/Saskatchewan"))

## merge the tcsset and surfset data sets on
joined <- merge(tcsset, surfset)

## Simplify some names
names(joined) <- c("Lake", "Date", "Temperature", "Conductivity", "Salinity",
                   "Depth", "pH", "Wind", "TIC", "DOY")

## Merge in elevd height data
joined <- merge(joined, subset(elevd, select = c("Abbreviation", "Elevation_m")),
                by.x = "Lake", by.y = "Abbreviation", sort = FALSE)
names(joined)[names(joined) == "Elevation_m"] <- "Elevation"

## select all lakes we're interested in
drop <- TRUE
if (drop) {
  joined <- subset(joined, Lake %in% c("B", "C", "WW", "D", "K", "L", "P"))
}

## choice of using or not using all superstations
regina <- TRUE

if (regina) {
  joined$Superstation <- "regina"
} else {
  joined <- merge(joined, superstations, sort = FALSE)
}

## merge joined with pressure and wind etc data
## add Year, Month and Day fields to aid the match
## FIXME: Still need to decide if and what to do with missing pressure data if using all stations
joined <- transform(joined, Year = as.numeric(format(Date, format = "%Y")),
                    Month = as.numeric(format(Date, format = "%m")),
                    Day = as.numeric(format(Date, format = "%d")))
joined <- merge(joined, pressuredata, sort = FALSE, all.x = TRUE)

# avoid ambiguity and change colname for met wind; add Date column
names(windsdata)[6] <- "metWind"

toDate <- function(year, month, day) {
  as.Date(paste(year, month, day, sep = "-"), tz="Canada/Saskatchewan") 
}

windsdata <- transform(windsdata, Date = toDate(Year, Month, Day))

## create appropriate averages for each instance (i.e. sampling date minus 14 days)
if(regina) {
  windsub <- subset(windsdata, Superstation == "regina")
  windgrab <- function(interval, WIND) {
    datalist <- WIND$metWind[WIND$Date <= interval[[2]][1] & WIND$Date >= interval[[2]][2]]
  }
} else {
  stop("to be developed") # i.e. case where, despite grabbing pressure data from regina, grabbing
  #   wind data from individual weather stations
}

if(regina) {
  windsub <- subset(windsdata, Superstation == "regina")
  spldf <- with(joined, split(joined, list(Lake), drop = TRUE)) # split into list by lakes
  sampledatelist <- lapply(spldf, "[", c("Lake", "Date")) # grab date and lake from each lake
  windraw <- list()
  windmeans <- list()
  winddf <- data.frame(Lake = character(0), sampleDate = character(0), meanWind = integer(0))
  for (i in seq_along(sampledatelist)) { # seq_along == 1:length(listObj)
      sampledates <- sampledatelist[[i]]
      interval <- list()
          for (i in seq_len(nrow(sampledates))) {
            interval[[i]] <- list(sampledates[['Lake']][1],
                                  seq(from = sampledates[i, 'Date'],
                                      by = "-13 days", length.out = 2))
          } # this gives a list the length of nrow(sampledates[[i]]) with lake id and date range
      windraw <- lapply(interval, windgrab, WIND = windsub) # gives list length of nrow(sampledates[[i]])
      #   with all 14 wind data (with by = -13 days)
      windmeans <- sapply(windraw, mean, na.rm = TRUE) # sapply unlists if possible!
      winddf <- rbind(winddf, data.frame(Lake = sampledates['Date'], sampleDate = sampledates['Lake'],
                                         meanWind = windmeans))
  }
} else {
  stop("not-regina to be developed") ## FIXME not-regina option null
}

joined <- merge(joined, winddf, by = c("Date", "Lake"), sort = FALSE, all.x = TRUE)

joined <- merge(joined, relhumdata, sort = FALSE, all.x = TRUE)

## need to change wind from km/h to m/s
joined <- transform(joined, WindMS = Wind * (1000/60/60))
joined <- transform(joined, meanWindMS = meanWind * (1000/60/60))

## need to change dic from mg/L to uM
joined <- transform(joined, TICumol = TIC / 0.012)
#   this is what Kerris' spreadsheet indicates for the unit conversion

## save joined as parameter data table
saveRDS(joined, "../data/private/params.rds")

## This is what was done before the creation of gasFluxFlex.R - redid using new joined but
##    essentially an archaic object now.
## Run the gas exchange equations on our data to create co2 flux
source("../functions/gasExchange.R")
joinedarchaic <- transform(joined,
                    co2Flux = gasExchange(temp = Temperature, cond = Conductivity, ph = pH, dic = TIC,
                                          pco2atm = 370, kpa = Pressure, wind = Wind, salt = Salinity))
#     archaic issue with naming in the database, dic = TIC is correct

## Save output for paper
saveRDS(joinedarchaic, "../data/private/gasFlux.rds")
