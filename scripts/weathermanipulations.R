## we have 1080 dfs in total; in met, they do not have the first 9 lines
## grab met from gasexchange_pressuredata.R

if (!file.exists("../data/dmet.rds")) {
  if (!file.exists("../data/met.rds")) {
    source("../scripts/gasexchange_pressuredata.R")}
  met <- readRDS("../data/met.rds")
  dmet <- do.call("rbind", met)
  saveRDS(dmet, "../data/dmet.rds")
} else {dmet <- readRDS("../data/dmet.rds")}

## get station details
if (!file.exists("../data/weatherstations.csv")) {
  stop("create station list for Qu'Appelle as ../data/weatherstations.csv before running this script")
}
stations <- read.csv("../data/weatherstations.csv")

## Combine stations
dmet$superstation[dmet$StationID == 3002 | dmet$StationID == 51441] <- "regina"
dmet$superstation[dmet$StationID == 3062 | dmet$StationID == 44203 |
                    dmet$StationID == 48977] <- "yorkton"
dmet$superstation[dmet$StationID == 3318] <- "outlook"
dmet$superstation[dmet$StationID == 2926 | dmet$StationID == 2925] <- "indhead"
unique(dmet$superstation) # yes all appear, no NAs

names(dmet)[grep("Date", names(dmet))] <- "datetime"

take <- c("StationID","datetime","Year","Month","Day","Time","superstation","Stn Press (kPa)",
          "Wind Spd (km/h)", "Temp (degC)", "Rel Hum (%)")
dmet <- dmet[take]
names(dmet)[names(dmet) == "Stn Press (kPa)"] <- "Pressure"
names(dmet)[names(dmet) == "Wind Spd (km/h)"] <- "Wind"
names(dmet)[names(dmet) == "Temp (degC)"] <- "Temperature"
names(dmet)[names(dmet) == "Rel Hum (%)"] <- "RelHum"

## Make datetime into Date object
dmet <- transform(dmet,
                  datetime = as.POSIXct(as.character(datetime), format = "%Y-%m-%d %H:%M",
                                        tz="Canada/Saskatchewan"),
                  superstation = as.factor(superstation))

## Delete a priori-known missing data - find station dates at
  ## http://climate.weather.gc.ca/advanceSearch/searchHistoricData_e.html

dmet$datetimepos <-as.POSIXlt(dmet$datetime) # note that year is expressed since 1900;
  # yday starts at 0 (so leap year max = 365)

leaptest <- function(year) {
  if (year %% 4 != 0) {
    print("not leap") # remaining contenders move on
  }
  else if (year %% 100 == 0 & year %% 400 != 0) { # exclude those that are divisible by both
    print("not leap")
  }
  else {print("is leap")}
} # --> http://disc.gsfc.nasa.gov/julian_calendar.shtml

dmet <- dmet[!(dmet$datetimepos$year == 113 & dmet$datetimepos$yday > 281 & dmet$StationID == 3002),]
  # excludes all from and including 10.10 i.e. runs to 09.10.
dmet <- dmet[!(dmet$datetimepos$year == 113 & dmet$datetimepos$yday < 282 & dmet$StationID == 51441),]
  # excludes all up to and including 09.10 i.e. runs from 10.10.

dmet <- dmet[!(dmet$datetimepos$year == 105 & dmet$datetimepos$yday > 276 & dmet$StationID == 3062),]
  # excludes all from and including 05.10 i.e. runs to 04.10.
dmet <- dmet[!((dmet$datetimepos$year == 105 & dmet$datetimepos$yday < 277 & dmet$StationID == 44203) |
                  (dmet$datetimepos$year == 111 & dmet$datetimepos$yday > 235 & dmet$StationID == 44203)),]
  # excludes all to and including 04.10 i.e. runs from 05.10;
  # excludes all from and including 25.08 i.e. runs to 24.08.
dmet <- dmet[!(dmet$datetimepos$year == 111 & dmet$datetimepos$yday < 236 & dmet$StationID == 48977),]
  # excludes all up to and including 24.08. i.e. runs from 25.08.

dmet <- dmet[!(dmet$datetimepos$year == 107 & dmet$datetimepos$yday < 168 & dmet$StationID == 3318),]
  # excludes all up to and including 17.06. i.e. runs from 18.06.
dmet <- dmet[!(dmet$datetimepos$year < 107 & dmet$StationID == 3318),]
  # no pressure data prior to this time

dmet <- dmet[!(dmet$datetimepos$year == 109 & dmet$datetimepos$yday < 234 & dmet$StationID == 2925),]
  # excludes all up to and including 22.08 i.e. runs from 23.08.
dmet <- dmet[!(dmet$datetimepos$year < 109 & dmet$StationID == 2925),]
  # no pressure data prior to this time
dmet <- dmet[!(dmet$StationID == 2926),]
  # realised that this station was not one of the hourly data ones - i.e. no kPa data

length(which(duplicated(dmet$datetime[dmet$superstation == "regina"])))
length(which(duplicated(dmet$datetime[dmet$superstation == "outlook"])))
length(which(duplicated(dmet$datetime[dmet$superstation == "yorkton"])))
length(which(duplicated(dmet$datetime[dmet$superstation == "indhead"])))

unique(diff(dmet$datetime[dmet$superstation == "regina"]))
unique(diff(dmet$datetime[dmet$superstation == "yorkton"]))
  dmet$datetime[diff(dmet$datetime[dmet$superstation == "yorkton"]) == 61]
  dmet$datetime[diff(dmet$datetime[dmet$superstation == "yorkton"]) == 59]
  # seems datetime-related not gap in data-related
unique(diff(dmet$datetime[dmet$superstation == "outlook"]))
unique(diff(dmet$datetime[dmet$superstation == "indhead"]))

## Split apart
spldmet <- with(dmet, split(dmet, list(superstation, Year, Month)))
spldmet2 <- with(dmet, split(dmet, list(superstation, Year, Month, Day)))

### Wrapper functions
meanPressure <- function(df) {
  meanP <- if (NROW(df) == 0) {
    NA
  } else {
    mean(df[["Pressure"]], na.rm = TRUE)
  }
  with(df, data.frame(Superstation = superstation[1],
                      StationID = StationID[1],
                      Year = Year[1],
                      Month = Month[1],
                      Pressure = meanP))
}
meanHumidPC <- function(df) {
  meanH <- if (NROW(df) == 0) {
    NA
  } else {
    mean(df[["RelHum"]], na.rm = TRUE)
  }
  with(df, data.frame(Superstation = superstation[1],
                      StationID = StationID[1],
                      Year = Year[1],
                      Month = Month[1],
                      RelHum = meanH))
}

meanWind <- function(df) {
  meanW <- if (NROW(df) == 0) {
    NA
  } else {
    mean(df[["Wind"]], na.rm = TRUE)
  }
  with(df, data.frame(Superstation = superstation[1],
                      StationID = StationID[1],
                      Year = Year[1],
                      Month = Month[1],
                      Day = Day[1],
                      Wind = meanW))
}
meanTemp <- function(df) {
  meanT <- if (NROW(df) == 0) {
    NA
  } else {
    mean(df[["Temperature"]], na.rm = TRUE)
  }
  with(df, data.frame(Superstation = superstation[1],
                      StationID = StationID[1],
                      Year = Year[1],
                      Month = Month[1],
                      Temperature = meanT))
}

pressure <- na.omit(do.call("rbind", lapply(spldmet, meanPressure)))
rownames(pressure) <- NULL

winds <- na.omit(do.call("rbind", lapply(spldmet2, meanWind)))
rownames(winds) <- NULL

temperature <- na.omit(do.call("rbind", lapply(spldmet, meanTemp)))
rownames(temperature) <- NULL

relhum <- na.omit(do.call("rbind", lapply(spldmet, meanHumidPC)))
rownames(relhum) <- NULL

head(pressure[order(pressure$Superstation),], n = 20) # just looking
head(winds[order(winds$Superstation),], n = 20) # just looking

## save for next script
saveRDS(pressure[order(pressure$Superstation),], "../data/pressuredata.rds")
saveRDS(winds[order(winds$Superstation),], "../data/windsdata.rds")
saveRDS(temperature[order(temperature$Superstation),], "../data/temperaturedata.rds")
saveRDS(relhum[order(relhum$Superstation),], "../data/relhumdata.rds")

