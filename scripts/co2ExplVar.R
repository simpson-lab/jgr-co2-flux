## This script creates rds's that contain all relevant data for running regressions against pCO2:
## pdo.rds, soi_nonstand.rds, soi_stand.rds, naoseasonal.rds, co2explained.rds
## 1. download weather pattern data from web: PDO, SOI, NAO
## 2. process routines sampling data tables
## 3. get flow et al data from other sources

extras <- FALSE # for running extra bits at end

## Part 1: all these data are complete
## ==================================================================================================

## download PDO url and create data frame structure with numerics
if (!file.exists("../data/pdo.txt")) {
  download.file("http://research.jisao.washington.edu/pdo/PDO.latest", "../data/pdo.txt")
}

file <- read.fwf(file = "../data/pdo.txt", skip = 32, n = 118,
                 widths = c(8, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7))

colnames(file) <- unlist(file[1,])
pdo <- file[-c(1, 2),] # all as factors and straight to numeric gives gibberish
pdo <- data.frame(lapply(pdo, as.character), stringsAsFactors=FALSE)
pdo$YEAR.... <- gsub("**", "", pdo$YEAR...., fixed = TRUE) # remove asterisks for numeric conversion
pdo <- data.frame(lapply(pdo, as.numeric)) 
#   these asterisks actually mean "Derived from OI.v2 SST fields" (erm...)
pdo[pdo == 99.9] <- NA
names(pdo) <- c("year", "jan", "feb", "mar", "apr","may", "jun", "jul", "aug", "sep", "oct", "nov",
                "dec") # unlist() earlier maintained spaces in the names...
# could've had month.abb here to c all months in a year

## create yearly means since I believe we want annual reso
pdo <- transform(pdo, mean = rowMeans(pdo[,-1], na.rm = TRUE))

## create rds for later
saveRDS(pdo, "../data/pdo.rds")

## download SOI url and create data frames with numeric entries
if (!file.exists("../data/soi.txt")) {
  download.file("http://www.cpc.noaa.gov/data/indices/soi", "../data/soi.txt")
  }

all_data <- read.fwf(file = "../data/soi.txt", skip = 3, widths = rep(6, 13))

## create separate df for the unstandardised data
table1 <- all_data[c(1:66),]
colnames(table1) <- unlist(table1[1,])
soinonst <- table1[-c(1),] # all as factors and straight to numeric gives gibberish
soinonst <- lapply(soinonst, as.character)
soinonst <- lapply(soinonst, as.numeric)
soinonst <- do.call(cbind, soinonst)
soinonst[soinonst == 99.9] <- NA
soinonst <- as.data.frame(soinonst)

## create rds for later
saveRDS(soinonst, "../data/soi_nonstand.rds")

## create separate df for the standardised data
table2 <- all_data[c(75:140),]
colnames(table2) <- unlist(table2[1,])
soist <- table2[-c(1),] # all as factors and straight to numeric gives gibberish
soist <- lapply(soist, as.character)
soist <- lapply(soist, as.numeric)
soist <- do.call(cbind, soist)
soist[soist == 99.9] <- NA
soist <- as.data.frame(soist)

## create rds for later
saveRDS(soist, "../data/soi_stand.rds")


## download NAO url and make data frame with numeric entries
if (!file.exists("../data/nao.txt")) {
  download.file("https://climatedataguide.ucar.edu/sites/default/files/nao_station_seasonal.txt",
                "../data/nao.txt")
}

naos <- readLines("../data/nao.txt")
naos <- naos[-c(1:2)]
naos <- lapply(naos, gsub, pattern = "    ", replacement = " ")
naos <- lapply(naos, gsub, pattern = "   ", replacement = " ")
naos <- lapply(naos, gsub, pattern = "  ", replacement = " ")

naosplit <- lapply(naos, strsplit, " ")
naovec <- unlist(naosplit)
naomat <- matrix(naovec, ncol=13, byrow=TRUE)
naoframe <- as.data.frame(naomat)
colnames(naoframe) <- c("year", "djf", "jfm", "fma", "mam", "amj", "mjj", "jja", "jas", "aso",
                      "son", "ond", "ndj")
naoframe <- lapply(naoframe, as.character)
naoframe <- lapply(naoframe, as.numeric)
naoframe <- do.call(cbind, naoframe)
naoframe <- as.data.frame(naoframe)
naoframe[naoframe <= -999] <- NA

## create rds for later
saveRDS(naoframe, "../data/naoseasonal.rds")

## download EMI (El Nino Modoki) url and create data frame structure with numerics
if (!file.exists("../data/emi.txt")) {
  download.file("http://www.jamstec.go.jp/frsgc/research/d1/iod/DATA/emi.monthly.txt", 
                "../data/emi.txt")
}

emi <- read.table(file = "../data/emi.txt", header = TRUE)

colnames(emi)[grep("EMI", colnames(emi))] <- "EMI"
emi$Date <- as.character(emi$Date)
emidates <- do.call(rbind, strsplit(emi$Date, ":"))
emi$Year <- as.integer(emidates[,1])
emi$Month <- as.integer(emidates[,2])

## create rds for later
saveRDS(emi, "../data/emi.rds")

## Part 2: 1039 rows, 415 rows of one or more NAs
## ==================================================================================================

## read in supporting monitoring data tables grabbed from database
## I created a Date2 column in OpenOffice to convert the current display of month as character to month
##    as numeric (former is how it came from database into excel)
if(any(!file.exists(c("../data/private/qpco2supportdata.csv","../data/private/Production.csv",
                      "../data/private/Chl.csv", "../data/private/qprodsupportdata.csv",
                      "../data/private/qprofilesoxtemp.csv","../data/private/qsecchietal.csv",
                      "../data/private/Lakes.csv","../data/private/qDICDOCupdate.csv")))) {
  stop("get all necessary private files from Emma")
}

routines <- read.csv("../data/private/qpco2supportdata.csv")
routines <- transform(routines, Date = as.Date(as.character(Date2), format = "%Y-%m-%d",
                                               tz="Canada/Saskatchewan"))
routines <- routines[with(routines, order(LAKE, Date)),]

produc <- read.csv("../data/private/Production.csv")
produc <- transform(produc, Date = as.Date(as.character(Date2), format = "%Y-%m-%d",
                                           tz="Canada/Saskatchewan"))

chl <- read.csv("../data/private/Chl.csv")
chl <- transform(chl, Date = as.Date(as.character(Date2), format = "%Y-%m-%d",
                                     tz="Canada/Saskatchewan"))

incub <- read.csv("../data/private/qprodsupportdata.csv")
incub <- transform(incub, Date = as.Date(as.character(Date), format = "%d-%m-%Y",
                                         tz="Canada/Saskatchewan"))

oxtemp <- read.csv("../data/private/qprofilesoxtemp.csv")
oxtemp <- rbind(oxtemp, read.csv("../data/private/qprofilesextrarecordsoxtemp.csv"))
oxtemp <- transform(oxtemp, Date = as.Date(as.character(Date), format = "%d-%m-%Y",
                                           tz="Canada/Saskatchewan"))


secchi <- read.csv("../data/private/qsecchietal.csv")
secchi <- transform(secchi, Date = as.Date(as.character(Date), format = "%Y-%m-%d",
                                              tz="Canada/Saskatchewan"))

lakes <- read.csv("../data/private/Lakes.csv")
colnames(lakes)[which(colnames(lakes) == "Abbreviation")] <- "LAKE"

if (file.exists("../data/private/qDICDOCupdate.csv")) {
  dicdocnew <- read.csv("../data/private/qDICDOCupdate.csv")
} else {
  stop("get 2013,2014 DIC data from Emma")
}
dicdocnew <- transform(dicdocnew, Date = as.Date(Date, format = "%d-%b-%y",
                                                              tz="Canada/Saskatchewan"))

## remove extra date column (originally retained in case wanna check that the date format
##    conversion worked in OpenOffice)
routines <- routines[,-which(colnames(routines)=="Date2")]
produc <- produc[,-which(colnames(produc)=="Date2")]
chl <- chl[,-which(colnames(chl)=="Date2")]

## !!!!!!! In database: lake names other than WW/WC *at least* in 2012 have been entered
##    with a trailing space which means that merges won't work!
chl <- transform(chl, LAKE = gsub(" ", "", LAKE))


## take out all chl data that isn't from an integrated sample, as that's what we're working with
##    (though worth thinking that surface could be played with since they're quite different!!)
chlsub <- subset(chl, TreatmentNewLabel == "Integrated")
chlsur <- subset(chl, TreatmentNewLabel == "Surface")

## split produc and chl for getting means for replicated production estimates
prodsplit <- with(produc, split(produc, list(LAKE, Date), drop = TRUE))
chlsplit <- with(chlsub, split(chlsub, list(LAKE, Date), drop = TRUE))
chlsursplit <- with(chlsur, split(chlsur, list(LAKE, Date), drop = TRUE))
## create wrapper that does colmeans for the columns we want and spits out df we can merge
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(LAKE = df['LAKE'][1,], Date = df['Date'][1,]), t(means))
}

## choose columns we want means for
chlmeans <- lapply(chlsplit, mycolMeans, cols = c("Chl_a_ug_L")) # , "Total_chl"
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL
chlmeans <- chlmeans[with(chlmeans, order(LAKE, Date)),]

chlsurmeans <- lapply(chlsursplit, mycolMeans, cols = c("Chl_a_ug_L")) # , "Total_chl"
chlsurmeans <- do.call(rbind, chlsurmeans)
rownames(chlsurmeans) <- NULL
chlsurmeans <- chlsurmeans[with(chlsurmeans, order(LAKE, Date)),]
names(chlsurmeans)[grep("Chl", names(chlsurmeans))] <- "Chl_a_ug_L_sur"

prodmeans <- lapply(prodsplit, mycolMeans, cols = c("LIGHT_O2_ppm", "DARK_O2_ppm", "LIGHT_Pois_O2_ppm",
                                                    "DARK_Pois_O2_ppm", "NetOxy_ppm", "RespOxy_ppm",
                                                    "Net_mgC_m3_h", "Resp_mgC_m3_H",
                                                    "LightAzidemgCm3h", "DarkAzidemgCm3h"))
prodmeans <- do.call(rbind, prodmeans)
rownames(prodmeans) <- NULL
prodmeans <- prodmeans[with(prodmeans, order(LAKE, Date)),]

## merge prodmeans with the other prod support data in another table
## using AlainsCode.R in private/ which nicole sent me, to estimate P's and R.
proddf <- merge(prodmeans, incub, by = c("LAKE", "Date"))
proddf <- transform(proddf, photO2ppm = LIGHT_O2_ppm - ProdAssayO2Start_ppm)
proddf <- transform(proddf, respO2ppm = DARK_O2_ppm - ProdAssayO2Start_ppm)
proddf <- transform(proddf, NPP_h = photO2ppm/Hours_incubation)
proddf <- transform(proddf, R_h = respO2ppm/Hours_incubation)
proddf <- transform(proddf, GPP_h = NPP_h - R_h)

## leave proddf for later if need to revisit calcs, select what we want
prodsub <- subset(proddf, select = c(LAKE, Date, GPP_h, NPP_h, R_h))

### nicole noticed an outlier and there it is: WW 2010-05-03 (GPP >15 cf mostly <1)
## changing this to NA for now
makena <- which(prodsub$GPP_h > 15)
prodsub[makena, c('GPP_h', 'NPP_h', 'R_h')] <- NA

## finlay et al lakewide estimate method doesn't make sense so doing this as alain and nicole have done
##    (finlay et al imply that secchi divided by lake volume, but this gives a strange proportion rather
##    than an upscaling to whole-lake NPP)
prodsub <- merge(prodsub, lakes[,c('LAKE', 'LakeArea_km2')])
prodsub <- merge(prodsub, secchi[c('LAKE', 'Date', 'Secchi_m')])
prodsub <- transform(prodsub, lakeGPP = GPP_h*(Secchi_m * (LakeArea_km2 * 100000)))
prodsub <- transform(prodsub, lakeNPP = NPP_h*(Secchi_m* (LakeArea_km2 * 100000)))
prodsub <- transform(prodsub, lakeR = R_h*(Secchi_m* (LakeArea_km2 * 100000)))

## merge all dfs by LAKE and Date
co2explained <- merge(chlmeans, prodsub, all = TRUE)
co2explained <- merge(co2explained, chlsurmeans, all=TRUE)
co2explained <- merge(co2explained, routines[,-which(names(routines) %in% 'RunNo')], all = TRUE)
# don't need this column which is cryptic anyway
co2expl <- merge(co2explained, oxtemp[,c('Date', 'LAKE', 'Temperature_deg_C','Oxygen_ppm')],
                 by = c('Date', 'LAKE'), all = TRUE)
## FIXME: minor discrepancies in nrow over merge steps... a few rows.. could check at some point

## airtemp, from env canada and weathermanipulations.R
if (!file.exists("../data/temperaturedata.rds")) {
  source("../scripts/weathermanipulations.R")
}
airtemp <- readRDS("../data/temperaturedata.rds")
## since we decided to use the regina station for all measurements, I will also here subset
##    to regina
airtemp <- subset(airtemp, Superstation == "regina", select = c(Year, Month, Temperature))
airtemp <- airtemp[order(airtemp$Year, airtemp$Month),]
rownames(airtemp) <- NULL #set new order for rows, otherwise in mixed order based on above command

## create yearmeans for temperature, too
airsplit <- with(airtemp, split(airtemp, list(Year)))
airmeans <- lapply(airsplit, colMeans, na.rm = TRUE)
airmeans <- do.call(rbind, airmeans)
colnames(airmeans)[which(colnames(airmeans) == "Temperature")] <- "AirTempAnnual"
airmeans <- airmeans[,-(which(colnames(airmeans) == "Month"))]
airmeans <- as.data.frame(airmeans)
airmeans$Year <- as.numeric(airmeans$Year)
airtemp <- merge(airtemp, airmeans)
colnames(airtemp)[which(colnames(airtemp) == "Temperature")] <- "AirTempMonthly"
colnames(airtemp)[which(colnames(airtemp) == "Year")] <- "YEAR"

## merge co2expl with monthly air temperature data....
co2expl <- transform(co2expl, Month = as.numeric(format(Date, format = "%m")))
co2expl <- merge(co2expl, airtemp)

## relhum, from env canada and weathermanipulations.R
if (!file.exists("../data/relhumdata.rds")) {
  source("../scripts/weathermanipulations.R")
}
relhum <- readRDS("../data/relhumdata.rds")
## since we decided to use the regina station for all measurements, I will also here subset
##    to regina
relhum <- subset(relhum, Superstation == "regina", select = c(Year, Month, RelHum))
relhum <- relhum[order(relhum$Year, relhum$Month),]
names(relhum)[which(names(relhum) == "Year")] <- "YEAR"
rownames(relhum) <- NULL #set new order for rows, otherwise in mixed order based on above command

## merge co2expl with monthly relative humidity data
co2expl <- merge(co2expl, relhum)

## insert DIC and DOC for 2013 and 2014
take <- dicdocnew[,c('LAKE', 'Date', 'TIC_mg_L', 'DOC_mg_L')]
take <- merge(co2expl[co2expl$YEAR >= 2013,], take, by.y = c("LAKE", "Date"), by.x = c("LAKE", "Date"), 
              all.x = TRUE) 
take$TIC_mg_L <- take$TIC_mg_L.y
take$DOC_mg_L <- take$DOC_mg_L.y
take <- take[,-which(names(take) %in% c('TIC_mg_L.x','TIC_mg_L.y', 'DOC_mg_L.x', 'DOC_mg_L.y'))]
take <- take[,c("YEAR", "Month", "Date", "LAKE", "Chl_a_ug_L", "GPP_h", "NPP_h", "R_h", "LakeArea_km2", 
                "Secchi_m", "lakeGPP", "lakeNPP", "lakeR", "Chl_a_ug_L_sur", "DOY", "pH_surface", "SRP_ug_L", "TDP_ug_L", 
                "TDN_ug_L", "NO2_ug_L", "NH4_ug_L", "DOC_mg_L", "TIC_mg_L", "SO4_mg_L", "Si_mg_L", 
                "NO3_ug_L", "Temperature_deg_C", "Oxygen_ppm", "AirTempMonthly", "AirTempAnnual", "RelHum")] 
# now in same co2expl order

takeout <- which(co2expl$YEAR >= 2013)
co2expl <- co2expl[-takeout,]

co2expl <- rbind(co2expl, take)

## save output for later
saveRDS(co2expl, "../data/private/co2explained.rds")

## bonus: let's see how many NAs in our data. (don't select Si and other optionals)
dataloss <- subset(co2expl, select = c("LAKE", "Date", "pH_surface", "Oxygen_ppm", "TIC_mg_L",
                                       "Temperature_deg_C", "SRP_ug_L", "TDN_ug_L", "NO3_ug_L",
                                       "NH4_ug_L", "DOC_mg_L", "Chl_a_ug_L", "GPP_h"))

nanumbers <- rowSums(is.na(dataloss))
dataloss <- dataloss[rowSums(is.na(dataloss)) > 0,]
nrow(dataloss)

## Part 3: ice-out, inflow, evaporation
## ==========================================================================================
if (extras) {
  ## kerri emailed me this (in PhDPapers/Data..) pCO2_vars_yearlyaverageswithclimate.csv, but
  ##    it has only a few years, so we need to supplement it.
  ## other files downloaded at this point are from Rich. Go to 2010/2011
  ## kerri used the "old" way of calculating annual flows, but nicole and I checking what data
  ##    we can get to bring us up to 2015. probably will be the "new way"
  ## "new way": monthly data set from rich (full csv in !git/fromrich, sensible csv in /private..)
  ##    filled in a few missing annual totals in openoffice before saving as csv
  ## "old way": annual data set from rich, full csv in !git/fromrich, sensible csv in /private..)
  ##  Interestingly: why does finlay et al 2009 say that inflow for crooked was regressed, since it
  ##    exists in the data set that rich sent that should correspond to what kerri did?
  
  ## ice-out
  iceout1 <- read.csv("../data/private/IceOut.csv") # pending data to 2015
  icedataloss <- iceout1[is.na(iceout1$ICEOUTDOY),]
  nrow(iceout1) # 120
  nrow(icedataloss) # 45
  
  ## inflow flows
  inflow1 <- read.csv("../data/private/Monthly_Inflow_Estimates2.csv") # pending data to 2015
  inflowold <- read.csv("../data/private/Lake_Inflows.csv") # 1994-2009
  
  othervars <- read.csv("../data/private/pCO2_vars_yearlyaverageswithclimate.csv")
  ## this doesn't however indicate which vars interpolated and which real, but has the separate inflow
  ##    sites and inflow data
  evap <- subset(othervars, select = "evap" %in% colnames(othervars))
}