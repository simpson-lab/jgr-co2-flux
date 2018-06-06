## received complete time series of data from kerri (.../fromkerri/): pco2foremmasept2015
##    and pco2foremmatidy - the latter has been used for sheet1, and sheet 2 grabbed from former,
##    saved as "data/private/kerri_co2flux_results_complete.csv" and
##    "data/private/kerri_co2flux_results_complete_2.csv"
## these files have some assumptions that differ from the original calcs in the Nature pub;
##    this script intends to reproduce exactly what was done before (rather than check
##    reproducibility of code as per co2data_comparisons.R) and also allow for comparison of
##    conclusions of regressions run with different options. Therefore necessary columns have been added
##    to params (wind and altitude) or embedded in gasExchangeFlex.R (salt = 0, missing DIC, pco2atm)
## altitude values from /fromkerri 2007 dic based on cond xls are identical to Elevation in params

## NOTE: *DIC values*: Kerri did DIC (mg/L) = 26.57 + 0.018 * Conductivity (see email 1st June 2015)
##    for missing DIC data points!
##       *CO2uM, pCO2*: kerri used the spreadsheet to calculate from DIC, but when DIC missing
##    (even with the conductivity relationship), she used the final regression of pH to CO2 and pH to pCO2
##    in order to fill in the missing numbers - so this should only have been evoked when DIC & cond missing:
##    log CO2 (uM) = 10.94 -1.124*pH  (r2 = 0.94)
##    log pCO2 (uatm) = 12.076-1.101*pH (r2 = 0.95)
## ergo since we using only existing data, the CO2 stuff is ok.

## for this iteration of outlier fixes we use the metadata FIXME sheet from 30th Nov 15

## get necessary data
if (!file.exists("../data/maunaloa.csv")) {
  source("../functions/getmaunaloa.R")
}
ml <- read.csv("../data/maunaloa.csv")

if (file.exists("../data/private/params.rds")) {
  params <- readRDS("../data/private/params.rds")
} else {
  source("pressuremanipulations.R")
  params <- readRDS("../data/private/params.rds")
  }

if (file.exists("../data/private/db_metadata_outliers.csv")) {
  checkdates <- read.csv("../data/private/db_metadata_outliers.csv")
} else {
  stop("get database metadata FIXME from Emma")
}

if (file.exists("../data/private/qDICDOCupdate.csv")) {
  dicdocnew <- read.csv("../data/private/qDICDOCupdate.csv")
} else {
  stop("get 2013,2014 DIC data from Emma")
}

### source functions that will run through the calculations
source("../functions/gasExchangeFlex.R")
source("../data/private/salcalc.R")


## Create necessary values for starting parameters
# insert pco2atm from Mauna Loa
mlsub <- subset(ml, select = c('Year', 'Month', 'pCO2'))
names(mlsub) <- c("Year", "Month", "pco2atm")
params <- merge(params, mlsub, by = c("Year", "Month"))

# insert kerri's wind values (finlayetal2009 p. 2556); Pasqua not done here but I calculated average
#   as 3.7. NOTE: if I calculate means of available data I have for all lakes, the values are different.
winds <- data.frame(Lake = c("B", "C", "D", "K", "L", "P", "WW"),
                    kerriWindMS = c(4.1, 4.2, 3.4, 3.3, 4.3, 3.7, 2.8))
params <- merge(params, winds, by = "Lake") # adds Kerri's wind values
replace <- which(names(params) %in% "WindMS")
names(params)[replace] <- "measuredWindMS" # wind of the *sampling date*

## create salcalc column to replace erratic salinity values
## assuming 1m depth adds .1 bar to air pressure, can take rough mean of all sites
##    and add that value for dbar; converted kpa to dbar (see also salinityrecalcproject.R)
kpa <- mean(params$Pressure)
dbar <- kpa/10 + 1
params <- transform(params, SalCalc = salcalc(Temperature, Conductivity, dbar))

checkdates <- transform(checkdates, Date = as.Date(Date, format = "%m/%d/%Y", 
                                                   tz="Canada/Saskatchewan"))
wrongs <- params[params$Date %in% c(checkdates$Date),c("Lake", "Date", "Year", "pH", 
                                                       "Conductivity","Salinity", "SalCalc")]
wrongrows <- which(params$Date %in% c(checkdates$Date))

## grab original row numbers of variable values that are wrong within the outlier subset
phtona <- as.numeric(rownames(wrongs[which(wrongs$pH < 7 | wrongs$pH > 11),]))
newcond <- as.numeric(rownames(wrongs[which(wrongs$Conductivity == 13.85),])) # typo cond
keepas <- as.numeric(rownames(wrongs[which((wrongs$Lake == "P" & wrongs$Year == 2006) |
                                             (wrongs$Lake == "K" & wrongs$Year == 1994) |
                                             (wrongs$Lake == "C" & wrongs$Year == 1996)),]))
# keepas are those with one date many lake and only one lake dodgy
all <- c(phtona, newcond, keepas)
salcondtona <- as.numeric(rownames(wrongs[-c(which(rownames(wrongs) %in% c(all))),]))
salcondtona2 <- subset(wrongs[c(which(rownames(wrongs) %in% c(phtona))),], Lake == "D" & Year == 2002)
salcondtona2 <- as.numeric(rownames(salcondtona2)) # one where not only pH but also these need NA

## replace the wrong values
params[phtona,'pH'] <- NA
params[newcond, 'Conductivity'] <- 1385
params[newcond, 'SalCalc'] <- with(params[newcond,], salcalc(Temperature, Conductivity, dbar))
params[c(salcondtona, salcondtona2),c('Conductivity', 'SalCalc')] <- NA

## found one more outlier for D in 1998, and in 2012
makena <- which(params$Conductivity >= 2700)
params[makena, c('Conductivity', 'SalCalc')] <- NA
makena <- which(params$TICumol > 7500)
params[makena, c('TIC', 'TICumol')] <- NA

## insert DIC for 2013 and 2014; convert new DICs into umol
dicdocnew <- transform(dicdocnew, Date = as.Date(Date, format = "%d-%b-%y", tz="Canada/Saskatchewan"))
take <- dicdocnew[,c('LAKE', 'Date', 'TIC_mg_L')]
take <- merge(params[params$Year >= 2013,], take, by.y = c("LAKE", "Date"), by.x = c("Lake", "Date"), all.x = TRUE) 
take$TIC <- take$TIC_mg_L
take <- take[,-which(names(take) == 'TIC_mg_L')]
take <- take[,c("Lake", "Year", "Month", "Superstation", "StationID", "Date", "Temperature", "Conductivity", "Salinity", 
        "Depth", "pH", "Wind", "TIC", "DOY", "Elevation", "Day", "Pressure", "meanWind", "RelHum", 
        "measuredWindMS", "meanWindMS", "TICumol", "pco2atm", "kerriWindMS", "SalCalc")] # now in same params order
take$TICumol <- take$TIC / 0.012

takeout <- which(params$Year >= 2013)
params <- params[-takeout,]

params <- rbind(params, take)

## create our two scenarios; parameters in function = temp, cond, ph, wind, kerri = FALSE, salt = NULL,
##    dic = NULL, alt = NULL, kpa = NULL, pco2atm = NULL, trace = FALSE
## CHECKED: ran co2data_comparisons.R, took joined object, plotted co2Flux from joined and params
##    (scenario = new) and the data are nearly reproducing!! (with the archaic flux rds, i.e. different
##    wind data)
scenario <- "new"
co2Flux <- switch(scenario,
                  new      = with(params, gasExchangeFlex(Temperature, Conductivity, pH, meanWindMS, 
                                                          kerri= FALSE,
                                             salt = SalCalc, dic = TICumol, kpa = Pressure,
                                             pco2atm = pco2atm)),
                  original = with(params, gasExchangeFlex(Temperature, Conductivity, pH, kerriWindMS, kerri= TRUE,
                                             dic = TICumol, kpa = Pressure,
                                             alt = Elevation)))
params$co2Flux <- co2Flux$fluxenh
params$lakepCO2 <- co2Flux$pco2

## save version of flux output
saveRDS(params, "../data/private/params-flux.rds")
