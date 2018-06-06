## regression model development for routines CO2 flux
## PNA data from Heather

## get necessary packages
library("ggplot2")
library("reshape2")

## get necessary explanatory data sets
co2expl <- readRDS("../data/private/co2explained.rds") # YEAR, Month, LAKE, each sampling date
precip <- readRDS("../data/precip.rds") # Year, Month, for snow and rain and both
nao <- readRDS("../data/naoseasonal.rds") # trimonth subsets for a whole year; all colnames lowercase
#   nao not in model right now
soistand <- readRDS("../data/soi_stand.rds") # monthly; all colnames uppercase
pdo <- readRDS("../data/pdo.rds") # monthly, and annual mean: all colnames lowercase
temp <- readRDS("../data/temperaturedata.rds")# Year, Month, Superstation...
pna <- read.csv("../data/dataforallteleconectionsuptoFeb2015.csv") # file from Heather 19th Jul 2016
emi <- readRDS("../data/emi.rds") # file from Heather 19th Jul 2016

## get CO2 flux estimates (alter calculation arguments in co2_scenarios.R)
fluxes <- readRDS("../data/private/params-flux.rds")

## choose params that we want in model from co2expl offerings
## at the mo, relative humidity is in as a proxy for evaporation effects.. Considered this better than using
##    just precipitation since that may not be directly linked due to advective-dominated precip
## NPP and R very neatly related; selecting GPP and R. Might add that the bias in
##    1:1 is generally in favour of NPP
regvars <- subset(co2expl,
                  select = c("YEAR", "Month", "Date", "DOY", "LAKE", "Chl_a_ug_L", "GPP_h",
                             "R_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "Chl_a_ug_L_sur",
                             "AirTempMonthly", "RelHum"))

## do whatever needs to be done with precip data. and March things
## ============================================================================
mycolSums <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  sum <- colSums(subdf, na.rm = TRUE)
  cbind(data.frame(Year = df['Year'][1,], df['Month'][1,], t(sum)))
}

# grab total snow fall in march
snowsplit <- with(precip, split(precip, list(Year, Month), drop = TRUE))
snowtotal <- lapply(snowsplit, mycolSums, cols = c("Total.Snow..cm."))
snowtotal <- do.call("rbind", snowtotal)
rownames(snowtotal) <- NULL
snowtotal <- subset(snowtotal, select = c(Year, Total.Snow..cm.), df..Month...1... == 3)
colnames(snowtotal)[which(colnames(snowtotal) == "Total.Snow..cm.")] <- "MarchSnowFallcm"


# grab total precip in march.. FIXME: snow-specific better?
monthsplit <- with(precip, split(precip, list(Year, Month), drop = TRUE))
monthtotal <- lapply(monthsplit, mycolSums, cols = c("Total.Precip..mm."))
monthtotal <- do.call("rbind", monthtotal)
rownames(monthtotal) <- NULL

# grab march values for temp
marchtemp <- subset(temp, select = c(Year, Temperature),
                                  subset = Month == 3 & Superstation == "regina")
colnames(marchtemp)[which(colnames(marchtemp) == "Temperature")] <- "MarchTemp"
## ============================================================================

## change colnames and designations to be compatible
colnames(regvars)[which(colnames(regvars) == "YEAR")] <- "Year"
colnames(regvars)[which(colnames(regvars) == "LAKE")] <- "Lake"

colnames(pdo)[which(colnames(pdo) == "year")] <- "Year"
pdostack <- stack(pdo, select = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                                  "sep", "oct", "nov", "dec"))
pdostack$Year <- rep(pdo$Year, times = 12)
colnames(pdostack)[which(colnames(pdostack) == "ind")] <- "Month"
colnames(pdostack)[which(colnames(pdostack) == "values")] <- "PDO"
pdostack <- transform(pdostack, Month = rep(c(1:12), each = 116))

names(soistand) <- c("Year", 1:12)
soistack <- stack(soistand, select = -Year)
soistack$Year <- rep(soistand$Year, times = 12)
colnames(soistack)[which(colnames(soistack) == "ind")] <- "Month"
colnames(soistack)[which(colnames(soistack) == "values")] <- "SOI"

## check NA situation
fluxna <- which(is.na(fluxes$co2Flux))
nrow(fluxes[fluxna,]) # due to missing DIC mostly
nrow(fluxes) # 362 NA out of 1141

nanumbers <- rowSums(is.na(regvars))
dataloss <- regvars[nanumbers > 0,]
nrow(dataloss) # 444 out of 1264

Nna <- which(is.na(regvars$TDN_ug_L)) #214
DOCna <- which(is.na(regvars$DOC_mg_L)) #148

## merge harmonised tables; start: 1264 rows in regvars
regvars <- merge(regvars, pdostack) # 1264
regvars <- merge(regvars, soistack) # 1264
regvars <- merge(regvars, marchtemp) # 1264
regvars <- merge(regvars, snowtotal) # 1264
regvars <- merge(regvars, pna[,c('year', 'month','PNA')], by.x = c('Year', 'Month'),
                 by.y = c('year','month'))
regvars <- merge(regvars, emi[,c('Year','Month','EMI')])

vardiff <- regvars$Date[regvars$Date %in% fluxes$Date == FALSE]
fluxdiff <- fluxes$Date[fluxes$Date %in% regvars$Date == FALSE]
vardiff[order(vardiff)] # 26 dates not in fluxes
fluxdiff[order(fluxdiff)]
# bottle production estimates only started in 1996.... so if we want -94 and -95 we can't use NPP, R.
regvars <- merge(regvars, fluxes[,c('Lake', 'Date', 'co2Flux', 'lakepCO2', 'meanWindMS','pH')], 
                 all = TRUE)

## subset to lakes we want
regvars <- subset(regvars, Lake %in% c("K", "L", "B", "C", "D", "WW", "P"))

## rename pH column for following scripts
names(regvars)[which(names(regvars) == 'pH')] <- 'pH_surface'

## deal with outliers and other data issues
## =============================================================================
## CO2 calculation variables already dealt with in co2_scenarios.R, and GPP high outlier
##    dealt with in co2ExplVar.R
fluxtona <- which(regvars$co2Flux > 500) # two points, one for L and one for WW
regvars$co2Flux[fluxtona] <- NA
regvars$lakepCO2[fluxtona] <- NA

## make grossly -ve DOC into NA (although I suspect this is a typo i.e. - sign added by accident)
lowdoc <- which(regvars$DOC_mg_L < -60)
regvars$DOC_mg_L[lowdoc] <- NA

## change NaN chlorophylls to NA
chlnan <- which(is.nan(regvars$Chl_a_ug_L))
regvars$Chl_a_ug_L[chlnan] <- NA
chlnan <- which(is.nan(regvars$Chl_a_ug_L_sur))
regvars$Chl_a_ug_L_sur[chlnan] <- NA

## any outliers left?
drops <- c("Month", "Year", "Lake", "DOY")
dropna <- which(is.na(regvars$co2Flux))
regmelt <- regvars[-dropna,]
regmelt <- regmelt[,!(names(regvars) %in% drops)]
regmelt <- melt(regmelt, id = "Date")

outplot <- ggplot(data = regmelt, aes(x= Date, y = value, group = variable)) +
  ylab("vars") +
  geom_point() +
  facet_wrap( "variable", scales = "free") +
  theme(legend.position = "top")
outplot

## save output for rmd and model development
saveRDS(regvars, "../data/private/regvars.rds")

## look at relationships etc.
with(regvars, plot(MarchTemp ~ MarchSnowFallcm))
with(regvars, plot(GPP_h ~ pH_surface))
with(regvars, plot(RelHum ~ PDO))
with(regvars, plot(pH_surface ~ MarchTemp))
with(regvars, plot(PDO ~ SOI))
with(regvars, plot(pH_surface ~ DOC_mg_L))
with(regvars[regvars$Lake == "K" | regvars$Lake == "P",], plot(pH_surface ~ format(Date, "%m")))
with(regvars, plot(pH_surface ~ AirTempMonthly))

## some lake diffs
co2plot <- ggplot(data = regvars, aes(x= Date, y = co2Flux, group = Lake)) +
  ylab("CO2 flux") +
  geom_point() +
  facet_wrap( "Lake" ) +
  theme(legend.position = "top")
co2plot
pdf("../data/private/CO2fluxroutines.pdf")
co2plot
dev.off()
co2plot <- ggplot(data = regvars, aes(x= format(Date, "%m"), y = co2Flux, group = Lake)) +
  ylab("CO2 flux") +
  geom_point() +
  facet_wrap( "Lake" ) +
  theme(legend.position = "top")
co2plot



