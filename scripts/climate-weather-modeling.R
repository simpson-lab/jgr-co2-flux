## honing in on appropriate inclusion of climate in the models
## also calculating the evaporative balance by month using the spei index
##    http://sac.csic.es/spei/database.html

## read in data and merge
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")
alldat <- merge(regvars, fluxes[,c('Lake', 'Date','Conductivity',"TIC")])

soi <- readRDS("../data/soi_stand.rds") # monthly; all colnames uppercase
pdo <- readRDS("../data/pdo.rds") # monthly, and annual mean: all colnames lowercase
precip <- readRDS("../data/precip.rds")
airtemp <- readRDS("../data/temperaturedata.rds")

speireg <- read.csv("../data/spei-regina.csv")
speicrook <- read.csv("../data/spei-crooked.csv")
speidief <- read.csv("../data/spei-diefenbaker.csv")
speikat <- read.csv("../data/spei-katepwa.csv")
speilast <- read.csv("../data/spei-lastmountain.csv")
speipas <- read.csv("../data/spei-pasqua.csv")
speibuff <- read.csv("../data/spei-buffalopound.csv")

## load packages
library('zoo')
library("mgcv")
library("ggplot2")
library("SPEI")

## ====================================================================================
## Create lagged data frames; needed for 1, 2, 5:11
## ====================================================================================
## create lagged SOI by three-month means (dec-feb, jan-mar, feb-apr)
colnames(soi)[grep("YEAR", colnames(soi))] <- "Year" 
soistack <- stack(soi, select = -Year)
soistack$Year <- rep(soi$Year, times = 12)
colnames(soistack)[which(colnames(soistack) == "ind")] <- "Month"
colnames(soistack)[which(colnames(soistack) == "values")] <- "SOI"
soistack <- transform(soistack, Month = rep(c(1:12), each = length(unique(soistack$Year))))
soistack$Date <- as.yearmon(paste(soistack$Year , soistack$Month , sep = "-" ))

soistack <- soistack[order(soistack$Date),]
rownames(soistack) <- NULL

allDecs <- which(soistack$Month==12)
djf <- lapply(allDecs[-length(allDecs)], seq, length.out=3) #remove final Dec with no next year
jfm <- lapply(djf, function(x) x+1) 
fma <- lapply(djf, function(x) x+2) 
mam <- lapply(djf, function(x) x+3) 
amj <- lapply(djf, function(x) x+4) 
mjj <- lapply(djf, function(x) x+5) 
jja <- lapply(djf, function(x) x+6) # til november next year to be applied

aso <- lapply(djf, function(x) x-4) # for the january next year
son <- lapply(djf, function(x) x-3) # for the feb 

rowsubs <- function(df, rowlist, colname) {
  mylist <- list() #lapply(1:length(rowlist), function(x) matrix(NA, nrow=3, ncol=4))
  for (i in 1:length(rowlist)) {
    mylist[[i]] <- df[rowlist[[i]],colname]
  }
  mylist
}
djfs <- rowsubs(soistack, djf, "SOI")
jfms <- rowsubs(soistack, jfm, "SOI")
fmas <- rowsubs(soistack, fma, "SOI")
mams <- rowsubs(soistack, mam, "SOI")
amjs <- rowsubs(soistack, amj, "SOI")
mjjs <- rowsubs(soistack, mjj, "SOI")
jjas <- rowsubs(soistack, jja, "SOI")

asos <- rowsubs(soistack, aso, "SOI")
sons <- rowsubs(soistack, son, "SOI")

meandjf <- as.data.frame(do.call(rbind, lapply(djfs, mean)))
meanjfm <- as.data.frame(do.call(rbind, lapply(jfms, mean)))
meanfma <- as.data.frame(do.call(rbind, lapply(fmas, mean)))
meanmam <- as.data.frame(do.call(rbind, lapply(mams, mean)))
meanamj <- as.data.frame(do.call(rbind, lapply(amjs, mean)))
meanmjj <- as.data.frame(do.call(rbind, lapply(mjjs, mean)))
meanjja <- as.data.frame(do.call(rbind, lapply(jjas, mean)))

meanaso <- as.data.frame(do.call(rbind, lapply(asos, mean)))
meanson <- as.data.frame(do.call(rbind, lapply(sons, mean)))

SOImeans <- rbind(meanaso, meanson, meandjf, meanjfm, meanfma, meanmam, meanamj, meanmjj, meanjja)
colnames(SOImeans) <- "SOImean"

allmeans <- SOImeans
allmeans$YearofDec <- rep(unique(soistack$Year)[-(nrow(meandjf) + 1)], times=9) # remove last year
allmeans$YearApplied <- rep(unique(soistack$Year)[-1], times = 9) 
allmeans$MonthApplied <- rep(c(1,2,5,6,7,8,9,10,11), each=nrow(meanfma))

## create lagged PDO by three-month means
colnames(pdo)[which(colnames(pdo) == "year")] <- "Year"
pdostack <- stack(pdo, select = c("jan", "feb", "mar", "apr", "may", "jun", "jul", "aug",
                                  "sep", "oct", "nov", "dec"))
pdostack$Year <- rep(pdo$Year, times = 12)
colnames(pdostack)[which(colnames(pdostack) == "ind")] <- "Month"
colnames(pdostack)[which(colnames(pdostack) == "values")] <- "PDO"
pdostack <- transform(pdostack, Month = rep(c(1:12), each = 116))
pdostack$Date <- as.yearmon(paste(pdostack$Year , pdostack$Month , sep = "-" ))

pdostack <- pdostack[order(pdostack$Date),]
rownames(pdostack) <- NULL

allDecs <- which(pdostack$Month==12)
djf <- lapply(allDecs[-length(allDecs)], seq, length.out=3) #remove last Dec with no next year
jfm <- lapply(djf, function(x) x+1) 
fma <- lapply(djf, function(x) x+2) 
mam <- lapply(djf, function(x) x+3) 
amj <- lapply(djf, function(x) x+4) 
mjj <- lapply(djf, function(x) x+5) 
jja <- lapply(djf, function(x) x+6) # til november next year to be applied

aso <- lapply(djf, function(x) x-4) # for the january next year
son <- lapply(djf, function(x) x-3) # for the feb 

djfs <- rowsubs(pdostack, djf, "PDO")
jfms <- rowsubs(pdostack, jfm, "PDO")
fmas <- rowsubs(pdostack, fma, "PDO")
mams <- rowsubs(pdostack, mam, "PDO")
amjs <- rowsubs(pdostack, amj, "PDO")
mjjs <- rowsubs(pdostack, mjj, "PDO")
jjas <- rowsubs(pdostack, jja, "PDO")

asos <- rowsubs(pdostack, aso, "PDO")
sons <- rowsubs(pdostack, son, "PDO")

meandjf <- as.data.frame(do.call(rbind, lapply(djfs, mean)))
meanjfm <- as.data.frame(do.call(rbind, lapply(jfms, mean)))
meanfma <- as.data.frame(do.call(rbind, lapply(fmas, mean)))
meanmam <- as.data.frame(do.call(rbind, lapply(mams, mean)))
meanamj <- as.data.frame(do.call(rbind, lapply(amjs, mean)))
meanmjj <- as.data.frame(do.call(rbind, lapply(mjjs, mean)))
meanjja <- as.data.frame(do.call(rbind, lapply(jjas, mean)))

meanaso <- as.data.frame(do.call(rbind, lapply(asos, mean)))
meanson <- as.data.frame(do.call(rbind, lapply(sons, mean)))

pdomeans <- rbind(meanaso, meanson, meandjf, meanjfm, meanfma, meanmam, meanamj, meanmjj, meanjja)
colnames(pdomeans) <- "PDOmean"
pdomeans$YearofDec <- rep(unique(pdostack$Year)[-(nrow(meandjf) + 1)], times=9) # remove last year
pdomeans$YearApplied <- rep(unique(pdostack$Year)[-1], times = 9) 
pdomeans$MonthApplied <- rep(c(1,2,5,6,7,8,9,10,11), each=nrow(meanfma))

allmeans <- merge(allmeans, pdomeans)

## ===============================================================================================
## merge winter averages with appropriate month and year in alldat
## ===============================================================================================
alldat <- merge(alldat, allmeans[,c("YearApplied","MonthApplied","SOImean","PDOmean")], 
                by.x=c("Year", "Month"), by.y = c("YearApplied", "MonthApplied"))

## ===============================================================================================
## calculate evapotranspiration stuff
## ===============================================================================================
## since we decided to use the regina station for all measurements, I will also here subset
##    to regina
airtemp <- subset(airtemp, Superstation == "regina", select = c(Year, Month, Temperature))
airtemp <- airtemp[order(airtemp$Year, airtemp$Month),]
rownames(airtemp) <- NULL #set new order for rows, otherwise in mixed order based on above command

mycolSums <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  sum <- colSums(subdf, na.rm = TRUE)
  cbind(data.frame(Year = df['Year'][1,], df['Month'][1,], t(sum)))
}

precipsplit <- with(precip, split(precip, list(Year, Month), drop = TRUE))
preciptotal <- do.call(rbind, lapply(precipsplit, mycolSums, 
                                     cols = "Total.Precip..mm.")) # assumes 1:10 ratio water:snow
colnames(preciptotal) <- c("Year", "Month", "Precipmm")
preciptotal <- preciptotal[order(preciptotal$Year, preciptotal$Month),]
rownames(preciptotal) <- NULL

weathers <- merge(preciptotal, airtemp)
weathers <- weathers[order(weathers$Year, weathers$Month),]
weathers$evtrans <- thornthwaite(weathers$Precipmm, lat=51.0) #approximate latitude
# monthly potential evapotranspiration (mm)
#speilakes <- spei(weathers$Precipmm-weathers$evtrans, 1)
## FIXME: this does not work for some reason and August is just NA (some other months
##    also depending on what lag and kernel I use....)
## Instead let's use spei available from website... http://sac.csic.es/spei/database.html
##    with Regina's coordinates

speireg$Month <- rep(1:12, times=nrow(speireg)/12)
speireg$Year <- rep(1901:2014, each=12) # select lags one and two
speireg$Lake <- rep("WW") # select lags one and two
speicrook$Month <- rep(1:12, times=nrow(speicrook)/12)
speicrook$Year <- rep(1901:2014, each=12) # select lags one and two
speicrook$Lake <- rep("C") # select lags one and two
speidief$Month <- rep(1:12, times=nrow(speidief)/12)
speidief$Year <- rep(1901:2014, each=12) # select lags one and two
speidief$Lake <- rep("D") # select lags one and two
speikat$Month <- rep(1:12, times=nrow(speikat)/12)
speikat$Year <- rep(1901:2014, each=12) # select lags one and two
speikat$Lake <- rep("K") # select lags one and two
speilast$Month <- rep(1:12, times=nrow(speilast)/12)
speilast$Year <- rep(1901:2014, each=12) # select lags one and two
speilast$Lake <- rep("L") # select lags one and two
speipas$Month <- rep(1:12, times=nrow(speipas)/12)
speipas$Year <- rep(1901:2014, each=12) # select lags one and two
speipas$Lake <- rep("P") # select lags one and two
speibuff$Month <- rep(1:12, times=nrow(speibuff)/12)
speibuff$Year <- rep(1901:2014, each=12) # select lags one and two
speibuff$Lake <- rep("B") # select lags one and two

speiall <- rbind(speireg, speibuff, speipas, speilast, speidief, speikat, speicrook)
speisub <- subset(speiall, Year >= 1994, select = c("SPEI01", "SPEI02", "Year", "Month", "Lake"))

weathers <- merge(weathers, speisub)
weathers <- merge(weathers, allmeans[,c("YearApplied","MonthApplied","SOImean","PDOmean")], 
                  by.x = c("Year", "Month"), 
                  by.y = c("YearApplied", "MonthApplied"))
weathers <- weathers[order(weathers$Year, weathers$Month),]

saveRDS(weathers, "../data/weathers.rds")

alldat <- merge(alldat, weathers)


## ===============================================================================================
## view some data
## ===============================================================================================
## subset and plot SOI
soisub <- subset(soistack, Year > 1993)
soisub$Month
with(soisub, plot(SOI ~ Date, col=ifelse(Month<5 | Month>8, 'blue', 'red')))
lines(soisub$SOI[order(soisub$Date)] ~ soisub$Date[order(soisub$Date)])

## subset and plot PDO
pdosub <- subset(pdostack, Year > 1993)
with(pdosub, plot(PDO ~ Date, col=ifelse(Month<5 | Month>8, 'blue', 'red')))
lines(pdosub$PDO[order(pdosub$Date)] ~ pdosub$Date[order(pdosub$Date)])

## run some models
condmod <- gam(Conductivity ~
                 ti(SOI) + ti(PDO) +
                 ti(PDO, SOI) +
                 s(Lake, bs = "re"), 
               data = alldat,
               select = TRUE, method = "REML", family = scat(),
               na.action = na.exclude,
               control = gam.control(nthreads = 3, trace = TRUE, 
                                     newton = list(maxHalf = 60)))
condmodlag <- gam(Conductivity ~ # ti model showed interaction required
                    te(SOImean, PDOmean) +
                    s(SPEI01, k=4) +
                    s(Lake, bs = "re"), 
                  data = alldat,
                  select = TRUE, method = "REML", family = scat(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE, 
                                        newton = list(maxHalf = 60)))
condmodlagpna <- gam(Conductivity ~ # ti model showed interaction required
                       s(PNA, k=4) +
                       s(SPEI01, k=4) +
                       s(Lake, bs = "re"), 
                     data = alldat,
                     select = TRUE, method = "REML", family = scat(),
                     na.action = na.exclude,
                     control = gam.control(nthreads = 3, trace = TRUE, 
                                           newton = list(maxHalf = 60)))
condmodplus <- gam(Conductivity ~
                     te(SOImean, PDOmean) +
                     s(Chl_a_ug_L) + # logging here didn't seem to improve resids
                     s(AirTempMonthly) +
                     s(Lake, bs = "re"), 
                   data = alldat,
                   select = TRUE, method = "REML", family = scat(),
                   na.action = na.exclude,
                   control = gam.control(nthreads = 3, trace = TRUE, 
                                         newton = list(maxHalf = 60)))

dicmod <- gam(TIC ~
                te(SOImean, PDOmean) + # ti model showed interaction required
                s(SPEI01, k=4) +
                s(Lake, bs = "re"), 
              data = alldat,
              select = TRUE, method = "REML", family = scat(),
              na.action = na.exclude,
              control = gam.control(nthreads = 3, trace = TRUE, 
                                    newton = list(maxHalf = 60)))
phmod <- gam(pH_surface ~
               te(SOImean, PDOmean) + # ti model showed interaction required
               s(SPEI01, k=4) +
               s(Lake, bs = "re"), 
             data = alldat,
             select = TRUE, method = "REML", family = scat(),
             na.action = na.exclude,
             control = gam.control(nthreads = 3, trace = TRUE, 
                                   newton = list(maxHalf = 60)))

plot(condmodlag, select = 3, scheme = 2)
with(alldat, text(x=PDOmean, y=SOImean, labels = Year))

fluxdicmod <- gam(co2Flux ~
                    s(TIC) +
                    s(TIC, by=Lake, m=1), 
                  data = alldat,
                  select = TRUE, method = "REML", family = scat(),
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE, 
                                        newton = list(maxHalf = 60)))
