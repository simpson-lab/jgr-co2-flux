## having done the CO2 flux ~ ... regression for high-reso routines data, this script
##    carries on and investigates what may be driving (covarying with) pH
## We know that the max that pH can achieve is related to the start-of-year pH, so
##    therefore we will also look at Year as a random effect!

## get necessary packages
library("ggplot2")
library("reshape2")
library("mgcv")

## load the data we need
if (!file.exists("../data/private/regvars.rds")) {
  source("../scripts/regression_routines.R")
}
regvars <- readRDS("../data/private/regvars.rds")

## create model; pH a function of production and climate removing effect of Year and Lake
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year))
phmod <- gam(pH_surface ~ s(Lake, Year, bs = "re") + s(Chl_a_ug_L) + s(GPP_h) +
               s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) + te(PDO, SOI), data = regvarf,
             select = TRUE, method = "REML", family = gaussian)
summary(phmod)
plot(phmod, pages = 1, pers = TRUE)
gam.check(phmod)

## create crude correlate of previous year's autumn productivity to see if signif
##    using chl a
## grab last two possible chl a data points
lastprod <- subset(regvars, Month == 8 | Month == 9, select = c(Lake, Year, Date, Chl_a_ug_L))
lastspl <- with(lastprod, split(lastprod, list(Lake, Year)))

## create wrappers that do colmeans for the rows we want and spits out df we can merge
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(Lake = df['Lake'][1,], Year = df['Year'][1,]), t(means))
}

grablast <- function(df) {
  grab <- c(nrow(df), nrow(df) -1)
  df <- df[grab,]
}

## apply wrappers
lastspl <- lapply(lastspl, grablast)
chlmeans <- lapply(lastspl, mycolMeans, cols = c("Chl_a_ug_L"))
chlmeans <- do.call(rbind, chlmeans)
rownames(chlmeans) <- NULL

## tidy chlmeans for merging: make mean apply to following year and use that to merge
chlmeans <- transform(chlmeans, PrevYear = Year + 1)
chlmeans <- chlmeans[,-which(names(chlmeans) == 'Date' | names(chlmeans) == 'Year')]
regvarl <- regvars

## merge by appropriate year
regvarl <- merge(regvarl, chlmeans, by.x = c('Lake', 'Year'), by.y = c('Lake', 'PrevYear'))
regvarl <- transform(regvarl, Year = as.factor(Year))

## create new regression to see if any sense in including it
phmod2 <- gam(pH_surface ~ s(Lake, bs = "re") + s(Year, bs = "re") + s(Chl_a_ug_L.x) +
                s(Chl_a_ug_L.y) +s(GPP_h) + s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) +
                te(PDO, SOI), data = regvarl,
             select = TRUE, method = "ML", family = scat())
summary(phmod2)
## nope, no sense

