## This document outlines the process used to develop the regression of fully 
##    resolved temporal series of all routines sites

## read in data set with all variables of interest
if (!file.exists("../data/private/regvars.rds")) {
  source("../scripts/regression_routines.R")
}
regvars <- readRDS("../data/private/regvars.rds")

## load necessary packages
library("car")
library("ggplot2")
library("reshape2")
library("mgcv")

## basic looking at data stuff
## ===========================
scatterplotMatrix(regvars[,-c(1:5)],pch=19,cex=.5,reg.line=F, 
                  spread=F,ellipse=F, col=c('gray60','#2957FF','#FF8000'),
                  col.axis='gray50')
scatterplotMatrix(scale(regvars[,-c(1:5)]),pch=19,cex=.5,reg.line=F, 
                  spread=F,ellipse=F, col=c('gray60','#2957FF','#FF8000'),
                  col.axis='gray50')
## can't see any strong correlations between variables chosen; PS TDP and TDN are rather tightly
##    correlated and TDP not in model... so TDN could be seen as a proxy for nutrient status,
##    rather than the effect of N uniquely.

## order lakes in prep for ggplot by water flow chain:
levels(regvars$Lake) 
# [1] "K"  "L"  "B"  "C"  "D"  "WW" "P"  "WC" "R"  "E"  "M"  "k" 
plotfactor <- data.frame(Lake = as.character(levels(regvars$Lake)), 
                         Level = c(5,3,2,6,1,7,4,8,9,10,11,12))
regvars <- merge(plotfactor, regvars)
regvars$Lake=factor(regvars$Lake, unique(regvars$Lake)[order(unique(regvars$Level))])

## melt regvars for plotting all predictors against CO2
drops <- c("Month", "Year", "Lake", "DOY", "Date")
dropna <- which(is.na(regvars$co2Flux))
regmelt <- regvars[-dropna,]
regmelt <- regmelt[,!(names(regvars) %in% drops)]
regmelt <- melt(regmelt, id = "co2Flux")

outplot <- ggplot(data = regmelt, aes(x= co2Flux, y = value, group = variable)) +
  ylab("vars") +
  geom_point() +
  facet_wrap( "variable", scales = "free") +
  geom_smooth(se=F, method='gam', formula=y~s(x), color='#2957FF') + # lulz
  theme(legend.position = "top")
outplot # what up pH

## one model for all lakes, lake not added as factor
gamall <- gam(co2Flux ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
                s(pH_surface, k = 20) + s(DOC_mg_L) + s(Oxygen_ppm) + 
                te(PDO, SOI), data = regvars,
              select = TRUE, method = "ML", family = scat(link = "identity"))

gam.check(gamall)
summary(gamall)

## one model for all lakes, lake added as factor
gamlake <- gam(co2Flux ~ Lake + s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
                 s(pH_surface, k = 20) + s(DOC_mg_L) + s(Oxygen_ppm) + 
                 te(PDO, SOI), data = regvars,
               select = TRUE, method = "ML", family = scat())
gam.check(gamlake)

## one model for all lakes, lake added as random effect
gamlake2 <- gam(co2Flux ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
                 s(pH_surface, k = 20) + s(DOC_mg_L) + s(Oxygen_ppm) + 
                 te(PDO, SOI) + s(Lake, bs = "re"), data = regvars,
               select = TRUE, method = "ML", family = scat())
gam.check(gamlake2)
summary(gamlake2)

## this is the version that Gavin helped create: first model pH, take residuals,
##    and model those based on the rest of the variables
## 1
phmod <- gam(co2Flux ~ s(pH_surface, k = 20), data = regvars,
             select = TRUE, method = "REML", family = scat(),
             na.action = na.exclude)
res <- resid(phmod, type = "pearson")
gam.check(phmod)
plot(phmod, pers = TRUE, pages = 1)
## FIXME: "Fitting terminated with step failure - check results carefully"

##2
mod <- gam(res ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
             s(DOC_mg_L) + s(Oxygen_ppm) + 
             te(PDO, SOI) + s(Lake, bs = "re"), data = regvars,
           select = TRUE, method = "REML", family = scat(),
           na.action = na.exclude)

plot(mod, pages = 1, pers = TRUE)
summary(mod)
gam.check(mod)
## note that gam stuff retains some form of mother index and NA row differences will 
##    not influence moving back and forth between models

## check residuals for temporal autocorrelation
regplus <- cbind(regvars, resid(mod))
names(regplus)[which(names(regplus) == "resid(mod)")] <- "resids"
ggplot(data = regplus, aes(x = format(Date, "%Y"), y = resids, group = Lake)) +
  ylab("residuals of 'mod'") +
  geom_point() + 
  #geom_line() + 
  facet_wrap( "Lake", scales = "free") +
  theme(legend.position = "top")
## doesn't seem to be any but could check formally!
## FIXME: best way to check autocorrelation formally for irregular spacing?
## +++ so yes, temporal autocorr may not be there because not uniquely varying
##    (e.g. pH may vary concordantly) but isn't that an autocorrelation structure
##    that should be dealt with in its own right?

## gam plots
plot(gamlake, pages=1, residuals=TRUE, pch=19, cex=0.25,
     scheme=1, col='#FF8000', shade=TRUE,shade.col='gray90')
plot(mod2, pages=1, residuals=TRUE, pch=19, cex=0.25,
     scheme=1, col='#FF8000', shade=TRUE,shade.col='gray90')

## evidence for DOC declining with the progression of summer?
ggplot(data = regvars, aes(x= Date, y = DOC_mg_L, group = Lake)) +
  ylab("DOC (mg/L)") +
  geom_point() + 
  facet_wrap( "Lake", scales = "free") +
  geom_smooth(se=F, method='gam', formula=y~s(x), color='#2957FF') + # lulz
  theme(legend.position = "top")

ggplot(data = regvars, aes(x= format(Date, "%m"), y = DOC_mg_L, group = Lake,
                           color = format(Date, "%Y"))) +
  ylab("DOC (mg/L)") +
  geom_point() +
  facet_wrap( "Lake", scales = "free") +
  theme(legend.position = "top")

