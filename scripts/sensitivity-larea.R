## sensitivity analysis adapted from sensitivity.R; here we use k600 values from 
##    Vachon and Prairie  dx.doi.org/10.1139/cjfas-2013-0241
## see https://cran.r-project.org/web/packages/pse/vignettes/pse_tutorial.pdf
##    https://cran.r-project.org/web/packages/sensitivity/sensitivity.pdf
## what indices to be used for expressing sensitivity?
## resources:
##  1. visual distribution cheat sheet http://www.itl.nist.gov/div898/handbook/eda/section3/eda366.htm
##  2. list of R shorts for distributions http://www.statmethods.net/advgraphs/probability.html
##  3. blurb on chisq: http://www.civil.uwaterloo.ca/brodland/EasyStats/EasyStats/Chi_squared_Distribution.html
##  4. fitdist issues: http://stats.stackexchange.com/questions/158163/why-does-this-data-throw-an-error-in-r-fitdistr
##  5. fitting distros: http://stats.stackexchange.com/questions/132652/how-to-determine-which-distribution-fits-my-data-best
##  6. chisq: https://www.youtube.com/watch?v=hcDb12fsbBU
##          http://www.civil.uwaterloo.ca/brodland/EasyStats/EasyStats/Chi_squared_Distribution.html

## source necessary packages:
library("mc2d")
library("pse")
library("fitdistrplus") # loads MASS too
library("reshape2")
library("extrafont")
library("ggplot2")
## source the function we want to decompose
source("../functions/gasExchangeSensitivityK600.R")

## set default for plot
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')

## source files with parameters and merge lake area data with 'fluxes'
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")
lakes <- read.csv("../data/private/Lakes.csv")

fluxes <- merge(fluxes, lakes[,c('Abbreviation','LakeArea_km2')], by.x="Lake", by.y="Abbreviation")
names(fluxes)[grep("Area", names(fluxes))] <- 'larea'

## define the parameters
factors <- c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", "TICumol", 
             "Pressure", "pco2atm", "larea")

## Get distributions for each variables
## ====================================
## basic distro plotting
dplott <- function(vect) {
  plot(density(vect, na.rm = TRUE), ylab = mean(vect, na.rm = TRUE))
}
opar <- par(mfrow = c(2,4))
apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure", "pco2atm")], 2, dplott) 
par(opar)

qqplott <- function(vect) {
  qqnorm(vect, ylab = mean(vect, na.rm = TRUE))
  qqline(vect)
}
opar <- par(mfrow = c(2,4))

apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure", "pco2atm")], 
      2, qqplott)  
par(opar)

apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                   "TICumol", "Pressure", "pco2atm")], 2, shapiro.test)   
# normal: NONE
## see this: http://stats.stackexchange.com/questions/16646/
##    what-is-a-good-index-of-the-degree-of-violation-of-normality-and-what-descriptive

## Use same values as in sensitivity.R, and add larea at the end with a uniform distro
##    going from our min to max
distro <- c("qweibull", "qnorm", "qlogis", "qnorm", "qweibull", "qnorm", "qnorm", "qunif",
            "qunif") 
# in order of 'factors'; retrieve args to list by calling the fit.[] object
props <- list( list(shape=4.02, scale=18.65), list(mean=1013, sd=499), 
               list(location=8.84,scale=0.31), list(mean=4.98, sd=0.83),
               list(shape=2.13, scale=0.68), list(mean=3821, sd=1068),
               list(mean=94.6, sd=0.17), list(min=356.9, max=402.2),
               list(min=0.5, max=500))

## Consider interrelationships between parameters
pairs(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure", "pco2atm")])
## WHAT T ~ meanWind!! (though a lot of noise)
## cond ~ SalCalc ~ TICumol --> these need a rank order combination

## make basic cube
latin <- LHS(gasExchangeK, factors = factors, N = 200, q = distro, q.arg = props, nboot = 50)

## introduce correlation structures
datacorr <- cor(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                          "TICumol", "Pressure", "pco2atm","larea")], method = "spearman",
                use = "complete.obs")
## there were instabilities, let's remove correlation with lake area and these things
##    to test
datacorr[9,] <- 0
datacorr[,9] <- 0

latincorr <- LHS(gasExchangeK, factors = factors, N = 500, q = distro, q.arg = props, 
                 nboot = 200, opts = list(COR = datacorr, eps = 0.1, maxIt=200))

## create an identical one for testing reproducibility (e.g. 200 had low sbma)
testcorr <- LHS(gasExchangeK, factors = factors, N = 500, q = distro, q.arg = props, 
                nboot = 200, opts = list(COR = datacorr, eps = 0.1, maxIt=200))
## test if N is large enough to produce reproducible results...
(testSbma <- sbma(latincorr, testcorr)) # > 90% agreement with eps 0.1, N=500
# NB the default is to use absolute values since -ve correlations get given almost
#   no weight by the method otherwise
targetLHS <- target.sbma(target=0.91, gasExchangeSens, factors, distro, props, 
                         opts = list(COR = datacorr, eps = 0.1, maxIt=200))

## look at diagnostic plots of objects
want <- latincorr

plotecdf(want)
abline(v=0)
plotscatter(want, ylim = c(-300,600), add.lm = FALSE)
plotprcc(want)


## now what happens when our DIC can explode to values like in Brian's data set?
distro <- c("qweibull", "qnorm", "qlogis", "qnorm", "qweibull", "qnorm", "qnorm", "qunif") 
# in order of 'factors'; retrieve args to list by calling the fit.[] object
propsdic <- list( list(shape=4.02, scale=18.65), list(mean=1013, sd=499), 
                  list(location=8.84,scale=0.31), list(mean=4.98, sd=0.83),
                  list(shape=2.13, scale=0.68), list(mean=8000, sd=1068),
                  list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))

diccube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distro, q.arg = propsdic, 
               nboot = 200)

plotscatter(diccube, ylim = c(-300,600))

### split by lake
lakesub <- fluxes[,c("Lake","Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                     "TICumol", "Pressure", "pco2atm")]
lakesplit <- with(lakesub, split(lakesub, list(Lake)))
lakesplit <- lapply(lakesplit, "[", -1) # remove lake before doing correlation matrix
lakesplit <- lakesplit[sapply(lakesplit, function(x) nrow(x) >= 4)] #remove unwanted lakes

lakecorr <- lapply(lakesplit, cor, method = "spearman",
                   use = "complete.obs")

## eyeball distros to guess at what to use now (probs same as for all lakes in most cases)
varnames <- c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
              "TICumol", "Pressure", "pco2atm")
#choose variable
varname <- varnames[1]

#plot chosen variable
ggplot(lakesub, aes_string(x=varname)) + geom_density(aes(colour=Lake, group=Lake)) +
  facet_wrap( "Lake", scales = "free") 

#write wrapper to look at probable distro
giveDistdf <- function(df) { 
  opar <- par(mfrow = c(2,4))
  apply(df[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
              "TICumol", "Pressure", "pco2atm")], 
        2, giveDist) 
  par(opar)
}

#apply to lake, plots in following order:
# "B"  "C"  "D"  "K"  "L"  "P"  "WW"
lapply(lakesplit, giveDistdf)

### play with this:
x <- lakesplit$C$Temperature
disttest <- "gamma"
## chisq may require this to run: start = list(df = 8) or some other number
## t requires same to run but behaves suspiciously (see plots with expected lines)
{if(any(is.na(x))) {
  remove <- which(is.na(x)) 
  fit.test <- fitdist(x[-remove], disttest)
  fit.norm <- fitdist(x[-remove], "norm")
}
  else { fit.test <- fitdist(x, disttest)
  fit.norm <- fitdist(x, "norm")}
}
plot(fit.norm)
plot(fit.test)

fit.norm$aic
fit.test$aic


## lake plots: 
want <- latincorr

plotecdf(want)
abline(v=0)
plotscatter(want)
plotprcc(want)

(testSbma <- sbma(wcube, dcube))
#apply wrapper, sort out rest of columns + colnames
prcclist <- lapply(prcclist, grabs)
prccdf <- do.call(rbind, c(prcclist, make.row.names= FALSE))
varnames <- c("Temperature","Conductivity","pH","meanWindMS","SalCalc","TICumol",
              "Pressure","pco2atm")
prccdf$var <- rep(varnames, times = 4)
prccdf$lake <- rep(c("all","B","D","L","W","P","C","K"), each = 8)
names(prccdf)[c(1:5)] <- c("original", "bias", "stderr", "minci", "maxci")

# plot the prcc's
pd <- position_dodge(0.7)
prccdf$initial <- substr(prccdf$lake, 1, 1)

prccplot <- ggplot(prccdf, aes(x=var, group=interaction(var,lake), colour = lake, 
                               label = initial)) + 
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1, position=pd) +
  geom_jitter(aes(y=original, label=initial), size=3, position=pd) +
  #geom_text(aes(y=original, label=initial), size=2, position=position_dodge(.3)) + 
  scale_color_brewer(palette="Paired") +
  ylim(-1,1) +
  ylab("PRCC")
prccplot
## FIXME: can't get this nice into r markdown...!

## plot only the main one
prcclist <- list(latincorr)
prcclist <- lapply(prcclist, grabs)
prccdf <- do.call(rbind, c(prcclist, make.row.names= FALSE))
varnames <- c("Temperature","Conductivity","pH","Wind","Salinity","DIC",
              "Pressure","Air pCO2")
prccdf$var <- varnames
names(prccdf)[c(1:5)] <- c("original", "bias", "stderr", "minci", "maxci")

prccdf$var  <- factor(prccdf$var, levels = (prccdf$var)[rev(order(abs(prccdf$original)))])

prccplot <- ggplot(prccdf, aes(x=var)) + 
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1) +
  geom_point(aes(y=original)) +
  #geom_text(aes(y=original, label=initial), size=2, position=position_dodge(.3)) + 
  #scale_color_brewer(palette="Paired") +
  theme_bw(base_size = 15) +
  theme(axis.text.x=element_text(angle=45, vjust=0.5), axis.title.x=element_blank()) +
  ylim(-1,1) +
  ylab("Correlation")
ggsave('../docs/private/prccplot.png', prccplot, width=8, height=4, units = 'in')

## general lake diffs in variables
melted <- melt(lakesub, id = "Lake")
melted$variable <- factor(melted$variable, levels = c("Temperature", "Conductivity", "pH", 
                                                      "SalCalc", "TICumol", "meanWindMS", 
                                                      "Pressure", "pco2atm", "larea"))

# create a list with strip labels
varnames <- list(
  'Temperature'=expression(paste("Temperature ("~degree*"C)")) ,
  'Conductivity'= expression(paste("Conductivity ("*mu*"S"~"cm"^{-1}*")")),
  'pH'="pH",
  'SalCalc' = "Salinity (ppt)",
  'TICumol' = expression(paste("DIC ("~mu*"mol"~"L"^{-1}*")")),
  'meanWindMS'= expression(paste("Mean Wind (m"~"s"^{-1}*")")),
  'Pressure' = 'Pressure (kPa)',
  'pco2atm' = expression(paste("Air"~italic(p)*"CO"[2]~"(ppm)"))
)

# Create a 'labeller' function, and push it into facet_grid call:
var_labeller <- function(variable,value){
  return(varnames[value])
}


## create ggplot of the simulated data
simdat <- cbind(latincorr$data, latincorr$res)
names(simdat)[names(simdat) == "cube$res"] <- "flux"
simmelt <- melt(simdat, id.vars = "flux")
simmelt$variable <- factor(simmelt$variable, levels = c("Temperature", "Conductivity", "pH", 
                                                        "SalCalc", "TICumol", "meanWindMS", 
                                                        "Pressure", "pco2atm"))
simplot <- ggplot(simmelt, aes(x=value,y=flux, group=variable)) +
  theme_bw(base_size = 12, base_family = 'Arial') +
  geom_point(size=0.5) +
  facet_wrap( "variable", scales = "free", labeller = var_labeller, ncol=2) +
  ylab(expression(paste(CO[2]~"flux (mmol"~"C "*m^{-2}*"d"^{-1}*')'))) +
  theme(axis.title.x = element_blank())

# save meltplots?
savemelt <- TRUE
if(savemelt) {
  ggsave("../docs/private/sensitivity-sims.pdf", simplot, width=12, height=15, units="cm")
}
## save LHS's?
saveLHS <- TRUE

if(saveLHS) {
  saveRDS(lareacube, file="../data/private/LHSall-lakes.rds")
  }

## save plots?
saveplots <- TRUE

if(saveplots) {
  saveRDS(prccplot, file = "../data/private/prccplot.rds")
  }
