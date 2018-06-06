## sensitivity analysis as to how much pH dominates the CO2 flux of our study sites
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

runextras <- FALSE

## source necessary packages:
library("mc2d")
library("pse")
library("fitdistrplus") # loads MASS too
library("reshape2")
library("extrafont")

## source the function we want to decompose
source("../functions/gasExchangeSensitivity.R")

## set default for plot
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## source files with parameters
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")

## define the parameters
factors <- c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", "TICumol", 
             "Pressure", "pco2atm")

## Get distributions for each variables
## ====================================
dplott <- function(vect) {
  plot(density(vect, na.rm = TRUE), ylab = mean(vect, na.rm = TRUE))
}
qqplott <- function(vect) {
  qqnorm(vect, ylab = mean(vect, na.rm = TRUE))
  qqline(vect)
}

if (runextras) {
  ## basic distro plotting
  opar <- par(mfrow = c(2,4))
  apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                  "TICumol", "Pressure", "pco2atm")], 2, dplott) 
  par(opar)
  
  opar <- par(mfrow = c(2,4))
  apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                  "TICumol", "Pressure", "pco2atm")], 
        2, qqplott)  
  par(opar)
  apply(fluxes[,c(7:8,11,13,17:18,22)], 2, shapiro.test)   
  # normal: NONE
  ## see this: http://stats.stackexchange.com/questions/16646/
  ##    what-is-a-good-index-of-the-degree-of-violation-of-normality-and-what-descriptive
}
## So let's try the packages on fitting distributions
giveDist <- function(vect) {
  if(any(is.na(vect))) {
    remove <- which(is.na(vect))
    descdist(vect[-remove], discrete = FALSE) }
  else {descdist(vect, discrete = FALSE)}
  legend("bottomleft", legend = head(vect))
}
if(runextras) {
  opar <- par(mfrow = c(2,4))
  apply(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                  "TICumol", "Pressure", "pco2atm")], 
        2, giveDist) 
  par(opar)
  ## --> normal:  dicumol, salcalc, (meanwind)
  ## --> unif: cond, pressure, pco2atm
  ## normal-logistic: pH               
  ## bleurgh: temperature
  
  ## Look at uncertain-distro variables and test with likely distro 
  ## "For the following named distributions, reasonable starting values will be computed 
  ##    if start is omitted (i.e. NULL) : "norm", "lnorm", "exp" and "pois", "cauchy", 
  ##    "gamma", "logis", "nbinom" (parametrized by mu and size), "geom", "beta", 
  ##    "weibull" from the stats package; "invgamma", "llogis", "invweibull", "pareto1", 
  ##    "pareto" from the actuar package. Note that these starting values may not be 
  ##    good enough if the fit is poor. The function uses a closed-form formula to fit the 
  ##    uniform distribution. If start is a list, then it should be a named list with the 
  ##    same names as in the d,p,q,r functions of the chosen distribution. 
  ##    If start is a function of data, then the function should return a named list with 
  ##    the same names as in the d,p,q,r functions of the chosen distribution. 
  x <- fluxes$Conductivity
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
  
  plot(fit.test)
  plot(fit.norm)
  
  ## lower = better, not in the absolute value sense
  fit.test$aic
  fit.norm$aic
  ## for temperature, weibull > normal > chisq | gamma --> weibull
  ## for pH, logis > normal > cauchy
  ## for cond, pressure: unif doesn't produce AIC but by eye, normal is better
  ## for DICumol, normal > logis > cauchy
  ## looking at dist plot, normal for windspeed
  ## for SalCalc, weibull > normal
  ## ngggg let's do unif for pco2atm
}
## Continue defining parameters for the chosen distros
distro <- c("qweibull", "qnorm", "qlogis", "qnorm", "qweibull", "qnorm", "qnorm", "qunif") 
# in order of 'factors'; retrieve args to list by calling the fit.[] object
props <- list( list(shape=4.02, scale=18.65), list(mean=1013, sd=499), 
               list(location=8.84,scale=0.31), list(mean=4.98, sd=0.83),
               list(shape=2.13, scale=0.68), list(mean=3821, sd=1068),
               list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))

## Consider interrelationships between parameters
pairs(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                "TICumol", "Pressure", "pco2atm")])
## WHAT T ~ meanWind!! (though a lot of noise)
## cond ~ SalCalc ~ TICumol --> these need a rank order combination

## make basic cube
if(runextras) {
  latin <- LHS(gasExchangeSens, factors = factors, N = 200, q = distro, q.arg = props, nboot = 50)
}
## introduce correlation structures
datacorr <- cor(fluxes[,c("Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                          "TICumol", "Pressure", "pco2atm")], method = "spearman",
                use = "complete.obs")

latincorr <- LHS(gasExchangeSens, factors = factors, N = 500, q = distro, q.arg = props, 
                 nboot = 200, opts = list(COR = datacorr, eps = 0.1, maxIt=200))
if(runextras) {
  ## create an identical one for testing reproducibility (e.g. 200 had low sbma)
  testcorr <- LHS(gasExchangeSens, factors = factors, N = 500, q = distro, q.arg = props, 
                  nboot = 200, opts = list(COR = datacorr, eps = 0.1, maxIt=200))
  ## test if N is large enough to produce reproducible results...
  (testSbma <- sbma(latincorr, testcorr)) # > 90% agreement with eps 0.1, N=500
  # NB the default is to use absolute values since -ve correlations get given almost
  #   no weight by the method otherwise
  targetLHS <- target.sbma(target=0.91, gasExchangeSens, factors, distro, props, 
                           opts = list(COR = datacorr, eps = 0.1, maxIt=200))
}
## look at diagnostic plots of objects
want <- latincorr

plotecdf(want)
abline(v=0)
plotscatter(want, ylim = c(-300,600), add.lm = FALSE)
plotprcc(want)

if(runextras){
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
}
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
if(runextras){
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
}

## lake-specific values
distroc <- c("qnorm", "qnorm", "qlogis", "qnorm", "qlogis", "qnorm", 
             "qnorm", "qunif") 
propsc <- list( list(mean=17.6, sd=4.93), list(mean=1182.2, sd=200), 
                list(location=8.79,scale=0.3), list(mean=4.97, sd=0.88),
                list(location=0.73, scale=0.05), list(mean=4391.4, sd=733.6),
                list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))


distrok <- c("qnorm", "qnorm", "qlogis", "qnorm", "qlogis", "qnorm", 
             "qnorm", "qunif") 
propsk <- list( list(mean=16.2, sd=5.67), list(mean=1241.3, sd=207), 
                list(location=8.94,scale=0.33), list(mean=4.96, sd=0.82),
                list(location=0.72, scale=0.05), list(mean=4171, sd=721.1),
                list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))


distrop <- c("qweibull", "qnorm", "qnorm", "qnorm", "qnorm", "qnorm", 
             "qnorm", "qnorm") 
propsp <- list( list(shape=3.93, scale=18.6), list(mean=1196.5, sd=272.7), 
                list(mean=8.84,sd=0.41), list(mean=4.98, sd=0.8),
                list(mean=0.72, sd=0.15), list(mean=4054.3, sd=608.6),
                list(mean=94.6, sd=0.16), list(mean=386.8, sd=7.09))

distrob <- c("qweibull", "qnorm", "qlogis", "qnorm", "qnorm", "qnorm", 
             "qnorm", "qunif") 
propsb <- list( list(shape=5.33, scale=19.6), list(mean=471.8, sd=83.2), 
                list(location=8.78,scale=0.32), list(mean=4.98, sd=0.83),
                list(mean=0.26, sd=0.04), list(mean=2719.1, sd=636.1),
                list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))


distrod <- c("qweibull", "qlogis", "qlogis", "qnorm", "qnorm", "qnorm", 
             "qnorm", "qunif") 
propsd <- list( list(shape=3.39, scale=17), list(location=359.6, scale=36.1), 
                list(location=8.74,scale=0.29), list(mean=4.98, sd=0.84),
                list(mean=0.23, sd=0.03), list(mean=2926.2, sd=720.5),
                list(mean=94.6, sd=0.17), list(min=356.9, max=402.2))

distrol <- c("qweibull", "qweibull", "qlogis", "qnorm", "qnorm", "qnorm", 
             "qnorm", "qunif") 
propsl <- list( list(shape=4.31, scale=18.34), list(shape=8.52, scale=1864), 
                list(location=8.85,scale=0.27), list(mean=4.93, sd=0.82),
                list(mean=1.07, sd=0.11), list(mean=4827.8, sd=930.6),
                list(mean=94.6, sd=0.16), list(min=358.9, max=402.2))


distrow <- c("qweibull", "qnorm", "qnorm", "qnorm", "qweibull", "qnorm", 
             "qnorm", "qunif") 
propsw <- list( list(shape=4.97, scale=19.58), list(mean=906, sd=295), 
                list(mean=9.01,sd=0.67), list(mean=5.07, sd=0.81),
                list(shape=3.54, scale=0.57), list(mean=3693, sd=1071),
                list(mean=94.6, sd=0.17), list(min=361.1, max=402.2))

## lake cubes:
bcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrob, q.arg = propsb, 
             nboot = 200, opts = list(COR = lakecorr$B, eps = 0.1, maxIt=200))
dcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrod, q.arg = propsd, 
             nboot = 200, opts = list(COR = lakecorr$D, eps = 0.1, maxIt=200))

lcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrol, q.arg = propsl, 
             nboot = 200, opts = list(COR = lakecorr$L, eps = 0.1, maxIt=200))
wcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrow, q.arg = propsw, 
             nboot = 200, opts = list(COR = lakecorr$WW, eps = 0.1, maxIt=200))
ccube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distroc, q.arg = propsc, 
             nboot = 200, opts = list(COR = lakecorr$C, eps = 0.1, maxIt=200))
pcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrop, q.arg = propsp, 
             nboot = 200, opts = list(COR = lakecorr$P, eps = 0.1, maxIt=200))
kcube <- LHS(gasExchangeSens, factors = factors, N = 500, q = distrok, q.arg = propsk, 
             nboot = 200, opts = list(COR = lakecorr$K, eps = 0.1, maxIt=200))

## lake plots: 
want <- latincorr

plotecdf(want)
abline(v=0)
plotscatter(want)
plotprcc(want)

(testSbma <- sbma(wcube, dcube))

## grab prcc's for a group plot:
prcclist <- list(latincorr, bcube, dcube, lcube, wcube, pcube, ccube, kcube)

# wrapper to grab the df from each object (which is a list)
grabs <- function(list) {
  grabdf <- as.data.frame(list$prcc[[1]]['PRCC'])
}

#apply wrapper, sort out rest of columns + colnames
prcclist <- lapply(prcclist, grabs)
prccdf <- do.call(rbind, c(prcclist, make.row.names= FALSE))
varnames <- c("Temperature","Conductivity","pH","meanWindMS","SalCalc","TICumol",
              "Pressure","pco2atm")
prccdf$var <- rep(varnames, times = 4)
prccdf$lake <- rep(c("all","B","D","L","W","P","C","K"), each = 8)
names(prccdf)[c(1:5)] <- c("original", "bias", "stderr", "minci", "maxci")

write.csv(prccdf,"../data/private/lhs-prcc.csv", row.names = FALSE)
# plot the prcc's
pd <- position_dodge(0.7)
prccdf$initial <- substr(prccdf$lake, 1, 1)

prccplot <- ggplot(prccdf, aes(x=var, group=interaction(var,lake), colour = lake)) + 
  geom_errorbar(aes(ymin=minci, ymax=maxci), width=.1, position=pd) +
  geom_jitter(aes(y=original), size=3, position=pd) +
  #geom_text(aes(y=original, label=initial), size=2, position=position_dodge(.3)) + 
  scale_color_brewer(palette="Paired") +
  ylim(-1,1) +
  ylab("PRCC")
prccplot

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

## general lake diffs in variables
# Create a 'labeller' function, and push it into facet_grid call:
var_labeller <- function(variable,value){
  return(varnames[value])
}


if(runextras) {
  melted <- melt(lakesub, id = "Lake")
  melted$variable <- factor(melted$variable, levels = c("Temperature", "Conductivity", "pH", 
                                                        "SalCalc", "TICumol", "meanWindMS", 
                                                        "Pressure", "pco2atm"))
  
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
  
  # run the ggplot
  meltplot <- ggplot(melted, aes(x=Lake,y=value, group=Lake)) +
    geom_boxplot(outlier.colour="black", outlier.shape=5,
                 outlier.size=1) +
    theme_bw(base_size = 14, base_family = 'Arial') +
    facet_wrap( "variable", scales = "free", labeller = var_labeller, ncol=2) +
    theme(axis.title = element_blank())
  meltplot
}
## FIXME: apparently var_labeller style labeling is deprecated but I can't get anything else to work
##    for the expressions so going with this for now... new version wants df's but df doesn't accept
##    expressions....
# vardf <- c(expression(paste("Temperature ("~degree*"C)")),
#                               expression(paste("Conductivity ("*mu*"S"~"cm"^{-1}*")")),
#                               "pH", expression(paste("Mean Wind (m"~"s"^{-1}*")")),
#                               "Salinity (ppt)", expression(paste("DIC ("~mu*"mol"~"L"^{-1}*")")),
#                               'Pressure (kPa)', expression(paste(italic(p)*"CO"[2]~"(ppm)")))
# as_labeller(vardf)  

## create ggplot of the simulated data
simdat <- cbind(latincorr$data, latincorr$res)
names(simdat)[names(simdat) == "latincorr$res"] <- "flux"
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

## save LHS's?
saveLHS <- TRUE

if(saveLHS) {
  saveRDS(latincorr, file="../data/private/LHSall-lakes.rds")
  saveRDS(bcube, file="../data/private/LHS-B.rds")
  saveRDS(dcube, file="../data/private/LHS-D.rds")
  saveRDS(lcube, file="../data/private/LHS-L.rds")
}

## save plots?
saveplots <- TRUE

if(saveplots) {
  ggsave("../docs/private/sensitivity-sims.pdf", simplot, width=12, height=15, units="cm")
  ggsave('../docs/private/prccplot.png', prccplot, width=8, height=4, units = 'in')
  saveRDS(prccplot, file = "../data/private/prccplot.rds")
}
