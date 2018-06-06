## following from regression_routines -pH and -models, this script creates lake-
##    specific models following meeting with Peter 25.01.2015
## final models chosen were co2mod and its residuals as resmodred for CO2 ~.
##    and egmodred.2 for pH ~ . These models are in this script run, the others
##    are set to not run
## UPDATE: archive-preparation and model changes: changed GPP to R since that works much better
##      with our hypotheses by incorporating at least one real respiration proxy
##    Also changing the final model, since after further source data changes oxygen is 
##      significant also with the climate indices
##    NOTE that TDN and DOC interfere with each other (check pairs() with logs); needs
##      discussing and therefore new TDN model is developed for paper

## set to whether or not to run the discarded models, and whether or not to plot output
##    from models that were kept
runextras <- FALSE
plotmods <- FALSE
runwasc <- FALSE

## load necessary packages
library("mgcv")
library("ggplot2")

## load necessary data
if (!file.exists("../data/private/regvars.rds")) {
  source("../scripts/regression_routines.R")
}
regvars <- readRDS("../data/private/regvars.rds")

if (!file.exists("../data/weathers.rds")) {
  source("../scripts/climate-weather-modeling.R")
}
weathers <- readRDS('../data/weathers.rds')

regvars <- merge(regvars, weathers)
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re

## In fact, good to log some of the more skewed distributions
## change -ve values to 0, then add 1
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))

## split regvarf by lake then remove all NAs from list
regsplit <- with(regvarf, split(regvarf, list(Lake)))
regsplit <- regsplit[lapply(regsplit,nrow) > 1]

## use final models created for CO2, adapted to see if any one Lake deviates
##    from the global trend
## =================================================================================
## 1. model most similar to ones done before
if(runextras) {
  egmod <- gam(pH_surface ~ 
                 s(Chl_a_ug_L) + s(Chl_a_ug_L, by = Lake, m = 1) +
                 # m = 1 means that penalty goes with 1st derivative... required
                 #    here because we're doing both s() and s(..,by = ...)
                 s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                 s(TDN_ug_L) + s(TDN_ug_L, by = Lake, m = 1) + 
                 s(DOC_mg_L) + s(DOC_mg_L, by = Lake, m = 1) +
                 s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                 te(PDO, SOI) + te(PDO, SOI, by = Lake, m = 1) +
                 s(Lake, Year, bs = "re"), # still need to explicitly keep Lake as re!!!
               data = regvarf,
               select = FALSE, method = "REML", family = gaussian,
               na.action = na.exclude,
               control = gam.control(nthreads = 3, trace = TRUE))
  
  ## save summary as txt document
  egmodsum <- summary(egmod)
  sink("../data/private/egmodsummary.txt")
  egmodsum
  sink()
}

## 2. for this model, we are removing some of the Lake-specificity where nothing signif
##    was found in previous model + using regvarf2
if(runextras) {
  egmod.red <- gam(pH_surface ~ 
                     s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                     s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                     s(log10(TDN_ug_L)) + 
                     s(log10(DOC_mg_L)) +
                     s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                     te(PDO, SOI) +
                     s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator),
                   data = regvarf2,
                   select = TRUE, method = "REML", family = gaussian,
                   na.action = na.exclude,
                   control = gam.control(nthreads = 3, trace = TRUE))
  # this setting makes the computation a bit more efficient 
  ## longer object length is not a multiple of shorter object length!!??
  
  egmodredsum <- summary(egmod.red)
  sink("../data/private/egmodredsummary.txt")
  egmodredsum
  sink()
}

# 3. model that removes further unnecessary elements to make fitting faster and simpler
#    AND introduces scat() again. Testing to see if the step failure can be averted
#    by setting maxHalf = 60
if(runextras){
  egmod.red2 <- gam(pH_surface ~
                      s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1,k=5) +
                      s(GPP_h) +
                      s(log10(TDN_ug_L)) +
                      s(log10(DOC_mg_L)) +
                      s(Oxygen_ppm) +
                      te(PDO, SOI) +
                      s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator
                    data = regvarf2,
                    select = TRUE, method = "REML", family = scat(),
                    na.action = na.exclude,
                    control = gam.control(nthreads = 3, trace = TRUE,
                                          newton = list(maxHalf = 60)))
  saveRDS(egmod.red2, "../data/private/egmodred2.rds")
  egmodred2sum <- summary(egmod.red2)
  
  sink("../docs/private/egmodred2summary.txt")
  egmodred2sum
  sink()
}
# 4. model that uses lagged SOI and PDO, and includes SPEI (see climate-weather-modeling.R)
#  Not here R instead (models with GPP and R still conclude the same -- it's not signif,
#   and other var effects don't change)
egmodlagged <- gam(pH_surface ~
                     s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1,k=5) +
                     s(R_h) +
                     s(log10(TDN_ug_L)) +
                     s(log10(DOC_mg_L)) +
                     s(Oxygen_ppm) +
                     s(SPEI02) +
                     te(PDOmean, SOImean) + # tested this and we need interaction
                     s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator
                   data = regvarf2,
                   select = TRUE, method = "REML", family = scat(),
                   na.action = na.exclude,
                   control = gam.control(nthreads = 3, trace = TRUE,
                                         newton = list(maxHalf = 60)))
#egmodlaggedti <- update(egmodlagged, . ~ . - te(PDOmean, SOImean) + ti(PDOmean) + ti(SOImean) 
#   + ti(PDOmean, SOImean))
#egmodlaggedtimarg <- update(egmodlagged, . ~ . - te(PDOmean, SOImean) 
#       + ti(PDOmean) + ti(SOImean))
#anova(egmodlaggedti, egmodlaggedtimarg, test = "LRT")
saveRDS(egmodlagged,"../data/private/egmodlagged.rds")

# take out insignificant vars
egmodlaggedsimp <- update(egmodlagged, . ~ . -s(R_h) - s(log10(TDN_ug_L)) -s(SPEI02) +
                            s(SPEI02, k=3)) #became squiggly
summary(egmodlaggedsimp)
saveRDS(egmodlaggedsimp,"../data/private/egmodlaggedsimp.rds")

# see mumdat below for why this created
egmodlaggedtdn <- update(egmodlagged, . ~ . -s(R_h) - s(log10(DOC_mg_L)) -s(SPEI02) +
                            s(SPEI02, k=3)) #became squiggly
summary(egmodlaggedtdn)
saveRDS(egmodlaggedtdn,"../data/private/egmodlaggedtdn.rds")

nullvars <- c("Chl_a_ug_L" , "Chl_a_ug_L", "R_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", 
              "SPEI02", "PDOmean", "SOImean", "Lake", "Year", "pH_surface")
mumdat <- regvarf2[complete.cases(regvarf2[,nullvars]),] # want to AIC compare, and NA diff
#   between vars
mummod <- update(egmodlagged, data=mumdat, na.action="na.fail")
mummodtdn <- update(mummod, .~.-s(log10(DOC_mg_L)))
mummoddoc <- update(mummod, .~.-s(log10(TDN_ug_L)))
AIC(mummod, mummoddoc, mummodtdn) # basically we can drop tdn with no added AIC,
#   but dropping DOC in favour of TDN makes it worse.

## testing to makes sure we really want SOI and PDO as te()
##  however couldn't get nointer to stabilise unless family was set to gaussian so
##      did this in order to compare.
if(runextras) {
  nointer <- gam(pH_surface ~ 
                   s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1,k=5) +
                   s(GPP_h) +
                   s(log10(TDN_ug_L)) +
                   s(log10(DOC_mg_L)) +
                   s(Oxygen_ppm) +
                   ti(PDO) + ti(SOI) +
                   s(Year, Lake, bs = "re", by = dummy), 
                 data = regvarf2,
                 select = TRUE, method = "REML", family = gaussian(),
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE,
                                       newton = list(maxHalf = 150)))
  
  inter <- gam(pH_surface ~ 
                 s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1,k=5) +
                 s(GPP_h) +
                 s(log10(TDN_ug_L)) +
                 s(log10(DOC_mg_L)) +
                 s(Oxygen_ppm) +
                 ti(PDO) + ti(SOI) + ti(PDO, SOI) +
                 s(Year, Lake, bs = "re", by = dummy), 
               data = regvarf2,
               select = TRUE, method = "REML", family = gaussian(),
               na.action = na.exclude,
               control = gam.control(nthreads = 3, trace = TRUE,
                                     newton = list(maxHalf = 60)))
  anova(nointer, inter, test = "LRT")
}
## plot output of model
if (plotmods) {
  pdf("../docs/private/egmod-reduced-log-covariates.pdf")
  op <- par(mar = c(4,4,1,1) + 0.1)
  plot(egmod.red2, pages = 4, scheme = 2)
  par(op)
  dev.off()
}

## use final models created for CO2, adapted to see if any one Lake deviates
##    from the global trend
## =================================================================================
## 1. can we use a blanket model for CO2 based on pH?
co2mod <- gam(co2Flux ~ 
                s(pH_surface) + s(pH_surface, by = Lake, m = 1) +
                # m = 1 means that penalty goes with 1st derivative... required
                #    here because we're doing both s() and s(..,by = ...)
                s(Lake, Year, bs = "re"), # still need to explicitly keep Lake as re!!!
              data = regvarf,
              select = TRUE, method = "REML", family = gaussian,
              na.action = na.exclude,
              control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))

co2modsum <- summary(co2mod)
if (plotmods) {
  sink("../data/private/co2modsummary.txt")
  co2modsum
  sink()
}

gam.check(co2mod)
saveRDS(co2mod, "../data/private/co2mod.rds")

## FIXME: the kindex and pvalues are still NA for this model too

if(runextras) {
  co2modnull <- gam(co2Flux ~ 
                      s(pH_surface) +
                      s(Lake, Year, bs = "re"), 
                    data = regvarf,
                    select = TRUE, method = "REML", family = gaussian,
                    na.action = na.exclude,
                    control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
  
  saveRDS(co2modnull, '../data/private/co2modnull.rds')
  anova(co2modnull, co2mod, test = "LRT")
  AIC(co2modnull, co2mod)
  # --> we need to take the more complex model
}

## plot output of model
if (plotmods) {
  pdf("../docs/co2mod.pdf")
  op <- par(mar = c(4,4,1,1) + 0.1)
  plot(co2mod, pages = 4, scheme = 2)
  par(op)
  dev.off()
}

## model residuals
res <- resid(co2mod, type = "pearson")

if(runextras) {
  resmod <- gam(res ~ 
                  s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                  s(GPP_h) + s(GPP_h, by = Lake, m = 1) +
                  s(log10(TDN_ug_L)) + s(log10(TDN_ug_L), by = Lake, m = 1) +
                  s(log10(DOC_mg_L)) + s(log10(DOC_mg_L), by = Lake, m = 1) +
                  s(Oxygen_ppm) + s(Oxygen_ppm, by = Lake, m = 1) +
                  te(PDO, SOI) +
                  s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator),
                data = regvarf2,
                select = TRUE, method = "REML", family = scat(),
                na.action = na.exclude,
                control = gam.control(nthreads = 3, trace = TRUE,
                                      newton = list(maxHalf = 60)))
  
  co2resmodsum <- summary(resmod)
  sink("../data/private/co2resmodsummary.txt")
  co2resmodsum
  sink()
  
  pdf("../docs/resmod-full-covariates.pdf")
  op <- par(mar = c(4,4,1,1) + 0.1)
  plot(resmod, pages = 5, scheme = 2)
  par(op)
  dev.off()
}

## some of the GPP lake-specific things look a bit far-fetched so took Lake away there.
if (runextras) {
resmodred <- gam(res ~ 
                   s(log10(Chl_a_ug_L)) + s(log10(Chl_a_ug_L), by = Lake, m = 1) +
                   s(GPP_h) +
                   s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) + 
                   s(Oxygen_ppm) + 
                   te(PDO, SOI) +
                   s(Lake, Year, bs = "re", by = dummy), # dummy is 0/1 indicator
                 data = regvarf2,
                 select = TRUE, method = "REML", family = scat(),
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE,
                                       newton = list(maxHalf = 60)))

  anova(resmodred, resmod, test = "LRT")
  ## the simpler model is just fine

saveRDS(resmodred, "../data/private/resmodred.rds")

co2resmodredsum <- summary(resmodred)

if (plotmods) {
  sink("../data/private/co2resmodredsummary.txt")
  co2resmodredsum
  ## s(log10(TDN_ug_L))          2.931e+00      9 25.756 9.97e-06 ***
  sink()
  
  pdf("../docs/resmodred.pdf")
  op <- par(mar = c(4,4,1,1) + 0.1)
  plot(resmodred, pages = 3, scheme = 2)
  par(op)
  dev.off()
  
  plot(resmodred, pages = 3, scheme = 2)
}}
## can we simplify even further?
if (runextras) {
  resmodredchl <- gam(res ~ 
                        s(log10(Chl_a_ug_L)) +
                        s(GPP_h) +
                        s(log10(TDN_ug_L)) + 
                        s(log10(DOC_mg_L)) + 
                        s(Oxygen_ppm) + 
                        te(PDO, SOI) +
                        s(Lake, Year, bs = "re"),
                      data = regvarf2,
                      select = TRUE, method = "REML", family = scat(),
                      na.action = na.exclude,
                      control = gam.control(nthreads = 3, trace = TRUE,
                                            newton = list(maxHalf = 60)))
  
  anova(resmodredchl, resmodred, test = "LRT")
  ## better to keep the more complicated model
}

## ==========================================================================
## RUN MODELS FOR WASCANA ONLY -- USE regvarf2 WITH <0 made to 0, +1
## ==========================================================================
if(runwasc) {
  ww <- subset(regvarf2, Lake == "WW") # 168 data points
  
  ## 1. Create models and output
  ## ====================================
  
  ## model CO2 on pH with Year as factor
  co2wwmod <- gam(co2Flux ~ 
                    s(pH_surface) +
                    s(Year, bs = "re"), 
                  data = ww,
                  select = TRUE, method = "REML", family = gaussian,
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, newton = list(maxHalf = 60), trace = TRUE))
  saveRDS(co2wwmod, "../data/private/co2wwmod.rds")
  
  ## save summary as txt document
  co2wwsum <- summary(co2wwmod)
  sink("../data/private/co2wwmodsummary.txt")
  co2wwsum
  sink()
  
  ## model with logged chla and TDN and Year as factor
  phwwmod <- gam(pH_surface ~ 
                   s(log10(Chl_a_ug_L)) + s(GPP_h) + s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) + s(Oxygen_ppm) + te(PDO, SOI) +
                   s(Year, bs = "re"), 
                 data = ww,
                 select = TRUE, method = "REML", family = gaussian,
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE))
  # a note regarding select = T | F: "select = FALSE just fits the model 
  #   without the second penalty on the linear part of the function. It seems to want to 
  #   keep some small amount of curvature but put a little shrinkage into the null space 
  #   just like the lasso would shrink a term away from the least squares fit." -GS
  saveRDS(phwwmod, "../data/private/phwwmod.rds")
  
  
  ## save summary as txt document
  phwwsum <- summary(phwwmod)
  sink("../data/private/phwwmodsummary.txt")
  phwwsum
  sink()
  
  ## one more model with co2 against everything
  dummymod <- gam(co2Flux ~ 
                    s(log10(Chl_a_ug_L)) + s(GPP_h) + s(log10(TDN_ug_L)) + 
                    s(log10(DOC_mg_L)) + s(Oxygen_ppm) + te(PDO, SOI) +
                    s(Year, bs = "re"), 
                  data = ww,
                  select = TRUE, method = "REML", family = gaussian,
                  na.action = na.exclude,
                  control = gam.control(nthreads = 3, trace = TRUE))
  
  saveRDS(dummymod, "../data/private/dummymod.rds")
  
  ## save summary as txt document
  dummysum <- summary(dummymod)
  sink("../data/private/dummysummary.txt")
  dummysum
  sink()
  
  ## AAANNND looking at SOI and PDO separately
  nointer <- gam(pH_surface ~ 
                   s(log10(Chl_a_ug_L)) + s(GPP_h) + s(log10(TDN_ug_L)) + 
                   s(log10(DOC_mg_L)) + s(Oxygen_ppm) +
                   ti(PDO) + ti(SOI) +
                   s(Year, bs = "re"), 
                 data = ww,
                 select = TRUE, method = "REML", family = gaussian,
                 na.action = na.exclude,
                 control = gam.control(nthreads = 3, trace = TRUE))
  
  inter <- gam(pH_surface ~ 
                 s(log10(Chl_a_ug_L)) + s(GPP_h) + s(log10(TDN_ug_L)) + 
                 s(log10(DOC_mg_L)) + s(Oxygen_ppm) + 
                 ti(PDO) + ti(SOI) + ti(PDO, SOI) +
                 s(Year, bs = "re"), 
               data = ww,
               select = TRUE, method = "REML", family = gaussian,
               na.action = na.exclude,
               control = gam.control(nthreads = 3, trace = TRUE))
  anova(nointer, inter, test = "LRT")
  
  # a note regarding select = T | F: "select = FALSE just fits the model 
  #   without the second penalty on the linear part of the function. It seems to want to 
  #   keep some small amount of curvature but put a little shrinkage into the null space 
  #   just like the lasso would shrink a term away from the least squares fit." -GS
  
  
  ## 2. Predict responses with most interesting terms and plot
  ## =========================================================
  
  ## predict CO2 based on pH: keep using iterms
  phwant <- with(ww, seq(min(`pH_surface`, na.rm = TRUE),
                         max(`pH_surface`, na.rm = TRUE),
                         length = N))
  
  ww.pdatc <- data.frame(pH_surface = phwant, Year = 2004) 
  ww.predc <- predict(co2wwmod, newdata = ww.pdatc, type = "iterms") 
  whichColsc <- grep("pH", colnames(ww.predc))
  
  ww.predcse <- predict(co2wwmod, newdata = ww.pdatc, type = "iterms", se.fit = TRUE)
  ww.predcse <- as.data.frame(ww.predcse$se.fit)
  whichColscse <- grep("pH", colnames(ww.predcse))
  
  ww.pdatc <- cbind(ww.pdatc, Fitted = ww.predc[, whichColsc], Fittedse = ww.predcse[,whichColscse])
  ww.pdatc <- with(ww.pdatc, transform(ww.pdatc, Fittedplus = Fitted + Fittedse))
  ww.pdatc <- with(ww.pdatc, transform(ww.pdatc, Fittedminus = Fitted - Fittedse))
  
  ## transform back into original pH value for clarity:
  # get mean/intercept pH used by the model
  shiftco2 <- attr(predict(co2wwmod, newdata = ww.pdatc, type = "iterms"), "constant")
  ww.pdatcnorm <- ww.pdatc
  ww.pdatcnorm <- with(ww.pdatcnorm, transform(ww.pdatcnorm, Fitted = Fitted + shiftco2, 
                                               Fittedplus = Fittedplus + shiftco2, 
                                               Fittedminus = Fittedminus + shiftco2))
  
  labdat1 <- data.frame(x = 7, y = -15, label = "mean flux: -23")
  # without this step, the resolution of the text is really off for some reason
  
  ## plot with all labels
  pdf("../data/private/wwco2mod.pdf", width = 10, height = 6.7)
  ggplot(ww.pdatcnorm, aes(x = pH_surface, y = Fitted)) +
    geom_line() + 
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25) +  
    geom_abline(slope = 0, intercept = shiftco2, linetype="dotted") +
    geom_text(data = labdat1, aes(label = label, x = x, y = y, size = 5), 
              show.legend = FALSE) +
    xlab('pH') + ylab('CO2 flux (mmolC/m2/d)')
  dev.off()
  
  ## plot with disabled labels and white background
  tiff("../data/private/wwco2mod-notext.tiff", width = 2.3, height = 1.5, 
       units = "in", res = 300)
  ggplot(ww.pdatcnorm, aes(x = pH_surface, y = Fitted)) +
    geom_line() + 
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25) +  
    theme_bw() +
    theme(axis.text = element_blank(), title = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  dev.off()
  
  ## extract predicted values based on chlorophyll for Matt's paper
  ## we can use raw, not logged, chl values, since the gam retains this information
  ##    and therefore we can feed it raw to predict() !!!
  N <- 200
  varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
  lakeWbar <- data.frame(t(colMeans(ww[, varWant], na.rm = TRUE)))
  lakeWbar$Year <- 2004
  
  ww.pdat <- with(droplevels(ww),
                  data.frame(`Chl_a_ug_L` = rep(seq(min(`Chl_a_ug_L`, na.rm = TRUE),
                                                    max(`Chl_a_ug_L`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
  ww.pdat <- merge(ww.pdat, lakeWbar)
  
  
  ## predict pH based on chlorophyll alone;
  ## iterms means that uncertainty in intercept is included; especially good for when
  ##   y axis units transformed back into original scale
  ww.pred <- predict(phwwmod, newdata = ww.pdat, type = "iterms")
  ww.predse <- predict(phwwmod, newdata = ww.pdat, type = "iterms", se.fit = TRUE)
  whichCols <- grep("Chl", colnames(ww.pred))
  ww.predse <- as.data.frame(ww.predse$se.fit)
  whichColsse <- grep("Chl", colnames(ww.predse))
  ww.pdat <- cbind(ww.pdat, Fitted = ww.pred[, whichCols], Fittedse = ww.predse[,whichColsse])
  ww.pdat <- with(ww.pdat, transform(ww.pdat, Fittedplus = Fitted + Fittedse))
  # simple subtraction addition works here because family=gaussian
  # if I had a log-transformed response, I'd have to change the se:
  ## "The standard errors are for the function on the log scale. You need to exp() the fitted 
  ##    values and the upper and lower confidence before you plot" -GS
  ww.pdat <- with(ww.pdat, transform(ww.pdat, Fittedminus = Fitted - Fittedse))
  
  ## transform back into original pH value for clarity:
  # get mean/intercept pH used by the model
  shiftph <- attr(predict(phwwmod, newdata = ww.pdat, type = "iterms"), "constant")
  ww.pdatnorm <- ww.pdat
  ww.pdatnorm <- with(ww.pdatnorm, transform(ww.pdatnorm, Fitted = Fitted + shiftph, 
                                             Fittedplus = Fittedplus + shiftph, 
                                             Fittedminus = Fittedminus + shiftph))
  
  labdat2 <- data.frame(x = 270, y = 9.1, label = "mean pH: 9.1")
  # without this step, the resolution of the text is really off for some reason
  
  ## plot with all labels
  pdf("../data/private/wwphmod.pdf", width = 10, height = 6.7)
  ggplot(ww.pdatnorm, aes(x = Chl_a_ug_L, y = Fitted)) +
    geom_line() + 
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25) +  
    geom_abline(slope = 0, intercept = shiftph, linetype="dotted") +
    geom_text(data = labdat2, aes(label = label, x = x, y = y, size = 5), 
              show.legend = FALSE) +
    xlab('Chlorophyll a (ug/L)') + ylab('pH')
  dev.off()
  
  ## plot with disabled labels and white background
  tiff("../data/private/wwphmod-notext.tiff", width = 2.3, height = 1.5, 
       units = "in", res = 300)
  ggplot(ww.pdatnorm, aes(x = Chl_a_ug_L, y = Fitted)) +
    geom_line() + 
    geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                alpha = 0.25) +  
    theme_bw() +
    theme(axis.text = element_blank(), title = element_blank(), 
          panel.grid.major = element_blank(), panel.grid.minor = element_blank())
  dev.off()
  
  if (runextras) {
    ## ===========================================================================
    ## Below are my initial apply applications before learning that I can introduce
    ##    by = ... into a gam call which makes it all better cause I am still using 
    ##    all the data that way and other obvious reasons
    ## ===========================================================================
    co2gams <- function(df) {
      if (nrow(df) < 3) { return(NA) } else {
        co2gam <- gam(co2Flux ~ s(pH_surface, k = 20), data = df,
                      select = TRUE, method = "ML", family = scat(),
                      na.action = na.exclude)
      }}
    
    resgams <- function(df, res) {
      if (nrow(df) < 3) { return(NA) } else {
        res <- as.vector(res)
        resmod <- gam(res ~ s(Chl_a_ug_L) + s(GPP_h) + s(TDN_ug_L) + 
                        s(DOC_mg_L) + s(Oxygen_ppm) + 
                        te(PDO, SOI) + s(Year, bs = "re"), data = df,
                      select = TRUE, method = "REML", family = scat(),
                      na.action = na.exclude) }
    }
    
    ## apply gam 
    regmods <- lapply(regsplit, co2gams)
    
    ## apply residual extraction to regmods
    allres <- lapply(regmods, resid, type = "pearson", na.action = na.exclude)
    
    ## apply resgams to the residuals
    ## FIXME: this not doing the right thing!
    resmods <- mapply(resgams, regsplit, allres)
    
    ## can lapply these to the above to check diagnostics however
    ##    over gam.check the plotting default means I don't know how to store objects
    lapply(regmods, gam.check)
    regsumms <- lapply(regmods, summary)
    ## deviance explained ranges from 68 - 80%
    lapply(regsumms, '[[', "dev.expl")
    
    lapply(regmods, plot, pers = TRUE, pages = 1)
    
    
    ## 2. pH ~ .
    ## ===========================================================================
    phgams <- function(df) {
      if (nrow(df) < 3) { return(NA) } else {
        gam(pH_surface ~ s(Lake, Year, bs = "re") + s(Chl_a_ug_L) + s(GPP_h) +
              s(TDN_ug_L) + s(DOC_mg_L) + s(Oxygen_ppm) + te(PDO, SOI), data = df,
            select = TRUE, method = "REML", family = gaussian, na.action = na.exclude)
      }}
    
    phmods <- lapply(regsplit, phgams)
    
    summary(phmod)
    plot(phmod, pages = 1, pers = TRUE)
    gam.check(phmod)
  }
}