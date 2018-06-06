### Get site specific data/predictions

## load packages
library('ggplot2')
library('viridis')
library('cowplot')
library('mgcv')
library('extrafont')
library('gridExtra')
library("reshape2")
library('grid')

## Set defaults
theme_set(theme_bw())

## load in data
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")
if (!file.exists("../data/weathers.rds")) {
  source("../scripts/climate-weather-modeling.R")
}
weathers <- readRDS('../data/weathers.rds')

## Load in gam models
co2modnull <- readRDS("../data/private/co2modnull.rds")
co2mod <- readRDS("../data/private/co2mod.rds")

egmod.red2 <- readRDS("../data/private/egmodred2.rds")
egmodlaggedsimp <- readRDS("../data/private/egmodlaggedsimp.rds")

source("../functions/geom_rug3.R")

## change to match what was used to create the models
regvars <- merge(regvars, weathers)
regvarf <- regvars
regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))

## subset flux df
lakesub <- fluxes[,c("Lake","Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                     "TICumol", "Pressure", "pco2atm")]

## what is mean pH? perhaps add a dotted line at this point in the plots
meanpH <- with(regvarf2, 
               mean(pH_surface, na.rm=TRUE)) 
meanco2 <- with(regvarf2, 
               mean(co2Flux, na.rm=TRUE))
meanO2 <- with(regvarf2, 
                mean(Oxygen_ppm, na.rm=TRUE))
meanchl <- with(regvarf2, 
                mean(Chl_a_ug_L, na.rm=TRUE))
meanTDN <- with(regvarf2, 
                mean(TDN_ug_L, na.rm=TRUE))
meanGPP <- with(regvarf2, 
                mean(GPP_h, na.rm=TRUE))
meanspei <- with(regvarf2, 
                 mean(SPEI02, na.rm=TRUE))
## limit prediction data frame to observed intervals
regsplit <- with(regvarf2, split(regvarf2, list(Lake)))
regsplit <- regsplit[sapply(regsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.
minmax <- function(df, colnames) {
  allmin <- as.data.frame(do.call(cbind, lapply(df[,colnames], min, na.rm = TRUE)))
  names(allmin) <- sapply(names(allmin), function(x) paste("min",x, sep = ""))
  allmax <- as.data.frame(do.call(cbind, lapply(df[,colnames], max, na.rm = TRUE)))
  names(allmax) <- sapply(names(allmax), function(x) paste("max",x, sep = ""))
  summ <- as.data.frame(cbind(allmin, allmax))
  summ <- cbind(summ, data.frame(Lake = df$Lake[1]))
  summ
}
minmaxes <- do.call(rbind, lapply(regsplit, minmax, colnames = c("GPP_h", "TDN_ug_L", "DOC_mg_L", 
                                      "Oxygen_ppm", "pH_surface", "Chl_a_ug_L")))
rownames(minmaxes) <- NULL

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top')
  

## generate predicted for co2 and ph for all signif variables
## co2 ~ ph: full model
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

co2.pdat <- with(droplevels(regvarf2),
                  data.frame(`pH_surface` = rep(seq(min(`pH_surface`, na.rm = TRUE),
                                                    max(`pH_surface`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
co2.pdat <- merge(co2.pdat, lakeXbar)
co2.pred <- predict(co2mod, newdata = co2.pdat, type = "terms", se.fit = TRUE)

whichCols <- grep("pH", colnames(co2.pred$fit))
whichColsSE <- grep("pH", colnames(co2.pred$se.fit))
co2.pdat <- cbind(co2.pdat, Fitted = rowSums(co2.pred$fit[, whichCols]), 
                  se.Fitted = rowSums(co2.pred$se.fit[, whichColsSE]))
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedplus = Fitted + se.Fitted))
co2.pdat <- with(co2.pdat, transform(co2.pdat, Fittedminus = Fitted - se.Fitted))

## make into original limits
shiftco2 <- attr(predict(co2mod, newdata = co2.pdat, type = "iterms"), "constant")
co2.pdatnorm <- co2.pdat
co2.pdatnorm <- with(co2.pdatnorm, transform(co2.pdatnorm, Fitted = Fitted + shiftco2))
co2.pdatnorm <- merge(co2.pdatnorm, minmaxes)
overs <- with(co2.pdatnorm, which(pH_surface < minpH_surface | pH_surface > maxpH_surface))
co2.pdatnorm <- co2.pdatnorm[-overs,]

labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

## add quantiles
phquants <- quantile(regvarf2$pH_surface, c(.05,.95), na.rm = TRUE)

co2plot <- ggplot(co2.pdatnorm, aes(x = pH_surface, y = Fitted, 
                    colour = ifelse(Lake == "WW", "Wascana", 
                                    ifelse(Lake == "D", "Diefenbaker",
                                           ifelse(Lake == 'K', "Katepwa", 
                                                  ifelse(Lake == 'P', "Pasqua", 
                                                         ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                ifelse(Lake=='L','Last Mountain',
                                                                       'Crooked')))))),
                    lty=ifelse(Lake == "WW", "Wascana", 
                               ifelse(Lake == "D", "Diefenbaker",
                                      ifelse(Lake == 'K', "Katepwa", 
                                             ifelse(Lake == 'P', "Pasqua", 
                                                    ifelse(Lake == 'B', 'Buffalo Pound', 
                                                           ifelse(Lake=='L','Last Mountain',
                                                                  'Crooked')))))))) +
  papertheme + 
  annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  #geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
   #           alpha = 0.25, fill='white') +  
  scale_linetype_manual(name='Lake', values = c("solid", "dotdash","longdash", "solid", "longdash", 
                                                "solid", "solid")) +
  scale_colour_manual(name="Lake", values = c("#b2abd2", "#5e3c99","#b2abd2", "#5e3c99", "#e66101",
                                              "#5e3c99", "#5e3c99"))+
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  geom_vline(xintercept = meanpH, linetype = 'dotted') +
  geom_rug3(aes(x=pH_surface, y=co2Flux), data = regvarf2, stat = "identity", position = "identity", 
           sides = "bl", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  theme(legend.position='top') +
  guides(colour=guide_legend(ncol=3, bycol =TRUE,title.position = 'left')) +
  xlab('pH') + ylab(expression(paste('CO'[2]*' (mmol m'^{-2}*d^{-1}*')')))

## oranges: #E6550D, #FD8D3C, blues: #A6BDDB, #67A9CF, #02818A
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)
## for the null model:
## for co2modnull
N <- 200
null.pdat <- with(droplevels(regvarf),
                  data.frame(`pH_surface` = rep(seq(7, 11, length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N))))
null.pred <- predict(co2modnull, newdata = null.pdat, type = "iterms") 
whichCols <- grep("pH", colnames(null.pred))

null.predse <- predict(co2modnull, newdata = null.pdat, type = "iterms", se.fit = TRUE)
null.predse <- as.data.frame(null.predse$se.fit)
whichColsse <- grep("pH", colnames(null.predse))

null.pdat <- cbind(null.pdat, Fitted = null.pred[, whichCols], Fittedse = null.predse[,whichColsse])
null.pdat <- with(null.pdat, transform(null.pdat, Fittedplus = Fitted + Fittedse))
null.pdat <- with(null.pdat, transform(null.pdat, Fittedminus = Fitted - Fittedse))

shiftnull <- attr(predict(co2modnull, newdata = null.pdat, type = "iterms"), "constant")
null.pdatnorm <- null.pdat
null.pdatnorm <- with(null.pdatnorm, transform(null.pdatnorm, Fitted = Fitted + shiftnull, 
                                               Fittedplus = Fittedplus + shiftnull, 
                                               Fittedminus = Fittedminus + shiftnull))
null.pdatnorm <- merge(null.pdatnorm, minmaxes)
overs <- with(null.pdatnorm, which(pH_surface < minpH_surface | pH_surface > maxpH_surface))
null.pdatnorm <- null.pdatnorm[-overs,]

labdatco2 <- data.frame(x = 7.5, y = meanco2 - 18, label = "mean flux")

nullplot <- ggplot(null.pdatnorm, aes(x = pH_surface, y = Fitted)) +
  geom_line() +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatco2, aes(label = label, x = x, y = y, size = 5), 
          #  show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanco2, linetype="dotted") +
  #geom_point(data=regvars, aes(x=pH_surface, y=co2Flux)) +
  geom_vline(xintercept = meanpH, linetype="dotted") +
  ylab(expression(paste(CO[2]~"flux (mmol"~"C "*m^{-2}*"d"^{-1}*')'))) + xlab('pH')

# plot something for demonstrating GAMs
lmplot <- ggplot(regvars, aes(x = pH_surface, y = co2Flux)) +
  geom_point() +
  theme_bw(base_size = 15) +
  stat_smooth(method = "lm", linetype = "dotted", col='black') +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())

nullplotnodots <- ggplot(null.pdatnorm, aes(x = pH_surface, y = Fitted)) +
  geom_line(linetype='dotted') +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_point(data=regvars, aes(x=pH_surface, y=co2Flux)) +
  theme(axis.text=element_blank(), axis.ticks=element_blank(), axis.title=element_blank())
gamdemo <- plot_grid(lmplot, nullplotnodots)
ggsave("../docs/private/gamdemo.png", gamdemo, width=10, height=4, units = 'in')

## for egmod.red2
## oxygen
N <- 200
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

oxy.pdat <- with(droplevels(regvarf2),
                  data.frame(`Oxygen_ppm` = rep(seq(min(`Oxygen_ppm`, na.rm = TRUE),
                                                    max(`Oxygen_ppm`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
oxy.pdat <- merge(oxy.pdat, lakeXbar)
oxy.pred <- predict(egmod.red2, newdata = oxy.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("Oxy", colnames(oxy.pred$fit))
whichColsSE <- grep("Oxy", colnames(oxy.pred$se.fit))
oxy.pdat <- cbind(oxy.pdat, Fitted = oxy.pred$fit[, whichCols], 
                  se.Fitted = oxy.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedplus = Fitted + se.Fitted))
oxy.pdat <- with(oxy.pdat, transform(oxy.pdat, Fittedminus = Fitted - se.Fitted))

shiftoxy <- attr(predict(egmod.red2, newdata = oxy.pdat, type = "iterms"), "constant")
oxy.pdatnorm <- oxy.pdat
oxy.pdatnorm <- with(oxy.pdatnorm, transform(oxy.pdatnorm, Fitted = Fitted + shiftoxy, 
                                             Fittedplus = Fittedplus + shiftoxy, 
                                             Fittedminus = Fittedminus + shiftoxy))
oxy.pdatnorm <- merge(oxy.pdatnorm, minmaxes)
overs <- with(oxy.pdatnorm, which(Oxygen_ppm < minOxygen_ppm | Oxygen_ppm > maxOxygen_ppm))
oxy.pdatnorm <- oxy.pdatnorm[-overs,]

labdatoxy <- data.frame(x = c(2.5, 11), y = c(meanpH + 0.03, 8.4), label = c("mean pH", 'mean oxygen'))

oxyquants <- quantile(regvarf2$Oxygen_ppm, c(.05,.95), na.rm = TRUE)

oxyplot <- ggplot(oxy.pdatnorm, aes(x = Oxygen_ppm, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=oxyquants[1], xmax=oxyquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanO2, linetype="dotted") +
  xlab(expression(paste(O[2]~"("*"mg"~L^{-1}*")"))) + ylab('pH')

## for GPP
N <- 200
varWant <- c("Oxygen_ppm", "TDN_ug_L", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

GPP.pdat <- with(droplevels(regvarf2),
                 data.frame(`GPP_h` = rep(seq(min(`GPP_h`, na.rm = TRUE),
                                                   max(`GPP_h`, na.rm = TRUE),
                                                   length = N),
                                               nlevels(Lake)),
                            Lake = rep(levels(Lake), each = N),
                            Year = rep(2004, prod(nlevels(Lake), N)),
                            dummy = rep(0, prod(nlevels(Lake), N))))
GPP.pdat <- merge(GPP.pdat, lakeXbar)
GPP.pred <- predict(egmod.red2, newdata = GPP.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("GPP", colnames(GPP.pred$fit))
whichColsSE <- grep("GPP", colnames(GPP.pred$se.fit))
GPP.pdat <- cbind(GPP.pdat, Fitted = GPP.pred$fit[, whichCols], 
                  se.Fitted = GPP.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
GPP.pdat <- with(GPP.pdat, transform(GPP.pdat, Fittedplus = Fitted + se.Fitted))
GPP.pdat <- with(GPP.pdat, transform(GPP.pdat, Fittedminus = Fitted - se.Fitted))

shiftGPP <- attr(predict(egmod.red2, newdata = GPP.pdat, type = "iterms"), "constant")
GPP.pdatnorm <- GPP.pdat
GPP.pdatnorm <- with(GPP.pdatnorm, transform(GPP.pdatnorm, Fitted = Fitted + shiftGPP, 
                                             Fittedplus = Fittedplus + shiftGPP, 
                                             Fittedminus = Fittedminus + shiftGPP))
GPP.pdatnorm <- merge(GPP.pdatnorm, minmaxes)
overs <- with(GPP.pdatnorm, which(GPP_h < minGPP_h | GPP_h > maxGPP_h))
GPP.pdatnorm <- GPP.pdatnorm[-overs,]

labdatGPP <- data.frame(x = 0, y = meanpH + 0.01, label = "mean pH")

GPPplot <- ggplot(GPP.pdatnorm, aes(x = GPP_h, y = Fitted)) +
  geom_line() +
  theme_bw(base_size = 15) +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  #geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanGPP, linetype='dotted')
  xlab(expression(paste('GPP ('~O[2]~h^{-1}*")"))) + ylab('pH')

## for TDN
N <- 200
varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

TDN.pdat <- with(droplevels(regvarf2),
                 data.frame(`TDN_ug_L` = rep(seq(min(`TDN_ug_L`, na.rm = TRUE),
                                              max(`TDN_ug_L`, na.rm = TRUE),
                                              length = N),
                                          nlevels(Lake)),
                            Lake = rep(levels(Lake), each = N),
                            Year = rep(2004, prod(nlevels(Lake), N)),
                            dummy = rep(0, prod(nlevels(Lake), N))))
TDN.pdat <- merge(TDN.pdat, lakeXbar)
TDN.pred <- predict(egmod.red2, newdata = TDN.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("TDN", colnames(TDN.pred$fit))
whichColsSE <- grep("TDN", colnames(TDN.pred$se.fit))
TDN.pdat <- cbind(TDN.pdat, Fitted = TDN.pred$fit[, whichCols], 
                  se.Fitted = TDN.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedplus = Fitted + se.Fitted))
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedminus = Fitted - se.Fitted))

shiftTDN <- attr(predict(egmod.red2, newdata = TDN.pdat, type = "iterms"), "constant")
TDN.pdatnorm <- TDN.pdat
TDN.pdatnorm <- with(TDN.pdatnorm, transform(TDN.pdatnorm, Fitted = Fitted + shiftTDN, 
                                             Fittedplus = Fittedplus + shiftTDN, 
                                             Fittedminus = Fittedminus + shiftTDN))
TDN.pdatnorm <- merge(TDN.pdatnorm, minmaxes)
overs <- with(TDN.pdatnorm, which(TDN_ug_L < minTDN_ug_L | TDN_ug_L > maxTDN_ug_L))
TDN.pdatnorm <- TDN.pdatnorm[-overs,]

labdatN <- data.frame(x = 3000, y = meanpH + 0.03, label = "mean pH")

TDNquants <- quantile(regvarf2$TDN_ug_L, c(.05,.95), na.rm = TRUE)

TDNplot <- ggplot(TDN.pdatnorm, aes(x = TDN_ug_L, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=TDNquants[1], xmax=TDNquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
                                     alpha = 0.25) +  
  #geom_text(data = labdatN, aes(label = label, x = x, y = y, size = 5), 
   #         show.legend = FALSE) +
  scale_x_log10(breaks = c(100,200,400,800,1600,3200,6400)) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanTDN, linetype = 'dotted') +
  xlab(expression(paste("TDN ("~mu*"g"~L^{-1}*")"))) + ylab('pH')

## for chl a
## Get site-specific data/predictions
N <- 1000
varWant <- c("GPP_h", "TDN_ug_L", "DOC_mg_L", "Oxygen_ppm", "PDO", "SOI")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

chla.pdat <- with(droplevels(regvarf2),
                  data.frame(`Chl_a_ug_L` = rep(seq(min(`Chl_a_ug_L`, na.rm = TRUE),
                                                    max(`Chl_a_ug_L`, na.rm = TRUE),
                                                    length = N),
                                                nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
chla.pdat <- merge(chla.pdat, lakeXbar)
chla.pred <- predict(egmod.red2, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

shiftchl <- attr(predict(egmod.red2, newdata = chla.pdat, type = "iterms"), "constant")
chl.pdatnorm <- chla.pdat
chl.pdatnorm <- with(chl.pdatnorm, transform(chl.pdatnorm, Fitted = Fitted + shiftchl))
chl.pdatnorm <- merge(chl.pdatnorm, minmaxes)
overs <- with(chl.pdatnorm, which(Chl_a_ug_L < minChl_a_ug_L | Chl_a_ug_L > maxChl_a_ug_L))
chl.pdatnorm <- chl.pdatnorm[-overs,]

labdatchl <- data.frame(x = 100, y = meanpH + 0.04, label = "mean pH")

chlquants <- quantile(regvarf2$Chl_a_ug_L, c(.05,.95), na.rm = TRUE)
chlaplot <- ggplot(chl.pdatnorm, aes(x = Chl_a_ug_L, y = Fitted, group = Lake, colour = 
                                       ifelse(Lake == "WW", "Wascana", 
                                                    ifelse(Lake == "D", "Diefenbaker",
                                                           ifelse(Lake == 'K', "Katepwa", 
                                                                  ifelse(Lake == 'P', "Pasqua", 
                                                                         ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                                ifelse(Lake=='L','Last Mountain',
                                                                                       'Crooked')))))),
                                           lty=ifelse(Lake == "WW", "Wascana", 
                                                      ifelse(Lake == "D", "Diefenbaker",
                                                             ifelse(Lake == 'K', "Katepwa", 
                                                                    ifelse(Lake == 'P', "Pasqua", 
                                                                           ifelse(Lake == 'B', 'Buffalo Pound', 
                                                                                  ifelse(Lake=='L','Last Mountain',
                                                                                         'Crooked')))))))) +
  papertheme + 
  annotate("rect", xmin=chlquants[1], xmax=chlquants[2], ymin=-Inf, ymax=Inf, alpha = 0.1, fill='gray60') +
  geom_line() + 
  #geom_text(data = labdatchl, aes(label = label, x = x, y = y, size = 5),
  #          show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanchl, linetype='dotted') +
  scale_linetype_manual(name='Lake', values = c("longdash", "solid","dotdash", "solid", "solid", 
                                                "solid", "solid")) +
  scale_colour_manual(name="Lake", values = c("#e66101", "#5e3c99","#5e3c99", "#5e3c99", "#5e3c99",
                                              "#5e3c99", "#b2abd2"))+
  theme(legend.position = 'top', legend.direction = "vertical", 
        axis.text.x = element_text(angle = 45)) +
  scale_x_log10(breaks = c(5,10,25,50,100,200,300)) +
  xlab(expression(paste("Chl"~italic(a)~"("~mu*"g"~"L"^{-1}*")"))) + 
  guides(colour=guide_legend(ncol=4,nrow=2,bycol =TRUE,title.position = 'left'),
         lty=guide_legend(ncol=4,nrow=2, bycol =TRUE,title.position = 'left')) +
  ylab('pH')
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## see http://stackoverflow.com/questions/11838278/
##    plot-with-conditional-colors-based-on-values-in-r
## make this bit into soi-pdo
N <- 100
simSOI <- c(-1, 0.3, 1.1)
SOIgroup <- factor(rep(c('-1', '0.3', '1.1'), each=100, times=7))  # 7*300
reptimes <- length(simSOI) #3
preddf <- data.frame(PDO = rep(seq(min(regvarf2$PDO, na.rm=TRUE),
                                      max(regvarf2$PDO, na.rm=TRUE),
                                      length = N), times=reptimes),
                     SOI = rep(simSOI, each = N)) #300
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

SOI.pdat <- with(droplevels(regvarf2),
                 data.frame(PDO = rep(preddf$PDO,
                                             nlevels(Lake)), # 300*7
                            SOI = rep(preddf$SOI,
                                        nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

SOI.pdat <- merge(SOI.pdat, lakeXbar)
SOI.pred <- predict(egmod.red2, newdata = SOI.pdat, type = "link")
SOI.pdat <- cbind(SOI.pdat, SOI.pred)
SOI.pdat$SOIgroup <- SOIgroup

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'pH'

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDO, SOI.pdat$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
SOI.pdat$pH[toofar] <- NA

labdatSOI <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat[SOI.pdat$Lake=='B',], aes(x = PDO, y = pH, group= SOI, col=SOIgroup, 
                                                     lty=SOIgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('E')~'      SOI')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('E')~'      SOI')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "SOI", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('PDO') + ylab('pH')
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## now slices for PDO
N <- 100
simPDO <- c(-1.5, 0, 1.6)
PDOgroup <- factor(rep(c('-1.5', '0', '1.6'), each=100, times=7))  # 7*300
reptimes <- length(simPDO) #3
preddf <- data.frame(SOI = rep(seq(min(regvarf2$SOI, na.rm=TRUE),
                                   max(regvarf2$SOI, na.rm=TRUE),
                                   length = N), times=reptimes),
                     PDO = rep(simPDO, each = N)) #300
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

PDO.pdat <- with(droplevels(regvarf2),
                 data.frame(SOI = rep(preddf$SOI,
                                      nlevels(Lake)), # 300*7
                            PDO = rep(preddf$PDO,
                                      nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

PDO.pdat <- merge(PDO.pdat, lakeXbar)
PDO.pred <- predict(egmod.red2, newdata = PDO.pdat, type = "link")
PDO.pdat <- cbind(PDO.pdat, PDO.pred)
PDO.pdat$PDOgroup <- PDOgroup

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'pH'

# need to take away places where PDO*PDO combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDO, PDO.pdat$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
PDO.pdat$pH[toofar] <- NA

labdatPDO <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-1.5', '0', '1.6'))

PDOplot <- ggplot(PDO.pdat[PDO.pdat$Lake=='B',], aes(x = SOI, y = pH, group= PDO, col=PDOgroup,
                                                     lty=PDOgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('F')~'      PDO')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('F')~ '      PDO')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "PDO", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH) +
  xlab('SOI') + ylab('pH')

## get expand grid object to plot a ggplot contour/heatmap
preddf <- with(regvarf2, expand.grid(PDO = seq(min(PDO), max(PDO), length = 100), 
                                   SOI = seq(min(SOI), max(SOI), length = 100)))
superN <- nrow(preddf)

varWant <- c("Oxygen_ppm", "GPP_h", "DOC_mg_L", "Chl_a_ug_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

comb.pdat <- with(droplevels(regvarf2),
                 data.frame(PDO = rep(preddf$PDO,
                                        nlevels(Lake)), 
                            SOI = rep(preddf$SOI,
                                        nlevels(Lake)),
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

comb.pdat <- merge(comb.pdat, lakeXbar, sort=FALSE)
comb.pred <- predict(egmod.red2, newdata = comb.pdat, type = "terms")

whichCols <- grep("SOI", colnames(comb.pred))
comb.pdat <- cbind(comb.pdat, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$PDO, comb.pdatnorm$SOI, regvarf2$PDO, regvarf2$SOI, dist=0.1)
comb.pdatnorm$pH <- comb.pdatnorm$Fitted
comb.pdatnorm$pH[toofar] <- NA

names(comb.pdat)[which(names(comb.pdat)=='SOI.pred')] <- 'pH'
comboplot <- ggplot(comb.pdatnorm, aes(x = SOI, y = PDO, z=pH)) + #, z=Fitted
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=pH)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent') +
  geom_point(data=regvarf2, aes(x=SOI, y=PDO, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  geom_vline(linetype='dotted', xintercept = c(-1,0.3, 1.1)) +
  geom_abline(slope = 0, intercept = c(-1.5, 0, 1.6), linetype='dashed')

## save objects:
plots <- list(TDNplot, oxyplot, GPPplot, chlaplot, SOIplot, PDOplot, comboplot, nullplot, co2plot)
plotnames <- c('TDNplot', 'oxyplot', 'GPPplot','chlaplot', 'SOIplot', 'PDOplot', 'comboplot',
               'nullplot','co2plot')

invisible( # this means I don't get the list [[1:3]] returned on screen
  lapply(
    seq_along(plots), 
    function(x) ggsave(filename=paste0("../docs/private/gam-plots-", plotnames[x], ".pdf"), 
                       plot=plots[[x]], scale=0.6) # width=7, height=5, units = 'in'
  ) )

## arrange plots
## ## FIXME: get relative sizes of plots!!
allgam <- plot_grid(co2plot, chlaplot, oxyplot, TDNplot, ncol = 2, rel_heights = c(2,1.5), 
                    labels='AUTO')
climgam <- grid.arrange(comboplot, SOIplot, PDOplot, ncol = 2, layout_matrix = cbind(c(1,1,2), c(1,1,3)))
#papergam <- grid.arrange(climgam, chlaplot, oxyplot, TDNplot, 
    #                     layout_matrix=rbind(c(2,1),c(2,1), c(2,1), c(3,1), c(3,1), c(4,1), c(4,1)))
    ## FIXME: how to get labels in grid.arrange??
splineplots <- plot_grid(chlaplot, oxyplot, TDNplot, ncol = 1, nrow=3, rel_heights = c(2, 1.5, 1.5), 
                     labels='AUTO')
vararrange <- plot_grid(splineplots, climgam, ncol = 2, labels=c('', "D"))
ggsave("../docs/private/ph-allgams.pdf", allgam, scale=0.77) #width=28, height=18, units = 'cm'
ggsave("../docs/private/climgam.pdf", climgam, scale=0.77, width = 7.5)
ggsave("../docs/private/paper-gams.pdf", vararrange, scale=1.2, width = 12)

## var vs time option for case where peter might want some temporal info
tempplot <- ggplot(regvarf2, aes(x = Year, y = pH_surface, group= Lake)) +
  papertheme +
  facet_wrap( "Lake", scales = "fixed", ncol=2) +
  geom_point() +
  theme(axis.text.x=element_text(angle=45, vjust=0.5), axis.title.x=element_blank()) +
  ylab('pH')

## or maybe this could be an approach
testing <- gam(data=regvars[regvars$Month > 4 & regvars$Month < 9,], pH_surface ~ s(Year) + 
                 s(Year, by=Lake, m=1, k=7) + 
                 s(DOY) + s(DOY, by=Lake, m=1) +
                 s(Lake, bs='re'), family='gaussian', na.action=na.exclude, select=TRUE)
N <- 500
varWant <- "DOY"
lakeXbar <- with(regvars, do.call(rbind, lapply(split(regvars[, varWant], droplevels(Lake)), 
                                                mean, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

Year.pdat <- with(droplevels(regvars),
                  data.frame(Year = rep(seq(min(`Year`, na.rm = TRUE),
                                            max(`Year`, na.rm = TRUE),
                                            length = N),
                                        nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N)
                  ))
Year.pdat <- merge(Year.pdat, lakeXbar)
names(Year.pdat)[which(names(Year.pdat)=='X_data')] <- 'DOY'

toofar <- c(which(Year.pdat$Lake == 'WW' & Year.pdat$Year < 1996),
            which(Year.pdat$Lake == 'P' & Year.pdat$Year < 2004))
Year.pdat <- Year.pdat[-toofar,]

Year.pred <- predict(testing, newdata = Year.pdat, type = "terms")
whichCols <- grep("Year", colnames(Year.pred))
Year.pdat <- cbind(Year.pdat, Fitted = rowSums(Year.pred[, whichCols]))

shiftYear <- attr(Year.pred, "constant")
Year.pdatnorm <- Year.pdat
Year.pdatnorm <- with(Year.pdatnorm, transform(Year.pdatnorm, Fitted = Fitted + shiftYear))


Yearplot <- ggplot(Year.pdatnorm, 
                   aes(x = Year, y = Fitted, group = Lake, 
                       colour = 
                         ifelse(Lake == "WW", "Wascana", 
                                ifelse(Lake == "D", "Diefenbaker",
                                       ifelse(Lake == 'K', "Katepwa", 
                                              ifelse(Lake == 'P', "Pasqua", 
                                                     ifelse(Lake == 'B', 'Buffalo Pound', 
                                                            ifelse(Lake=='L','Last Mountain',
                                                                   'Crooked')))))),
                       lty=ifelse(Lake == "WW", "Wascana", 
                                  ifelse(Lake == "D", "Diefenbaker",
                                         ifelse(Lake == 'K', "Katepwa", 
                                                ifelse(Lake == 'P', "Pasqua", 
                                                       ifelse(Lake == 'B', 'Buffalo Pound', 
                                                              ifelse(Lake=='L','Last Mountain',
                                                                     'Crooked')))))))) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() + 
  xlab("Year") + 
  guides(colour=guide_legend(ncol=2,bycol =TRUE,title.position = 'left')) + 
  scale_linetype_manual(name='Lake', values = c("solid", "longdash","dotdash", "solid", "solid", 
                                                "solid", "longdash")) +
  scale_colour_manual(name="Lake", values = c("#5e3c99", "#b2abd2","#5e3c99", "#5e3c99", "#b2abd2",
                                              "#e66101", "#e66101"))+
  #c("black", "#02818A","black", "black", "#02818A",
  #"#A6BDDB", "#A6BDDB"))+
  ylab('pH') + #02818A", "#A6BDDB" 
  #theme(legend.position = "none")+
  #guides(colour=guide_legend(nrow=1,byrow=TRUE))+
  scale_x_continuous(breaks=c(1995, 2000,2005, 2010, 2014), labels = c(1995, 2000,2005, 2010, 2014))

## general lake diffs in variables -- remove those deriving from same weather station
melted <- melt(lakesub[,c('Lake', 'Temperature', 'Conductivity', 'pH', 'SalCalc', 'TICumol')], id = "Lake")

# create a list with strip labels
varnames <- list(
  'Temperature'=expression(paste("Water temperature ("~degree*"C)")) ,
  'Conductivity'= expression(paste("Conductivity ("*mu*"S"~"cm"^{-1}*")")),
  'pH'="pH",
  #'meanWindMS'= expression(paste("Mean Wind (m"~"s"^{-1}*")")),
  'SalCalc' = "Salinity (ppt)",
  'TICumol' = expression(paste("DIC ("~mu*"mol"~"L"^{-1}*")"))#,
  #'Pressure' = 'Pressure (kPa)',
  #'pco2atm' = expression(paste("Air"~italic(p)*"CO"[2]~"(ppm)"))
)

# Create a 'labeller' function, and push it into facet_grid call:
var_labeller <- function(variable,value){
  return(varnames[value])
}

# run the ggplot
melted$Lake <- as.character(melted$Lake)
melted$Lake[melted$Lake == 'WW'] <- "W"
melted$Lake <- factor(melted$Lake)
meltplot <- ggplot(melted, aes(x=Lake,y=value, group=Lake)) +
  geom_boxplot(outlier.colour="black", outlier.shape=5,
               outlier.size=1) +
  theme_bw(base_size = 14, base_family = 'Arial') +
  facet_wrap( "variable", scales = "free", labeller = var_labeller, ncol=2, strip.position = "left") +
  theme(axis.title = element_blank())

boxgrid <- plot_grid(meltplot, Yearplot, ncol = 1, nrow=2, 
                     rel_heights = c(3,1.5), labels = c('a.','b.'))
ggsave("../docs/private/summaryfig.pdf", boxgrid, width=20, height=40, units = "cm")

## ==================================================================================
## Testing effect by time plots
## ==================================================================================
testing1 <- predict(egmod.red2, type = 'terms')
testing <- as.data.frame(testing1)
tosum <- grep("Chl", colnames(testing))
chleffect <- rowSums(testing[,tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$Chl <- chleffect
names(testing) <- c("GPP", "TDN", "DOC", "Oxy", "PDO-SOI", 'LakeYear', 'Chl')
testing$Date <- regvarf2$Date #is this ok? assuming order is preserved
testing$Year <- regvarf2$Year
testing$Lake <- regvarf2$Lake
testing$DOY <- regvarf2$DOY
ggplot(testing, aes(x = DOY, y = Chl, group=Lake, colour = Lake)) +
  geom_point() +
  theme_bw(base_size = 15) +
  facet_wrap('Year') +  
  #geom_text(data = labdatGPP, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  xlab(expression(paste('GPP ('~O[2]~h^{-1}*")"))) + ylab('pH')
## could do relative influence of each term by time, or absolute... by month, year,
##    Lake..... a lot of options

