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

source("../functions/geom_rug3.R")

## load in data
regvars <- readRDS("../data/private/regvars.rds")
fluxes <- readRDS("../data/private/params-flux.rds")
if (!file.exists("../data/weathers.rds")) {
  source("../scripts/climate-weather-modeling.R")
}
weathers <- readRDS('../data/weathers.rds')

## Load in gam models
egmodlagged <- readRDS("../data/private/egmodlaggedsimp.rds")
egmodtdn <- readRDS("../data/private/egmodlaggedtdn.rds")

co2mod <- readRDS("../data/private/co2mod.rds")

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
meanR <- with(regvarf2, 
              mean(R_h, na.rm=TRUE))
meanspei <- with(regvarf2, 
                 mean(SPEI02, na.rm=TRUE))
meanox <- with(regvarf2, 
               mean(Oxygen_ppm, na.rm=TRUE)) 
meanDOC <- with(regvarf2, 
                mean(DOC_mg_L, na.rm=TRUE))

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
minmaxes <- do.call(rbind, lapply(regsplit, minmax, colnames = c("R_h", "TDN_ug_L", "DOC_mg_L", 
                                                                 "Oxygen_ppm", "pH_surface", "Chl_a_ug_L")))
rownames(minmaxes) <- NULL

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=18, base_family = 'Arial') +
  theme(legend.position='top')

## summary plot for all vars
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
## add yearly means as rug
regmod <- regvars
regmod <- ddply(regmod, .(Year, Lake), dplyr::summarise, meanpH=mean(pH_surface, na.rm=TRUE))

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
  geom_rug3(inherit.aes=FALSE, data=regmod, aes(y=meanpH), stat = "identity", 
            position = "identity", 
            sides = "l", na.rm = FALSE, show.legend = NA, alpha=0.3) +
  #c("black", "#02818A","black", "black", "#02818A",
  #"#A6BDDB", "#A6BDDB"))+
  ylab('pH') + #02818A", "#A6BDDB" 
  #theme(legend.position = "none")+
  #guides(colour=guide_legend(nrow=1,byrow=TRUE))+
  scale_x_continuous(breaks=c(1995, 2000,2005, 2010, 2014), labels = c(1995, 2000,2005, 2010, 2014))

## general lake diffs in variables -- remove those deriving from same weather station
lakesub <- fluxes[,c("Lake","Temperature", "Conductivity", "pH", "meanWindMS", "SalCalc", 
                     "TICumol", "Pressure", "pco2atm")]
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

## CO2-pH model
## generate predicted for co2 and ph 
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
co2.pdatnorm$realname <- with(co2.pdatnorm,ifelse(Lake == "WW", "Wascana", 
                             ifelse(Lake == "D", "Diefenbaker",
                                    ifelse(Lake == 'K', "Katepwa", 
                                           ifelse(Lake == 'P', "Pasqua", 
                                                  ifelse(Lake == 'B', 'Buffalo Pound', 
                                                         ifelse(Lake=='L','Last Mountain',
                                                                'Crooked')))))))
co2plot <- ggplot(co2.pdatnorm, aes(x = pH_surface, y = Fitted, 
                                    colour = realname, lty=realname)) +
  papertheme + 
  annotate("rect", xmin=phquants[1], xmax=phquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
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
  #theme(legend.position='top') +
  guides(colour=guide_legend(nrow=3, ncol=3, byrow =TRUE,title.position = 'left',
                             override.aes=list(size=0.7)),
         linetype = guide_legend(nrow=3, ncol=3, byrow=TRUE, keywidth = unit(3,'line'))) +
   xlab('pH') + ylab(expression(paste('CO'[2]*' (mmol m'^{-2}*d^{-1}*')')))

## for SPEI
N <- 200
varWant <- c("TDN_ug_L", "Chl_a_ug_L", "PDOmean", "SOImean", "Oxygen_ppm", "R_h", "DOC_mg_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

spei.pdat <- with(droplevels(regvarf2),
                  data.frame(`SPEI02` = rep(seq(min(`SPEI02`, na.rm = TRUE),
                                                max(`SPEI02`, na.rm = TRUE),
                                                length = N),
                                            nlevels(Lake)),
                             Lake = rep(levels(Lake), each = N),
                             Year = rep(2004, prod(nlevels(Lake), N)),
                             dummy = rep(0, prod(nlevels(Lake), N))))
spei.pdat <- merge(spei.pdat, lakeXbar)
spei.pred <- predict(egmodlagged, newdata = spei.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("SPEI", colnames(spei.pred$fit))
whichColsSE <- grep("SPEI", colnames(spei.pred$se.fit))
spei.pdat <- cbind(spei.pdat, Fitted = spei.pred$fit[, whichCols], 
                   se.Fitted = spei.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
spei.pdat <- with(spei.pdat, transform(spei.pdat, Fittedplus = Fitted + se.Fitted))
spei.pdat <- with(spei.pdat, transform(spei.pdat, Fittedminus = Fitted - se.Fitted))

shiftspei <- attr(predict(egmodlagged, newdata = spei.pdat, type = "iterms"), "constant")
spei.pdatnorm <- spei.pdat
spei.pdatnorm <- with(spei.pdatnorm, transform(spei.pdatnorm, Fitted = Fitted + shiftspei, 
                                               Fittedplus = Fittedplus + shiftspei, 
                                               Fittedminus = Fittedminus + shiftspei))

speiquants <- quantile(regvarf2$SPEI02, c(.05,.95), na.rm = TRUE)

speiplot <- ggplot(spei.pdatnorm, aes(x = SPEI02, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=speiquants[1], xmax=speiquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_rug3(aes(x=SPEI02), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  #geom_text(data = labdatoxy, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanspei, linetype="dotted") +
  xlab("SPEI index (dry to wet)") + ylab('pH') +
  ylim(c(8.7, 9.2)) +
  ggtitle('d.') + theme(plot.title = element_text(hjust = -0.1, size=14, face="bold"))

# for TDN
N <- 200
varWant <- c("Chl_a_ug_L", "PDOmean", "SOImean", "SPEI02", "Oxygen_ppm", "R_h", "DOC_mg_L")
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
TDN.pred <- predict(egmodtdn, newdata = TDN.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("TDN", colnames(TDN.pred$fit))
whichColsSE <- grep("TDN", colnames(TDN.pred$se.fit))
TDN.pdat <- cbind(TDN.pdat, Fitted = TDN.pred$fit[, whichCols], 
                  se.Fitted = TDN.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

# make into original limits
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedplus = Fitted + se.Fitted))
TDN.pdat <- with(TDN.pdat, transform(TDN.pdat, Fittedminus = Fitted - se.Fitted))

shiftTDN <- attr(predict(egmodlagged, newdata = TDN.pdat, type = "iterms"), "constant")
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
  theme_bw(base_size=12, base_family = 'Arial') +
  theme(legend.position='top') +
  annotate("rect", xmin=TDNquants[1], xmax=TDNquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_rug3(aes(x=TDN_ug_L), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  #geom_text(data = labdatN, aes(label = label, x = x, y = y, size = 5), 
  #         show.legend = FALSE) +
  scale_x_log10(breaks = c(100,200,400,800,1600,3200,6400)) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanTDN, linetype = 'dotted') +
  xlab(expression(paste("TDN ("~mu*"g"~L^{-1}*")"))) + ylab('pH')

## for Oxygen_ppm
N <- 200
varWant <- c("Chl_a_ug_L", "PDOmean", "SOImean", "SPEI02","R_h", "DOC_mg_L", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

Oxygen.pdat <- with(droplevels(regvarf2),
                    data.frame(`Oxygen_ppm` = rep(seq(min(`Oxygen_ppm`, na.rm = TRUE),
                                                      max(`Oxygen_ppm`, na.rm = TRUE),
                                                      length = N),
                                                  nlevels(Lake)),
                               Lake = rep(levels(Lake), each = N),
                               Year = rep(2004, prod(nlevels(Lake), N)),
                               dummy = rep(0, prod(nlevels(Lake), N))))
Oxygen.pdat <- merge(Oxygen.pdat, lakeXbar)
Oxygen.pred <- predict(egmodlagged, newdata = Oxygen.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("Oxygen", colnames(Oxygen.pred$fit))
whichColsSE <- grep("Oxygen", colnames(Oxygen.pred$se.fit))
Oxygen.pdat <- cbind(Oxygen.pdat, Fitted = Oxygen.pred$fit[, whichCols], 
                     se.Fitted = Oxygen.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

## make into original limits
Oxygen.pdat <- with(Oxygen.pdat, transform(Oxygen.pdat, Fittedplus = Fitted + se.Fitted))
Oxygen.pdat <- with(Oxygen.pdat, transform(Oxygen.pdat, Fittedminus = Fitted - se.Fitted))

shiftOxygen <- attr(predict(egmodlagged, newdata = Oxygen.pdat, type = "iterms"), "constant")
Oxygen.pdatnorm <- Oxygen.pdat
Oxygen.pdatnorm <- with(Oxygen.pdatnorm, transform(Oxygen.pdatnorm, Fitted = Fitted + shiftOxygen, 
                                                   Fittedplus = Fittedplus + shiftOxygen, 
                                                   Fittedminus = Fittedminus + shiftOxygen))
Oxygen.pdatnorm <- merge(Oxygen.pdatnorm, minmaxes)
overs <- with(Oxygen.pdatnorm, which(Oxygen_ppm < minOxygen_ppm | Oxygen_ppm > maxOxygen_ppm))
Oxygen.pdatnorm <- Oxygen.pdatnorm[-overs,]

Oxygenquants <- quantile(regvarf2$Oxygen_ppm, c(.05,.95), na.rm = TRUE)

Oxygenplot <- ggplot(Oxygen.pdatnorm, aes(x = Oxygen_ppm, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=Oxygenquants[1], xmax=Oxygenquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_rug3(aes(x=Oxygen_ppm), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanox, linetype = 'dotted') +
  xlab(expression(paste("Oxygen ("~m*"g"~L^{-1}*")"))) + ylab('pH')+
  ggtitle('b.') + theme(plot.title = element_text(hjust = -0.1, size=14, face="bold"))

## for DOC
N <- 200
varWant <- c("Chl_a_ug_L", "PDOmean", "SOImean", "SPEI02", "Oxygen_ppm", "R_h", "TDN_ug_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar)))

DOC.pdat <- with(droplevels(regvarf2),
                 data.frame(`DOC_mg_L` = rep(seq(min(`DOC_mg_L`, na.rm = TRUE),
                                                 max(`DOC_mg_L`, na.rm = TRUE),
                                                 length = N),
                                             nlevels(Lake)),
                            Lake = rep(levels(Lake), each = N),
                            Year = rep(2004, prod(nlevels(Lake), N)),
                            dummy = rep(0, prod(nlevels(Lake), N))))
DOC.pdat <- merge(DOC.pdat, lakeXbar)
DOC.pred <- predict(egmodlagged, newdata = DOC.pdat, type = "terms", se.fit = TRUE)
whichCols <- grep("DOC", colnames(DOC.pred$fit))
whichColsSE <- grep("DOC", colnames(DOC.pred$se.fit))
DOC.pdat <- cbind(DOC.pdat, Fitted = DOC.pred$fit[, whichCols], 
                  se.Fitted = DOC.pred$se.fit[, whichColsSE])
limits <- aes(ymax = Fitted + se.Fitted, ymin= Fitted - se.Fitted)

# make into original limits
DOC.pdat <- with(DOC.pdat, transform(DOC.pdat, Fittedplus = Fitted + se.Fitted))
DOC.pdat <- with(DOC.pdat, transform(DOC.pdat, Fittedminus = Fitted - se.Fitted))

shiftDOC <- attr(predict(egmodlagged, newdata = DOC.pdat, type = "iterms"), "constant")
DOC.pdatnorm <- DOC.pdat
DOC.pdatnorm <- with(DOC.pdatnorm, transform(DOC.pdatnorm, Fitted = Fitted + shiftDOC, 
                                             Fittedplus = Fittedplus + shiftDOC, 
                                             Fittedminus = Fittedminus + shiftDOC))
DOC.pdatnorm <- merge(DOC.pdatnorm, minmaxes)
overs <- with(DOC.pdatnorm, which(DOC_mg_L < minDOC_mg_L | DOC_mg_L > maxDOC_mg_L))
DOC.pdatnorm <- DOC.pdatnorm[-overs,]

labdatDOC <- data.frame(x = 100, y = meanpH + 0.03, label = "mean pH")

DOCquants <- quantile(regvarf2$DOC_mg_L, c(.05,.95), na.rm = TRUE)

DOCplot <- ggplot(DOC.pdatnorm, aes(x = DOC_mg_L, y = Fitted)) +
  papertheme +
  annotate("rect", xmin=DOCquants[1], xmax=DOCquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
  geom_line() +
  geom_ribbon(aes(ymin = Fittedminus, ymax = Fittedplus), 
              alpha = 0.25) +  
  geom_rug3(aes(x=DOC_mg_L), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  scale_x_log10(breaks = c(10,20,40,60,80)) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanDOC, linetype = 'dotted') +
  xlab(expression(paste("DOC ("~m*"g"~L^{-1}*")"))) + ylab('pH') +
  ggtitle('c.') + theme(plot.title = element_text(hjust = -0.1, size=14, face="bold"))

## for chl a
## Get site-specific data/predictions
N <- 1000
varWant <- c("TDN_ug_L", "PDOmean", "SOImean", "SPEI02", "R_h", "DOC_mg_L", "Oxygen_ppm")
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
chla.pred <- predict(egmodlagged, newdata = chla.pdat, type = "terms")
whichCols <- grep("Chl", colnames(chla.pred))
chla.pdat <- cbind(chla.pdat, Fitted = rowSums(chla.pred[, whichCols]))

shiftchl <- attr(predict(egmodlagged, newdata = chla.pdat, type = "iterms"), "constant")
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
  annotate("rect", xmin=chlquants[1], xmax=chlquants[2], ymin=-Inf, ymax=Inf, alpha = 0.3, fill='gray60') +
  geom_line() + 
  geom_rug3(aes(x=Chl_a_ug_L), data = regvarf2, stat = "identity", position = "identity", 
            sides = "b", na.rm = FALSE, show.legend = NA, inherit.aes = FALSE, alpha=0.3) +
  #geom_text(data = labdatchl, aes(label = label, x = x, y = y, size = 5),
  #          show.legend = FALSE, inherit.aes = FALSE) +
  geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  geom_vline(xintercept = meanchl, linetype='dotted') +
  scale_linetype_manual(name='Lake', values = c("dotdash", "solid","solid", "solid", "solid", 
                                                "solid", "longdash")) +
  scale_colour_manual(name="Lake", values = c("#5e3c99","#b2abd2", "#b2abd2","#b2abd2", "#b2abd2",
                                              "#b2abd2", "#e66101"))+
  theme(legend.position = 'top', legend.direction = "vertical", 
        axis.text.x = element_text(angle = 45)) +
  scale_x_log10(breaks = c(5,10,25,50,100,200,300)) +
  xlab(expression(paste("Chl"~italic(a)~"("~mu*"g"~"L"^{-1}*")"))) + 
  guides(colour=guide_legend(ncol=3,nrow=4,bycol =TRUE,title.position = 'left'),
         lty=guide_legend(ncol=3,nrow=4, bycol =TRUE,title.position = 'left')) +
  ylab('pH') +
  ggtitle('a.') + theme(plot.title = element_text(hjust = -0.1, size=14, face="bold"))
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## see http://stackoverflow.com/questions/11838278/
##    plot-with-conditional-colors-based-on-values-in-r
## make this bit into soi-pdo
N <- 100
simSOI <- c(-1.1, 0.3, 1.1)
SOIgroup <- factor(rep(c('-1.1', '0.3', '1.1'), each=100, times=7))  # 7*300
reptimes <- length(simSOI) #3
preddf <- data.frame(PDOmean = rep(seq(min(regvarf2$PDOmean, na.rm=TRUE),
                                       max(regvarf2$PDOmean, na.rm=TRUE),
                                       length = N), times=reptimes),
                     SOImean = rep(simSOI, each = N)) #300
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L", "R_h", "Oxygen_ppm","DOC_mg_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

SOI.pdat <- with(droplevels(regvarf2),
                 data.frame(PDOmean = rep(preddf$PDOmean,
                                          nlevels(Lake)), # 300*7
                            SOImean = rep(preddf$SOImean,
                                          nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

SOI.pdat <- merge(SOI.pdat, lakeXbar)
SOI.pred <- predict(egmodlagged, newdata = SOI.pdat, type = "link")
SOI.pdat <- cbind(SOI.pdat, SOI.pred)
SOI.pdat$SOIgroup <- SOIgroup

names(SOI.pdat)[which(names(SOI.pdat)=='SOI.pred')] <- 'pH'

# need to take away places where PDO*SOI combo has not occurred in the data!!
toofar <- exclude.too.far(SOI.pdat$PDOmean, SOI.pdat$SOImean, 
                          regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
SOI.pdat$pH[toofar] <- NA

labdatSOI <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
SOI.pdat$SOIgroup <- factor(SOI.pdat$SOIgroup, levels = c('-1.1', '0.3', '1.1'))

SOIplot <- ggplot(SOI.pdat[SOI.pdat$Lake=='B',], aes(x = PDOmean, y = pH, group= SOImean, col=SOIgroup, 
                                                     lty=SOIgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('b.')~'      SOI')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('b.')~'      SOI')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "SOI", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH, linetype="dotted") +
  xlab('PDO') + ylab('pH')
## '#e66101','#fdb863','#b2abd2','#5e3c99' (order red:light:dark)

## now slices for PDO
N <- 100
simPDO <- c(-0.5, 0.5, 1.5)
PDOgroup <- factor(rep(c('-0.5', '0.5', '1.5'), each=100, times=7))  # 7*300
reptimes <- length(simPDO) #3
preddf <- data.frame(SOImean = rep(seq(min(regvarf2$SOImean, na.rm=TRUE),
                                       max(regvarf2$SOImean, na.rm=TRUE),
                                       length = N), times=reptimes),
                     PDOmean = rep(simPDO, each = N)) #300
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L","Oxygen_ppm","R_h","DOC_mg_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

PDO.pdat <- with(droplevels(regvarf2),
                 data.frame(SOImean = rep(preddf$SOImean,
                                          nlevels(Lake)), # 300*7
                            PDOmean = rep(preddf$PDOmean,
                                          nlevels(Lake)), # 300*7
                            Lake = rep(levels(Lake), each = superN),
                            Year = rep(2004, prod(nlevels(Lake), superN)),
                            dummy = rep(0, prod(nlevels(Lake), superN))))

PDO.pdat <- merge(PDO.pdat, lakeXbar)
PDO.pred <- predict(egmodlagged, newdata = PDO.pdat, type = "link")
PDO.pdat <- cbind(PDO.pdat, PDO.pred)
PDO.pdat$PDOgroup <- PDOgroup

names(PDO.pdat)[which(names(PDO.pdat)=='PDO.pred')] <- 'pH'

# need to take away places where PDO*PDO combo has not occurred in the data!!
toofar <- exclude.too.far(PDO.pdat$PDOmean, PDO.pdat$SOImean, regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
PDO.pdat$pH[toofar] <- NA

labdatPDO <- data.frame(x = 2.5, y = meanpH + 0.08, label = "mean pH")

#reorder factors
PDO.pdat$PDOgroup <- factor(PDO.pdat$PDOgroup, levels = c('-0.5', '0.5', '1.5'))

PDOplot <- ggplot(PDO.pdat[PDO.pdat$Lake=='B',], aes(x = SOImean, y = pH, group= PDOmean, col=PDOgroup,
                                                     lty=PDOgroup)) +
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_line() +
  scale_color_manual(name=expression(paste(bold('c.')~'      PDO')), values = c("#5e3c99", "#b2abd2", "#e66101"))+
  scale_linetype_manual(name=expression(paste(bold('c.')~ '      PDO')), values = c("solid", "solid","longdash")) +
  #scale_colour_brewer(name = "PDO", type = 'qual', palette = 'Dark2', direction=1) +
  #geom_abline(slope = 0, intercept = meanpH) +
  xlab('SOI') + ylab('pH')

## get expand grid object to plot a ggplot contour/heatmap
preddf <- with(regvarf2, expand.grid(PDOmean = seq(min(PDOmean), max(PDOmean), length = 100), 
                                     SOImean = seq(min(SOImean), max(SOImean), length = 100)))
superN <- nrow(preddf)

varWant <- c("SPEI02", "Chl_a_ug_L", "TDN_ug_L","Oxygen_ppm", "R_h", "DOC_mg_L")
lakeXbar <- with(regvarf2, do.call("rbind",
                                   lapply(split(regvarf2[, varWant], droplevels(Lake)), 
                                          colMeans, na.rm = TRUE)))
lakeXbar <- transform(lakeXbar, Lake = factor(rownames(lakeXbar))) # 7 lakes, 7 rows

comb.pdat <- with(droplevels(regvarf2),
                  data.frame(PDOmean = rep(preddf$PDOmean,
                                           nlevels(Lake)), 
                             SOImean = rep(preddf$SOImean,
                                           nlevels(Lake)),
                             Lake = rep(levels(Lake), each = superN),
                             Year = rep(2004, prod(nlevels(Lake), superN)),
                             dummy = rep(0, prod(nlevels(Lake), superN))))

comb.pdat <- merge(comb.pdat, lakeXbar, sort=FALSE)
comb.pred <- predict(egmodlagged, newdata = comb.pdat, type = "terms")

whichCols <- grep("SOI", colnames(comb.pred))
comb.pdat <- cbind(comb.pdat, Fitted = comb.pred[, whichCols])

shiftcomb <- attr(comb.pred, "constant")
comb.pdatnorm <- comb.pdat
comb.pdatnorm <- with(comb.pdatnorm, transform(comb.pdatnorm, Fitted = Fitted + shiftcomb))

#exclude things too far from real
toofar <- exclude.too.far(comb.pdatnorm$PDOmean, comb.pdatnorm$SOImean, regvarf2$PDOmean, regvarf2$SOImean, dist=0.1)
comb.pdatnorm$pH <- comb.pdatnorm$Fitted
comb.pdatnorm$pH[toofar] <- NA

names(comb.pdat)[which(names(comb.pdat)=='SOI.pred')] <- 'pH'
comboplot <- ggplot(comb.pdatnorm, aes(x = SOImean, y = PDOmean, z=pH)) + #, z=Fitted
  theme_bw(base_size=14, base_family = 'Arial') +
  theme(legend.position='top') +
  geom_raster(aes(fill=pH)) + # change to turn grey background into nothing
  scale_fill_viridis(na.value='transparent', option= "plasma", 
                     name=expression(paste(bold('a.')~'      pH'))) +
  geom_point(data=regvarf2, aes(x=SOImean, y=PDOmean, z=NULL)) +
  theme(legend.key.width=unit(2,"cm")) +
  geom_vline(linetype='dashed', xintercept = c(-1.1, 0.3, 1.1)) +
  geom_abline(slope = 0, intercept = c(-0.5, 0.5, 1.5), linetype='dashed') +
  ylab("PDO") + xlab("SOI")

## save objects:
plots <- list(Oxygenplot, chlaplot, DOCplot, SOIplot, PDOplot, comboplot, speiplot, TDNplot)
plotnames <- c('Oxygenplot', 'chlaplot', 'DOCplot', 'SOIplot', 'PDOplot', 'comboplot', 'SPEIplot', 
               'TDNplot')

invisible( # this means I don't get the list [[1:3]] returned on screen
  lapply(
    seq_along(plots), 
    function(x) ggsave(filename=paste0("../docs/private/gam-plots-lagged-", plotnames[x], ".pdf"), 
                       plot=plots[[x]], scale=0.6) # width=7, height=5, units = 'in'
  ) )

## arrange plots
climgam <- grid.arrange(comboplot, SOIplot, PDOplot, speiplot, ncol = 2, 
                        layout_matrix = cbind(c(1,1,2,4), c(1,1,3,4)))
metabplots <- grid.arrange(chlaplot, Oxygenplot, DOCplot, ncol = 2,
                           layout_matrix=cbind(c(1,1), c(2,3)))
ggsave("../docs/private/ph-metab-lagged.pdf", metabplots, scale=1.1, width=10) #width=28, height=18, units = 'cm'
ggsave("../docs/private/ph-clim-lagged.pdf", climgam, scale=0.9, width = 8, height=13)
ggsave("../docs/private/ph-TDN-lagged.pdf", TDNplot, scale=0.9, width = 6, height=4)
ggsave("../docs/private/gam-plots-co2plot.pdf", co2plot, scale=1.1, width = 7, height=5.5)

## time plots
## ==================================================================================
## Testing effect by time plots
## ==================================================================================
testing1 <- predict(egmodlagged, type = 'terms')
testing <- as.data.frame(testing1)
tosum <- grep("Chl", colnames(testing))
chleffect <- rowSums(testing[,tosum], na.rm = TRUE)
testing <- testing[,-tosum]
testing$Chl <- chleffect
names(testing) <- c("DOC", "Oxygen", "PDO-SOI", 'LakeYear','SPEI', 'Chl')
testing$Date <- regvarf2$Date #is this ok? assuming order is preserved
testing$Year <- regvarf2$Year
testing$Lake <- as.character(regvarf2$Lake)
testing$Lake[which(testing$Lake == "WW")] <- "W"
testing$DOY <- regvarf2$DOY
testing$Month <- regvarf2$Month

## messy... so much variability and outliers.. perhaps do a mean effect by lake and month?
testsplit <- with(testing, split(testing, list(Lake, Year, Month)))
testsplit <- testsplit[sapply(testsplit, function(x) dim(x)[1]) > 0] #remove empties for lakes R, E etc.
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(Lake = df['Lake'][1,], Year = df['Year'][1,],Month = df['Month'][1,]), t(means))
}

varWant = c("Oxygen", "SPEI", "LakeYear", "Chl", "PDO-SOI", "DOC")
testmeans <- do.call(rbind, lapply(testsplit, mycolMeans, cols=varWant))
rownames(testmeans) <- NULL

oxes <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
               aes(x=Month, y=Oxygen, group = Month)) + #fill=Lake, alpha=0.2
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap('Lake', nrow=2) +
  #ylim(c(-0.25, 0.5))+
  ylab("Oxygen effect") +
  theme(axis.title.x=element_blank())
chls <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
               aes(x=Month, y=Chl, group = Month)) + #, fill=Lake, alpha=0.2
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap('Lake', nrow=2) +
  #ylim(c(-0.25, 0.5))+
  ylab(expression(paste("Chl"~italic('a')~"effect"))) +
  theme(axis.title.x=element_blank())
docs <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
               aes(x=Month, y=DOC, group = Month)) + #, fill=Lake, alpha=0.2
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap('Lake', nrow=2) +
  #ylim(c(-0.25, 0.5))+
  ylab("DOC effect") +
  theme(axis.title.x=element_blank())

speis <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
                aes(x=Month, y=SPEI, group = Month, alpha=0.2)) + #fill=Lake, 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  #facet_wrap('Lake', nrow=2) +
  #ylim(c(-0.25, 0.5))+
  ylab("SPEI effect") +
  theme(axis.title.x=element_blank())

climate <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
                  aes(x=Month, y=`PDO-SOI`, group = Month, alpha=0.2)) + #fill=Lake, 
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  #facet_wrap('Lake', nrow=2) +
  #ylim(c(-0.25, 0.5))+
  ylab("PDO*SOI effect") +
  theme(axis.title.x=element_blank())
lakeyear <- ggplot(testmeans[testmeans$Month < 10 & testmeans$Month > 3,], 
                   aes(x=Month, y=LakeYear, group = Month, fill=Lake, alpha=0.2)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-0.05, ymax=0.05, alpha=0.3, fill="grey60") +
  geom_boxplot(show.legend = FALSE) +
  facet_wrap('Lake', nrow=2)+
  #ylim(c(-0.25, 0.5)) +
  theme(axis.title.x=element_blank())
alltimers <- plot_grid(oxes, chls, docs,speis, climate, ncol = 2, rel_heights= c(2, 1), labels='AUTO')
alltimers <- grid.arrange(oxes, chls, docs,speis, climate, ncol=2, 
                          layout_matrix=cbind(c(1,1,2,2,3,3,4), c(1,1,2,2,3,3,5)),
                          top="Effect of predictors by month (5:9)")
ggsave("../docs/private/alltimers.pdf", alltimers, width=15, height=25, units="cm")
## try melting
testmelt <- melt(testmeans, id.vars = c('Lake', 'Year', 'Month'), variable.name = "Variable")

ggplot(testmelt, aes(x=factor(Month), y= value, fill=Lake)) +
  geom_boxplot(aes(y=value, x=factor(Month))) +
  facet_wrap('Variable')

ggplot(testmelt, aes(x=factor(Month), y= value)) +
  geom_boxplot(aes(y=value, x=factor(Month))) +
  geom_smooth(method="gam",  se=FALSE, na.rm = TRUE, aes(group=Lake, col=Lake))  +
  facet_wrap('Variable')

