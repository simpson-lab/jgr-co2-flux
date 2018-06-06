## code for appendix figures + reviewer response plots following revision

library('ggplot2')
library("mgcv")
library("plyr")
library("dplyr")
library('reshape2')
library('gridExtra')
library('extrafont')

## read in regvars and co2expl and models
regvars <- readRDS("../data/private/regvars.rds")
expl <- readRDS("../data/private/co2explained.rds")

weathers <- readRDS('../data/weathers.rds')

egmodlagged <- readRDS("../data/private/egmodlaggedsimp.rds")
egmodtdn <- readRDS("../data/private/egmodlaggedtdn.rds")

## create a theme to save linespace in plots
papertheme <- theme_bw(base_size=10, base_family = 'Arial') +
  theme(legend.position='top')

## create function for shared legends
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## my lazy labeler fix
mylabel_parsed <- function (labels, multi_line = FALSE) 
{
  labels <- label_value(labels, multi_line = multi_line)
  if (multi_line) {
    lapply(unname(labels), lapply, function(values) {
      c(parse(text = as.character(values)))
    })
  }
  else {
    lapply(labels, function(values) {
      values <- paste0("list(", values, ")")
      lapply(values, function(expr) c(parse(text = expr)))
    })
  }
}
## put data together for facetted plot
nrow(expl)
nrow(regvars)

names(expl)[names(expl) %in% c("YEAR","LAKE")] <- c("Year","Lake")
names(expl)[names(expl) %in% names(regvars)] # some overlap

extras <- names(expl)[!names(expl) %in% names(regvars)]
regvars <- merge(regvars, expl[,c("Year","Lake","Month","Date",extras)])

subdf <- regvars[,c("Date","Lake","Year","Month","DOY","Chl_a_ug_L", "TDN_ug_L","DOC_mg_L","TDP_ug_L",
                     "TIC_mg_L"  , "pH_surface")]
monthmeans <- aggregate(cbind(Chl_a_ug_L, TDN_ug_L,DOC_mg_L,TDP_ug_L,
          TIC_mg_L  , pH_surface)~Month + Lake + Year, data=subdf, mean, na.rm=TRUE, na.action=NULL)

monthmeans <- with(subdf, split(subdf, list(Month, Lake, Year), drop=TRUE)) 
monthmeans <- monthmeans[sapply(monthmeans, function(x) dim(x)[1]) > 0] 
mycolMeans <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  means <- colMeans(subdf, na.rm = TRUE)
  cbind(data.frame(Lake = df['Lake'][1,], Month = df['Month'][1,], Year = df['Year'][1,], t(means)))
}

monthmeans <- do.call(rbind, lapply(monthmeans, mycolMeans, cols=c("Chl_a_ug_L", "TDN_ug_L","DOC_mg_L",
                                                                   "TDP_ug_L","TIC_mg_L"  , "pH_surface")))
monthmelt <- melt(monthmeans, id.vars = c("Lake","Year","Month"))
monthmelt$realname <- factor(monthmelt$variable,
                        labels=c('plain(Chl)~italic(a)~mu~L^{-1}', 'plain(TDN)~mu~N~L^{-1}', 
                                 'plain(DOC)~plain(mg)~C~L^{-1}', 'plain(TDP)~mu~P~L^{-1}', 
                                 'plain(DIC)~mu~C~L^{-1}',"pH"))
monthmelt <- droplevels(monthmelt)
monthmelt$shortlake <- monthmelt$Lake
monthmelt$Lake <- factor(monthmelt$Lake, labels=c("Katepwa","Last Mountain","Buffalo Pound",
                                                      "Crooked","Diefenbaker","Wascana","Pasqua"))                         

## plot..
appendixplot <- ggplot(monthmelt[monthmelt$Month > 4 & monthmelt$Month <9,], 
                       aes(y=value, x=Month,group=Lake)) +
  papertheme +
  geom_point(aes(col=Lake), alpha=0.8, size=.8) +
  geom_line(aes(col=Lake)) +
  facet_wrap(realname~Year, scales = "free_y",labeller = mylabel_parsed, dir="h")+ 
  theme(strip.background = element_rect(fill="white", colour = "white")) +
  guides(col=guide_legend(nrow=1, override.aes = list(alpha=1, size=2))) +
  scale_color_manual(values=c('#8c510a','#d8b365','#fc8d59','black','#c7eae5','#5ab4ac','#01665e')) 

ggsave("../docs/private/flux-appendix.pdf", appendixplot, width=18, height=10)

## plot diagnostics for the flux model
# check tdn mod if want
#gam.check(egmodtdn)
#summary(egmodtdn)

type <- "deviance"  ## "pearson" & "response" are other valid choices
resid <- residuals(egmodlagged, type = type)
linpred <- napredict(egmodlagged$na.action, egmodlagged$linear.predictors)
observed.y <- napredict(egmodlagged$na.action, egmodlagged$y)

## change qq plot to capitalised axis labels
myQQ <- function (object, rep = 0, level = 0.9, s.rep = 10, 
                  type = c("deviance", 
                           "pearson", "response"), pch = ".", rl.col = 2, rep.col = "gray80", ...) {
  type <- match.arg(type)
  ylab <- "Deviance residuals"
  if (inherits(object, c("glm", "gam"))) {
    if (is.null(object$sig2)) 
      object$sig2 <- summary(object)$dispersion
  }
  else stop("object is not a glm or gam")
  object$na.action <- NULL
  D <- residuals(object, type = type)
  if (object$method %in% c("PQL", "lme.ML", "lme.REML", "lmer.REML", 
                           "lmer.ML", "glmer.ML")) {
    qqnorm(D, ylab = ylab, pch = pch, ...)
    return()
  }
  lim <- Dq <- NULL
  if (rep == 0) {
    fam <- fix.family.qf(object$family)
    if (is.null(fam$qf)) 
      rep <- 50
    level <- 0
  }
  n <- length(D)
  if (rep > 0) {
    fam <- fix.family.rd(object$family)
    if (!is.null(fam$rd)) {
      dm <- matrix(0, n, rep)
      for (i in 1:rep) {
        yr <- fam$rd(object$fitted.values, object$prior.weights, 
                     object$sig2)
        object$y <- yr
        dm[, i] <- sort(residuals(object, type = type))
      }
      Dq <- quantile(as.numeric(dm), (1:n - 0.5)/n)
      alpha <- (1 - level)/2
      if (alpha > 0.5 || alpha < 0) 
        alpha <- 0.05
      if (level > 0 && level < 1) 
        lim <- apply(dm, 1, FUN = quantile, p = c(alpha, 
                                                  1 - alpha))
      else if (level >= 1) 
        lim <- level
    }
  }
  else {
    U <- (1:n - 0.5)/n
    if (!is.null(fam$qf)) {
      dm <- matrix(0, n, s.rep)
      for (i in 1:s.rep) {
        U <- sample(U, n)
        q0 <- fam$qf(U, object$fitted.values, object$prior.weights, 
                     object$sig2)
        object$y <- q0
        dm[, i] <- sort(residuals(object, type = type))
      }
      Dq <- sort(rowMeans(dm))
    }
  }
  if (!is.null(Dq)) {
    qqplot(Dq, D, ylab = ylab, xlab = "Theoretical quantiles", 
           ylim = range(c(lim, D)), pch = pch, ...)
    abline(0, 1, col = rl.col)
    if (!is.null(lim)) {
      if (level >= 1) 
        for (i in 1:rep) lines(Dq, dm[, i], col = rep.col)
      else {
        n <- length(Dq)
        polygon(c(Dq, Dq[n:1], Dq[1]), c(lim[1, ], lim[2, 
                                                       n:1], lim[1, 1]), col = rep.col, border = NA)
      }
      abline(0, 1, col = rl.col)
    }
    points(Dq, sort(D), pch = pch, ...)
    return(invisible(Dq))
  }
  else qqnorm(D, ylab = ylab, pch = pch, ...)
}
pdf("../docs/private/appendix2.pdf", onefile=TRUE)
op <- par(mfrow = c(2,2))
myQQ(egmodlagged, rep = 0, level = 0.9, type = type, rl.col = 2, 
       rep.col = "gray80", main="QQ plot")
hist(resid, xlab = "Residuals", main = "Histogram of residuals")

plot(linpred, resid, main = "Residuals vs linear predictor", 
     xlab = "Linear predictor", ylab = "Residuals")
plot(fitted(egmodlagged), observed.y, xlab = "Fitted Values", 
     ylab = "Response", main = "Response vs Fitted Values")
par(op)
dev.off()

## save text of model output
## save summary as txt document
egmodsum <- summary(egmodlagged)
sink("../docs/private/egmodsummary.txt")
egmodsum
sink()

## hypothetical prediction:
newdata <-  data.frame(Chl_a_ug_L = c(200, 5),DOC_mg_L = c(50,5), Oxygen_ppm = c(8,3), 
                       SPEI02 = c(-1.5, 0.2), PDOmean=c(1.5, 0.2),
                       SOImean = c(-1.1,-1.1), Year=c(2009, 2012), Lake=c("WW","B"), dummy=c(1,1))
newdata <-  data.frame(Chl_a_ug_L = 40,DOC_mg_L = 12, Oxygen_ppm = 9, 
                       SPEI02 = 0.2, PDOmean=0.01,
                       SOImean = 0.21, Year=1995, Lake=c("WW","B","D","L","K","C"), dummy=1)

predict(egmodlagged, newdata=newdata, se.fit = TRUE)
#summary(predict(egmodlagged))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#  7.765   8.648   8.861   8.881   9.108  10.057     169 
#summary(egmodlagged$y)
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#7.050   8.590   8.820   8.878   9.200  10.940 

## plot of predicted vs observed pH over time
regvarf <- regvars
regvarf <- merge(regvarf, weathers)

regvarf <- transform(regvarf, Year = as.factor(Year)) # make Year into factor for re
regvarf2 <- regvarf
regvarf2$`Chl_a_ug_L`[regvarf2$`Chl_a_ug_L` <=0 ] <- 0
regvarf2$`DOC_mg_L`[regvarf2$`DOC_mg_L` <=0 ] <- 0
regvarf2$`Chl_a_ug_L` <- regvarf2$`Chl_a_ug_L` + 1
regvarf2$`DOC_mg_L` <- regvarf2$`DOC_mg_L` + 1
regvarf2 <- transform(regvarf2, dummy = rep(1, nrow(regvarf2)))
regvarf2$preds <- predict(egmodlagged, newdata=regvarf2, type = "response")

realnames <- data.frame(Lake = c(as.character(unique(regvarf2$Lake))), 
                        LakeName = c("BuffaloPound","Crooked","Diefenbaker","Katepwa",
                                     "LastMountain","Pasqua","Wascana"))
regvarf2 <- merge(regvarf2, realnames)

regvarf2$error <- regvarf2$pH_surface - regvarf2$preds

summaries <- ddply(regvarf2, .(LakeName, Year, Month), dplyr::summarise, meanpH = mean(pH_surface, na.rm=TRUE),
                   meanPred = mean(preds, na.rm = TRUE), meanError = mean(error, na.rm = TRUE))

testdf <- reshape2::melt(regvarf2, id.vars = c('LakeName', 'Date', 'Month', 'Year'), 
                         measure.vars = c('preds', 'pH_surface','error'))
testdf$Class <- ifelse(testdf$variable == "preds", "Predicted", 
                       ifelse(testdf$variable == "error", "Error", "Measured"))

testdf2 <- melt(summaries, measure.vars = c('meanpH',"meanPred","meanError"))
testdf2$Class <- ifelse(testdf2$variable == "meanPred", "Predicted", 
                       ifelse(testdf2$variable == "meanError", "Error", "Measured"))
testdf2 <- testdf2[-which(testdf2$Year %in% c(1994, 1996, 2002)),]
testdf2$cohort <- ifelse(testdf2$Year %in% c(2006:2014), "one",'two')

firstyears <- ggplot(testdf2[testdf2$cohort == 'two' & testdf2$Class!="Error",], 
                     aes(y=value, x=Month, group=LakeName, col=LakeName)) +
  papertheme +
  geom_point(size=1, alpha=0.8) +
  facet_grid(Class~Year) + 
  scale_color_manual(values=c('#8c510a','#d8b365','#fc8d59','black','#c7eae5','#5ab4ac','#01665e')) +
  theme(axis.text.x = element_text(angle=45)) +
  geom_smooth(aes(y=value), se=FALSE, size=0.5, alpha=1) + ylab('pH')

secondyears <-
ggplot(testdf2[testdf2$cohort == 'one' & testdf2$Class!="Error" & testdf2$Month > 4 & testdf2$Month <10,], aes(y=value, x=Month, group=LakeName, col=LakeName)) +
  papertheme +
  geom_point(size=1, alpha=0.8) +
  facet_grid(Class~Year) + 
  scale_color_manual(values=c('#8c510a','#d8b365','#fc8d59','black','#c7eae5','#5ab4ac','#01665e')) +
  theme(axis.text.x = element_text(angle=45)) +
  geom_smooth(aes(y=value), se=FALSE, size=0.5, alpha=1) +ylab('pH')

allyears <- ggplot(testdf2[testdf2$Class=="Measured" & testdf2$Month > 4 & testdf2$Month <10,], 
       aes(y=value, x=Month, group=LakeName, col=LakeName)) +
  papertheme +
  geom_point(size=1, alpha=0.8) +
  facet_grid(LakeName~Year) + 
  scale_color_manual("Lake", values=c('#8c510a','#d8b365','#fc8d59','black','#542788','#5ab4ac','#01665e')) +
  theme(axis.text.x = element_text(angle=45),
        legend.position = 'none', strip.text = element_text(size=6)) +
  geom_smooth(data=testdf2[testdf2$Class=="Predicted" & testdf2$Month > 4 & testdf2$Month <10,], 
              aes(y=value), se=FALSE, size=0.5, alpha=1) +ylab('pH: predicted (line) & measured (points)')

#allyears <- grid.arrange(firstyears, secondyears, ncol=1)
#ggsave(plot = allyears, filename = "../docs/private/truevspred.pdf", height=7, width=10)
ggsave(plot = allyears, filename = "../docs/private/truevspred-allyears.pdf", height=6, width=10)

## =================================================================================================
## indication of intra-annual variation between May and Sep for each lake all years
myRange <- function(dat) {range=max(dat, na.rm=TRUE)-min(dat,na.rm = TRUE)}

vardf <- ddply(regvars[regvars$Month>4 & regvars$Month <10,], .(Lake, Year), summarise, chl=myRange(Chl_a_ug_L),
               tdn=myRange(TDN_ug_L), doc=myRange(DOC_mg_L), pH=myRange(pH_surface),
               tdp=myRange(TDP_ug_L), resp=myRange(R_h))
is.na(vardf) <- sapply(vardf, is.infinite)


vardfmad <- ddply(regvars[regvars$Month>4 & regvars$Month <10,], .(Lake, Year), summarise, chl=mad(Chl_a_ug_L, na.rm = TRUE),
               tdn=mad(TDN_ug_L, na.rm = TRUE), doc=mad(DOC_mg_L, na.rm = TRUE), pH=mad(pH_surface, na.rm = TRUE),
               tdp=mad(TDP_ug_L, na.rm = TRUE), resp=mad(R_h, na.rm = TRUE))
is.na(vardfmad) <- sapply(vardfmad, is.infinite)

deg <- 29.531
vardf <- droplevels(vardf)
vardfmad <- droplevels(vardfmad)
lakelist <- levels(vardfmad$Lake)

# need to create the difference in values to create a legible rose plot close to origin
switch <- vardfmad
varmelt <- melt(switch, id.vars = c("Lake","Year"))
varmelt$Lake <- as.character(varmelt$Lake)
varmelt <- by(varmelt, list(varmelt$Year, varmelt$variable), function(x) {x <- x[order(x$value),]
x$order <- 1:nrow(x)
return(x)})
varmelt <- do.call(rbind, varmelt)

vartest <- by(varmelt, list(varmelt$Year, varmelt$variable), function(x) {
  x$diff <- c(x$value[1], diff(x$value))
  return(x)
})
vartest <- do.call(rbind,vartest)
vartest$value[vartest$variable=="tdp"& vartest$Lake=="WW"&vartest$value > 1000] <- NA

## order data so that geom_col with identity position plots everything is the correct order
ordered_data <- vartest[order(-vartest$value), ]
ordered_data$Lake <- factor(ordered_data$Lake)
ordered_data$Lake <- mapvalues(ordered_data$Lake, from = c(levels(ordered_data$Lake)), 
                               to = c("Buffalo Pound","Crooked","Diefenbaker","Katepwa","Last Mountain",
                                      "Pasqua","Wascana"))

## create function to plot all that i want
plotall <- function(df, varname) {
  ymax <- max(df$value[df$variable==varname], na.rm = TRUE)
  yinter <- ymax/8
  segments <- data.frame(x=seq(1,21,1), xend=seq(1,21,1), y=rep(200), yend=rep(ymax))
  
  plott <- ggplot(df[df$variable==varname,], aes(x=factor(Year), y=value)) +#, fill=Lake
    theme_bw(base_size=9, base_family = 'Arial') +
    geom_hline(yintercept = seq(0,ymax,yinter), size=0.3, col="grey60",lty=3) +
    geom_vline(xintercept=seq(1,21,1), size=0.3, col='grey30', lty=2) +
    geom_col(position='identity', width=1, size=0.5, aes(fill=Lake)) +
    scale_x_discrete() +
    coord_polar(start=-0.15) +
    scale_fill_manual(name="Lake", 
                      values=c('#1b9e77','#d95f02','#7570b3','#e7298a','#66a61e','#e6ab02','#a6761d')) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          axis.title.x = element_blank(),axis.title.y=element_blank(), plot.margin = unit(c(17,0,0,0),'pt'),
          plot.title = element_text(hjust = 0.5, size=9)) }

## specify what I want
vartitles <- c(expression("Chl"~a~mu*g~L^{-1}),expression("TDN"~mu*g~L^{-1}),
               expression("DOC"~mg~L^{-1}),"pH",expression("TDP"~mu*g~L^{-1}),
               expression("R"~O[2]~h^{-1}))
iwant <- c(levels(ordered_data$variable))

plotlist <- lapply(iwant, plotall, df=ordered_data)
plotlist <- Map(function(x,y) {x <- x+
  ggtitle(y)}, plotlist, vartitles)
nullplot <- plotlist[[1]] + theme(legend.position="bottom", legend.direction = "horizontal") +
  guides(fill=guide_legend(nrow=1, ncol = 7)) 
plotlist <- lapply(plotlist, function(x) {x + theme(legend.position = "none") + 
    cowplot::panel_border(remove=TRUE)})
plotlist <- lapply(plotlist, ggplotGrob)
plottest <- lapply(plotlist, function(x) {x$widths <- plotlist[[1]]$widths
x$heights <- plotlist[[1]]$heights
return(x)})

## create plots without legend
justplots <- 
  cowplot::plot_grid(plotlist = plottest, ncol=3, labels = c(letters[1:length(plotlist)]), 
                     label_fontface = "plain", hjust=-1.3, vjust=2.2) #vjust = 0.9, , hjust=-0.5
## create 'plot' of legend
legendplot <- cowplot::get_legend(nullplot)

## plot the result
p <- 
cowplot::plot_grid( justplots, legendplot, ncol = 1, rel_heights = c(1, .3))
ggsave("../docs/private/rangeplot.pdf", width = 10,height=7)


## ================================================================================================
## what is the variability in pH when set against DIC?

dicsum <- ddply(regvars, .(Lake,Year), summarise, 
                pHvar = max(pH_surface, na.rm = TRUE)-min(pH_surface, na.rm=TRUE),
                meanDIC=mean(TIC_mg_L, na.rm = TRUE), meanpH = mean(pH_surface, na.rm = TRUE))

ggplot(dicsum, aes(y=meanDIC, x=pHvar, group=Lake)) +
  papertheme +
  geom_point() +
  facet_wrap(~Lake) + 
  ylab("Mean annual DIC (mg/L)") + xlab("Annual range in pH")

ggplot(dicsum, aes(y=meanDIC, x=meanpH, group=Lake)) +
  papertheme +
  geom_point() +
  facet_wrap(~Lake) + 
  ylab("Mean annual DIC (mg/L)") + xlab("Mean annual pH")
