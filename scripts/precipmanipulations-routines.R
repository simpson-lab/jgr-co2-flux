## archaically used gasexchange_precipdata.R to source regina info for routines analysis,
##    now updated to 24,5,2016 script fromthebottomoftheheap since met office updated web things
## Make info suitable to combine with other variables

## source necessary functions:
source("../functions/getData.R")
source("../functions/genURLS.R")

## my station data frame: for regina gilmour as per heather's searches
stations <- data.frame(StationID = c(3007), start = c(1994), end = c(2014))

## create precipitation data frame
if (!dir.exists('../data/precipdata/routines')) {dir.create('../data/precipdata/routines',
                                                            recursive = TRUE)}
precip <- getData(stations, folder = "../data/precipdata/routines", timeframe='daily')
precips <- do.call("rbind", precip)
rownames(precips) <- NULL

names(precips)[grep("Date", names(precips))] <- "Date"

## make Date into R date object
precips <- transform(precips, Date= as.Date(as.character(Date), format = "%Y-%m-%d", 
                                            tz="Canada/Saskatchewan"))

saveRDS(precips, "../data/precip.rds") # save for future

## hmmm... annual snow total? wrapper for which col want
mycolSums <- function(df, cols) {
  df <- as.data.frame(df)
  subdf <- subset(df, select = cols)
  sum <- colSums(subdf, na.rm = TRUE)
  cbind(data.frame(Year = df['Year'][1,], t(sum)))
}

precipsplit <- with(precips, split(precips, list(Year), drop = TRUE))
snowtotal <- lapply(precipsplit, mycolSums, cols = c("Total.Snow..cm."))
snowtotal <- do.call("rbind", snowtotal)
rownames(snowtotal) <- NULL
