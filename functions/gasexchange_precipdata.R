## adapted scripts from bottomoftheheap/...canadian climate data
## we want precipitation data, but not avaiable for hourly csv downloads so need separate scripts
##    from the hourly wind et al
## note that character denotations in columns are: 
## A = Accumulated
## C = Precipitation occurred, amount uncertain
## E = Estimated
## F = Accumulated and estimated
## L = Precipitation may or may not have occurred
## M = Missing
## N = Temperature missing but known to be > 0
## S = More than one occurrence
## T = Trace
## Y = Temperature missing but known to be   < 0
## [empty] = No data available
## ^ = The value displayed is based on incomplete data
## † = Data for this day has undergone only preliminary quality checking
## ‡ = Partner data that is not subject to review by the National Climate Archives

## genURLS function
genURLS <- function(id, start, end) {
  years <- seq(start, end, by = 1)
  nyears <- length(years)
  years <- rep(years, each = 1) #12
  URLS <- paste0("http://climate.weather.gc.ca/climateData/bulkdata_e.html?timeframe=2&Prov=SK&StationID=",
                 id,
                 "&dlyRange=1970-05-01|2015-11-22&cmdB1=Go",
                 "&Year=",
                 years,
                 "&Month=11", # ps hourly code had months object for target but tested one year's URL
                 #    and the month set to 11 still brought in the whole year.. also timeframe = 1 vs 2
                 #    seems to have no effect on what happens..!?
                 "&Day=27",
                 "&format=csv",
                 "&stationID=",
                 id)
  list(urls = URLS, ids = rep(id, nyears), years = years) 
}
## getData function

getData <- function(stations, folder, verbose = TRUE, delete = TRUE) {
  ## form URLS
  urls <- lapply(seq_len(NROW(stations)),
                 function(i, stations) {
                   genURLS(stations$StationID[i],
                           stations$start[i],
                           stations$end[i])
                 }, stations = stations)
  
  ## check the folder exists and try to create it if not
  if (!file.exists(folder)) {
    warning(paste("Directory:", folder,
                  "doesn't exist. Will create it"))
    fc <- try(dir.create(folder))
    if (inherits(fc, "try-error")) {
      stop("Failed to create directory '", folder,
           "'. Check path and permissions.", sep = "")
    }
  }
  
  ## Extract the data from the URLs generation
  URLS <- unlist(lapply(urls, '[[', "urls"))
  sites <- unlist(lapply(urls, '[[', "ids"))
  years <- unlist(lapply(urls, '[[', "years"))
  
  ## filenames to use to save the data
  fnames <- paste(sites, years, "data.csv", sep = "-") 
  fnames <- file.path(folder, fnames)
  
  nfiles <- length(fnames)
  
  ## set up a progress bar if being verbose
  if (isTRUE(verbose)) {
    pb <- txtProgressBar(min = 0, max = nfiles, style = 3)
    on.exit(close(pb))
  }
  
  out <- vector(mode = "list", length = nfiles)
  cnames <- c("Date/Time", "Year", "Month","Day", "Data Quality", "Max Temp (degC)",
              "Max Temp Flag", "Min Temp (degC)", "Min Temp Flag", 
              "Mean Temp (degC)", "Mean Temp Flag", "Heat Deg Days (degC)",
              "Heat Deg Days Flag", "Cool Deg Days (degC)",
              "Cool Deg Days Flag", "Total Rain (mm)", "Total Rain Flag",
              "Total Snow (cm)", "Total Snow Flag", "Total Precip (mm)",
              "Total Precip Flag", "Snow on Grnd (cm)", "Snow on Grnd Flag",
              "Dir of Max Gust (10s deg)", "Dir of Max Gust Flag", 
              "Spd of Max Gust (km/h)", "Spd of Max Gust Flag")
  
  for (i in seq_len(nfiles)) {
    curfile <- fnames[i]
    
    ## Have we downloaded the file before?
    if (!file.exists(curfile)) {    # No: download it
      dload <- try(download.file(URLS[i], destfile = curfile, quiet = TRUE))
      if (inherits(dload, "try-error")) { # If problem, store failed URL...
        out[[i]] <- URLS[i]
        if (isTRUE(verbose)) {
          setTxtProgressBar(pb, value = i) # update progress bar...
        }
        next                             # bail out of current iteration
      }
    }
    
    ## Must have downloaded, try to read file
    ## skip first 16 rows of header stuff
    ## encoding must be latin1 or will fail - may still be problems with character set
    cdata <- try(read.csv(curfile, skip = 25, encoding = "latin1"), silent = TRUE)
    
    ## Did we have a problem reading the data?
    if (inherits(cdata, "try-error")) { # yes handle read problem
      ## try to fix the problem with dodgy characters
      cdata <- readLines(curfile) # read all lines in file
      cdata <- gsub("\x87", "x", cdata) # remove the dodgy symbol for partner data in Data Quality
      cdata <- gsub("\xb0", "deg", cdata) # remove the dodgy degree symbol in column names
      cdata <- gsub("\x86", "NA", cdata) # remove x86 causing trouble for daily data!!
      writeLines(cdata, curfile)          # write the data back to the file
      ## try to read the file again, if still an error, bail out
      cdata <- try(read.csv(curfile, skip = 25, encoding = "latin1"), silent = TRUE)
      if (inherits(cdata, "try-error")) { # yes, still!, handle read problem
        if (delete) {
          file.remove(curfile) # remove file if a problem & deleting
        }
        out[[i]] <- URLS[i]    # record failed URL...
        if (isTRUE(verbose)) {
          setTxtProgressBar(pb, value = i) # update progress bar...
        }
        next                  # bail out of current iteration
      }
    }
    
    ## Must have (eventually) read file OK, add station data
    cdata <- cbind.data.frame(StationID = rep(sites[i], NROW(cdata)),
                              cdata)
    names(cdata)[-1] <- cnames
    out[[i]] <- cdata
    
    if (isTRUE(verbose)) { # Update the progress bar
      setTxtProgressBar(pb, value = i)
    }
  }
  
  out                                 # return
}



