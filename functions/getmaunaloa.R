## get mauna loa CO2

download.file("ftp://aftp.cmdl.noaa.gov/data/trace_gases/co2/flask/surface/co2_mlo_surface-flask_1_ccgg_month.txt", 
              destfile = "../data/maunaloa.csv")
ml <- read.csv("../data/maunaloa.csv", skip = 70, sep = "", encoding = "latin1", header = FALSE,
               col.names = c("Site", "Year", "Month", "pCO2"))
write.csv(ml, "../data/maunaloa.csv", row.names = FALSE)

