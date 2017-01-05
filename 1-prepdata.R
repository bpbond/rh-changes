# Prep data
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "1-prepdata.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# 1. Get SRDB data, filter for XXXXX

# 2. Match with climate data

# 3. Match with GPP data


crufile <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.tmp.dat.nc.gz"

extract_cru_data <- function(crufile, lon, lat, midyear, nyears, file_startyear = 1901) {
  
  assert_that(length(lon) == length(lat))
  
  compressed <- grepl("gz$", crufile)
  if(compressed) {
    printlog("Decompressing", crufile)
    ncfile <- R.utils::gunzip(crufile, remove = FALSE, overwrite = TRUE)
  } else {
    ncfile <- crufile
  }
  
  library(raster) # 2.5.8
  nc <- brick(ncfile)

  # Find nearest neighbors for all lon/lat pairs
  x <- rep(NA_real_, length(lon))
  for(i in seq_along(lon)) {
    sp <- SpatialPoints(cbind(lon[i], lat[i]))
    midyear_layer = (midyear - file_startyear) * 12
    start_layer = midyear_layer - nyears / 2 * 12
    
    nc %>%
      extract(sp, layer = start_layer, nl = nyears * 12) %>%
      mean(na.rm = TRUE) ->
      x
  }
  
  # Clean up
  if(compressed) {
    printlog("Removing", ncfile)
    file.remove(ncfile)
  }
  
  x
}



printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
