# Prep the ESA-CCI data (per Referee 2 suggestion)
# Ben Bond-Lamberty July 2017
#
# Unzip the folders for each year; concatenate the daily data; compute mean and var
# Uses CDO version 1.7.1rc1 for processing netcdf files
# This script uses a lot (tens of GBs) of disk space for temporary files


source("0-functions.R")

SCRIPTNAME  	<- "1-cci.R"
PROBLEM       <- FALSE
SKIP_EXISTING <- TRUE

# Downloaded 6 June 2017 from http://data.ceda.ac.uk/neodc/esacci/soil_moisture/data/daily_files/COMBINED/v02.2/
CCI_DATA <- "~/Data/ESA-CCI/"

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

tmpdir <- "~/Desktop/tmp/"
if(!dir.exists(tmpdir)) dir.create(tmpdir)

for(f in list.files(CCI_DATA, "*.zip", full.names = TRUE)) {
  printlog(SEPARATOR)
  
  yrpos <- regexpr("[0-9]{4}", f)[1]
  year <- substr(f, yrpos, yrpos + 3)
  yearfile <- file.path(tmpdir, paste0(year, ".nc"))
  meanfile <- file.path(CCI_DATA, paste0(year, "-mean.nc"))
  stdfile <- file.path(CCI_DATA, paste0(year, "-std.nc"))
  
  if(file.exists(meanfile) && file.exists(stdfile) && SKIP_EXISTING) {
    printlog(year, "files exist - skipping")
    next
  }
  
  printlog("Decompressing", basename(f))
  xfiles <- unzip(f, exdir = tmpdir)
  exdir <- dirname(xfiles[1])
  
  # Concatenate daily values into a single file
  printlog("Concatenating daily data...")
  system2("cdo", c("-s",   # silent
                   "select,name=sm", file.path(exdir, "*.nc"), yearfile),
          stderr = NULL)
  # Compute mean and variance
  printlog("Calculating", meanfile)
  system2("cdo", c("yearmean", yearfile, meanfile))
  printlog("Calculating", stdfile)
  system2("cdo", c("yearstd", yearfile, stdfile))
  
  list.files(tmpdir, pattern = "*.nc", recursive = TRUE, full.names = TRUE) %>% 
    file.remove()
  
}

# Now we want to concatenate the annual means and std values into a single file
printlog("Creating final mean and stdfiles...")
all_meanfiles <- list.files(CCI_DATA, "-mean.nc", full.names = TRUE)
final_meanfile <- file.path(CCI_DATA, "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-fv02.2.nc_means.nc")
system2("cdo", c("-s",   # silent
                 "select,name=sm", all_meanfiles, final_meanfile))
if(file.exists(final_meanfile)) {
  printlog("Final meanfile created OK; removing annual files")
  file.remove(all_meanfiles)
}

all_stdfiles <- list.files(CCI_DATA, "-std.nc", full.names = TRUE)
final_stdfile <- file.path(CCI_DATA, "ESACCI-SOILMOISTURE-L3S-SSMV-COMBINED-fv02.2.nc_stds.nc")
system2("cdo", c("-s",   # silent
                 "select,name=sm", all_stdfiles, final_stdfile))
if(file.exists(final_stdfile)) {
  printlog("Final stdfile created OK; removing annual files")
  file.remove(all_stdfiles)
}


printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
