# Prep the FLUXNET data
#
# Ben Bond-Lamberty January 2017

source("0-functions.R")

SCRIPTNAME  	<- "1-fluxnet.R"
PROBLEM       <- FALSE

FLUXNET_DATA <- "~/Data/FLUXNET2015/"

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# Extract the *_SUBSET_YY_* files (for GPP) and *_SUBSET_HH_* (for 
# nighttime Reco) from the FLUXNET zip files and save
# Downloaded 30 Jan 2017 from http://fluxnet.fluxdata.org (ftp.fluxdata.org/.fluxnet_downloads_86523/)
td <- tempdir()
td <- "~/Desktop/test/"
files <- list.files(FLUXNET_DATA, pattern = "zip$", full.names = TRUE)
stopifnot(length(files) > 0)
d <- list()
for(f in files) {
  printlog("Unzipping", basename(f))
  zf <- utils::unzip(f, list = TRUE)
  annual_file <- utils::unzip(f, files = zf$Name[grep("SUBSET_YY", zf$Name)], exdir = td)
  hourly_file <- utils::unzip(f, files = zf$Name[grep("SUBSET_HH", zf$Name)], exdir = td)
  
  # Read in the extracted annual file 
  stopifnot(length(annual_file) == 1)
  printlog("Reading", basename(annual_file))
  readr::read_csv(annual_file, na = "-9999") %>%
    select(TIMESTAMP, TA_F, P_F, NEE_VUT_REF_QC, GPP_DT_VUT_REF, GPP_NT_VUT_REF, RECO_DT_VUT_REF, RECO_NT_VUT_REF) %>%
    mutate(filename = annual_file) ->
    d_annual
  file.remove(annual_file)
  
  # Read in the extracted hourly files and compute annual nighttime Reco
  if(length(hourly_file) == 1) {
    printlog("Reading", basename(hourly_file))
    readr::read_csv(hourly_file, na = "-9999") %>%
      filter(NIGHT == 1) %>%
      select(TIMESTAMP_START, TIMESTAMP_END, NEE_VUT_REF) %>%
      mutate(YEAR_START = as.numeric(substr(TIMESTAMP_START, 1, 4))) %>%
      group_by(YEAR_START) %>%
      summarise(NEE_VUT_REF_NIGHT = mean(as.numeric(NEE_VUT_REF)) *
                  # convert from µmol/m2/s to gC/m2/yr. Here n() is number of half-hours
                  12 / 10^6 * 60 * 30 * n()) %>%
      left_join(d_annual, by = c("YEAR_START" = "TIMESTAMP")) %>%
      rename(TIMESTAMP = YEAR_START) ->
      d[[f]]
    file.remove(hourly_file)
  } else {
    d[[f]] <- d_annual
  }
}


# Combine with site data (in particular lon/lat information)
printlog(SEPARATOR)
sitedata <- read_csv("ancillary/fluxdata_sites.csv", col_types = "ccdddcdi")

printlog("Combining flux data and merging with site data...")
bind_rows(d) %>%
  separate(filename, into = c("FLX", "SITE_ID"), extra = "drop", sep = "_") %>%
  select(-FLX) %>%
  left_join(sitedata, by = "SITE_ID") ->
  fluxnet

save_data(fluxnet, scriptfolder = FALSE)

p <- qplot(RECO_DT_VUT_REF, NEE_VUT_REF_NIGHT, data = fluxnet, na.rm = TRUE) + geom_abline()
print(p)
save_plot("Reco_night")

print(lm(NEE_VUT_REF_NIGHT ~ RECO_DT_VUT_REF, data=fluxnet))

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
