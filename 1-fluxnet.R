# Prep the FLUXNET data
#
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "X-fluxnet.R"
PROBLEM       <- FALSE

FLUXNET_DATA <- "~/Data/FLUXNET/"

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# Extract the *_SUBSET_YY_* files from the FLUXNET zip files and save
# Downloaded 30 Jan 2017 from http://fluxnet.fluxdata.org (ftp.fluxdata.org/.fluxnet_downloads_86523/)
td <- tempdir()
td <- "~/Desktop/test/"
files <- list.files(FLUXNET_DATA, pattern = "zip$", full.names = TRUE)
for(f in files) {
  printlog("Unzipping", basename(f))
  zf <- utils::unzip(f, list = TRUE)
  utils::unzip(f, files = zf$Name[grep("SUBSET_YY", zf$Name)], exdir = td)
}

# Read in the extracted files into a list
printlog(SEPARATOR)
d <- list()
xfiles <- list.files(td)
for(f in xfiles) {
  printlog("Reading", f)
  readr::read_csv(file.path(td, f), na = "-9999") %>%
    select(TIMESTAMP, TA_F, P_F, NEE_VUT_REF, NEE_VUT_REF_QC, RECO_NT_VUT_REF, GPP_NT_VUT_REF) %>%
    mutate(filename = f) ->
    d[[f]]
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

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
