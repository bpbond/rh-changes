# Find longest time series in the SRDB database
#
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "X-ameriflux.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

readxl::read_excel("/Users/d3x290/Documents/Work/Data-ongoing/Soil respiration database/srdb-data.xlsx") %>%
  mutate(lat = round(Latitude, 1), 
         lon = round(Longitude, 1),
         Site_name = substr(Site_name, 1, 10)) %>%
  filter(!is.na(lat), Manipulation == "None") %>%
  filter(is.na(Duplicate_record)) %>%
  filter(!is.na(Rs_annual) | !is.na(Rh_annual)) %>%
  group_by(lon, lat) %>% 
  summarise(years = n_distinct(Study_midyear),
            FLUXNET = paste(unique(FLUXNET_SITE_ID), collapse = " "),
            minyear = min(Study_midyear), 
            maxyear = max(Study_midyear),
            sites = paste(unique(Site_name), collapse = " ")) %>% 
  arrange(desc(years)) %>%
  filter(years >= 6) -> 
  x

print(x)

# Load in Fluxnet data and quickly find nearest station?


printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
