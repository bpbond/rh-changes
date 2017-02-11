# Prep data for heterotrophic respiration analysis
#
# Load SRDB; filter for 'good' data (unmanipulated ecosystems, IRGA/GC only, etc);
# spatially match with CRU climate, Max Planck GPP, and MODIS GPP datasets
#
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "X-ancillary.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# 1. ORNL DAAC litter data (Post)

# Downloaded January 11, 2017 from https://daac.ornl.gov/VEGETATION/guides/Global_Litter_Carbon_Nutrients.html
readr::read_csv("/Users/d3x290/Data/ORNL/GLOBAL_LITTER_CARBON_NUTRIENTS_1244/data/Global_Litterfall.csv", 
         skip = 9, na = "-99",
         col_names = c("SITE_ID", "PLOTNAME", "LATITUDE", "LONGITUDE", "ELEVATION", "MEANTEMP", "BIOTEMP", "PRECIP", "OLSON", "LIFEZONE", "MATTHEWS", "VEGETATION", "MANAGED", "AGE", "REFERENCE1", "REFERENCE2", "LFLLW", "SWFLLW", "FINFLLW", "RPRFLLW", "LTRFLLW", "LWFLLW", "LLEAFW", "LSWW", "LFINW", "LRPRW", "LLTRW", "FLEAFW", "FSWW", "FFINW", "FRPRW", "FLTRW", "TLEAFW", "TSWW", "TFINW", "TRPRW", "TLTRW", "LWW", "LFLLN", "SWFLLN", "FINFLLN", "RPRFLLN", "LTRFLLN", "LWFLLN", "LLEAFN", "LSWN", "LFINN", "LRPRN", "LLTRN", "FLEAFN", "FSWN", "FFINN", "FRPRN", "FLTRN", "TLEAFN", "TSWN", "TFINN", "TRPRN", "TLTRN", "LWN", "LFLLP", "SWFLLP", "FINFLLP", "RPRFLLP", "LTRFLLP", "LWFLLP", "LLEAFP", "LSWP", "LFINP", "LRPRP", "LLTRP", "FLEAFP", "FSWP", "FFINP", "FRPRP", "FLTRP", "TLEAFP", "TSWP", "TFINP", "TRPRP", "TLTRP", "LWP", "LFLLCA", "SWFLLCA", "FINFLLCA", "RPRFLLCA", "LTRFLLCA", "LWFLLCA", "LLEAFCA", "LSWCA", "LFINCA", "LRPRCA", "LLTRCA", "FLEAFCA", "FSWCA", "FFINCA", "FRPRCA", "FLTRCA", "TLEAFCA", "TSWCA", "TFINCA", "TRPRCA", "TLTRCA", "LWCA", "LFLLMG", "SWFLLMG", "FINFLLMG", "RPRFLLMG", "LTRFLLMG", "LWFLLMG", "LLEAFMG", "LSWMG", "LFINMG", "LRPRMG", "LLTRMG", "FLEAFMG", "FSWMG", "FFINMG", "FRPRMG", "FLTRMG", "TLEAFMG", "TSWMG", "TFINMG", "TRPRMG", "TLTRMG", "LWMG", "LFLLK", "SWFLLK", "FINFLLK", "RPRFLLK", "LTRFLLK", "LWFLLK", "LLEAFK", "LSWK", "LFINK", "LRPRK", "LLTRK", "FLEAFK", "FSWK", "FFINK", "FRPRK", "FLTRK", "TLEAFK", "TSWK", "TFINK", "TRPRK", "TLTRK", "LWK", "LFLLNA", "SWFLLNA", "FINFLLNA", "RPRFLLNA", "LTRFLLNA", "LWFLLNA", "LLEAFNA", "LSWNA", "LFINNA", "LRPRNA", "LLTRNA", "FLEAFNA", "FSWNA", "FFINNA", "FRPRNA", "FLTRNA", "TLEAFNA", "TSWNA", "TFINNA", "TRPRNA", "TLTRNA", "LWNA", "LFLLMN", "SWFLLMN", "FINFLLMN", "RPRFLLMN", "LTRFLLMN", "LWFLLMN", "LLEAFMN", "LSWMN", "LFINMN", "LRPRMN", "LLTRMN", "FLEAFMN", "FSWMN", "FFINMN", "FRPRMN", "FLTRMN", "TLEAFMN", "TSWMN", "TFINMN", "TRPRMN", "TLTRMN", "LWMN")) %>%
  dplyr::select(LATITUDE, LONGITUDE, MEANTEMP, PRECIP,
         LTRFLLW, LWFLLW,
         TLTRW, LWW,
         REFERENCE1) ->
  post

# Extract a (publication) year
pattern <- "19[0-9]{2}"
post1 <- post[grep(pattern, post$REFERENCE1),]
re <- regexec(pattern, post1$REFERENCE1)
post1$year <- as.numeric(unlist(regmatches(post1$REFERENCE1, re)))

post1 %>%
  # Compute (steady state) k as flux divided by pool
  mutate(k = LTRFLLW / TLTRW) %>%
  # Filter out a few real oddballs
  filter(k < 10,  year > 1970) ->  
  post1

p <- qplot(year, k, data = post1, color = MEANTEMP, size = TLTRW)
print(p)

m <- lm(k ~ MEANTEMP * PRECIP + year, data = post1)
print(summary(m))

m1 <- MASS::stepAIC(m, direction = "both")
print(summary(m1))

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
