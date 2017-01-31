# Heterotrophic respiration - main analysis script
#
# Ben Bond-Lamberty January 2016

source("0-functions.R")

SCRIPTNAME  	<- "4-analysis.R"
PROBLEM       <- FALSE

MAX_FLUXNET_DIST <- 1.0  # km
MIN_NEE_QC <- 0.5

library(broom)  # 0.4.1

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

srdb <- read_csv("outputs/srdb_filtered.csv") 

print_dims(srdb)


# ----------------------------- FLUXNET data analysis

printlog("Filtering to MAX_FLUXNET_DIST =", MAX_FLUXNET_DIST)
x <- subset(x, FLUXNET_DIST <= MAX_FLUXNET_DIST)
print_dims(x)
printlog("Filtering to MIN_NEE_QC =", MIN_NEE_QC)
x <- subset(x, NEE_VUT_REF_QC >= MIN_NEE_QC)
print_dims(x)

x %>%
  mutate(tmp_trend_label = if_else(tmp_trend > 0, "Warming", "Cooling"),
         pre_trend_label = if_else(pre_trend > 0, "Wetter", "Drier")) ->
  x

x %>%
  group_by(tmp_trend_label) %>%
  do(mod = lm(Rs_annual/GPP_NT_VUT_REF ~ Year, data = .)) %>%
  tidy(mod) %>%
  print
  

p <- ggplot(x, aes(Year, Rs_annual / GPP_NT_VUT_REF, group = FLUXNET_SITE_ID)) + 
  geom_point(aes(color = FLUXNET_SITE_ID)) + 
  facet_grid(tmp_trend_label ~ pre_trend_label, scales = "free") + 
  geom_smooth(method = "lm", color = 'grey', fill = NA) + 
  geom_smooth(method = "lm", group = 1)

print(p)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
