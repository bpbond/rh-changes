# Test the variability of ISIMIP RH trends (per Referee 1)
#
# Ben Bond-Lamberty November 2017

source("0-functions.R")

SCRIPTNAME  	<- "X-ref1_isimip.R"
PROBLEM       <- FALSE

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

# Start by plotting global RH versus the Hashimoto data
read_csv("ancillary/isimip-rh/all_annual_global_rh.csv") %>% 
  gather(model, rh, -year) ->
  isimip_global_rh

library(ncdf4)
# Downloaded August 25, 2017 from http://cse.ffpri.affrc.go.jp/shojih/data/index.html
# Processed using `cdo fldmean` to produce annual values
# These annual data start in 1901; extract 1971-2010
printlog("Reading Hashimoto data (again) and plotting versus ISIMIP...")
nc <- nc_open("~/Data/Hashimoto/RH_yr_Hashimoto2015_global.nc")
hashimoto <- ncvar_get(nc, "co2", start = c(1, 1, 1, 70), count = c(-1, -1, 1, 40))
global_area <- 1.293606e+14 # total land area (m2) in Hashimoto file
hashimoto <- tibble(year = 1971:2010, rh = hashimoto * global_area)          

p <- ggplot(isimip_global_rh, aes(year, rh / 1e15, color = model)) + 
  geom_point(alpha = 0.5) + geom_smooth(method = "lm")
p <- p + geom_point(data = hashimoto, color = "black", alpha = 0.5) +
  geom_smooth(data = hashimoto, method = "lm", color = "black", size = 3)
p <- p + ggtitle("ISIMIP versus observations (in black)") + ylab("Global RH (Pg C)")
print(p)
save_plot("isimip_vs_hashimoto")

# Now look at all the site-specific RH values from the ISIMIP models
# How variable are they?
read_csv("ancillary/isimip-rh/all_annual_site_rh.csv") %>% 
  gather(year, rh, matches("[0-9]{4}")) %>% 
  mutate(year = as.integer(year)) ->
  isimip_site_rh

library(broom)
printlog("Computing slopes of modeled RH for each model and site...")
isimip_site_rh %>% 
  filter(!is.na(year), !is.na(rh)) %>% 
  group_by(model, site) %>% 
  do(sitemod = lm(rh ~ year, data = .)) %>%
  tidy(sitemod) %>%
  filter(term == "year", !is.na(p.value)) %>% 
  select(model, site, estimate, p.value) ->
  model_site_slopes

print_dims(model_site_slopes)
model_site_slopes %>% 
  mutate(slope_sign = sign(estimate), significant = p.value < 0.05) %>% 
  group_by(slope_sign, significant) %>% 
  summarise(n = n(), percent = round(n / nrow(model_site_slopes) * 100, 0)) %>% 
  printlog

p <- ggplot(model_site_slopes, aes(model, estimate)) + 
  geom_boxplot() +
  geom_jitter(aes(color = p.value < 0.05), alpha = 0.25) + 
  ylab("RH slope over time (gC/m2/yr)")
save_plot("isimip_site_slopes")

save_data(model_site_slopes, scriptfolder = FALSE)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
