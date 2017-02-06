# Heterotrophic respiration - main analysis script
#
# Ben Bond-Lamberty February 2017

source("0-functions.R")

SCRIPTNAME  	<- "4-analysis.R"
PROBLEM       <- FALSE

MAX_FLUXNET_DIST <- 1.0  # km
MIN_NEE_QC <- 0.5
SRDB_MINYEAR <- 1989
MAX_FLUX_TO_GPP <- 5   # Exclude ratios above this value; chosen based on distribution

library(broom)  # 0.4.1
library(Kendall) # 2.2
library(MASS) # 7.3.45


# Save a 2x2 grid plot of linear model diagnostics
save_model_diagnostics <- function(m, modelname = deparse(substitute(m))) {
  old.par <- par()
  pdf(file.path(outputdir(), paste0(modelname, ".pdf")))
  par(mfrow = c(2, 2))
  plot(m)
  dev.off()
  par(old.par)
}

# --------------------- Main --------------------------- 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

read_csv(SRDB_FILTERED_FILE) %>%
  print_dims() %>%
  filter(Study_midyear >= SRDB_MINYEAR) ->
  srdb

printlog("Filtered for studies after", SRDB_MINYEAR)


# -------------- 1. SRDB Rh:Rs analysis ------------------- 

printlog(SEPARATOR)
printlog("SRDB Rh:Rs analysis")

s1 <- subset(srdb, !is.na(Stage) & !is.na(Rh_annual) & !is.na(Rs_annual))
m1_rh_rs <- lm(Rh_annual/Rs_annual ~ Study_midyear * Stage + mat_hadcrut4 * map_hadcrut4, data = s1)
m1_rh_rs <- MASS::stepAIC(m1_rh_rs, direction = "both")
print(anova(m1_rh_rs))

m1_rh_rs_trend <- summary(m1_rh_rs)$coefficients["Study_midyear", "Pr(>|t|)"]
save_model_diagnostics(m1_rh_rs)

printlog("Mann-Kendall trend test:")
mk1_rh_rs <- MannKendall(s1$Rh_annual / s1$Rs_annual)
print(mk1_rh_rs)

srdb %>%
  filter(Year >= 1989, !is.na(Rs_annual), !is.na(Rh_annual)) %>%
  mutate(Year = as.integer(Year),
         yeargroup = cut(Year, breaks = c(1989, 1994, 1999, 2004, 2009, 2014),
                         labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014"))) %>%
  group_by(yeargroup) %>% 
  mutate(group_midyear = mean(Year),
         group = paste0(yeargroup, " (N = ", n(), ")")) -> 
  s3

p1_rh_rs <- ggplot(s3, aes(Rs_annual, Rh_annual, color = group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_color_grey("Year", start = 0.8, end = 0.2) +
  annotation_logticks() +
  xlab(expression(R[S]~(g~C~m^-2~yr^-1))) +
  ylab(expression(R[H]~(g~C~m^-2~yr^-1)))

p1_rh_rs <- p1_rh_rs + coord_cartesian(xlim=c(80, 3300), ylim=c(70, 2000))
printlog("NOTE we are plotting this graph with one point cut off:")
printlog(s3[which.min(s3$Rs_annual), c("Rs_annual", "Rh_annual")])

print(p1_rh_rs)
save_plot("1-srdb-rh-rs")


# ------------- 2. SRDB Rh:climate analysis --------------- 

printlog(SEPARATOR)
printlog("SRDB Rh climate analysis")

# Compute anomalies from the HadCRUT baseline
srdb$tmp_anom <- srdb$tmp_hadcrut4 - srdb$mat_hadcrut4
srdb$pre_anom <- srdb$pre_hadcrut4 - srdb$map_hadcrut4
srdb$pet_anom <- srdb$pet - srdb$pet_norm

m2_rh_climate <- lm(sqrt(Rh_annual) ~ mat_hadcrut4 + tmp_anom + map_hadcrut4 + pre_anom + pet_norm + pet_anom + Study_midyear * Stage, 
         data = srdb)
m2_rh_climate <- stepAIC(m2_rh_climate, direction = "both")
print(summary(m2_rh_climate))
save_model_diagnostics(m2_rh_climate)


# --------------- 3. FLUXNET analysis --------------------- 

printlog(SEPARATOR)
printlog("FLUXNET data analysis")
printlog("Filtering to MAX_FLUXNET_DIST =", MAX_FLUXNET_DIST)
printlog("Filtering to MIN_NEE_QC =", MIN_NEE_QC)

s3 %>%
  filter(FLUXNET_DIST <= MAX_FLUXNET_DIST) %>%
  filter(NEE_VUT_REF_QC >= MIN_NEE_QC) %>%
  print_dims %>%
  mutate(tmp_trend_label = if_else(tmp_trend > 0, "Warming", "Cooling"),
         pre_trend_label = if_else(pre_trend > 0, "Wetter", "Drier")) ->
  s3

# There are almost no FLUXNET sites with cooling trends, so we just 
# Looking at precip effect
s3 %>%
  group_by(pre_trend_label) %>%
  do(mod = lm(Rs_annual/gpp_fluxnet ~ Year, weights = YearsOfData, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Year") %>%
  print ->
  s3_models

printlog("Mann-Kendall trend test, drier-trend sites:")
s3_dry <- filter(s3, pre_trend_label == "Drier")
mk3_fluxnet_dry <- MannKendall(s3_dry$Rs_annual / s3_dry$gpp_fluxnet)
print(mk3_fluxnet_dry)
printlog("Mann-Kendall trend test, wetter-trend sites:")
s3_wet <- filter(s3, pre_trend_label == "Wetter")
mk3_fluxnet_wet <- MannKendall(s3_wet$Rs_annual / s3_wet$gpp_fluxnet)
print(mk3_fluxnet_wet)

p_fluxnet <- ggplot(s3, aes(Year, Rs_annual / gpp_fluxnet, group = FLUXNET_SITE_ID)) + 
  geom_point(aes(color = IGBP), na.rm = TRUE) + 
  facet_grid( ~ pre_trend_label, scales = "free") + 
  geom_smooth(method = "lm", color = "grey", fill = NA, na.rm = TRUE) + 
  geom_smooth(method = "lm", group = 1, na.rm = TRUE) +
  ylab(expression(R[S]:GPP[fluxnet]))

print(p_fluxnet)
save_plot("gpp_fluxnet")


# ----------- 4. Remotely sensed GPP analysis -------------- 

printlog(SEPARATOR)
printlog("Remote sensing analysis")

srdb %>%
  dplyr::select(Study_midyear, Biome, Leaf_habit, mat_hadcrut4, map_hadcrut4,
                gpp_beer, gpp_modis, Rs_annual, Rh_annual) %>%
  rename(`Beer et al.` = gpp_beer,
         MODIS = gpp_modis) %>%
  gather(Flux, fluxvalue, Rs_annual, Rh_annual) %>%
  gather(GPP, gppvalue, `Beer et al.`, MODIS) %>%
  filter(gppvalue > 0, !is.na(Leaf_habit)) -> 
  s_gpp

s_gpp %>%
  filter(fluxvalue / gppvalue >= MAX_FLUX_TO_GPP) ->
s_gpp_excluded
s_gpp %>%
  filter(fluxvalue / gppvalue < MAX_FLUX_TO_GPP) ->
  s_gpp_included

p_gpp_remotesensing <- ggplot(s_gpp_included, aes(Study_midyear, fluxvalue / gppvalue, color = Leaf_habit)) +
  geom_point() +
  geom_smooth(data = subset(s_gpp_included, Leaf_habit %in% c("Deciduous", "Evergreen")), method = "lm") +
  facet_grid(Flux ~ GPP, scales = "free") +
  scale_color_discrete("Leaf habit") +
  xlab("Year") +
  ylab("Respiration:GPP") +
  coord_cartesian(ylim = c(0, 2))
  
print(p_gpp_remotesensing)
save_plot("gpp_remotesensing")

printlog("Rs:MODIS GPP trend tests")
s_gpp_modis_rs <- subset(s_gpp, GPP == "MODIS" & Flux == "Rs_annual")
mk_gpp_modis_rs <- MannKendall(s_gpp_modis_rs$fluxvalue / s_gpp_modis_rs$gppvalue)
print(mk_gpp_modis_rs)

m_gpp_modis_rs <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + Study_midyear * Leaf_habit, 
                     data = s_gpp_modis_rs)
m_gpp_modis_rs <- stepAIC(m_gpp_modis_rs, direction = "both")
print(summary(m_gpp_modis_rs))
save_model_diagnostics(m_gpp_modis_rs)

printlog("Rh:MODIS GPP trend tests")
s_gpp_modis_rh <- subset(s_gpp, GPP == "MODIS" & Flux == "Rh_annual")
mk_gpp_modis_rh <- MannKendall(s_gpp_modis_rh$fluxvalue / s_gpp_modis_rh$gppvalue)
print(mk_gpp_modis_rh)
printlog("Rh:MODIS GPP trend tests (dbf)")
s_gpp_modis_rh_dbf <- subset(s_gpp_modis_rh, Leaf_habit == "Deciduous")
mk_gpp_modis_rh_dbf <- MannKendall(s_gpp_modis_rh_dbf$fluxvalue / s_gpp_modis_rh_dbf$gppvalue)
print(mk_gpp_modis_rh_dbf)
printlog("Rh:MODIS GPP trend tests (enf)")
s_gpp_modis_rh_enf <- subset(s_gpp_modis_rh, Leaf_habit == "Evergreen")
mk_gpp_modis_rh_enf <- MannKendall(s_gpp_modis_rh_enf$fluxvalue / s_gpp_modis_rh_enf$gppvalue)
print(mk_gpp_modis_rh_enf)

m_gpp_modis_rh <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + Study_midyear * Leaf_habit, 
         data = s_gpp_modis_rh)
m_gpp_modis_rh <- stepAIC(m_gpp_modis_rh, direction = "both")
print(summary(m_gpp_modis_rh))
save_model_diagnostics(m_gpp_modis_rh)


printlog("Rs:Beer GPP trend tests")
s_gpp_beer_rs <- subset(s_gpp, GPP == "Beer et al." & Flux == "Rs_annual")
mk_gpp_beer_rs <- MannKendall(s_gpp_beer_rs$fluxvalue / s_gpp_beer_rs$gppvalue)
print(mk_gpp_beer_rs)

m_gpp_beer_rs <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + Study_midyear * Leaf_habit, 
                     data = s_gpp_beer_rs)
m_gpp_beer_rs <- stepAIC(m_gpp_beer_rs, direction = "both")
print(summary(m_gpp_beer_rs))
save_model_diagnostics(m_gpp_beer_rs)


printlog("Rh:Beer GPP trend tests")
s_gpp_beer_rh <- subset(s_gpp, GPP == "Beer et al." & Flux == "Rh_annual")
mk_gpp_beer_rh <- MannKendall(s_gpp_beer_rh$fluxvalue / s_gpp_beer_rh$gppvalue)
print(mk_gpp_beer_rh)
printlog("Rh:Beer GPP trend tests (dbf)")
s_gpp_beer_rh_dbf <- subset(s_gpp_beer_rh, Leaf_habit == "Deciduous")
mk_gpp_beer_rh_dbf <- MannKendall(s_gpp_beer_rh_dbf$fluxvalue / s_gpp_beer_rh_dbf$gppvalue)
print(mk_gpp_beer_rh_dbf)
printlog("Rh:Beer GPP trend tests (enf)")
s_gpp_beer_rh_enf <- subset(s_gpp_beer_rh, Leaf_habit == "Evergreen")
mk_gpp_beer_rh_enf <- MannKendall(s_gpp_beer_rh_enf$fluxvalue / s_gpp_beer_rh_enf$gppvalue)
print(mk_gpp_beer_rh_enf)

m_gpp_beer_rh <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + Study_midyear * Leaf_habit, 
                     data = s_gpp_beer_rh)
m_gpp_beer_rh <- stepAIC(m_gpp_beer_rh, direction = "both")
print(summary(m_gpp_beer_rh))
save_model_diagnostics(m_gpp_beer_rh)


# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
