# Heterotrophic respiration - main analysis script
#
# Ben Bond-Lamberty February 2017

source("0-functions.R")

SCRIPTNAME  	<- "4-analysis.R"
PROBLEM       <- FALSE

MAX_FLUXNET_DIST <- 5  # km
MIN_NEE_QC <- 0.5
SRDB_MINYEAR <- 1989
MAX_FLUX_TO_GPP <- 5   # Exclude ratios above this value; chosen based on distribution

library(broom)  # 0.4.1
library(Kendall) # 2.2
library(MASS) # 7.3.45


# Save a 2x2 grid plot of linear model diagnostics
save_model_diagnostics <- function(m, modelname = deparse(substitute(m))) {
  pdf(file.path(outputdir(), paste0(modelname, ".pdf")))
  old.par <- par(mfrow = c(2, 2))
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

srdb %>%
  filter(!is.na(Stage), !is.na(Leaf_habit), !is.na(Rs_annual)) ->
  s_rh_rs
m1_rh_rs <- lm(Rh_annual/Rs_annual ~ Study_midyear * Stage + 
                 Study_midyear * Partition_method +
                 Study_midyear * Leaf_habit +
                 mat_hadcrut4 * map_hadcrut4, data = s_rh_rs)
m1_rh_rs <- MASS::stepAIC(m1_rh_rs, direction = "both")
print(anova(m1_rh_rs))

m1_rh_rs_signif <- anova(m1_rh_rs)["Study_midyear", "Pr(>F)"]
save_model_diagnostics(m1_rh_rs)

printlog("Mann-Kendall trend test:")
mk1_rh_rs <- MannKendall(s_rh_rs$Rh_annual / s_rh_rs$Rs_annual)
print(mk1_rh_rs)

srdb %>%
  filter(Year >= 1989, !is.na(Rs_annual), !is.na(Rh_annual)) %>%
  mutate(Year = as.integer(Year),
         yeargroup = cut(Year, breaks = c(1989, 1994, 1999, 2004, 2009, 2014),
                         labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014"))) %>%
  group_by(yeargroup) %>% 
  mutate(group_midyear = mean(Year),
         group = paste0(yeargroup, " (N = ", n(), ")")) -> 
  s_rh_rs

# Compute summary statistics
s_rh_rs %>%
  group_by(yeargroup) %>%
  summarise(rh_rs_mean = pn(mean(Rh_annual / Rs_annual), 2),
            rh_rs_sd = pn(sd(Rh_annual / Rs_annual), 2), 
            n = n()) ->
  rh_rs_summary

# Make Figure 1
p1_rh_rs <- ggplot(s_rh_rs, aes(Rs_annual, Rh_annual, color = group)) +
  scale_x_log10() + scale_y_log10() +
  scale_color_grey("Year", start = 0.8, end = 0.2) +
  annotation_logticks() +
  xlab(expression(R[S]~(g~C~m^-2~yr^-1))) +
  ylab(expression(R[H]~(g~C~m^-2~yr^-1))) + 
  coord_cartesian(xlim=c(80, 3300), ylim=c(70, 2000))

p_inset <- ggplot(s_rh_rs, aes(Rh_annual / Rs_annual, color = yeargroup, fill = yeargroup)) + 
  geom_density(alpha = 0.5) + 
  xlab(expression(R[H]:R[S])) + ylab("") +
  scale_fill_grey(start = 0.8, end = 0.2, guide = FALSE) +
  scale_color_grey(start = 0.8, end = 0.2, guide = FALSE) +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 6), axis.title.x = element_text(size = 8))

p1_rh_rs <- p1_rh_rs + 
  annotation_custom(grob = ggplotGrob(p_inset), xmin = log10(60), xmax = log10(800), ymin = log10(600), ymax = log10(2300)) 

p1_rh_rs <- p1_rh_rs + geom_point() + geom_smooth(method = "lm", se = FALSE)

printlog("NOTE we are plotting this graph with one point cut off:")
printlog(s_rh_rs[which.min(s_rh_rs$Rs_annual), c("Rs_annual", "Rh_annual")])

print(p1_rh_rs)
save_plot("1-srdb-rh-rs", ptype = ".png")


# ------------- 2. SRDB Rh:climate analysis --------------- 

printlog(SEPARATOR)
printlog("SRDB Rh climate analysis")

# Compute anomalies from the HadCRUT baseline
srdb$tmp_anom <- srdb$tmp_hadcrut4 - srdb$mat_hadcrut4
srdb$pre_anom <- srdb$pre_hadcrut4 - srdb$map_hadcrut4
srdb$pet_anom <- srdb$pet - srdb$pet_norm

m2_rh_climate <- lm(sqrt(Rh_annual) ~ mat_hadcrut4 + tmp_anom + map_hadcrut4 + pre_anom + pet_norm + pet_anom + Stage * Leaf_habit, 
                    data = srdb)
m2_rh_climate <- stepAIC(m2_rh_climate, direction = "both")
print(anova(m2_rh_climate))
m2_rh_climate_map_signif <- anova(m2_rh_climate)["map_hadcrut4", "Pr(>F)"]
m2_rh_climate_pet_signif <- anova(m2_rh_climate)["pet_norm", "Pr(>F)"]
m2_rh_climate_stage_signif <- anova(m2_rh_climate)["Stage", "Pr(>F)"]
m2_rh_climate_tanom_signif <- anova(m2_rh_climate)["tmp_anom", "Pr(>F)"]
save_model_diagnostics(m2_rh_climate)


# --------------- 3. FLUXNET analysis --------------------- 

printlog(SEPARATOR)
printlog("FLUXNET data analysis")
printlog("Filtering to MAX_FLUXNET_DIST =", MAX_FLUXNET_DIST)
printlog("Filtering to MIN_NEE_QC =", MIN_NEE_QC)

# Side analysis: check how choice of MAX_FLUXNET_DIST affects our N
test <- seq(0.1, 10, by = 0.1)
result <- rep(NA, length(test))
for(i in seq_along(test)) {
  filter(srdb, !is.na(Rs_annual), 
         FLUXNET_DIST <= test[i], 
         NEE_VUT_REF_QC >= MIN_NEE_QC) -> x
  result[i] <- nrow(x)
}
df <- tibble(x = test, y = result)
p <- qplot(x, y, data = df, geom="line") + xlab("MAX_FLUXNET_DIST") + ylab("N")
print(p)
save_plot("fluxnet_dist_n")

# OK now for real, here we go: filter SRDB for available data,
# distance to FLUXNET tower, and NEE quality
srdb %>%
  filter(!is.na(Rs_annual),
         FLUXNET_DIST <= MAX_FLUXNET_DIST,
         NEE_VUT_REF_QC >= MIN_NEE_QC) %>%
  # We only allow one observation per site per year
  # Otherwise Harvard Forest swamps everything!
  group_by(FLUXNET_SITE_ID, IGBP, tmp_trend, pre_trend, Stage, Year) %>%
  summarise(Rs_annual = mean(Rs_annual),
            NEE_VUT_REF = mean(NEE_VUT_REF),
            RECO_NT_VUT_REF = mean(RECO_NT_VUT_REF),
            gpp_fluxnet = mean(gpp_fluxnet),
            YearsOfData = mean(YearsOfData),
            mat_hadcrut4 = mean(mat_hadcrut4),
            map_hadcrut4 = mean(map_hadcrut4)) %>%
  mutate(tmp_trend_label = if_else(tmp_trend > 0, "Warming", "Cooling"),
         pre_trend_label = if_else(pre_trend > 0, "Wetter", "Drier")) ->
  s_fluxnet

s_fluxnet_dry <- subset(s_fluxnet, pre_trend_label == "Drier")
m_fluxnet_dry <- lm(Rs_annual/gpp_fluxnet ~ Year * Stage + mat_hadcrut4 * map_hadcrut4, 
                    data = s_fluxnet_dry, weights = YearsOfData)
m_fluxnet_dry <- MASS::stepAIC(m_fluxnet_dry, direction = "both")
print(summary(m_fluxnet_dry))

s_fluxnet_wet <- subset(s_fluxnet, pre_trend_label == "Wetter")
m_fluxnet_wet <- lm(Rs_annual/gpp_fluxnet ~ Year * Stage + mat_hadcrut4 * map_hadcrut4, 
                    data = s_fluxnet_wet,
                    weights = YearsOfData)
m_fluxnet_wet <- MASS::stepAIC(m_fluxnet_wet, direction = "both")
print(summary(m_fluxnet_wet))


# Experimented with fitting a mixed-effects model (with site as random effect)
# Doesn't seem to add/change much
# s_fluxnet_wet <- subset(s_fluxnet, pre_trend_label == "Wetter")
# m_fluxnet_wet <- lme(Rs_annual/gpp_fluxnet ~ Year * Stage + mat_hadcrut4 * map_hadcrut4, 
#                     data = s_fluxnet_wet,
#                     random = ~ 1 | FLUXNET_SITE_ID,
#                     weights = YearsOfData)
# m_fluxnet_wet <- MASS::stepAIC(m_fluxnet_wet, direction = "both")
# print(summary(m_fluxnet_wet))


s_fluxnet_wet_nohf <- subset(s_fluxnet_wet, FLUXNET_SITE_ID != "US-Ha1")
m_fluxnet_wet_nohf <- lm(Rs_annual/gpp_fluxnet ~ Year * Stage + mat_hadcrut4 * map_hadcrut4, 
                         data = s_fluxnet_wet_nohf, weights = YearsOfData)
m_fluxnet_wet_nohf <- MASS::stepAIC(m_fluxnet_wet_nohf, direction = "both")
print(summary(m_fluxnet_wet_nohf))


printlog("Mann-Kendall trend test, drier-trend sites:")
mk3_fluxnet_dry <- MannKendall(s_fluxnet_dry$Rs_annual / s_fluxnet_dry$gpp_fluxnet)
print(mk3_fluxnet_dry)
printlog("Mann-Kendall trend test, wetter-trend sites:")
mk3_fluxnet_wet <- MannKendall(s_fluxnet_wet$Rs_annual / s_fluxnet_wet$gpp_fluxnet)
print(mk3_fluxnet_wet)
printlog("Mann-Kendall trend test, wetter-trend sites w/o HF:")
mk3_fluxnet_wet_nohf <- MannKendall(s_fluxnet_wet_nohf$Rs_annual / s_fluxnet_wet_nohf$gpp_fluxnet)
print(mk3_fluxnet_wet_nohf)

s_fluxnet %>%
  group_by(FLUXNET_SITE_ID) %>%
  summarise(Year = min(Year), 
            ratio = (Rs_annual / gpp_fluxnet)[which.min(Year)],
            IGBP = unique(IGBP),
            pre_trend_label = unique(pre_trend_label)) ->
  s_fluxnet_labels

p_fluxnet <- ggplot(s_fluxnet, aes(Year, Rs_annual / gpp_fluxnet, group = FLUXNET_SITE_ID)) + 
  geom_point(aes(color = IGBP), na.rm = TRUE) + 
  facet_grid( ~ pre_trend_label, scales = "free") + 
  geom_smooth(method = "lm", color = "grey", fill = NA, na.rm = TRUE) + 
  geom_smooth(method = "lm", group = 1, na.rm = TRUE) +
  ylab(expression(R[S]:GPP[fluxnet]))

p_fluxnet <- p_fluxnet + geom_text(data = s_fluxnet_labels, 
                                   aes(y = ratio, label = FLUXNET_SITE_ID, color = IGBP), 
                                   size = 2, nudge_y = 0.05,
                                   alpha = 0.75)

print(p_fluxnet)

save_plot("gpp_fluxnet")


# ----------- 4. Remotely sensed GPP analysis -------------- 

printlog(SEPARATOR)
printlog("Remote sensing analysis")

srdb %>%
  dplyr::select(Study_midyear, Biome, Leaf_habit, Partition_method,
                mat_hadcrut4, map_hadcrut4,
                gpp_beer, gpp_modis, Rs_annual, Rh_annual) %>%
  rename(`Beer` = gpp_beer,
         MODIS = gpp_modis) %>%
  gather(Flux, fluxvalue, Rs_annual, Rh_annual) %>%
  gather(GPP, gppvalue, Beer, MODIS) %>%
  filter(gppvalue > 0, !is.na(Leaf_habit), !is.na(fluxvalue)) -> 
  s_gpp

# Make pretty facet labels
s_gpp$Prettyflux <- "R[H]"
s_gpp$Prettyflux[s_gpp$Flux == "Rs_annual"] <- "R[S]"

s_gpp %>%
  filter(fluxvalue / gppvalue >= MAX_FLUX_TO_GPP) ->
  s_gpp_excluded
s_gpp %>%
  filter(fluxvalue / gppvalue < MAX_FLUX_TO_GPP) ->
  s_gpp_included

p_gpp_remotesensing <- ggplot(s_gpp_included, aes(Study_midyear, fluxvalue / gppvalue, color = Leaf_habit)) +
  geom_point() +
  geom_smooth(data = subset(s_gpp_included, Leaf_habit %in% c("Deciduous", "Evergreen")), method = "lm", show.legend = FALSE) +
  facet_grid(Prettyflux ~ GPP, scales = "free", labeller = label_parsed) +
  scale_color_discrete("Leaf habit") +
  xlab("Year") +
  ylab("Respiration:GPP") +
  coord_cartesian(ylim = c(0, 2))

print(p_gpp_remotesensing)
save_plot("gpp_remotesensing", ptype = ".png")

printlog("Rs:MODIS GPP trend tests")
s_gpp_modis_rs <- subset(s_gpp, GPP == "MODIS" & Flux == "Rs_annual")
mk_gpp_modis_rs <- MannKendall(s_gpp_modis_rs$fluxvalue / s_gpp_modis_rs$gppvalue)
print(mk_gpp_modis_rs)

m_gpp_modis_rs <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + 
                       Study_midyear * Leaf_habit + Study_midyear * Partition_method, 
                     data = s_gpp_modis_rs)
m_gpp_modis_rs <- stepAIC(m_gpp_modis_rs, direction = "both")
print(summary(m_gpp_modis_rs))
m_gpp_modis_rs_signif <- anova(m_gpp_modis_rs)["Study_midyear", "Pr(>F)"]
m_gpp_modis_rs_leaf_signif <- anova(m_gpp_modis_rs)["Leaf_habit", "Pr(>F)"]
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

m_gpp_modis_rh <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + 
                       Study_midyear * Leaf_habit + Study_midyear * Partition_method,
                     data = s_gpp_modis_rh)
m_gpp_modis_rh <- stepAIC(m_gpp_modis_rh, direction = "both")
print(summary(m_gpp_modis_rh))
m_gpp_modis_rh_signif <- anova(m_gpp_modis_rh)["Study_midyear", "Pr(>F)"]
m_gpp_modis_rh_leaf_signif <- anova(m_gpp_modis_rh)["Leaf_habit", "Pr(>F)"]
save_model_diagnostics(m_gpp_modis_rh)


printlog("Rs:Beer GPP trend tests")
s_gpp_beer_rs <- subset(s_gpp, GPP == "Beer" & Flux == "Rs_annual")
mk_gpp_beer_rs <- MannKendall(s_gpp_beer_rs$fluxvalue / s_gpp_beer_rs$gppvalue)
print(mk_gpp_beer_rs)

m_gpp_beer_rs <- lm(fluxvalue / gppvalue ~ mat_hadcrut4 + map_hadcrut4 + Study_midyear * Leaf_habit, 
                    data = s_gpp_beer_rs)
m_gpp_beer_rs <- stepAIC(m_gpp_beer_rs, direction = "both")
print(summary(m_gpp_beer_rs))
m_gpp_beer_rs_signif <- anova(m_gpp_beer_rs)["Study_midyear", "Pr(>F)"]
m_gpp_beer_rs_leaf_signif <- anova(m_gpp_beer_rs)["Leaf_habit", "Pr(>F)"]
save_model_diagnostics(m_gpp_beer_rs)


printlog("Rh:Beer GPP trend tests")
s_gpp_beer_rh <- subset(s_gpp, GPP == "Beer" & Flux == "Rh_annual")
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
m_gpp_beer_rh_signif <- anova(m_gpp_beer_rh)["Study_midyear", "Pr(>F)"]
m_gpp_beer_rh_leaf_signif <- anova(m_gpp_beer_rh)["Leaf_habit", "Pr(>F)"]
save_model_diagnostics(m_gpp_beer_rh)


# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
