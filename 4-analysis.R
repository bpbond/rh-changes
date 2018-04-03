# Heterotrophic respiration - main analysis script
# Ben Bond-Lamberty February 2017 (last update April 2018)
#
# A lot goes on here, but there are four main sections below:
# 1. The Rh:Rs analysis: is the ratio changing over time?
# 2. Is Rh responding to climate anomalies? (Minor; parallels 2010 paper)
# 3. Is Rs:FLUXNET GPP rising over time?
# 4. For FLUXNET only, is nighttime NEE:GPP rising?
# 5. Remote sensing: are Rh:GPP, Rs:GPP, Rh:SIF, and Rs:SIF changing over time?

# Then there are three more sections added to address Referee 1 concerns:
# 6. Does RH:ISIMIP GPP rise?
# 7. How sensitive are the results to choice of start/end date?
# 8. What about site-specific (longitudinal) trends?

# Note all the figures can be searched as "Fig1", "Fig2", "EDF1", "EDF2", etc.
# (Some of the EDF figures are generated in other scripts.)

source("0-functions.R")

SCRIPTNAME  	<- "4-analysis.R"
PROBLEM       <- FALSE

MAX_FLUXNET_DIST <- 5  # km
MIN_NEE_QC <- 0.5
SRDB_MINYEAR <- 1989
MAX_FLUX_TO_GPP <- 5   # Exclude ratios above this value; chosen based on distribution

library(broom)  # 0.4.4
library(Kendall) # 2.2
library(MASS) # 7.3-47
library(mblm) # 0.12
library(scales) # 0.5.0

# Save a 2x2 grid plot of linear model diagnostics
save_model_diagnostics <- function(m, modelname = deparse(substitute(m))) {
  pdf(file.path(outputdir(), paste0(modelname, ".pdf")))
  old.par <- par(mfrow = c(2, 2))
  plot(m)
  dev.off()
  par(old.par)
}

# Add a result to summary table
overall_results <- list()
add_result <- function(testname, x, y, n, p, p_landcover, p_disturbance) {
  printlog("Fitting Theil-Sen for", testname, "and adding to results...")
  m_ts <- mblm::mblm(y ~ x) %>% summary  # Theil-Sen
  overall_results[[testname]] <<- tibble(test = testname,
                                         n = n,
                                         p = p,
                                         p_theil.sen = m_ts$coefficients[2, 4],
                                         p_land.cover = p_landcover,
                                         p_disturbance = p_disturbance)
}

# --------------------- Main --------------------------- 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

read_csv(SRDB_FILTERED_FILE, guess_max = 1e6) %>%
  print_dims() %>%
  
  # Per Referee 2, we look at land cover, splitting trees (by far the largest Ecosystem_type) into 
  # deciduous and coniferous categories
  mutate(Land_cover = if_else(Ecosystem_type == "Forest", paste(Leaf_habit, Ecosystem_type), Ecosystem_type)) ->
  srdb_complete

# Group most land cover categories (with less than 100 observations) into "Other"
srdb_complete %>% group_by(Land_cover) %>% summarise(n = n()) %>% filter(n >= 100) -> over100
srdb_complete$Land_cover[! srdb_complete$Land_cover %in% over100$Land_cover] <- "Other"
lvls <- unique(srdb_complete$Land_cover)
lvls <- c(lvls[-which(lvls == "Other")], "Other")
srdb_complete$Land_cover <- factor(srdb_complete$Land_cover, levels = lvls)
srdb_complete %>% group_by(Land_cover) %>% summarise(n = n()) 

srdb_complete %>%
  filter(Study_midyear >= SRDB_MINYEAR,
         Ecosystem_state != "Managed",
         Ecosystem_type != "Agriculture") ->
  srdb

printlog("Filtered for studies after", SRDB_MINYEAR)

printlog("MAT range is", paste(round(range(srdb$mat_hadcrut4, na.rm = TRUE), 1), collapse = ", "))
printlog("MAP range is", paste(round(range(srdb$map_hadcrut4, na.rm = TRUE), 0), collapse = ", "))


# -------------- 1. SRDB Rh:Rs analysis ------------------- 

printlog(SEPARATOR)
printlog("SRDB Rh:Rs analysis")

srdb %>% group_by(Biome, Ecosystem_type) %>% 
  summarise(SOC_biome = median(SOC, na.rm = TRUE)) -> 
  biome_medians

srdb %>%
  dplyr::select(Year, Biome, Ecosystem_type, Rh_annual, Rs_annual, Study_midyear, Stage, Partition_method, 
                Land_cover, mat_hadcrut4, map_hadcrut4, SOC) %>%
  filter(!is.na(Stage),
         !is.na(Rs_annual), !is.na(Rh_annual), 
         Year >= SRDB_MINYEAR) %>%
  mutate(Year = as.integer(Year),
         yeargroup = cut(Year, breaks = c(1989, 1998, 2006, 2014),
                         labels = c("1990-1998", "1999-2006", "2007-2014"))) %>%
  group_by(yeargroup) %>% 
  mutate(group_midyear = mean(Year),
         group = paste0(yeargroup, " (N = ", n(), ")")) %>%
  # There are a few studies with NA for Partition_method; assume they're exclusions (by far the most common)
  replace_na(list(Partition_method = "Exclusion")) %>%
  # A few studies are missing SOC from the soilsgrid1km database; replace with biome median 
  left_join(biome_medians, by = c("Biome", "Ecosystem_type")) %>%
  mutate(SOC = if_else(is.na(SOC), SOC_biome, SOC)) -> 
  s_rh_rs

s_rh_rs %>%
  group_by(yeargroup) %>%
  summarise(`Rh:Rs` = mean(Rh_annual / Rs_annual) %>% round(2),
            stdev = sd(Rh_annual / Rs_annual) %>% round(2),
            n = n()) %>%
  print()


m1_rh_rs <- lm(Rh_annual/Rs_annual ~ Study_midyear * Stage + 
                 Study_midyear * Partition_method +
                 Study_midyear * Land_cover +
                 mat_hadcrut4 * map_hadcrut4 ^ 2 + 
                 Study_midyear * SOC, 
               data = s_rh_rs)
m1_rh_rs <- MASS::stepAIC(m1_rh_rs, direction = "both")
print(anova(m1_rh_rs))

m1_rh_rs_signif <- anova(m1_rh_rs)["Study_midyear", "Pr(>F)"]
save_model_diagnostics(m1_rh_rs)
add_result("Rh:Rs ~ time", s_rh_rs$Study_midyear, s_rh_rs$Rh_annual / s_rh_rs$Rs_annual, 
           n = length(m1_rh_rs$residuals),
           p = anova(m1_rh_rs)["Study_midyear", "Pr(>F)"], 
           p_landcover = anova(m1_rh_rs)["Land_cover", "Pr(>F)"],
           p_disturbance = anova(m1_rh_rs)["Stage", "Pr(>F)"])

# Compute summary statistics
s_rh_rs %>%
  group_by(yeargroup) %>%
  summarise(year_min = min(Study_midyear),
            year_max = max(Study_midyear),
            rh_rs_mean = pn(mean(Rh_annual / Rs_annual), 2),
            rh_rs_sd = pn(sd(Rh_annual / Rs_annual), 2), 
            n = n()) ->
  rh_rs_summary

# Make Figure 1
p1_rh_rs <- ggplot(s_rh_rs, aes(Rs_annual, Rh_annual, color = group)) +
  scale_x_log10(label = comma) + scale_y_log10(label = comma) +
  annotation_logticks() +
  xlab(expression(R[S]~(g~C~m^-2~yr^-1))) +
  ylab(expression(R[H]~(g~C~m^-2~yr^-1))) + 
  coord_cartesian(xlim = c(80, 3300), ylim = c(70, 2000))
p1_rh_rs_clr <- p1_rh_rs + scale_color_discrete("Year")
p1_rh_rs_clr2 <- p1_rh_rs + scale_color_brewer("Year")

p_inset <- ggplot(s_rh_rs, aes(Rh_annual / Rs_annual, color = yeargroup, fill = yeargroup)) + 
  geom_density(alpha = 0.5) + 
  xlab(expression(R[H]:R[S])) + ylab("") +
  theme(axis.ticks.y = element_blank(), axis.text.y  = element_blank(),
        axis.text.x = element_text(size = 6), axis.title.x = element_text(size = 8))
p_inset_clr <- p_inset +  scale_fill_discrete(guide = FALSE) +
  scale_color_discrete(guide = FALSE)
p_inset_clr2 <- p_inset +  scale_fill_brewer(guide = FALSE) +
  scale_color_brewer(guide = FALSE)

p1_rh_rs_clr <- p1_rh_rs_clr + 
  annotation_custom(grob = ggplotGrob(p_inset_clr), xmin = log10(60), xmax = log10(800), ymin = log10(600), ymax = log10(2300)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE)
p1_rh_rs_clr2 <- p1_rh_rs_clr2 + 
  annotation_custom(grob = ggplotGrob(p_inset_clr2), xmin = log10(60), xmax = log10(800), ymin = log10(600), ymax = log10(2300)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE)

printlog("NOTE we are plotting this graph with one point cut off:")
printlog(s_rh_rs[which.min(s_rh_rs$Rs_annual), c("Rs_annual", "Rh_annual")])

print(p1_rh_rs_clr)
save_plot("Fig1-srdb-rh-rs-clr", ptype = ".png", width = 9, height = 8)
print(p1_rh_rs_clr2)
save_plot("Fig1-srdb-rh-rs-clr2", ptype = ".png", width = 9, height = 8)

p_rh_rs_time <- ggplot(s_rh_rs, aes(Study_midyear, Rh_annual/Rs_annual, color = Biome)) +
  geom_point() + geom_smooth(method = "lm", aes(group = 1)) +
  xlab("Year") + ylab(expression(R[H]:R[S]))
print(p_rh_rs_time)
save_plot("Fig1-rh_rs_time", width = 6, height = 4)

# Simple regression for Figure 1 - appears in Extended Data Table 1
printlog("Formula for lines in Figure 1:")
mx <- lm(log(Rh_annual) ~ log(Rs_annual) * yeargroup, data = s_rh_rs)
mx <- MASS::stepAIC(mx, direction="both")
print(summary(mx))
print(anova(mx))


# ------------- 2. SRDB Rh:climate analysis --------------- 

printlog(SEPARATOR)
printlog("SRDB Rh climate analysis")

# Compute anomalies from the HadCRUT baseline
srdb$tmp_anom <- srdb$tmp_hadcrut4 - srdb$mat_hadcrut4
srdb$pre_anom <- srdb$pre_hadcrut4 - srdb$map_hadcrut4
srdb$pet_anom <- srdb$pet - srdb$pet_norm

srdb %>%
  filter(!is.na(mat_hadcrut4), !is.na(tmp_anom), 
         !is.na(mat_hadcrut4), !is.na(pre_anom),
         !is.na(pet_norm), !is.na(pet_anom),
         !is.na(Rh_annual), !is.na(Stage)) ->
  s_rh_climate
m2_rh_climate <- lm(sqrt(Rh_annual) ~ tmp_hadcrut4 * pre_hadcrut4 ^ 2 * pet + 
                      tmp_hadcrut4 * Stage + pre_hadcrut4 * Stage + 
                      tmp_hadcrut4 * Leaf_habit + pre_hadcrut4 * Leaf_habit, 
                    data = s_rh_climate)
m2_rh_climate <- stepAIC(m2_rh_climate, direction = "both")
print(summary(m2_rh_climate))
print(anova(m2_rh_climate))

m2_rh_climate_pre_signif <- anova(m2_rh_climate)["pre_hadcrut4", "Pr(>F)"]
m2_rh_climate_pet_signif <- anova(m2_rh_climate)["pet", "Pr(>F)"]
m2_rh_climate_stage_signif <- anova(m2_rh_climate)["Stage", "Pr(>F)"]
save_model_diagnostics(m2_rh_climate)
add_result("Rh ~ time, climate", s_rh_climate$Study_midyear, s_rh_climate$Rh_annual,
           n = length(m2_rh_climate$residuals), 
           p = NA, p_landcover = NA, p_disturbance = m2_rh_climate_stage_signif)

# Global flux computation

printlog("Global flux computation...")
# This is a very simple prediction, based on model above, and deriving
# variance from treating the whole world as evergreen, deciduous, aggrading,
# and mature

# Global grid of climate data
read_csv("outputs/crudata_annual.csv.gz", guess_max = 1e6) %>%
  rename(pre_hadcrut4 = pre, tmp_hadcrut4 = tmp) %>%
  mutate(Leaf_habit = "Deciduous", Stage = "Mature") ->
  globalclim

globalclim$predict1 <- predict(m2_rh_climate, globalclim)
globalclim$Leaf_habit = "Evergreen"
globalclim$predict2 <- predict(m2_rh_climate, globalclim)
globalclim$Stage = "Aggrading"
globalclim$predict3 <- predict(m2_rh_climate, globalclim)
globalclim$Leaf_habit = "Deciduous"
globalclim$predict4 <- predict(m2_rh_climate, globalclim)

# Compute area-weighted fluxes
globalclim %>%
  dplyr::select(-lon, -lat, -Leaf_habit, -Stage) %>%
  gather(case, sqrt_rh_gCm2, predict1, predict2, predict3, predict4) %>%
  mutate(rh_PgC = sqrt_rh_gCm2 ^ 2 * area_km2 * 1000 * 1000 / 1e15) %>%
  group_by(case, year) %>%
  summarise(rh_PgC = sum(rh_PgC, na.rm = TRUE),
            tmp_hadcrut4 = weighted.mean(tmp_hadcrut4, area_km2, na.rm = TRUE)) ->
  gp

# Fit linear models to the predictions to get smoothed values for Q10, etc.
gp %>%
  group_by(case) %>% 
  do(data.frame(year = .$year, 
                rh_fit = predict(lm(rh_PgC ~ year, data = .)),
                tmp_fit = predict(lm(tmp_hadcrut4 ~ year, data = .)))) %>%
  right_join(gp, by = c("case", "year")) -> 
  gp

# Summarise: compute first and last values and Q10
gp %>%
  group_by(case) %>%
  summarise(first_rh = first(rh_fit), last_rh = last(rh_fit),
            first_tmp = first(tmp_fit), last_tmp = last(tmp_fit),
            q10 = (last_rh / first_rh) ^ (10 / (last_tmp - first_tmp))) %>%
  print ->
  gp_summary

p_prediction <- qplot(year, rh_PgC, data = gp, geom = "line", group = case) + 
  geom_smooth(method = "lm", group = 1)
print(p_prediction)
save_plot("global_rh_prediction")

save_data(gp_summary, fname = "global_prediction")


# --------------- 3. FLUXNET analysis --------------------- 

printlog(SEPARATOR)
printlog("FLUXNET data analysis")

printlog("Non-NA IGBP =", sum(!is.na(srdb$IGBP)))
s <- subset(srdb, !is.na(IGBP))
printlog("... & ecosystem match =", sum(s$FLUXNET_ECOSYSTEM_MATCH))
s <- subset(s, FLUXNET_ECOSYSTEM_MATCH)
printlog("... & dist =", sum(s$FLUXNET_DIST <= MAX_FLUXNET_DIST))
s <- subset(s, FLUXNET_DIST <= MAX_FLUXNET_DIST)
printlog("... & QC =", sum(!is.na(s$NEE_VUT_REF_QC) & s$NEE_VUT_REF_QC >= MIN_NEE_QC))

printlog("Filtering to MAX_FLUXNET_DIST =", MAX_FLUXNET_DIST)
printlog("Filtering to MIN_NEE_QC =", MIN_NEE_QC)

# Filter SRDB for available data, distance to FLUXNET tower, and NEE quality
srdb %>%
  filter(!is.na(Rs_annual),
         FLUXNET_ECOSYSTEM_MATCH,
         FLUXNET_DIST <= MAX_FLUXNET_DIST,
         NEE_VUT_REF_QC >= MIN_NEE_QC) %>%
  # We only allow one observation per site per year
  group_by(FLUXNET_SITE_ID, IGBP, tmp_trend, pre_trend, Land_cover, Stage, Year) %>%
  summarise(Rs_annual = mean(Rs_annual),
            Rh_annual = mean(Rh_annual),
            RECO_NT_VUT_REF = mean(RECO_NT_VUT_REF),
            gpp_fluxnet = mean(gpp_fluxnet),
            YearsOfData = mean(YearsOfData),
            pre_hadcrut4 = mean(pre_hadcrut4),
            mat_hadcrut4 = mean(mat_hadcrut4),
            map_hadcrut4 = mean(map_hadcrut4)) %>%
  mutate(tmp_trend_label = if_else(tmp_trend > 0, "Warming", "Cooling"),
         pre_trend_label = if_else(pre_trend > 0, "Wetter", "Drier")) ->
  s_fluxnet

save_data(s_fluxnet, scriptfolder = FALSE)

s_fluxnet %>%
  gather(flux, fluxvalue, gpp_fluxnet, Rs_annual) %>%
  ggplot(aes(Year, fluxvalue, color = flux)) + 
  geom_line() +
  facet_wrap(~FLUXNET_SITE_ID, scales = "free")
save_plot("fluxnet_site_diagnostic")

m_fluxnet <- lm(Rs_annual/gpp_fluxnet ~ Year * Land_cover + mat_hadcrut4 * map_hadcrut4 ^ 2, 
                data = s_fluxnet, weights = YearsOfData)
m_fluxnet <- MASS::stepAIC(m_fluxnet, direction = "both")
print(anova(m_fluxnet))
add_result("Rs/GPPfluxnet ~ time", s_fluxnet$Year, s_fluxnet$Rs_annual / s_fluxnet$gpp_fluxnet,
           n = length(m_fluxnet$residuals), p = anova(m_fluxnet)["Year", "Pr(>F)"], 
           p_landcover = anova(m_fluxnet)["Land_cover", "Pr(>F)"],
           p_disturbance = anova(m_fluxnet)["Stage", "Pr(>F)"])

p <- qplot(Year, Rs_annual/gpp_fluxnet, color = pre_trend_label, data = s_fluxnet) + 
  geom_smooth(method = "lm", na.rm = TRUE)
print(p)
save_plot("fluxnet_basic")

s_fluxnet_nohf <- subset(s_fluxnet, FLUXNET_SITE_ID != "US-Ha1")
m_fluxnet_nohf <- lm(Rs_annual/gpp_fluxnet ~ Year * Land_cover + mat_hadcrut4 * map_hadcrut4, 
                     data = s_fluxnet_nohf, weights = YearsOfData)
m_fluxnet_nohf <- MASS::stepAIC(m_fluxnet_nohf, direction = "both")
print(anova(m_fluxnet_nohf))
print(p %+% s_fluxnet_nohf)
save_plot("fluxnet_nohf_basic")

# Make plot
s_fluxnet %>%
  group_by(FLUXNET_SITE_ID, IGBP) %>%
  summarise(Year = min(Year), 
            ratio = (Rs_annual / gpp_fluxnet)[which.min(Year)],
            pre_trend_label = paste(pre_trend_label, collapse = " ")) ->
  s_fluxnet_labels

p_fluxnet <- ggplot(s_fluxnet, aes(Year, Rs_annual / gpp_fluxnet, group = FLUXNET_SITE_ID)) + 
  geom_point(aes(color = IGBP), na.rm = TRUE) + 
  #  geom_line(aes(color = IGBP), na.rm = TRUE) +
  geom_smooth(method = "lm", color = "grey", fill = NA, na.rm = TRUE) + 
  geom_smooth(method = "lm", group = 1, na.rm = TRUE) +
  ylab(expression(R[S]:GPP[fluxnet]))

p_fluxnet <- p_fluxnet + geom_text(data = s_fluxnet_labels, 
                                   aes(y = ratio, label = FLUXNET_SITE_ID, color = IGBP), 
                                   size = 2, nudge_y = 0.05,
                                   alpha = 0.75)
print(p_fluxnet)

save_plot("gpp_fluxnet")

# How many _sites_ have positive trends, versus not?
s_fluxnet %>%
  group_by(FLUXNET_SITE_ID) %>%
  do(sitemod = lm(Rs_annual / gpp_fluxnet ~ Year, data = .)) %>%
  tidy(sitemod) %>%
  filter(term == "Year") %>%
  mutate(trend = if_else(sign(estimate) > 0, "Rising", "Not rising")) %>%
  group_by(trend) %>%
  summarise(n = n()) ->
  x
fluxnet_Rs_GPP_sitetrends = x$n
names(fluxnet_Rs_GPP_sitetrends) <- x$trend
print(fluxnet_Rs_GPP_sitetrends)

# --------------- 4. FLUXNET-only analysis --------------------- 

# Rodrigo's suggestion: look at ratio of annual nighttime NEE to
# annual GPP in the FLUXNET data. Because nighttime NEE is presumably
# dominated by soil respiration, we would expect, if Rh is increasing
# relative to GPP, that NEEnight would also increase relative to GPP.
printlog(SEPARATOR)
readr::read_csv("outputs/fluxnet.csv") %>%
  filter(NEE_VUT_REF_NIGHT > 0, GPP_DT_VUT_REF > 0,
         NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF < 1,
         !is.na(TA_F), !is.na(P_F), !is.na(GPP_DT_VUT_REF)) ->
  s_fluxnet_only

# Record lengths of FLUXNET sites
s_fluxnet_only %>% 
  group_by(SITE_ID) %>% 
  summarise(nyears = n()) %>% 
  dplyr::select(-SITE_ID) %>% 
  summarise_all(funs(mean, median, sd, min, max)) -> 
  s_fluxnet_nyears

# Make table for acknowledging FLUXNET publications
s_fluxnet_only %>%
  group_by(SITE_ID) %>%
  summarise(Years = paste0(min(Year), "-", max(Year)), Ref = "") %>%
  save_data(fname = "SupplementaryTable1.csv")

s_fluxnet_only %>%
  group_by(SITE_ID) %>%
  summarise(nyears = n()) %>%
  filter(nyears > 4) %>%
  left_join(s_fluxnet_only, by = "SITE_ID") %>%
  group_by(IGBP, SITE_ID) %>%
  do(mod = lm(NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF ~ Year, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Year") %>%
  mutate(rising = if_else(estimate > 0, "Rising", "Not\nrising")) -> site_trends 
printlog("Sites with positive NEE night/GPP trends =", sum(site_trends$estimate > 0), "of", nrow(site_trends))
printlog("Sites with SIGNIFICANT positive NEE night/GPP trends =", sum(site_trends$p.value < 0.05 & site_trends$estimate > 0), "of", nrow(site_trends))

p_fluxnet_sitecounts <- ggplot(site_trends, aes(x = rising)) + 
  geom_bar() + facet_wrap(~IGBP) +
  xlab(expression(NEE[night]:GPP[fluxnet]))

s_fluxnet_only %>%
  group_by(IGBP) %>%
  summarise(nyears = n()) %>%
  filter(nyears > 2) %>%
  left_join(s_fluxnet_only, by = "IGBP") %>%
  group_by(IGBP) %>%
  do(mod = lm(NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF ~ Year, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Year") -> igbp_trends 
printlog("IGBPs with positive NEE night/GPP trends =", sum(igbp_trends$estimate > 0), "of", nrow(igbp_trends))
printlog("IGBPs with SIGNIFICANT positive NEE night/GPP trends =", sum(igbp_trends$p.value < 0.05 & igbp_trends$estimate > 0), "of", nrow(igbp_trends))

p_fluxnet_only <- ggplot(s_fluxnet_only, aes(Year, NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF, group = SITE_ID)) + 
  geom_point(aes(color = IGBP), na.rm = TRUE) + 
  geom_smooth(method = "lm", color = "grey", linetype = 2, size = 0.2, fill = NA, na.rm = TRUE) + 
  geom_smooth(method = "lm", group = 1, na.rm = TRUE) +
  ylab(expression(NEE[night]:GPP[fluxnet])) +
  facet_wrap(~ IGBP) +
  coord_cartesian(xlim = c(1990, 2016)) +
  ylim(c(0, 1))

print(p_fluxnet_only)
save_plot("EDF6-fluxnet_only")
save_plot("EDF6-fluxnet_only", ptype = ".png", height = 6, width = 8)

printlog("Fitting temporal model...")
m_fluxnet_only <- lm(NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF ~ Year * TA_F + Year * P_F + Year * IGBP, data = s_fluxnet_only)
m_fluxnet_only <- MASS::stepAIC(m_fluxnet_only, direction = "both")
print(anova(m_fluxnet_only))
m_fluxnet_only_trend_signif <- anova(m_fluxnet_only)["Year", "Pr(>F)"]
m_fluxnet_only_precip_trend_signif <- anova(m_fluxnet_only)["Year:p_trend", "Pr(>F)"]
save_model_diagnostics(m_fluxnet_only)
add_result("FLUXNET NEEnight / GPP ~ time", s_fluxnet_only$Year, s_fluxnet_only$NEE_VUT_REF_NIGHT / s_fluxnet_only$GPP_DT_VUT_REF,
           n = length(m_fluxnet_only$residuals), p = m_fluxnet_only_trend_signif, 
           p_landcover = anova(m_fluxnet_only)["IGBP", "Pr(>F)"],
           p_disturbance = NA)

# ----------- 5. GPP and SIF analysis -------------- 

printlog(SEPARATOR)
printlog("Remote sensing analysis")

srdb %>%
  dplyr::select(Study_midyear, Biome, Land_cover, Partition_method, Stage,
                mat_hadcrut4, map_hadcrut4,
                gpp_beer, gpp_modis, 
                SCIA_SIF, GOME2_SIF,
                Rs_annual, Rh_annual, SOC) %>%
  rename(`GPP[MTE]` = gpp_beer,
         `GPP[MODIS]` = gpp_modis) %>%
  # since we're going to divide, can't have SIF crossing 1!
  # rescale to GPP range for convenience 
  mutate(`SIF[SCIAMACHY]` = rescale(SCIA_SIF, 0, 3500)) %>%
  mutate(`SIF[GOME2]` = rescale(GOME2_SIF, 0, 3500)) %>%
  gather(Flux, fluxvalue, Rs_annual, Rh_annual) %>%
  mutate(Prettyflux = if_else(Flux == "Rs_annual", "R[S]", "R[H]")) %>%
  gather(GPPSIF, gppsifvalue, `GPP[MTE]`, `GPP[MODIS]`, `SIF[SCIAMACHY]`, `SIF[GOME2]`) %>%
  mutate(GPPSIF = factor(GPPSIF, levels = c("GPP[MTE]", "GPP[MODIS]", "SIF[SCIAMACHY]", "SIF[GOME2]"))) %>%
  filter(!is.na(Land_cover), !is.na(gppsifvalue), !is.na(SOC)) ->
  s_gppsif

s_gppsif %>%
  filter(fluxvalue / gppsifvalue <= MAX_FLUX_TO_GPP) ->
  s_gppsif_included

save_data(s_gppsif_included, fname = basename(SRDB_GPPSIF_FILE), scriptfolder = FALSE)

p_gppsif_base <- ggplot(s_gppsif_included, aes(Study_midyear, fluxvalue / gppsifvalue, color = Land_cover)) +
  geom_point(alpha = I(0.75), size = 0.5) +
  facet_grid(GPPSIF ~ Prettyflux, scales = "free", labeller = label_parsed) +
  scale_color_discrete("Land cover") +
  xlab("Year") +
  ylab("Respiration:(GPP or SIF)") +
  coord_cartesian(ylim = c(0, 2)) + 
  scale_x_continuous(breaks = c(1990, 1995, 2000, 2005, 2010, 2015),
                     labels = c(1990, 1995, 2000, 2005, 2010,""))
p_gppsif <- p_gppsif_base + 
  geom_smooth(data = filter(s_gppsif_included, 
                            ! Land_cover %in% c("Other"), 
                            (GPPSIF != "SIF[GOME2]" | Flux == "Rs_annual")),
              method = "lm", show.legend = FALSE)
print(p_gppsif )
save_plot("Fig2-gppsif", ptype = ".png", height = 8, width = 7)

printlog("Trend tests")
s_gppsif_resids <- list()
for(dataset in unique(s_gppsif_included$GPPSIF)) {
  for(f in unique(s_gppsif_included$Flux)) {
    d <- filter(s_gppsif_included, Flux == f, GPPSIF == dataset, !is.na(mat_hadcrut4), !is.na(map_hadcrut4))
    dname <- make.names(paste("s", dataset, f, sep = "_"))
    assign(dname, d)
    save_data(d, fname = dname, scriptfolder = FALSE)
    
    # Linear model without time--for a residual plot (per Referee 3)
    m <- lm(fluxvalue / gppsifvalue ~ mat_hadcrut4 + map_hadcrut4 ^ 2 + Land_cover + Stage + SOC, 
            data = d)
    m <- stepAIC(m, direction = "both", trace = 0)
    d$residual <- m$residuals
    s_gppsif_resids[[paste(dataset, f)]] <- d   # save this for plotting later
    
    # Linear model
    m <- lm(fluxvalue / gppsifvalue ~ mat_hadcrut4 + map_hadcrut4 ^ 2 + 
              Study_midyear * Land_cover + 
              Study_midyear * Stage + SOC, 
            data = d)
    m <- stepAIC(m, direction = "both", trace = 0)
    print(summary(m))
    mn <- make.names(paste("m", dataset, f, sep = "_"))
    printlog(mn)
    save_model_diagnostics(m, modelname = mn)
    assign(mn, m)
    
    add_result(paste(f, "/", dataset, "~ time"), d$Study_midyear, d$fluxvalue / d$gppsifvalue,
               n = length(m$residuals), p = anova(m)["Study_midyear", "Pr(>F)"], 
               p_landcover = anova(m)["Land_cover", "Pr(>F)"], 
               p_disturbance = anova(m)["Stage", "Pr(>F)"])
  }
}

printlog("Making residual plot...")
s_gppsif_resids <- bind_rows(s_gppsif_resids)
save_data(s_gppsif_resids)
p_gppsif_resid <- ggplot(s_gppsif_resids, aes(Study_midyear, residual)) +
  geom_point(alpha = I(0.75), size = 0.5) +
  facet_grid(GPPSIF ~ Prettyflux, scales = "free", labeller = label_parsed) +
  xlab("Year") + ylab("Model residual") + geom_smooth(method = "lm")
print(p_gppsif_resid)
save_plot("EDF3-gppsif-residuals", ptype = ".png", height = 8, width = 7)


# ----------- 6. ISIMIP analysis (per Referee 1) -------------- 

printlog("Testing ratio of RH to ISIMIP GPP...")

srdb %>%
  dplyr::select(Longitude, Latitude, Rh_annual, mat_hadcrut4, map_hadcrut4, Study_midyear, gpp_isimip,
                Land_cover, Stage, SOC) %>%
  na.omit ->
  srdb_isimip

# Linear model
m_isimip <- lm(Rh_annual / gpp_isimip ~ mat_hadcrut4 + map_hadcrut4 ^ 2 + 
                 Study_midyear * Land_cover + 
                 Study_midyear * Stage + SOC, 
               data = srdb_isimip)
m_isimip <- stepAIC(m_isimip, direction = "both", trace = 0)
print(summary(m_isimip))

add_result("Rh / GPPisimip ~ time", srdb_isimip$Study_midyear, srdb_isimip$Rh_annual / srdb_isimip$gpp_isimip,
           n = length(m_isimip$residuals), p = anova(m_isimip)["Study_midyear", "Pr(>F)"], 
           p_landcover = anova(m_isimip)["Land_cover", "Pr(>F)"],
           p_disturbance = anova(m_isimip)["Stage", "Pr(>F)"])

# Plot
p_isimip <- ggplot(srdb, aes(Study_midyear, Rh_annual / gpp_isimip)) + 
  geom_point(aes(color = Land_cover)) + geom_smooth(method = "lm") +
  xlab("Year") + ylab(expression(R[H]:GPP[isimip])) + 
  theme(legend.position = c(0.25, 0.85)) + scale_color_discrete("")
print(p_isimip)
save_plot("isimip", width = 7, height = 4)

# New test per Referee 1: how well explained by climate is ISIMIP RH?
printlog(SEPARATOR)
printlog("How well is ISIMIP RH explained by climate?")
read_csv("ancillary/isimip-rh/all_annual_site_rh.csv") %>% 
  gather(year, rh, matches("[0-9]{4}")) %>% 
  mutate(year = as.integer(year), lat = round(lat, 2), lon = round(lon, 2)) %>% 
  # Compute a mean ISIMIP RH for each lon/lat/year
  group_by(lat, lon, year) %>% 
  summarise(rh_isimip = mean(rh)) ->
  isimip_site_rh

# Merge the ISIMIP GPP data with our observational data
srdb_isimip %>% 
  mutate(Study_midyear = round(Study_midyear, 0),
         Longitude = round(Longitude, 2), 
         Latitude = round(Latitude, 2)) %>% 
  left_join(isimip_site_rh, by = c("Longitude" = "lon", "Latitude" = "lat", "Study_midyear" = "year")) ->
  srdb_isimip

# Now fit our basic linear model
m_rh_isimip <- lm(rh_isimip ~ Study_midyear * Stage + 
                    Study_midyear * Land_cover +
                    mat_hadcrut4 * map_hadcrut4 ^ 2 + 
                    Study_midyear * SOC, 
                  data = srdb_isimip)
m_rh_isimip <- MASS::stepAIC(m_rh_isimip, direction = "both")
print(summary(m_rh_isimip))


# ----------- 7. Sensitivity to start/end date (per Referee 1) -------------- 

printlog(SEPARATOR)
printlog("Checking sensitivity of results to start date (note this generates some errors in red)...")

startdate_results <- list()
for(minyr in 1981:(max(srdb_complete$Study_midyear) - 5)) {
  printlog(minyr)
  for(maxyr in (minyr + 5):max(srdb_complete$Study_midyear)) {
    # Repeatedly fit our basic Rh/Rs time model to shorter and shorter time periods
    d <- srdb_complete %>%
      filter(Study_midyear >= minyr, Study_midyear <= maxyr,
             Ecosystem_state != "Managed") %>%
      dplyr::select(Rh_annual, Rs_annual, Study_midyear, Stage, Partition_method, Land_cover,
                    map_hadcrut4, mat_hadcrut4, SOC)
    d <- d[complete.cases(d),]
    try({
      m <- lm(Rh_annual/Rs_annual ~ Study_midyear * Stage + 
                Study_midyear * Partition_method +
                Study_midyear * Land_cover +
                mat_hadcrut4 * map_hadcrut4 ^ 2 + 
                Study_midyear * SOC, 
              data = d)
      m <- MASS::stepAIC(m, direction = "both", trace = 0)
      startdate_results[[paste(minyr, maxyr)]] <- tibble(min_year = minyr, max_year = maxyr,
                                                         p = round(anova(m)["Study_midyear", "Pr(>F)"], 3),
                                                         n = length(m$residuals))
    })
    
  }
}
startdate_results <- bind_rows(startdate_results)


p_startdate <- ggplot(startdate_results, aes(min_year, p, color = (p <= 0.05))) + 
  geom_point() +
  geom_hline(yintercept = 0.05, linetype = 2) +
  geom_hline(yintercept = 0.10, linetype = 2) +
  geom_vline(xintercept = 1989, linetype = 2) +
  guides(color = FALSE) +
  geom_line(aes(y = n / max(n)), color = "black") +
  xlab("First year in dataset") + ylab("P-value of Rh:Rs trend")
print(p_startdate)
save_plot("startdate", width = 7, height = 4)

p_startdate_2d <- ggplot(startdate_results, aes(min_year, max_year)) + 
  geom_tile(aes(fill = p < 0.05)) + 
  xlab("First year in dataset") + ylab("Last year in dataset") +
  annotate(geom = "point", x = 1989, y = 2014, pch = 5, size = 3)
print(p_startdate_2d)
save_plot("EDF1-startdate_2d", width = 7, height = 4)

# ----- 8. Site-specific trends, both managed and unmanaged (per Referee 1) -------------- 

printlog(SEPARATOR)
printlog("Doing Referee 1 site-level analysis...")

R1_MIN_TIMESPAN <- 8
R1_MIN_OBS <- 3
srdb_complete %>% 
  mutate(Longitude = round(Longitude, 1), Latitude = round(Latitude, 1)) ->
  srdb_complete_rounded

srdb_complete_rounded %>% 
  group_by(Longitude, Latitude, Ecosystem_state, Land_cover) %>% 
  summarise(nyears = length(unique(Study_midyear)), 
            pre_trend = mean(pre_trend),
            tmp_trend = mean(tmp_trend),
            year_range = round(max(Study_midyear) - min(Study_midyear), 0),
            Site_name = substr(paste(unique(Site_name), collapse = " "), 1, 20)) %>% 
  filter(year_range >= R1_MIN_TIMESPAN, nyears >= R1_MIN_OBS) %>%
  arrange(desc(year_range)) -> 
  longterm_sites_list

printlog("There are", nrow(longterm_sites_list), "sites with at least", R1_MIN_OBS, "observations over at least", R1_MIN_TIMESPAN, "years")
print(longterm_sites_list)
save_data(longterm_sites_list)

# Filter the complete database to these sites, run regressions
srdb_complete_rounded %>%
  dplyr::select(Study_midyear, Longitude, Latitude, Ecosystem_state, Land_cover,
                Rs_annual, Rh_annual, tmp_trend, pre_trend) %>%
  semi_join(longterm_sites_list, by = c("Longitude", "Latitude")) ->
  srdb_complete_rounded

srdb_complete_rounded %>%
  filter(!is.na(Rs_annual)) %>%
  group_by(Longitude, Latitude, Ecosystem_state, Land_cover, tmp_trend, pre_trend) %>% 
  do(rsmod = lm(Rs_annual ~ Study_midyear, data = .)) %>%
  tidy(rsmod) %>%
  filter(term == "Study_midyear") %>%
  mutate(trend = if_else(sign(estimate) > 0, "Rising", "Not rising")) ->
  rs_site_individual

p <- ggplot(rs_site_individual, aes(tmp_trend, pre_trend, color = trend)) + 
  geom_point() + ggtitle("Rs longitudinal trends") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
print(p)
save_plot("Rs_longitudinal_trends_climate")

rs_site_individual %>%
  group_by(Land_cover, tmp_trend, pre_trend, trend) %>% 
  summarise(n = n()) %>%
  complete(nesting(Land_cover, trend), fill = list(n = 0)) ->
  rs_site_summary

printlog(rs_site_summary %>% spread(trend, n))

p_rs_sites <- ggplot(rs_site_summary, aes(trend, n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Land_cover) +
  xlab(expression(R[S]~trend)) + ylab("Number of sites")
print(p_rs_sites)
save_plot("rs_all_sites", height = 6, width = 8)

# ...and again for Rh
srdb_complete_rounded %>%
  filter(!is.na(Rh_annual)) %>%
  group_by(Longitude, Latitude, Ecosystem_state, Land_cover, tmp_trend, pre_trend) %>% 
  do(rhmod = lm(Rh_annual ~ Study_midyear, data = .)) %>%
  tidy(rhmod) %>%
  filter(term == "Study_midyear") %>%
  mutate(trend = if_else(sign(estimate) > 0, "Rising", "Not\nrising")) ->
  rh_site_individual

p <- ggplot(rh_site_individual, aes(tmp_trend, pre_trend, color = trend)) + 
  geom_point(size = 3) + ggtitle("Rh longitudinal trends") +
  xlab("Temperature trend (C/yr/yr)") + ylab("Precipitation trend (mm/yr/yr)") +
  geom_hline(yintercept = 0) + geom_vline(xintercept = 0)
print(p)
save_plot("Rh_longitudinal_trends_climate")

rh_site_individual %>%
  group_by(Land_cover, trend) %>% 
  summarise(n = n()) ->
  rh_site_summary

printlog(rh_site_summary %>% spread(trend, n))

p_rh_sites <- ggplot(rh_site_summary, aes(trend, n)) + 
  geom_bar(stat = "identity") + 
  facet_wrap(~Land_cover) +
  xlab(expression(R[H]~trend)) + ylab("Number of sites")
print(p_rh_sites)
save_plot("rh_all_sites", height = 3, width = 8)

p_rh_site_climatespace <- ggplot(srdb, aes(tmp_trend, pre_trend)) +
  geom_vline(xintercept = 0) + geom_hline(yintercept=0) + 
  geom_jitter(size = 2, alpha = 0.2, color = "darkgrey") + 
  geom_point(data = longterm_sites_list, aes(color = Land_cover), size = 3) +
  xlab("1991-2010 temp trend, degC/yr") + ylab("1991-2010 precip trend, mm/yr")
print(p_rh_site_climatespace)
save_plot("EDF9-rh_site_climatespace", height = 5, width = 8)

printlog("Climate space of main dataset:")
print(summary(dplyr::select(srdb, tmp_trend, pre_trend)))
printlog("Climate space of longterm sites:")
print(summary(dplyr::select(longterm_sites_list, tmp_trend, pre_trend)))

# Late Referee 1 suggestion - look at trends across all these longitudinal sites

srdb_complete_rounded %>% 
  semi_join(dplyr::select(ungroup(rh_site_individual), Longitude, Latitude)) -> 
  srdb_13sites

printlog("Models across all the longitudinal sites:")
m_longitudinal_rh <- lm(Rh_annual ~ Study_midyear + tmp_trend * pre_trend, data = srdb_13sites)
print(summary(m_longitudinal_rh))

# ----- 9. Some SRDB stats (per Referee 3) -------------- 



# ----------------- Multipanel figure 2 --------------------- 

printlog("Doing new multipanel Figure 2...")
library(cowplot)  # 0.9.2
theme_set(theme_bw())
xts <- 7  # x axis text size
sts <- 8  # strip (facet) text size
p_gppsif1 <- p_gppsif + theme(axis.text = element_text(size = xts),
                              strip.text = element_text(size = sts))
p_isimip1 <- p_isimip + theme(axis.text = element_text(size = xts),
                              legend.text = element_text(size = 5),
                              legend.background = element_rect(fill = "transparent", colour = NA),
                              legend.key.height = unit(0.25, "cm"))
p_rh_sites1 <- p_rh_sites + theme(axis.text = element_text(size = xts),
                                  strip.text = element_text(size = 5))
p_fluxnet_sitecounts1 <- p_fluxnet_sitecounts + theme(axis.text = element_text(size = xts),
                                                      strip.text = element_text(size = 6))

fig2_multipanel <- ggdraw() +
  draw_plot(p_gppsif1 + guides(colour = FALSE), 0, 0, 0.6, 1) +
  # draw_plot(p_gppsif1 + guides(colour = guide_legend(NULL, nrow = 2)) + 
  #             theme(legend.position = "bottom"), 0, 0, 0.6, 1) +
  draw_plot(p_isimip1, 0.6, 0.5, 0.39, 0.5) +
  #  draw_plot(p_fluxnet_sitecounts1, 0.6, 0, 0.4, 0.5) +
  draw_plot(p_rh_sites1, 0.6, 0, 0.4, 0.5) +
  draw_plot_label(c("a", "b", "c"), c(0, 0.6, 0.6), c(1, 1, 0.5), size = 15)

cowplot::save_plot(file.path(outputdir(), "Fig2-multipanel.pdf"), fig2_multipanel, base_aspect_ratio = 1.5)


# ----------------------- Clean up ------------------------- 

overall_results_table <- bind_rows(overall_results)
overall_results_table %>%
  mutate(p = pclean(p), 
         p_theil.sen = pclean(p_theil.sen),
         p_land.cover = pclean(p_land.cover),
         p_disturbance = pclean(p_disturbance)) ->
  overall_results_table_clean
printlog(overall_results_table_clean)
save_data(overall_results_table)
save_data(overall_results_table_clean)

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
