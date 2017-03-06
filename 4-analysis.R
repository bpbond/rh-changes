# Heterotrophic respiration - main analysis script
# Ben Bond-Lamberty February 2017
#
# A lot goes on here, but there are four main sections below:
# 1. The Rh:Rs analysis: is the ratio changing over time?
# 2. Is Rh responding to climate anomalies? (Minor; parallels 2010 paper)
# 3. Is Rs:FLUXNET GPP rising over time?
# 4. For FLUXNET only, is nighttime NEE:GPP rising?
# 5. Remote sensing: are Rh:GPP, Rs:GPP, Rh:SIF, and Rs:SIF changing over time?
# This script saves plots, but also leaves a bunch of stuff in-memory
# for the RMarkdown script to access.

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
  filter(!is.na(Stage), !is.na(Leaf_habit), 
         !is.na(Rs_annual), !is.na(Rh_annual), 
         Year >= 1989) %>%
  mutate(Year = as.integer(Year),
         yeargroup = cut(Year, breaks = c(1989, 1994, 1999, 2004, 2009, 2014),
                         labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014"))) %>%
  group_by(yeargroup) %>% 
  mutate(group_midyear = mean(Year),
         group = paste0(yeargroup, " (N = ", n(), ")")) -> 
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

srdb %>%
  filter(!is.na(mat_hadcrut4), !is.na(tmp_anom), 
         !is.na(mat_hadcrut4), !is.na(pre_anom),
         !is.na(pet_norm), !is.na(pet_anom),
         !is.na(Rh_annual), !is.na(Leaf_habit), !is.na(Stage)) ->
  s_rh_climate
m2_rh_climate <- lm(sqrt(Rh_annual) ~ tmp_hadcrut4 * pre_hadcrut4 * pet + Stage * Leaf_habit, 
     data = s_rh_climate)
m2_rh_climate <- stepAIC(m2_rh_climate, direction = "both")
print(anova(m2_rh_climate))

m2_rh_climate_pre_signif <- anova(m2_rh_climate)["pre_hadcrut4", "Pr(>F)"]
m2_rh_climate_pet_signif <- anova(m2_rh_climate)["pet", "Pr(>F)"]
m2_rh_climate_stage_signif <- anova(m2_rh_climate)["Stage", "Pr(>F)"]
save_model_diagnostics(m2_rh_climate)

# Global flux computation

read_csv("outputs/crudata_annual.csv.gz") %>%
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
globalclim %>%
  dplyr::select(-lon, -lat, -Leaf_habit, -Stage) %>%
  gather(case, sqrt_rh_gCm2, predict1, predict2, predict3, predict4) %>%
  mutate(rh_PgC = sqrt_rh_gCm2 ^ 2 * area_km2 * 1000 * 1000 / 1e15) %>%
  group_by(year, case) %>%
  summarise(rh_PgC = sum(rh_PgC, na.rm = TRUE),
            tmp_hadcrut4 = weighted.mean(tmp_hadcrut4, area_km2, na.rm = TRUE)) %>%
  summarise(rh_PgC_sd = sd(rh_PgC), 
            rh_PgC = mean(rh_PgC), 
            tmp_hadcrut4 = mean(tmp_hadcrut4)) ->
  gp

p_prediction <- qplot(year, rh_PgC, data = gp, geom = "line") + 
  geom_ribbon(aes(ymin=rh_PgC - rh_PgC_sd, ymax = rh_PgC + rh_PgC_sd), alpha = 0.25) +
  geom_smooth(method = "lm")
print(p_prediction)
save_plot("global_rh_prediction")

slope_model <- lm(rh_PgC ~ year, data = gp)
gp$predict <- predict(slope_model)
global_rh_begin <- gp$rh_PgC[1]
global_rh_end <- gp$rh_PgC[nrow(gp)]
global_rh_end_sd <- gp$rh_PgC_sd[nrow(gp)]
global_q10 <- (global_rh_end / global_rh_begin) ^ (10 / (gp$tmp_hadcrut4[nrow(gp)] - gp$tmp_hadcrut4[1]))

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
  group_by(FLUXNET_SITE_ID, IGBP, tmp_trend, pre_trend, Leaf_habit, Stage, Year) %>%
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

printlog("Mann-Kendall trend test:")
mk3_fluxnet <- MannKendall(s_fluxnet$Rs_annual / s_fluxnet$gpp_fluxnet)
print(mk3_fluxnet)

m_fluxnet <- lm(Rs_annual/gpp_fluxnet ~ Year * Leaf_habit + mat_hadcrut4 * map_hadcrut4, 
                data = s_fluxnet, weights = YearsOfData)
m_fluxnet <- MASS::stepAIC(m_fluxnet, direction = "both")
print(anova(m_fluxnet))
p <- qplot(Year, Rs_annual/gpp_fluxnet, color=pre_trend_label, data=s_fluxnet) + 
  geom_smooth(method = "lm", na.rm = TRUE)
print(p)
save_plot("fluxnet_basic")

s_fluxnet_nohf <- subset(s_fluxnet, FLUXNET_SITE_ID != "US-Ha1")
m_fluxnet_nohf <- lm(Rs_annual/gpp_fluxnet ~ Year * Leaf_habit + mat_hadcrut4 * map_hadcrut4, 
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
            pre_trend_label = unique(pre_trend_label)) ->
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

printlog("Computing warming/cooling/drying/wetting trends from tower data")
s_fluxnet_only %>%
  group_by(SITE_ID) %>% 
  summarise(nyears = n()) %>%
  filter(nyears > 2) %>%
  left_join(s_fluxnet_only, by = "SITE_ID") ->
  s_fluxnet_only

s_fluxnet_only %>%
  group_by(SITE_ID) %>%
  do(modt = lm(TA_F ~ Year, data = .),
     modp = lm(P_F ~ Year, data = .)) ->
  fluxnet_only_mods

fluxnet_only_mods %>%
  tidy(modt) %>%
  filter(term == "Year") %>%
  mutate(t_trend = if_else(sign(estimate) > 0, "Warmer", "Cooler")) %>%
  dplyr::select(SITE_ID, t_trend) ->
  t_trends
fluxnet_only_mods %>%
  tidy(modp) %>%
  filter(term == "Year") %>%
  mutate(p_trend = if_else(sign(estimate) > 0, "Wetter", "Drier")) %>%
  dplyr::select(SITE_ID, p_trend) %>%
  right_join(s_fluxnet_only, by = "SITE_ID") %>%
  left_join(t_trends, by = "SITE_ID") %>%
  filter(!is.na(t_trend), !is.na(p_trend)) ->
  s_fluxnet_only

p_fluxnet_only <- ggplot(s_fluxnet_only, aes(Year, NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF, group = SITE_ID)) + 
  geom_point(aes(color = IGBP), na.rm = TRUE) + 
  geom_smooth(method = "lm", color = "grey", fill = NA, na.rm = TRUE) + 
  geom_smooth(method = "lm", group = 1, na.rm = TRUE) +
  ylab(expression(NEE[night]:GPP[fluxnet])) +
  facet_grid(t_trend ~ p_trend) +
  ylim(c(0, 1))

print(p_fluxnet_only)
save_plot("fluxnet_only")

printlog("Fitting temporal model...")
m_fluxnet_only <- lm(NEE_VUT_REF_NIGHT / GPP_DT_VUT_REF ~ Year * t_trend * p_trend, data = s_fluxnet_only)
m_fluxnet_only <- MASS::stepAIC(m_fluxnet_only, direction = "both")
print(anova(m_fluxnet_only))
m_MODIS.GPP_Rs_annual_precip_trend_signif <- anova(m_fluxnet_only)["Year:p_trend", "Pr(>F)"]
save_model_diagnostics(m_fluxnet_only)


# ----------- 5. GPP and SIF analysis -------------- 

printlog(SEPARATOR)
printlog("Remote sensing analysis")

rescale <- function(x, a, b) {
  ((b - a) * (x - min(x, na.rm = TRUE))) / 
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + a
}

srdb %>%
  dplyr::select(Study_midyear, Biome, Leaf_habit, Partition_method, Stage,
                mat_hadcrut4, map_hadcrut4,
                gpp_beer, gpp_modis, 
                SCIA_SIF, GOME2_SIF,
                Rs_annual, Rh_annual) %>%
  rename(`Beer~GPP` = gpp_beer,
         `MODIS~GPP` = gpp_modis) %>%
  # since we're going to divide, can't have SIF crossing 1!
  # rescale to GPP range for convenience 
  mutate(`SCIAMACHY~SIF` = rescale(SCIA_SIF, 0, 3500)) %>%
  mutate(`GOME2~SIF` = rescale(GOME2_SIF, 0, 3500)) %>%
  gather(Flux, fluxvalue, Rs_annual, Rh_annual) %>%
  mutate(Prettyflux = if_else(Flux == "Rs_annual", "R[S]", "R[H]")) %>%
  gather(GPPSIF, gppsifvalue, `Beer~GPP`, `MODIS~GPP`, `SCIAMACHY~SIF`, `GOME2~SIF`) %>%
  mutate(GPPSIF = factor(GPPSIF, levels = c("Beer~GPP", "MODIS~GPP", "SCIAMACHY~SIF", "GOME2~SIF"))) %>%
  filter(!is.na(Leaf_habit), !is.na(gppsifvalue)) ->
  s_gppsif

s_gppsif %>%
  filter(fluxvalue / gppsifvalue <= MAX_FLUX_TO_GPP) ->
  s_gppsif_included
s_gppsif %>%
  filter(fluxvalue / gppsifvalue > MAX_FLUX_TO_GPP) ->
  s_gppsif_excluded

p_gppsif_base <- ggplot(s_gppsif_included, aes(Study_midyear, fluxvalue / gppsifvalue, color = Leaf_habit)) +
  geom_point() +
  facet_grid(GPPSIF ~ Prettyflux, scales = "free", labeller = label_parsed) +
  scale_color_discrete("Leaf habit") +
  xlab("Year") +
  ylab("Respiration:(GPP or SIF)") +
  coord_cartesian(ylim = c(0, 2))
p_gppsif <- p_gppsif_base + 
  geom_smooth(data = subset(s_gppsif_included, Leaf_habit %in% c("Deciduous", "Evergreen")), method = "lm", show.legend = FALSE)
print(p_gppsif )
save_plot("2-gppsif", ptype = ".png")

s_gppsif1 <- subset(s_gppsif_included, GPPSIF %in% c("Beer~GPP", "MODIS~GPP", "SCIAMACHY~SIF"))
p_gppsif1 <- p_gppsif_base %+% s_gppsif1 +
  geom_smooth(data = subset(s_gppsif1, Leaf_habit %in% c("Deciduous", "Evergreen")), method = "lm", show.legend = FALSE)
print(p_gppsif1 )
save_plot("2-gppsif_scia", ptype = ".png")

printlog("Trend tests")
results <- list()
for(dataset in unique(s_gppsif_included$GPPSIF)) {
  for(f in unique(s_gppsif_included$Flux)) {
    d <- filter(s_gppsif_included, Flux == f, GPPSIF == dataset, !is.na(mat_hadcrut4), !is.na(map_hadcrut4))
    dname <- make.names(paste("s", dataset, f, sep = "_"))
    assign(dname, d)
    save_data(d, fname = dname, scriptfolder = FALSE)
    
    printlog("Trend tests for", f, dataset)
    
    # Mann-Kendall
    mk <- MannKendall(d$fluxvalue / d$gppsifvalue)
    print(mk)
    assign(make.names(paste("mk", dataset, f, sep = "_")), mk)
    
    # Linear model
    m <- lm(fluxvalue / gppsifvalue ~ mat_hadcrut4 + map_hadcrut4 + 
              Study_midyear * Leaf_habit + 
              Study_midyear * Stage, 
            data = d)
    m <- stepAIC(m, direction = "both", trace = 0)
    print(summary(m))
    mn <- make.names(paste("m", dataset, f, sep = "_"))
    printlog(mn)
    save_model_diagnostics(m, modelname = mn)
    assign(mn, m)
    
    signif <- anova(m)["Study_midyear", "Pr(>F)"]
    assign(paste0(mn, "_signif"), signif)
    leaf_signif <- anova(m)["Leaf_habit", "Pr(>F)"]
    assign(paste0(mn, "_leaf_signif"), leaf_signif)
    leaf_trend_signif <- anova(m)["Study_midyear:Leaf_habit", "Pr(>F)"]
    assign(paste0(mn, "_leaf_trend_signif"), leaf_trend_signif)
    results[[paste(f, dataset)]] <- tibble(flux = f, 
                                           dataset = make.names(dataset),
                                           n = nrow(d),
                                           `m-k` = pclean(mk$sl),
                                           time_signif = pclean(signif),
                                           leaf_signif = pclean(leaf_trend_signif))
  }
}
rs_results <- bind_rows(results)

# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
