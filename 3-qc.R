# Heterotrophic respiration - quality control
# Ben Bond-Lamberty February 2017
#
# Make a bunch of plots so we can check whether the various datasets
# (Fluxnet, Beer GPP, MODIS GPP, SRDB) agree with each other, or whether
# I've screwed up something up.

source("0-functions.R")

SCRIPTNAME  	<- "3-qc.R"
PROBLEM       <- FALSE


# --------------------- Main --------------------------- 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

read_csv("outputs/srdb_filtered.csv") %>%
  print_dims() ->
  srdb

# Check to see if Rs_annual, Ra_annual, Rh_annual are within 5%

lbls <- subset(srdb, abs(Rs_annual - (Ra_annual+Rh_annual)) > Rs_annual * 0.05)
qplot(Rs_annual, Ra_annual + Rh_annual, data=srdb) + geom_text(data=lbls, aes(label=Record_number), nudge_y = 100)
save_plot("Rs_partitioning")

# Check to see if climate data matches/makes sense
# Note that the FLUXNET data is from the nearest tower,
# which can be a looooong way away
# (That's fine--we filter this to <1 km before analyzing)

srdb %>%
  dplyr::select(Biome, mat_srdb, mat_fluxnet, mat_hadcrut4) %>%
  gather(dataset, value, -Biome, -mat_hadcrut4) %>%
  qplot(mat_hadcrut4, value, color = Biome, data = ., na.rm = TRUE) + 
  facet_grid(dataset ~ .) + geom_abline() + 
  geom_smooth(method = "lm", aes(group = 1), na.rm = TRUE)
save_plot("mat_comparison")

srdb %>%
  dplyr::select(Biome, map_srdb, map_fluxnet, map_hadcrut4) %>%
  gather(dataset, value, -Biome, -map_hadcrut4) %>%
  qplot(map_hadcrut4, value, color = Biome, data = ., na.rm = TRUE) + 
  facet_grid(dataset ~ .) + geom_abline() + 
  geom_smooth(method = "lm", aes(group = 1), na.rm = TRUE)
save_plot("map_comparison")

# Check to see if GPP matches/makes sense

srdb %>%
  dplyr::select(Biome, GOME2_SIF, SCIA_SIF, gpp_fluxnet, gpp_beer, gpp_modis, gpp_srdb) %>%
  gather(dataset, value, -Biome, -gpp_modis) %>%
  qplot(gpp_modis, value, color = Biome, data = ., na.rm = TRUE) + 
  facet_grid(dataset ~ ., scales = "free") + geom_abline() + 
  geom_smooth(method = "lm", aes(group = 1), na.rm = TRUE)
save_plot("gpp_comparison")

srdb %>%
  dplyr::select(Biome, Study_midyear, GOME2_SIF, SCIA_SIF, gpp_fluxnet, gpp_beer, gpp_modis, gpp_srdb) %>%
  gather(dataset, value, -Study_midyear, -Biome) %>%
  qplot(Study_midyear, value, color = Biome, data = ., na.rm = TRUE) + 
  facet_grid(dataset ~ ., scales = "free")
save_plot("gpp_time")

qplot(FLUXNET_DIST, gpp_fluxnet - gpp_modis, data = srdb, log = "x", color = Biome, na.rm = TRUE)
save_plot("gpp_distance")

qplot(RECO_NT_VUT_REF, ER, data = srdb, color = Biome, na.rm = TRUE) +
  geom_abline() + geom_smooth(method = "lm", aes(group = 1), na.rm = TRUE)
save_plot("reco_comparison")

# Climate space distibution

read_csv("outputs/crudata_period.csv.gz") %>%
  print_dims() %>%
  # bin by temp and precip
  mutate(tmp_round = round(tmp / 2, 0) * 2, 
         pre_round = round(pre / 300, 0) * 300) %>%
  group_by(tmp_round, pre_round) %>%
  summarise(area_km2 = sum(area_km2)) ->
  crudata_period

p <- ggplot(crudata_period, aes(tmp_round, pre_round)) + 
  geom_tile(aes(fill = area_km2)) + 
  scale_fill_continuous(low = "lightgrey", high = "black", guide = FALSE) +
  geom_point(data = srdb, aes(mat_hadcrut4, map_hadcrut4, color = Study_midyear), alpha = I(0.5)) + 
  scale_color_continuous(guide = FALSE) +
  xlab(expression(MAT~(degree*C))) + ylab("MAP (mm)")
print(p)
save_plot("climate_space", ptype = ".png")


# How does choice of max fluxnet distance affect our N?

test <- seq(0.1, 10, by = 0.1)  # km
result <- rep(NA, length(test))
for(i in seq_along(test)) {
  filter(srdb, !is.na(Rs_annual), 
         FLUXNET_DIST <= test[i], 
         NEE_VUT_REF_QC >= 0.5) -> x
  result[i] <- nrow(x)
}
df <- tibble(x = test, y = result)
p <- qplot(x, y, data = df, geom = "line", na.rm = TRUE) + 
  xlab("MAX_FLUXNET_DIST") + ylab("N")
print(p)
save_plot("fluxnet_dist_n")


# Latitude/geographic biases over time

srdb %>% 
  qplot(Study_midyear, Latitude, data = ., color = Biome) + 
  xlim(c(1989, 2015)) + 
  geom_smooth(method="lm", color = "blue", group = 1)
save_plot("latitude")

srdb %>% 
  qplot(Study_midyear, Longitude, data = ., color = Biome) + 
  xlim(c(1989, 2015)) + 
  geom_smooth(method="lm", color = "blue", group = 1)
save_plot("longitude")

if(require(maps) & require(mapdata)) {
  world <- map_data("world")
  srdb$long <- srdb$Longitude
  srdb$lat <- srdb$Latitude
  p_base <- ggplot(subset(srdb, Study_midyear >= 1989), aes(x = long, y = lat)) + 
    geom_path(data = world, aes(group = group)) +
    scale_y_continuous(breaks = (-2:2) * 30) +
    scale_x_continuous(breaks = (-4:4) * 45) +
    coord_fixed(xlim = c(-180, 180), ylim = c(-90, 90)) 
  
  print(p_base + geom_point(aes(color = Study_midyear)))
  save_plot("worldmap")
  print(p_base + geom_point(aes(alpha = Study_midyear)) + scale_alpha_continuous(range = c(1, 0)))
  save_plot("worldmap-alpha")
  
}

# Vanessa's request:
# histogram of biome, ecosystem type, and stage, by decade 
# (the combination of all 3, biome + type, and each separately). 
# BUT…I’d like to see the histograms by decade, but then duplicate
# the plots to shift the decadal start point by 5 years.
srdb %>%
  filter(Study_midyear >= 1989) %>%
  dplyr::select(Study_midyear, Biome, Ecosystem_type, Stage, map_hadcrut4, mat_hadcrut4) %>%
  gather(thing, value, -Study_midyear, -Biome, -map_hadcrut4, -mat_hadcrut4) ->
  srdb1
srdb1$yearbin <- as.character(cut(srdb1$Study_midyear, breaks = 3))
srdb2 <- srdb1
srdb2$yearbin <- "All"
srdb_combined <- bind_rows(srdb1, srdb2)
srdb_combined %>%
  qplot(value, data = .) + facet_grid(yearbin ~ thing, scales="free") +
  theme(axis.text.x = element_text(angle = 90))
save_plot("distributions")

srdb_combined %>% 
  qplot(mat_hadcrut4, map_hadcrut4, color=Biome, data=.) + facet_grid(yearbin~.)
save_plot("climate_space_time")

srdb1$yearbin <- as.character(cut(srdb1$Study_midyear, breaks = 4))
srdb_combined <- bind_rows(srdb1, srdb2)
print(last_plot() %+% srdb_combined)
save_plot("distributions-alt")


# A figure for my talk at Stanford, showing growth in SRDB
# between 2010 Nature paper and now
srdb %>%
  filter(Rs_annual < 4000) %>%
  select(Study_midyear, Record_number, Rs_annual, Rh_annual, Biome) %>%
  gather(respiration, value, Rs_annual, Rh_annual) ->
  srdb_stan

p_stan <- ggplot(srdb_stan, aes(Study_midyear, value)) +
  geom_point(aes(color = Record_number > 3400), alpha = 0.5) +
  facet_grid(respiration ~ ., scales = "free_y") +
  scale_color_discrete(guide = FALSE) +
  xlab("Year") + ylab(expression(R~(g~C~m^{-2}~yr^{-1})))
print(p_stan)
save_plot("database_growth")


# Responding to Reviewer 1, see how well remote-sensing indices track tower GPP over time

rescale <- function(x, a, b) {
  ((b - a) * (x - min(x, na.rm = TRUE))) / 
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + a
}

print(SEPARATOR)
read_csv("outputs/fluxnet_remotesensing_comparison.csv") %>%
  arrange(Year) %>%
  print_dims ->
  fluxnet

# Data from Min Chen, 5 June 2017
read_csv("ancillary/fluxsite_SIF_annual-GOME2.csv") %>%
  mutate(sensor = "GOME2") ->
  gome2
read_csv("ancillary/fluxsite_SIF_annual-SCIMACHY.csv") %>%
  mutate(sensor = "SCIMACHY") %>%
  bind_rows(gome2) %>%
  gather(Year, sif_value, -Site, -sensor) %>%
  mutate(sif_value = rescale(sif_value, 0, 3500)) %>%
  spread(sensor, sif_value) %>%
  mutate(Year = as.integer(Year)) ->
  sif

fluxnet <- left_join(fluxnet, sif, by = c("SITE_ID" = "Site", "Year"))
library(Kendall)
print(summary(lm(gpp_modis - GPP_DT_VUT_REF ~ Year, data = fluxnet)))
printlog("Mann-Kendall trend test:")
print(MannKendall(fluxnet$gpp_modis - fluxnet$GPP_DT_VUT_REF))

print(summary(lm(gpp_mte - GPP_DT_VUT_REF ~ Year, data = fluxnet)))
printlog("Mann-Kendall trend test:")
print(MannKendall(fluxnet$gpp_mte - fluxnet$GPP_DT_VUT_REF))

print(summary(lm(GOME2 - GPP_DT_VUT_REF ~ Year, data = fluxnet)))
printlog("Mann-Kendall trend test:")
print(MannKendall(fluxnet$GOME2 - fluxnet$GPP_DT_VUT_REF))

print(summary(lm(SCIMACHY - GPP_DT_VUT_REF ~ Year, data = fluxnet)))
printlog("Mann-Kendall trend test:")
print(MannKendall(fluxnet$SCIMACHY - fluxnet$GPP_DT_VUT_REF))

fluxnet %>%
  mutate(gpp_modis_diff = gpp_modis - GPP_DT_VUT_REF,
         gpp_mte_diff = gpp_mte - GPP_DT_VUT_REF,
         gome2_diff = GOME2 - GPP_DT_VUT_REF,
         scimachy_diff = SCIMACHY - GPP_DT_VUT_REF) %>%
  gather(diffvar, value, gpp_modis_diff, gpp_mte_diff, gome2_diff, scimachy_diff) ->
  fluxnet_plot

p <- ggplot(fluxnet_plot, aes(Year, value, linetype = diffvar != "gpp_mte_diff")) + geom_jitter() + 
  facet_grid(diffvar ~ ., scales = "free_y") + 
  geom_smooth(method = "lm") +
  guides(linetype = FALSE) +
  ylab("Satellite - tower GPP difference (gC)")
print(p)
save_plot("fluxnet_comparison", height = 6, width = 6)


# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
