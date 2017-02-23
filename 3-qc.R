# Heterotrophic respiration - quality control
# Make a bunch of plots so we can check whether the various datasets
# (Fluxnet, Beer GPP, MODIS GPP, SRDB) agree with each other
# Ben Bond-Lamberty February 2017

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

qplot(mat_hadcrut4, map_hadcrut4, data = srdb, color = Biome)
save_plot("climate_space")

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



# How choice of max fluxnet distance affect our N?

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


# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
