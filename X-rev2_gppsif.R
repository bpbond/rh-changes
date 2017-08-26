# Test: 
#
# Method: 
#
# Ben Bond-Lamberty August 2017

source("0-functions.R")

SCRIPTNAME  	<- "X-rev2_gppsif.R"
PROBLEM       <- FALSE

library(MASS) # 7.3.45

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

read_csv(SRDB_GPPSIF_FILE) %>% 
  mutate(year_round = round(Study_midyear, 0)) ->
  s_gppsif

# First, construct our adjustment table, that sets how much we'll adjust
# the observed GPP values by

yrs <- 1990:2014

rs_adjustment_25yr <- 0.1  # remote sensing missing 10% (half)

results <- list()
for(rs_adjustment_25yr in seq(0, 0.9, by = 0.1)) {
  printlog(SEPARATOR)
  
  base <- tibble(year = yrs, gpp = rep(1, length.out = length(yrs)))
  
  base %>%
    mutate(GPPSIF = "GPP[MODIS]",
           adjust = seq(1, 1 / (1 - rs_adjustment_25yr), length.out = length(yrs)),
           gpp_adj = gpp * adjust) ->
    modis
  
  base %>%
    mutate(GPPSIF = "GPP[MTE]",
           adjust = seq(1 / (1 - rs_adjustment_25yr), 1, length.out = length(yrs)),
           gpp_adj = gpp * adjust) ->
    mte
  
  base %>%
    mutate(GPPSIF = "SIF[SCIAMACHY]",
           adjust = seq(1 / (1 + (rs_adjustment_25yr / 2)), 1 / (1 - (rs_adjustment_25yr / 2)), length.out = length(yrs)),
           gpp_adj = gpp * adjust) %>% 
    filter(year %in% 1996:2005) %>%
    mutate(year = 2005:2014) ->
    sif
  
  bind_rows(modis, mte, sif, mutate(sif, GPPSIF = "SIF[GOME2]")) %>%
    gather(which, value, gpp, gpp_adj) ->
    adjustment_data
  
  p <- ggplot(adjustment_data, aes(year, value, linetype = which)) + geom_line() + facet_wrap(~GPPSIF) +
    ggtitle(paste0("Assuming r.s. misses ", rs_adjustment_25yr * 100, "% of gpp increase over 25 years")) +
    scale_linetype_discrete(labels = c("Observed by satellite", "'Real' (adjusted)"))
  print(p)
  save_plot(paste("gpp_adjustment", rs_adjustment_25yr, sep = "_"))
  
  adjustment_data %>% 
    dplyr::select(year, GPPSIF, adjust) %>% 
    right_join(s_gppsif, by = c("GPPSIF", "year" = "year_round")) ->
    s_gppsif_adj
  
  
  printlog("Trend tests")
  for(dataset in unique(s_gppsif_adj$GPPSIF)) {
    for(f in unique(s_gppsif_adj$Flux)) {
      d <- filter(s_gppsif_adj, Flux == f, GPPSIF == dataset, !is.na(mat_hadcrut4), !is.na(map_hadcrut4))
      dname <- make.names(paste(dataset, f, rs_adjustment_25yr))
      
      # The 'adjust' values apply to the *increase*, not the absolute GPP
      # Compute trend
      trend <- coefficients(lm(gppsifvalue ~ Study_midyear, data = d))["Study_midyear"]
      # Compute an adjusted GPP, which is trend * years * adjustment added to original GPP
      d$gppsif_adjusted <- d$gppsifvalue + (d$Study_midyear - min(d$Study_midyear)) * trend * d$adjust
      
      # Linear model
      m <- lm(fluxvalue / gppsif_adjusted ~ Study_midyear, data = d)
#      m <- stepAIC(m, direction = "both", trace = 0)
      print(summary(m))
      
      results[[dname]] <- tibble(dataset = dataset,
                                 flux = f,
                                 rs_adjustment_25yr = rs_adjustment_25yr,
                                 n = length(m$residuals),
                                 trend = m$coefficients["Study_midyear"],
                                 p = anova(m)["Study_midyear", "Pr(>F)"])
    }
  }
  
}  # end of rs_adjustment_25yr loop

printlog(SEPARATOR)
results <- bind_rows(results)

p <- ggplot(filter(results, trend >= 0), aes(rs_adjustment_25yr, p, color = dataset)) +
  geom_line() + facet_grid(flux ~ ., scales = "free") +
  xlab("Fraction of GPP missed by satellites") + 
  ylab("p-value for trend of R/GPPadjusted") +
  geom_hline(yintercept = 0.05, linetype = 2)
print(p)
save_plot("satellite_missing_gpp")

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")

