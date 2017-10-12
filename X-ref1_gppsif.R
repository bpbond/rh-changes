# Test whether the R:GPP/SIF trends are robust or not (per Reviewer 1)
#
# Ben Bond-Lamberty September 2017

source("0-functions.R")

SCRIPTNAME  	<- "X-ref1_gppsif.R"
PROBLEM       <- FALSE

# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

read_csv(SRDB_GPPSIF_FILE) %>% 
  mutate(year_round = round(Study_midyear, 0)) ->
  s_gppsif

all_data <- list()
results <- list()

# We assume a real-world GPP increase of 0.5%/yr. This is very conservative,
# i.e. it's at the high end of estimates (Anav et al. 2015, Ito et al. 2017)
# and so produces the largest effect to test robustness of R:GPP trends
real_gpp_increase <- 0.5   # %/yr, i.e. 1.005 each year

# Assume satellites are seeing anywhere from ALL of the GPP increase above
# (missing_percent_per_year = 0) to NONE of it (missing_percent_per_year = 100)
for(missing_percent_per_year in seq(0, 100, by = 10)) {
  printlog(SEPARATOR) 
  printlog(missing_percent_per_year)
  
  for(dataset in unique(s_gppsif$GPPSIF)) {   # MODIS, MTE, GPP, SIF
    for(f in unique(s_gppsif$Flux)) {         # Rh or Rs
      d <- filter(s_gppsif, Flux == f, GPPSIF == dataset, !is.na(mat_hadcrut4), !is.na(map_hadcrut4))
      dname <- make.names(paste("gppadjust", dataset, f, missing_percent_per_year))
      
      m <- lm(fluxvalue / gppsifvalue ~ Study_midyear, data = d)
      trend <- coefficients(m)["Study_midyear"]
      p <- anova(m)["Study_midyear", "Pr(>F)"]
      
      # We assume that every year satellites are missing some percentage of GPP increase
      # so compute an adjusted value. For example, if satellites are missing 0.1%/yr, and we
      # have a value of 1000 after 10 years into record, the adjusted value is
      # 1000 * (1.001 ^ 10) = 1010. 
      # (This assumes that every record is 'good' (no error) at beginning, and accumulates
      # error with time. This doesn't quite match reality, in which different records are calibrated
      # to different time periods and then extended, but seems a reasonable simplification.)
      unaccounted_gpp <- real_gpp_increase * missing_percent_per_year/100 # in %/yr
      d$gppsif_adjusted <- d$gppsifvalue * (1 + unaccounted_gpp/100) ^ (d$Study_midyear - min(d$Study_midyear))
      m_adj <- lm(fluxvalue / gppsif_adjusted ~ Study_midyear, data = d)
      trend_adj <- coefficients(m_adj)["Study_midyear"]
      p_adj <- anova(m_adj)["Study_midyear", "Pr(>F)"]
      
      cat(missing_percent_per_year, dataset, f, trend, p, trend_adj, p_adj, "\n")
      d$p_adj <- p_adj
      d$sat_missing <- missing_percent_per_year
      
      results[[dname]] <- tibble(dataset = dataset,
                                 flux = f,
                                 missing_percent_per_year = missing_percent_per_year,
                                 n = length(m_adj$residuals),
                                 trend = trend,
                                 trend_adj = trend_adj,
                                 p = p, p_adj = p_adj)
      # Save for later use
      all_data[[dname]] <- d
    }
  }
}  # end of missing_percent_per_year loop


printlog(SEPARATOR)
results <- bind_rows(results)

# This is the response to reviewers figure, and Extended Data Figure 4
p <- ggplot(results, aes(missing_percent_per_year, p_adj, color = dataset)) +
  geom_line() + facet_grid(flux ~ ., scales = "free") +
  xlab("Percent of GPP increase missed by satellites") + 
  ylab("p-value for trend of R/GPPadjusted") +
  geom_hline(yintercept = 0.05, linetype = 2) + ylim(c(0, 0.1))
print(p)
save_plot("satellite_missing_gpp")

# Also make some diagnostic figures
all_data <- bind_rows(all_data)
all_data$sat_missing_bin <- cut(all_data$sat_missing, breaks = 4)

for(dataset in unique(s_gppsif$GPPSIF)) {
  for(f in unique(s_gppsif$Flux)) {
    d <- filter(all_data, GPPSIF == dataset, Flux== f)
    dname <- make.names(paste("gppadjust", dataset, f))
    
    p <- qplot(Study_midyear, fluxvalue / gppsif_adjusted, data = d) + 
      geom_smooth(method = "lm", color = "red") + facet_wrap(~sat_missing_bin) +
      ggtitle(dname) + ylim(c(0, 2))
    print(p)
    save_plot(dname)
  }
}


printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
