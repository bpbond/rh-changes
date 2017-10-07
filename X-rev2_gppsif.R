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

yrs <- 1990:2014

all_data <- list()
results <- list()
for(missing_gpp_per_year in seq(0, 30, by = 2)) {
  printlog(SEPARATOR) 
  
  printlog("Trend tests")
  for(dataset in unique(s_gppsif$GPPSIF)) {   # MODIS, MTE, GPP, SIF
    for(f in unique(s_gppsif$Flux)) {         # Rh or Rs
      d <- filter(s_gppsif, Flux == f, GPPSIF == dataset, !is.na(mat_hadcrut4), !is.na(map_hadcrut4))
      dname <- make.names(paste("gppadjust", dataset, f, missing_gpp_per_year))
      
      trend <- coefficients(lm(gppsifvalue ~ Study_midyear, data = d))["Study_midyear"]      # Compute an adjusted GPP, which is GPP + years * adjustment/yr
      d$gppsif_adjusted <- d$gppsifvalue + (d$Study_midyear - min(d$Study_midyear)) * missing_gpp_per_year
      d$missing_gpp_per_year <- missing_gpp_per_year
      
      # Linear model
      m <- lm(fluxvalue / gppsif_adjusted ~ Study_midyear, data = d)
      #      m <- stepAIC(m, direction = "both", trace = 0)
      print(summary(m))
      p <- anova(m)["Study_midyear", "Pr(>F)"]
      sat_missing <- missing_gpp_per_year / (missing_gpp_per_year + abs(trend)) * 100
      d$p <- p
      d$sat_missing <- sat_missing
      results[[dname]] <- tibble(dataset = dataset,
                                 flux = f,
                                 missing_gpp_per_year = missing_gpp_per_year,
                                 n = length(m$residuals),
                                 trend = trend,
                                 newtrend = m$coefficients["Study_midyear"],
                                 sat_missing = sat_missing,
                                 p = p)
      # Save for later use
      all_data[[dname]] <- d
    }
  }
  
}  # end of rs_adjustment_25yr loop


printlog(SEPARATOR)
results <- bind_rows(results)

results$signed_p_value <- results$p * sign(results$trend)

p <- ggplot(results, aes(sat_missing, p, color = dataset)) +
  geom_line() + facet_grid(flux ~ ., scales = "free") +
  xlab("Percent of GPP increase missed by satellites") + 
  ylab("p-value for trend of R/GPPadjusted") +
  geom_hline(yintercept = 0.05, linetype = 2) +
  coord_cartesian(xlim = c(0, 50)) # there's no way satellites are missing more than 50%!
print(p)
save_plot("satellite_missing_gpp")

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

