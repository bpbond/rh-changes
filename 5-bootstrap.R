# Heterotrophic respiration - bootstrap analysis script
#
# Ben Bond-Lamberty March 2017

source("0-functions.R")

SCRIPTNAME  	<- "5-bootstrap.R"
PROBLEM       <- FALSE

BOOTSTRAP_N <- 1000
set.seed(12345)

RS_DATASET <- "outputs/s_MODIS.GPP_Rs_annual"
library(broom)  # 0.4.1

# --------------------- Main --------------------------- 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

s_fluxnet <- read_csv("outputs/s_fluxnet.csv")
read_csv(RS_DATASET) %>%
  filter(Flux == "Rs_annual") %>%
  dplyr::select(Study_midyear, fluxvalue, gppsifvalue, Leaf_habit) %>%
  mutate(rs_to_gpp = fluxvalue / gppsifvalue) ->
  s_rs

# --------------------- FLUXNET --------------------------- 

# What fraction of sites exhibit positive trends?
printlog("Computing fraction of FLUXNET sites with positive Rs/GPP trends")
s_fluxnet %>%
  group_by(FLUXNET_SITE_ID) %>%
  do(mod = lm(Rs_annual / gpp_fluxnet ~ Year, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Year") %>%
  mutate(trend = sign(estimate), signif = p.value < 0.05) %>%
  group_by(trend) %>%
  summarise(n = n()) %>%
  mutate(percent = round(n / sum(n) * 100, 0)) %>%
  print ->
  fluxnet_trends

fluxnet_positive_trends <- fluxnet_trends$percent[which(fluxnet_trends$trend == 1)]

# --------------------- Remote sensing bootstrap 1 (fake data) --------------------------- 

# Progressively replace 5%, 10%, ...95% of 'real' MODIS dataset
# with no-change fake data. Compute trends for each group, using
# 1000 trials each.
printlog("Bootstrap 1: progressively replace data with no-trend fake data")
results <- list()
for(fdf in seq(0.05, 0.95, by = 0.1)) {
  printlog("Fake data fraction =", fdf)
  for(i in 1:BOOTSTRAP_N) {
    df <- s_rs
    df$fdf <- fdf
    df$bootstrap <- i
    df$rs_to_gpp[sample(1:nrow(s_rs), size = fdf * nrow(s_rs))] <- median(s_rs$rs_to_gpp)
    results[[paste(fdf, i)]] <- df
  }
}

printlog("Computing trends over all bootstrap samples...")
bind_rows(results) %>%
  group_by(fdf, bootstrap) %>%
  do(mod = lm(rs_to_gpp ~ Study_midyear, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Study_midyear") %>%
  mutate(trend = sign(estimate), signif = p.value < 0.05) ->
  bs1_results

p1 <- ggplot(bs1_results, aes(fdf, p.value, group = fdf)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot() + 
  scale_y_log10() +
  geom_hline(yintercept = 0.05, linetype = 2) + 
  xlab("Fake (no trend) data fraction") +
  ylab("Significance of Rs:GPP trend over time") +
  ggtitle(paste(basename(RS_DATASET), "N =", BOOTSTRAP_N))

print(p1)
save_plot("bootstrap1")

# --------------------- Remote sensing bootstrap 2 (data size) --------------------------- 

# Progressively use smaller subsets of the MODIS data.
# Compute trends for each group, using 1000 trials each.
printlog("Bootstrap 2: progressively smaller dataset size")
results <- list()
for(dfrac in seq(0.05, 0.95, by = 0.1)) {
  printlog("Data fraction =", fdf)
  for(i in 1:BOOTSTRAP_N) {
    s_rs %>%
      mutate(dfrac = dfrac, bootstrap = i) %>%
      sample_n(dfrac * nrow(s_rs)) ->
      results[[paste(dfrac, i)]]
  }
}

printlog("Computing trends over all bootstrap samples...")
bind_rows(results) %>%
  group_by(dfrac, bootstrap) %>%
  do(mod = lm(rs_to_gpp ~ Study_midyear, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Study_midyear") %>%
  mutate(trend = sign(estimate), signif = p.value < 0.05) ->
  bs2_results

p2 <- ggplot(bs2_results, aes(dfrac, p.value, group = dfrac)) +
  geom_jitter(alpha = 0.1) +
  geom_boxplot() + 
  scale_y_log10() +
  geom_hline(yintercept = 0.05, linetype = 2) + 
  geom_vline(xintercept = nrow(s_fluxnet) / nrow(s_rs), linetype = 2) +
  xlab("Fraction of original data used") +
  ylab("Significance of Rs:GPP trend over time") +
  ggtitle(paste(basename(RS_DATASET), "N =", BOOTSTRAP_N))

print(p2)
save_plot("bootstrap2")

# ----------------------- Clean up ------------------------- 

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
