# Heterotrophic respiration - main analysis script
#
# Ben Bond-Lamberty January 2017

source("0-functions.R")

SCRIPTNAME  	<- "4-analysis.R"
PROBLEM       <- FALSE

MAX_FLUXNET_DIST <- 1.0  # km
MIN_NEE_QC <- 0.5

library(broom)  # 0.4.1
library(Kendall) # 2.2

# ==============================================================================
# -------------- Main ------------------- 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

srdb <- read_csv("outputs/srdb_filtered.csv") 

print_dims(srdb)

# -------------- SRDB Rh:Rs analysis ------------------- 

srdb %>%
  filter(Year >= 1989, !is.na(Rs_annual), !is.na(Rh_annual)) %>%
  mutate(Year = as.integer(Year),
         yeargroup = cut(Year, breaks = c(1989, 1994, 1999, 2004, 2009, 2014),
                         labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014"))) %>%
  group_by(yeargroup) %>% 
  mutate(group = paste0(yeargroup, " (N = ", n(), ")")) -> 
  x


p1 <- ggplot(x, aes(Rs_annual, Rh_annual, color = group)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10() +
#  scale_color_grey("Study year") +
  annotation_logticks() +
  xlab(expression(R[S]~(g~C~m^-2~yr^-1))) +
  ylab(expression(R[H]~(g~C~m^-2~yr^-1)))

p1 <- p1 + coord_cartesian(xlim=c(80, 3200), ylim=c(70, 2000))
printlog("NOTE we are plotting this graph with one point cut off:")
printlog(x[which.min(x$Rs_annual), c("Rs_annual", "Rh_annual")])

print(p1)
save_plot("srdb-rh-rs")


# --------------- FLUXNET analysis --------------------- 

printlog(SEPARATOR)
printlog("FLUXNET data analysis")
printlog("Filtering to MAX_FLUXNET_DIST =", MAX_FLUXNET_DIST)
x <- subset(srdb, FLUXNET_DIST <= MAX_FLUXNET_DIST)
print_dims(x)
printlog("Filtering to MIN_NEE_QC =", MIN_NEE_QC)
x <- subset(x, NEE_VUT_REF_QC >= MIN_NEE_QC)
print_dims(x)

x %>%
  mutate(tmp_trend_label = if_else(tmp_trend > 0, "Warming", "Cooling"),
         pre_trend_label = if_else(pre_trend > 0, "Wetter", "Drier")) ->
  x

# There are almost no FLUXNET sites with cooling trends, so we just 
# Looking at precip effect
x %>%
  group_by(pre_trend_label) %>%
  do(mod = lm(Rs_annual/GPP_NT_VUT_REF ~ Year, weights = YearsOfData, data = .)) %>%
  tidy(mod) %>%
  filter(term == "Year") %>%
  print

printlog("Mann-Kendall trend test, drier-trend sites:")
x1 <- filter(x, pre_trend_label == "Drier")
print(MannKendall(x1$Rs_annual / x1$GPP_NT_VUT_REF))
printlog("Mann-Kendall trend test, wetter-trend sites:")
x2 <- filter(x, pre_trend_label == "Wetter")
print(MannKendall(x2$Rs_annual / x2$GPP_NT_VUT_REF))


p <- ggplot(x, aes(Year, Rs_annual / GPP_NT_VUT_REF, group = FLUXNET_SITE_ID)) + 
  geom_point(aes(color = IGBP)) + 
  facet_grid( ~ pre_trend_label, scales = "free") + 
  geom_smooth(method = "lm", color = "grey", fill = NA) + 
  geom_smooth(method = "lm", group = 1) +
  ylab(expression(R[S]:GPP[fluxnet]))

print(p)
save_plot("fluxnet")

printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
