# Test: how many grid cells show positive, and SIGNIFICANT positive, RH trends?
# This is to address Reviewer 1's concern that our site- and FLUXNET-specific
# show inconsistent results that we don't fully explain.
#
# Method: load Hashimoto (2015) dataset; apply a linear model to the 1990-2012 data 
# for each global grid cell; extract slope and significance information and report them.
#
# Ben Bond-Lamberty August 2017

source("0-functions.R")

SCRIPTNAME  	<- "X-rev2_hashimoto.R"
PROBLEM       <- FALSE


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)

library(ncdf4)
# Downloaded August 25, 2017 from http://cse.ffpri.affrc.go.jp/shojih/data/index.html
ncfiles <- c("~/Data/Hashimoto/RH_yr_Hashimoto2015.nc",
             "~/Data/Hashimoto/RS_yr_Hashimoto2015.nc")

# Fit a linear model (i.e. over time) to each grid cell
f <- function(rh) { 
  df <- data.frame(x = seq_along(rh), y = rh)
  m <- tryCatch(lm(y ~ x, data = df), 
                error = function(e) NA)
}

do_fitting <- function(co2) {
  printlog("Fitting linear model to each grid cell (this is slow)...")
  mods <- apply(co2, c(1, 2), FUN = f)  # slow
  
  printlog("Extracting slopes...")
  slopes <- apply(mods, c(1, 2), FUN = function(x) 
    if(!is.na(x)) x[[1]]$coefficients[["x"]] else NA)
  slopes <- matrix(slopes, nrow = nrow(mods), ncol = ncol(mods))
  
  printlog("Extracting slope p-values...")
  signif <- apply(mods, c(1, 2), FUN = function(x) 
    if(!is.na(x)) summary(x[[1]])$coefficients["x", "Pr(>|t|)"] else NA)
  signif <- matrix(signif, nrow = nrow(mods), ncol = ncol(mods))
  
  image(signif)
  hist(signif)
  
  ncells <- sum(!is.na(slopes))
  printlog("Total cells =", ncells)
  pos_slope <- sum(slopes > 0, na.rm = TRUE)
  printlog("Cells with positive slope =", pos_slope,
           "or", round(pos_slope / ncells * 100, 0), "%")
  signif_pos_slope <- sum(slopes > 0 & signif < 0.05, na.rm = TRUE)
  printlog("Cells with significant positive slope =", signif_pos_slope, 
           "or", round(signif_pos_slope / ncells * 100, 0), "%")
}

for(ncfile in ncfiles) {
  printlog(SEPARATOR)
  printlog("Reading", ncfile)
  nc <- nc_open(ncfile)
  
  #co2 <- ncvar_get(nc, "co2", start = c(199, 299, 1, 90), count = c(3, 3, 1, 22))
  #co2[2,2,] <- NA # punch a hole for testing
  
  # These annual data start in 1901; extract 1990-2012
  co2 <- ncvar_get(nc, "co2", start = c(1, 1, 1, 90), count = c(-1, -1, 1, 23))
  printlog("Dimensions =", paste(dim(co2), collapse = " "))
  
  do_fitting(co2)
  
  printlog(SEPARATOR)
  printlog("Fitting for last 10 years only...")
  do_fitting(co2[,,14:23])
}


# MTE GPP dataset
# These monthly data start in 1982; extract 1990-2011
printlog(SEPARATOR)
ncfile <- "~/Data/MaxPlanck/201715151429EnsembleGPP_GL.nc"
printlog("Reading", ncfile)
nc <- nc_open(ncfile)
gpp <- ncvar_get(nc, "gpp", start = c(1, 1, 96), count = c(-1, -1, -1))
do_fitting(gpp)

coprintlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")

