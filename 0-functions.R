# Template for R analysis script
# Ben Bond-Lamberty September 2015

# This is my starting point for most analysis or data processing scripts
# Most critically, it provides lightweight logging, with sessionInfo() 
# written at the bottom of every log; easy ggplot and data saving; 
# logged csv[.gz|zip] read/write; and a few other handy things.

SCRIPTNAME		<- "Rscript.R"
PROBLEM			<- FALSE
OUTPUT_DIR		<- "outputs/"
SEPARATOR		<- "-------------------"

SRDB_FILTERED_FILE <- file.path(OUTPUT_DIR, "srdb_filtered.csv")
SRDB_GPPSIF_FILE <- file.path(OUTPUT_DIR, "gpp_sif_filtered.csv")

# -----------------------------------------------------------------------------
# Print dimensions of data frame
print_dims <- function(d, dname = deparse(substitute(d))) {
  stopifnot(is.data.frame(d) | is.matrix(d))
  printlog(dname, "rows =", nrow(d), "cols =", ncol(d))
  invisible(d)
} # print_dims

# -----------------------------------------------------------------------------
# Return matrix of memory consumption
object_sizes <- function() {
  rev(sort(sapply(ls(envir = .GlobalEnv), function(object.name) 
    object.size(get(object.name)))))
} # object_sizes

# -----------------------------------------------------------------------------
# Return output directory (perhaps inside a script-specific folder)
# If caller specifies `scriptfolder=FALSE`, return OUTPUT_DIR
# If caller specifies `scriptfolder=TRUE` (default), return OUTPUT_DIR/SCRIPTNAME
outputdir <- function(scriptfolder = TRUE) {
  output_dir <- OUTPUT_DIR
  if(scriptfolder) output_dir <- file.path(output_dir, sub(".R$", "", SCRIPTNAME))
  if(!file.exists(output_dir)) dir.create(output_dir)
  output_dir
} # outputdir

# -----------------------------------------------------------------------------
# Save a ggplot figure
save_plot <- function(pname, p = last_plot(), ptype = ".pdf", scriptfolder = TRUE, ...) {
  fn <- file.path(outputdir(scriptfolder), paste0(pname, ptype))
  printlog("Saving", fn)
  ggsave(fn, p, ...)
} # save_plot

# -----------------------------------------------------------------------------
# Save a data frame
save_data <- function(df, fname = paste0(deparse(substitute(df)), ".csv"), scriptfolder = TRUE, gzip = FALSE, ...) {
  fn <- file.path(outputdir(scriptfolder), fname)
  printlog("Saving", fn)   
  write_csv(df, fn, ...)
  if(gzip & require(R.utils)) {
    R.utils::gzip(fn, overwrite = TRUE)
  }
} # save_data

# -----------------------------------------------------------------------------
# Open a netCDF file and return handle (using ncdf4 package)
open_ncdf <- function(fn, datadir = ".") {
  if(is.null(datadir)) {  # NULL signifies absolute path
    fqfn <- fn 
  } else {
    fqfn <- file.path(datadir, fn)      
  }
  printlog("Opening", fqfn)
  nc_open(fqfn)
} # open_ncdf

# -----------------------------------------------------------------------------
# Open a (possibly compressed) csv file and return data
read_csv <- function(fn, datadir = ".", ...) {
  if(is.null(datadir)) {  # NULL signifies absolute path
    fqfn <- fn 
  } else {
    fqfn <- file.path(datadir, fn)      
  }
  printlog("Opening", fqfn)
  invisible(readr::read_csv(fqfn, progress = FALSE, ...))
} # read_csv

# -----------------------------------------------------------------------------
# Read data from the clipboard
paste_data <- function(header=TRUE) {
  read.table(pipe("pbpaste"), header = header)
} # paste_data

# -----------------------------------------------------------------------------
is_outlier <- function(x, devs = 3.2) {
  # See: Davies, P.L. and Gather, U. (1993).
  # "The identification of multiple outliers" (with discussion)
  # J. Amer. Statist. Assoc., 88, 782-801.
  
  x <- na.omit(x)
  lims <- median(x) + c(-1, 1) * devs * mad(x, constant = 1)
  x < lims[ 1 ] | x > lims[2]
} # is_outlier

# -----------------------------------------------------------------------------
# 'Pretty n' function to round a numeric value and print that # of digits
pn <- function(x, n) {
  formatC(round(unlist(x), n), digits = n, format = "f")
} # pn

# -----------------------------------------------------------------------------
# 'Clean p value' function to pretty-print p value(s), specifically
pclean <- function(x, digits = 3, printP = TRUE) {
  x <- as.vector(x)
  ltstring <- paste0("< 0.", paste(rep("0", digits - 1), collapse = ""), "1")
  valstring <- ifelse(x < 10 ^ -digits, ltstring, pn(x, digits))
  if(printP) {
    paste("P", ifelse(x < 10 ^ -digits, valstring, paste("=", valstring)))
  } else {
    valstring
  }
} # pclean

# -----------------------------------------------------------------------------
# Rescale vector x to the range [a,b]
rescale <- function(x, a, b) {
  ((b - a) * (x - min(x, na.rm = TRUE))) / 
    (max(x, na.rm = TRUE) - min(x, na.rm = TRUE)) + a
}


if(!file.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

library(dplyr)    # 0.7.4
library(tidyr)    # 0.7.1
library(readr)    # 1.1.1
library(ggplot2)  # 2.2.2
theme_set(theme_bw() + theme(panel.grid.minor = element_blank()))
library(luzlogr)  # 0.2.0
library(R.utils)
library(assertthat)
