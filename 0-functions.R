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


# -----------------------------------------------------------------------------
# Print dimensions of data frame
print_dims <- function(d, dname = deparse(substitute(d))) {
  stopifnot(is.data.frame(d) | is.matrix(d))
  printlog(dname, "rows =", nrow(d), "cols =", ncol(d))
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
  invisible(readr::read_csv(fqfn, ...))
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

if(!file.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR)
}

# if(require(checkpoint) & exists("CHECKPOINTDATE"))
#   try(checkpoint(CHECKPOINTDATE)) # 'try' b/c errors w/o network (issue #171)
library(ggplot2)
theme_set(theme_bw())
library(luzlogr)  # 0.2.0
library(readr)
library(R.utils)
library(assertthat)
library(tidyr)
library(dplyr)
