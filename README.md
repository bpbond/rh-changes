# rh-changes

Data and code for analyses in (Bond-Lamberty et al., "[XXX]")[https://doi.org/10.1038/s41586-018-0358-x].

Our [2010 study](http://www.nature.com/nature/journal/v464/n7288/full/nature08930.html) concluded by saying yes, we think soil respiration (Rs) is increasing, but we don't know if this is simply an acceleration of the C cycle, or constitutes a climate feedback. (Or both.) A reasonable _a priori_ expectation is that an overall warming climate will result in losses of soil organic carbon (SOC) and increased heterotrophic respiration (Rh), above and beyond any climate-driven increases in GPP and/or the total soil-to-atmosphere Rs flux.  **Can we test this from observations?**

Potential lines of evidence might include:
- Rh/Rs temporal trend - this was the original idea :)
- Rh response to climate anomalies
- Rs and Rh relationships to [FLUXNET (MTE) GPP](http://dx.doi.org/10.1029/2010JG001566), [MODIS GPP](http://dx.doi.org/10.1016/j.rse.2004.12.011), and [SIF](http://dx.doi.org/10.1111/j.1365-2486.2009.01908.x) - increasing would imply losses not fueled by increased GPP
- Site-specific changes over time, for example at [FLUXNET](http://fluxnet.fluxdata.org) sites

What we did:
- Updated the global [Soil Respiration Database](https://github.com/bpbond/srdb) (see also [publication](http://www.biogeosciences.net/7/1915/2010/)) with data through 2015
- Matched temporally- and spatially-resolved respiration data with the ancillary datasets
- Used linear models to examine whether Rh, Rh:Rs, Rh:GPP, Rs:GPP, etc., change over a 25-year period
- Examined [FLUXNET](http://fluxnet.fluxdata.org//data/fluxnet2015-dataset/) data for changes in respiration
- Thanks to great reviewer feedback, added analyses looking at [ISIMIP](http://dx.doi.org/10.1088/1748-9326/12/1/010301) and longitudinal (site-specific) Rh changes, as well as _many_ tests (h/t Referee 1...) for robustness of the results

To re-run our analysis:
- All scripts used are included in this repository
- Unless you want to rebuild everything from the underlying SRDB, Hadley, MODIS, MTE, etc., datasets (which are not included in this repo, although URLs are given for all), it's simplest to use a pre-processed dataset and start with the main analysis script `4-analysis.R`
- Copy `reproducibility/srdb-filtered.csv` and `reproducibility/fluxnet.csv` into `outputs/` (which you may need to create first)
- Run the main analysis script. It uses a bunch of R packages, the names and version numbers of which are listed in the script and in `0-functions.R`. A quick one-liner to install the necessary packages: `install.packages(c("dplyr", "broom", "Kendall", "MASS", "mblm", "scales”, ”tidyr", "readr", "ggplot2", "luzlogr", "R.utils", "assertthat", "cowplot"))`
- Note that the original script logs, including R session information details, are archived in `reproducibility/`

Things in this root directory:

File/folder | Description
----------- | -------------
0-functions.R | Utility functions and shared settings
1-cci.R | Prep script: process ESA-CCI soil moisture data
1-fluxnet.R | Prep script: process FLUXNET2015 data
2-prepdata.R | Prep script: match SRDB data with all the various ancillary datasets
3-qc.R | QC script: make sure we haven't screwed something up 
**4-analysis.R** | **Main analysis script**
5-bootstrap.R | Sensitivity/test script: look at FLUXNET data
5-ref1_gppsif.R | Sensitivity/test script: examine whether 'missed' satellite data are inducing a false trend
5-ref1_hashimoto.R | Sensitivity/test script: look at probability of 'seeing' a significant change in Rh at the site level
5-ref1_ismip.R | Sensitivity/test script: test the variability of ISIMIP Rh data
ancillary/ | Ancillary data (FLUXNET and ISIMIP)
inputs/ | SRDB and SIF data
LICENSE | License
README.md | This file...
**reproducibility/** | **Archived log files as well as data to run `4-analysis.R`**
rh-changes.Rproj | [RStudio](https://www.rstudio.com) project file
