# rh-changes

Our [2010 study](http://www.nature.com/nature/journal/v464/n7288/full/nature08930.html) concluded by saying yes, we think soil respiration (Rs) is increasing, but we don't know if this is simply an acceleration of the C cycle, or constitutes a climate feedback. (Or both.) A reasonable _a priori_ expectation is that an overall warming climate will result in losses of soil organic carbon (SOC) and increased heterotrophic respiration (Rh), above and beyond any climate-driven increases in GPP and/or the total soil-to-atmosphere Rs flux.  **Can we test this from observations?**

Potential lines of evidence might include:
- Rh/Rs temporal trend - this was the original idea :)
- Rh response to climate anomalies
- Rs and Rh relationships to [FLUXNET (MTE) GPP](http://dx.doi.org/10.1029/2010JG001566), [MODIS GPP](http://dx.doi.org/10.1016/j.rse.2004.12.011), and [SIF](http://dx.doi.org/10.1111/j.1365-2486.2009.01908.x) - increasing would imply losses not fueled by increased GPP

What we did:
- Updated the global [Soil Respiration Database](http://www.biogeosciences.net/7/1915/2010/) with data through 2015
- Matched temporally- and spatially-resolved respiration data with the ancillary datasets
- Used linear models to examine whether Rh, Rh:Rs, Rh:GPP, Rs:GPP, etc., change over a 25-year period
- Examined [FLUXNET](http://fluxnet.fluxdata.org//data/fluxnet2015-dataset/) data for changes in respiration
- Thanks to great referee feedback, added analyses looking at [ISIMIP](http://dx.doi.org/10.1088/1748-9326/12/1/010301) and longitudinal (site-specific) Rh changes, as well as many tests for robustness of the results

To re-run our analysis:
- All scripts used are included in this repository
- Unless you want to rebuild everything from the underlying SRDB, Hadley, MODIS, MTE, etc., datasets, it's simplest to use a pre-processed dataset and start with the main analysis script `4-analysis.R`
- Copy `reproducibility/srdb-filtered.csv` folder to `outputs/srdb-filtered.csv`
- Run the main analysis script. It uses a bunch of R packages, the names and version numbers of which are listed in the script and in `0-functions.R`
- Note that the original script logs, including R session information details, are archived in `reproducibility/`
