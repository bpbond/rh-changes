# rh-changes

The [2010 Nature paper](http://www.nature.com/nature/journal/v464/n7288/full/nature08930.html) concluded by saying yes, we think Rs is increasing, but we don't know if this is an acceleration of the C cycle, or a climate feedback (or both). Our _a priori_ expectation is that an overall warming climate will result in  **Can we test this from observations?**

Potential lines of evidence might include:
- Rh/Rs temporal trend
- Rh response to climate anomalies
- Rs relationship to Fluxnet GPP
- Rs relationship to MODIS and Beer GPP
- Rs relationship to SIF
- Rs rises over time

What we did:
- Updated the global [Soil Respiration Database](http://www.biogeosciences.net/7/1915/2010/) with data through 2015
- Matched temporally- and spatially-resolved respiration data with the ancillary datasets
- Used linear models to examine whether Rh, Rs, and Rh:Rs change over a 25-year period

To re-run our analysis:
- All scripts used are included in this repository
- Unless you want to rebuild everything from the underlying SRDB, Hadley, MODIS, MTE, etc., datasets, it's simplest to use a pre-processed dataset and start with the main analysis script `4-analysis.R`
- Copy the `srdb-filtered.csv` file from the `reproducibility/` folder to `outputs/` (you may need to create the latter if no scripts have been run yet)
- Run the analysis script. Note that the original script logs, including R session information details, are archived in `reproducibility/`
