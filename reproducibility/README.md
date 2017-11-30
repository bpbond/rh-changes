# reproducibility

Ecological and climate change science has to be reproducible and open.

http://dx.doi.org/10.1038/sdata.2017.114
http://dx.doi.org/10.1126/science.1197962
http://dx.doi.org/10.1111/j.1365-2486.2012.02693.x

This folder contains the `srdb_filtered.csv` file that is the output of the `2-prepdata.R` script (which is expensive to run and uses lots of large datasets not included in this repository). This file can be copied into `outputs/srdb_filtered.csv` and used to run the main analysis script, `4-analysis.R`.

This folder also contains the log files produced by every script. These contain the analytical results, but also (and perhaps more importantly) the R package and session information used.
