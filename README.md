# rh-changes

The 2010 Nature paper concluded by saying yeah, Rs is increasing, but we don't know if this is an acceleration of the C cycle, or a feedback.

With eight (if I enter 2014-2016) more years' data, **can we test this**?

1. We think we know, roughly, how much Rh should appear at a given Rs level (2004 GCB paper). Has this ratio changed (increased) over time? Do we have enough data to test 1990s, 2000s, 2010s?
	- naively graphing Rh/Rs over time = increasing trend
	- try fitting same log-log model as 2004 paper by decade or 5 years
	
2. Rs is correlated with GPP (2010 Biogeosci paper). Has this changed/increased over time? One test would be the GPP numbers in the database; another would be to use a temporal GPP product: Max Planck (https://www.bgc-jena.mpg.de/bgi/index.php/Services/Overview) ("gpp") or MODIS (http://www.ntsg.umt.edu/project/mod17).
	- naively graphing Rs/GPP over time = flat
	- naively graphing Rh/GPP over time = negative
	- naively graphing Rh/gpp over time = increasing trend
	- naively graphing Rh/gpp over time = increasing trend

3. Following the 2010 methodology, test climate effects on Rh (N=206 in 2010; presumably 50% more now?).
  - quick linear model shows temporal trend after tmp_norm, pre_norm, pet_norm, gpp

For this will need climate data (CRU); GPP data; updated SRDB.
SRDB - prioritize ones with RH
