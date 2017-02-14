# Prep data for heterotrophic respiration analysis
#
# Load SRDB; filter for 'good' data (unmanipulated ecosystems, IRGA/GC only, etc);
# spatially match with CRU climate, Max Planck GPP, MODIS GPP, FLUXNET, and SoilGrids1km datasets
#
# Ben Bond-Lamberty January 2017

source("0-functions.R")

SCRIPTNAME  	<- "2-prepdata.R"
PROBLEM       <- FALSE

APPEND_ONLY <- FALSE

library(raster) # 2.5.8


# -----------------------------------------------------------------------------
# For every SRDB record, 'expand' so there's one year per row
expand_datayears <- function(fd) {
  x <- list()
  for(i in seq_len(nrow(fd))) {
    x[[i]] <- tibble(Record_number = fd$Record_number[i],
                     FLUXNET_SITE_ID = fd$FLUXNET_SITE_ID[i],
                     Longitude = fd$Longitude[i],
                     Latitude = fd$Latitude[i],
                     Year = seq(ceiling(fd$Study_midyear[i] - 0.5 - fd$YearsOfData[i] / 2), 
                                floor(fd$Study_midyear[i] - 0.5 + fd$YearsOfData[i] / 2)))
  }
  bind_rows(x)  
}


# -----------------------------------------------------------------------------
# For every SRDB record, find distance and ID of nearest FLUXNET tower 
match_fluxnet <- function(d, fluxnet) {
  library(fossil)  # 0.3.7
  
  printlog("Starting FLUXNET nearest-neighbor matching...")
  x <- d[c("Longitude", "Latitude")]
  y <- fluxnet[c("LOCATION_LONG", "LOCATION_LAT")]
  names(y) <- names(x)
  z <- rbind(x, y)
  coordinates(z) <- ~ Longitude + Latitude
  dists <- fossil::earth.dist(z, dist = FALSE)
  dists <- dists[(nrow(x)+1):nrow(dists), 1:nrow(x)]
  
  d$FLUXNET_DIST <- NA_real_
  d$FLUXNET_SITE_ID <- NA_character_
  for(i in seq_len(nrow(d))) {
    d$FLUXNET_DIST[i] <- dists[,i][which.min(dists[,i])]
    d$FLUXNET_SITE_ID[i] <- fluxnet$SITE_ID[which.min(dists[,i])]
  }
  d  
}


# -----------------------------------------------------------------------------
# Extract data from a raster brick or raster stack, given vectors of lon/lat/time info
# This is general-purpose and called by both extract_ncdf_data and extract_geotiff_data below
extract_data <- function(rasterstack, varname, lon, lat, midyear, nyears, 
                         file_startyear, file_layers,
                         baseline = c(1961, 1990), 
                         trendline = c(1991, 2010),
                         print_every = 100) {
  
  printlog(SEPARATOR)
  printlog("Starting extraction for varname:", varname)
  
  # Results vectors: variable, variable normal (1961-1990 by default),
  # trend (1991-2010) and trend significance
  x <- normx <- trend <- trend_p <- rep(NA_real_, length(lon))
  
  # Find nearest neighbors for all lon/lat pairs
  for(i in seq_along(lon)) {
    sp <- SpatialPoints(cbind(lon[i], lat[i]))
    
    if(is.null(file_startyear)) {
      start_layer <- nlayers <- 1  # no time dimension
    } else {
      midyear_layer <- (midyear[i] - file_startyear) * 12
      start_layer <- midyear_layer - nyears[i] / 2 * 12
      nlayers <- nyears[i] * 12
    }
    
    printit <- print_every & i %% print_every == 0
    if(printit) {
      printlog("Extracting", i, lon[i], lat[i], midyear[i], nyears[i], "- layers", start_layer, 
               "to", start_layer + nlayers - 1, "...")
    }
    
    # Weirdly, raster::extract does not throw an error if we pass it a negative start layer
    # It just rolls merrily along, returning from the beginning of the file
    if(start_layer < 0 | start_layer + nlayers - 1 > file_layers) {
      printlog(varname, "layer out of bounds #", i)
      next
    }
    
    # Extract the information for this point, over as many years as needed, then average
    rasterstack %>%
      raster::extract(sp, layer = start_layer, nl = nlayers) %>%
      mean(na.rm = TRUE) ->
      x[i]
    
    if(printit) {
      printlog(varname, "value:", x[i])
    }
    
    # Calculate baseline (normal) data
    if(!is.null(baseline)) {
      start_layer <- max(1, (baseline[1] - file_startyear) * 12 + 1)
      nlayers <- (baseline[2] - baseline[1] + 1) * 12
      rasterstack %>%
        raster::extract(sp, layer = start_layer, nl = nlayers) %>%
        mean(na.rm = TRUE) ->
        normx[i]
      if(printit) {
        printlog(varname, "normal:", normx[i])
      }
    }
    
    # Calculate trend
    if(!is.null(trendline)) {
      start_layer <- max(1, (trendline[1] - file_startyear) * 12 + 1)
      nlayers <- (trendline[2] - trendline[1] + 1) * 12
      rasterstack %>%
        raster::extract(sp, layer = start_layer, nl = nlayers) ->
        vals
      tibble(year = rep(trendline[1]:trendline[2], each = 12), 
             x = as.numeric(vals)) %>%
        group_by(year) %>%
        summarise(x = mean(x, na.rm = TRUE)) ->
        vals
      
      try({
        m <- lm(x ~ year, data = vals)
        trend[i] <- m$coefficients[2]
        trend_p[i] <- summary(m)$coefficients[2, 4]
      }, silent = TRUE)
      
      if(printit) {
        printlog(varname, "trend:", trend[i], "p =", trend_p[i])
      }
    }
  }
  
  # Assemble output data set and return  
  out <- tibble(x = x)
  names(out) <- c(varname)
  
  if(!is.null(baseline)) {
    out <- bind_cols(out, tibble(normx = normx))
    names(out)[2] <- paste0(varname, "_norm")
  }
  if(!is.null(trendline)) {
    out <- bind_cols(out, tibble(trend = trend, trend_p = trend_p))
    names(out)[3:4] <- paste0(varname, c("_trend", "_trend_p"))
  } 
  out
}

# -----------------------------------------------------------------------------
# Extract MODIS NPP data given vectors of lon/lat/time info
extract_geotiff_data <- function(directory, varname, lon, lat, midyear, nyears, file_startyear,
                                 print_every = 100) {
  
  # Decompress if necessary
  zipfiles <- list.files(directory, pattern = "*.tif.gz$", full.names = TRUE)
  for(f in zipfiles) {
    printlog("Decompressing", basename(f))
    R.utils::gunzip(f, remove = FALSE, overwrite = TRUE)
  }
  
  files <- list.files(directory, pattern = "*.tif$", full.names = TRUE)
  printlog("Creating raster stack from", length(files), "files in", directory)
  
  # NB: extracting data from a a stack is **much** slower than from a brick(), as
  # used below. I don't anticipate running this program often enough to care, but FYI.
  nc <- stack(as.list(files))
  
  out <- extract_data(nc, varname, lon, lat, midyear, nyears, 
                      file_startyear = file_startyear, file_layers = length(files), 
                      baseline = NULL, trendline = NULL, print_every)
  
  # Clean up if we decompressed anything
  for(f in zipfiles) {
    printlog("Removing", f)
    file.remove(gsub(".gz$", "", f))  # remove the unzipped file
  }
  
  out
}

# -----------------------------------------------------------------------------
# Extract CRU and Max Planck data given vectors of lon/lat/time info
extract_ncdf_data <- function(filename, lon, lat, midyear, nyears, file_startyear,
                              baseline = c(1961, 1990),
                              trendline = c(1991, 2010),
                              print_every = 100) {
  
  assert_that(length(lon) == length(lat))
  assert_that(length(lon) == length(midyear))
  assert_that(length(midyear) == length(nyears))
  
  # Load file, decompressing first if necessary
  compressed <- grepl("gz$", filename)
  if(compressed) {
    printlog("Decompressing", filename)
    ncfile <- R.utils::gunzip(filename, remove = FALSE, overwrite = TRUE)
  } else {
    ncfile <- filename
  }
  nc <- brick(ncfile)
  varname <- attributes(nc@data)$zvar
  
  out <- extract_data(nc, varname, lon, lat, midyear, nyears, 
                      file_startyear = file_startyear, file_layers = nc@data@nlayers, 
                      baseline, trendline, print_every)
  
  # Clean up
  if(compressed) {
    printlog("Removing", ncfile)
    file.remove(ncfile)
  }
  out
}


# ==============================================================================
# Main 

openlog(file.path(outputdir(), paste0(SCRIPTNAME, ".log.txt")), sink = TRUE)
printlog("Welcome to", SCRIPTNAME)
all_data <- list()

# -------------- 1. Get SRDB data and filter ------------------- 

srdb <- read_csv("inputs/srdb-data.csv", col_types = "dcicicccccdddddccddccccccccddcdddddcddcddddididdddddddddddcccccddddddddcddddddcdcddddddddddddddddddddddc")
print_dims(srdb)

printlog("Filtering...")
srdb %>%
  filter(!is.na(Longitude), !is.na(Latitude), 
         !is.na(Rs_annual) | !is.na(Rh_annual),
         !is.na(Study_midyear), !is.na(YearsOfData),
         is.na(Duplicate_record),
         Ecosystem_state != "Managed", 
         Manipulation == "None",
         Meas_method %in% c("IRGA", "Gas chromatography")) %>%
  dplyr::select(Record_number, Quality_flag, Study_midyear, 
                YearsOfData, Longitude, Latitude, 
                Biome, Ecosystem_type, Leaf_habit, Stage, Soil_drainage,
                MAT, MAP, Study_temp, Study_precip, Partition_method,
                Rs_annual, Rh_annual, Ra_annual, RC_annual, GPP, ER) ->
  srdb
print_dims(srdb)

stopifnot(!any(duplicated(srdb$Record_number)))

if(file.exists(SRDB_FILTERED_FILE)) {
  old_data <- read_csv(SRDB_FILTERED_FILE)
  if(APPEND_ONLY) {
    printlog("Filtering pre-calculated data...")
    srdb <- subset(srdb, !(Record_number %in% old_data$Record_number))
  }
} else {
  old_data <- tibble()
}

if(!nrow(srdb)) {
  closelog()
  stop("No rows of data--nothing to do!")
}


# -------------- 2. SIF ------------------- 

printlog("Joining with SIF data...")
read_csv("inputs/SIF.csv", col_types = "iddidd") %>%
  dplyr::select(Record_number, GOME2_SIF, SCIA_SIF) %>%
  group_by(Record_number) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  right_join(srdb, by = "Record_number") ->
  srdb

stopifnot(!any(duplicated(srdb$Record_number)))

# -------------- 3. FLUXNET ------------------- 

# Start by finding the nearest Fluxnet station, and its distance in km
fluxnet <- read_csv("outputs/fluxnet.csv", col_types = "iddddddccdddcdi")
srdb <- match_fluxnet(srdb, fluxnet)

# Expand the srdb data so that we have an entry for every integer year;
# merge with the Fluxnet data; and put back together
printlog("Building merge data by expanding SRDB years...")
srdb %>%
  dplyr::select(Record_number, Study_midyear, YearsOfData, Longitude, Latitude, FLUXNET_SITE_ID) %>%
  expand_datayears ->
  srdb_expanded

save_data(srdb_expanded)

printlog("Computing FLUXNET means as necessary and merging back in...")
srdb_expanded %>%
  dplyr::select(-Longitude, -Latitude) %>%
  left_join(fluxnet, by = c("FLUXNET_SITE_ID" = "SITE_ID", "Year" = "TIMESTAMP")) %>%
  dplyr::select(-SITE_NAME) %>%
  rename(mat_fluxnet = MAT, map_fluxnet = MAP) %>%
  group_by(FLUXNET_SITE_ID, Record_number, IGBP) %>%
  summarise_all(mean) %>%
  right_join(srdb, by = c("Record_number", "FLUXNET_SITE_ID")) ->
  srdb

printlog("Checking for ecosystem type match between FLUXNET and SRDB")
fem <- rep(FALSE, nrow(srdb))  # 'fluxnet ecosystem match'
for(i in seq_len(nrow(srdb))) {
  igbp <- srdb$IGBP[i]
  et <- srdb$Ecosystem_type[i]
  lh <- srdb$Leaf_habit[i]
  
  if(is.na(igbp)) {
    fem[i] <- FALSE
  } else if(igbp == "CRO") { 
    fem[i] <- et == "Agriculture"
  } else if(igbp %in% c("CSH", "OSH")) {
    fem[i] <- et == "Shrubland"
  } else if(igbp %in% c("DBF", "DNF")) {
    fem[i] <- et == "Forest" & lh %in% c("Deciduous", "Mixed")
  } else if(igbp %in% c("EBF", "ENF")) {
    fem[i] <- et == "Forest" & lh %in% c("Evergreen", "Mixed")
  } else if(igbp == "GRA") {
    fem[i] <- et == "Grassland"
  } else if (igbp == "MF") {
    fem[i] <- et == "Forest"
  } else if (igbp %in% c("SAV", "WSA")) {
    fem[i] <- et == "Savanna"
  } else if (igbp %in% "WET") {
    fem[i] <- et == "Wetland"
  } else {
    stop("Don't know ", igbp)
  }
}
fem[is.na(fem)] <- FALSE
srdb$FLUXNET_ECOSYSTEM_MATCH <- fem

all_data[["srdb"]] <- srdb


# -------------- 3. Match with CRU climate data ------------------- 

fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.tmp.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/tmp/cru_ts3.24.1901.2015.tmp.dat.nc.gz
all_data[["tmp"]] <- extract_ncdf_data(fn, srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = 1901)
fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.pre.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/pre/cru_ts3.24.1901.2015.pre.dat.nc.gz
pre <- extract_ncdf_data(fn, srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = 1901)
all_data[["pre"]] <- pre * 10 # to mm
fn <- "/Users/d3x290/Data/CRU/cru_ts3.24.1901.2015.pet.dat.nc.gz"
# Downloaded 5 Jan 2017 from https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_3.24/cruts.1609301803.v3.24/pet/cru_ts3.24.1901.2015.pet.dat.nc.gz
all_data[["pet"]] <- extract_ncdf_data(fn, srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = 1901)


# -------------- 4. Match with Max Planck GPP data ------------------- 

fn <- "/Users/d3x290/Data/MaxPlanck/201715151429EnsembleGPP_GL.nc.gz"
# Downloaded 5 Jan 2017 from https://www.bgc-jena.mpg.de/geodb/tmpdnld/201715151429EnsembleGPP_GL.nc
# See https://www.bgc-jena.mpg.de/bgi/index.php/Services/Overview
gpp <- extract_ncdf_data(fn, srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, baseline = NULL, trendline = NULL, file_startyear = 1982)
all_data[["gpp"]] <- gpp * 1000 * 60 * 60 * 24 * 365  # Convert from kgC/m2/s to gC/m2/yr


# -------------- 5. Match with MODIS GPP data ------------------- 

dir <- "/Users/d3x290/Data/MODIS_GPP/"
# Downloaded 6 Jan 2017 from http://www.ntsg.umt.edu/project/mod17
modisgpp <- extract_geotiff_data(dir, "modisgpp", srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = 2000)
modisgpp <- modisgpp * 0.1 # scale factor, per README file; results in gC/m2
modisgpp <- modisgpp * 12 # from mean monthly value to annual sum
# There are some crazy (>10,000 gC/m2) values in MODIS GPP. Remove those
modisgpp$modisgpp[modisgpp$modisgpp > 10000] <- NA
all_data[["modisgpp"]] <- modisgpp


# -------------- 6. Match with SoilGrids1km data ------------------- 

# Downloaded 9 Jan 2017 from ftp://ftp.soilgrids.org/data/archive/12.Apr.2014/
dir <- "/Users/d3x290/Data/soilgrids1km/BLD/"
bd <- extract_geotiff_data(dir, "BD", srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = NULL)
dir <- "/Users/d3x290/Data/soilgrids1km/ORCDRC/"
orc <- extract_geotiff_data(dir, "ORC", srdb$Longitude, srdb$Latitude, srdb$Study_midyear, srdb$YearsOfData, file_startyear = NULL)

all_data[["soc"]] <- tibble(SOC = bd$BD * orc$ORC / 1000)  # kg C in top 1 m

modisgpp <- tibble(modisgpp = 1:nrow(srdb))


# -------------- Done!  ------------------- 

# Combine the various spatial data with the SRDB data and save
bind_cols(all_data) %>%
  rename(gpp_modis = modisgpp, 
         gpp_beer = gpp, 
         gpp_srdb = GPP,
         gpp_fluxnet = GPP_NT_VUT_REF,
         tmp_hadcrut4 = tmp,
         pre_hadcrut4 = pre,
         mat_hadcrut4 = tmp_norm,
         map_hadcrut4 = pre_norm,
         mat_srdb = MAT,
         map_srdb = MAP) %>%
  bind_rows(old_data) ->
  srdb_filtered


save_data(srdb_filtered, scriptfolder = FALSE, fname = basename(SRDB_FILTERED_FILE))


printlog("All done with", SCRIPTNAME)
closelog()

if(PROBLEM) warning("There was a problem - see log")
