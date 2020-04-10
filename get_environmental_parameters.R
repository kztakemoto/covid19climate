args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
    stop("dataset type needed (global or US)", call.=FALSE)
} else if (length(args)==1) {
    if(args[1] != "global" && args[1] != "US"){
        stop("invalid dataset type", call.=FALSE)
    } else {
        dataset_type <- args[1]
    }
} else {
    stop("too many auguments", call.=FALSE)
}

# load libraries
library(raster)
library(incidence)

# path to COVID-19 directory (https://github.com/CSSEGISandData/COVID-19)
dir <- "./"
if(dataset_type == "global"){
    # first index for count data
    first_idx_count <- 5
    # generate filename
    filename <- paste(dir, "COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv",sep="")
} else {
    # first index for count data
    first_idx_count <- 12
    # generate filename
    filename <- paste(dir, "COVID-19/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_US.csv",sep="")
}

# path to directory containing environmental data
dir_env <- "~/workspace/dataset/"
# The radius of a buffer around each point from which to extract cell values
if(dataset_type == "global"){
    extract_buffer <- 50000
} else {
    extract_buffer <- 10000
}

# fix observational period (if NA, the period between the first date and peak date considered)
fixed_observational_period <- NA
# threshold of confirmed cases
case_count_threshold <- 30

# load data
d <-read.csv(filename, check.names = FALSE)
if(dataset_type == "US"){
    names(d)[[10]] <- "Long"
}
# check data table dimension
size <- dim(d)
# get date
date_reference <- as.Date(names(d)[first_idx_count:size[[2]]], tryFormats = c("%m/%d/%y"))

## load basis WorldClim data ##########################################
# Elevation
# https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_elev.zip
raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_elev.tif", sep="")
# load raster file
r_elevation <- raster(raster_filename)

## Bioclimatic variables Current conditions
# https://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_bio.zip
## annual mean temperature
raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_1.tif", sep="")
# load raster file
r_bio1_current <- raster(raster_filename)
# compute slope (spatial gradient) of current temperature
r_slope_bio1_current <- terrain(r_bio1_current, opt='slope', neighbors=4)

## Temperature Seasonality
# load raster file
raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_4.tif", sep="")
# load raster file
r_bio4_current <- raster(raster_filename)

## Precipitation Seasonality
# load raster file
raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_bio/wc2.1_2.5m_bio_15.tif", sep="")
# load raster file
r_bio15_current <- raster(raster_filename)

## Bioclimatic variables Past conditions
# http://biogeo.ucdavis.edu/data/climate/worldclim/1_4/grid/pst/21k/wc_2_5m_CCSM_21k_bio.zip
# mean temperature
raster_filename <- paste(dir_env, "world_clim/wc_past_climate_2_5m/wc_2_5m_CCSM_21k_bio_1.bil", sep="")
# load raster file
r_bio1_past <- raster(raster_filename)

## load SEDAC data ##############################################
# Global Human Footprint (Geographic), v2 (1995 – 2004)
# https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic/data-download
raster_filename <- paste(dir_env, "hfp-global-geo-grid/hfp_global_geo_grid/hf_v2geo/w001001.adf", sep="")
# load raster file
r_human_footprint <- raster(raster_filename)

# Population Density, v4.11 (2000, 2005, 2010, 2015, 2020)
# https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11/data-download
raster_filename <- paste(dir_env, "population_density/gpw_v4_population_density_rev11_2020_2pt5_min.tif", sep="")
# load raster file
r_population_density <- raster(raster_filename)

## Global Aridity and PET Data ####################################
# https://figshare.com/articles/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/3
# generate filename
raster_filename <- paste(dir_env, "Global_Aridity_and_PET/ai_et0/ai_et0.tif", sep="")
# load raster file
r_aridity <- raster(raster_filename)

## Gridded global datasets for Gross Domestic Product and Human Development Index over 1990-2015 #####
# https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0
# GDP per capita
# generate filename
raster_filename <- paste(dir_env, "GDP_HDI/GDP_per_capita_PPP_1990_2015_v2.nc", sep="")
# load raster file
r_GDP <- raster(raster_filename)

# Human Development Index
# generate filename
raster_filename <- paste(dir_env, "GDP_HDI/HDI_1990_2015_v2.nc", sep="")
# load raster file
r_HDI <- raster(raster_filename)

environmental_data <- c()
# computation for each sample
for(sample_id in 1:size[[1]]){
    # get cumulative count data
    cumulative_count <- c(0, as.numeric(d[sample_id, first_idx_count:size[[2]]]))
    # convert into count data
    count <- diff(cumulative_count)
    count <- ifelse(count < 0, 0, count)
    # add date information
    names(count) <- date_reference

    # generate metadata
    if(dataset_type == "global"){
        metadata <- c(paste("global_",sample_id,sep=""))
    } else {
        metadata <- c(paste("US_",sample_id,sep=""))
    }
    tmp_metadata <- d[sample_id, 1:(first_idx_count - 1)]
    for(meta in tmp_metadata ){
        if(is.factor(meta)){
            meta <- as.character(meta)
        }
        metadata <- c(metadata, meta)
    }

    # check whether count data available 
    if(length(which(cumulative_count >= case_count_threshold)) > 0){
        # get the first date when confirmed the cases and its index
        first_idx <- as.integer(which(cumulative_count >= case_count_threshold)[[1]] - 1)
        first_date <- date_reference[[first_idx]]

        # generate onset date data
        onset_date_data <- c()
        for(i in which(date_reference >= first_date)){
            onset_date_data <- c(onset_date_data, rep(as.character(date_reference[[i]]), count[[i]]))
        }
        onset_date_data <- as.Date(onset_date_data)
        # conver to incidence format
        i.1 <- incidence(onset_date_data, interval = 1)

        # get end date for the analysis and its index
        # if not considered fix observational period 
        if(is.na(fixed_observational_period)){
            # find peak
            peak_data <- estimate_peak(i.1)
            # set target end date to peak date
            target_end_date <- peak_data$estimated
        } else {
            target_end_date <- first_date + (fixed_observational_period - 1)
        }
        period_idx <- which(date_reference <= target_end_date)
        end_idx <- as.integer(period_idx[[length(period_idx)]])
        end_date <- date_reference[[end_idx]]
        # actual observational period (days) between the first date and end date
        observational_period <- as.integer(round(end_date - first_date) + 1)

        # get median date
        median_date <- first_date + observational_period / 2
        # get month
        month <- as.integer(format(median_date, "%m"))
    } else {
        # get month
        month <- as.integer(format(date_reference[[length(date_reference)]], "%m"))
    }

    ## WorldClim month database ##################################################
    ## minimum temperature (°C)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_tmin.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_tmin/wc2.1_2.5m_tmin_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    temp.min <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## maximum temperature (°C)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_tmax.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_tmax/wc2.1_2.5m_tmax_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    temp.max <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## average temperature (°C)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_tavg.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_tavg/wc2.1_2.5m_tavg_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    temp.ave.set <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer)[[1]]
    temp.ave <- mean(temp.ave.set[is.finite(temp.ave.set)])

    ## diurnal temperature range
    d.temp.range <- temp.max - temp.min

    ## temperature seasonality
    temp.seasonality <- extract(r_bio4_current, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)
    
    ## precipitation (mm)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_prec.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_prec/wc2.1_2.5m_prec_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    precipitation <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## precipitation seasonality
    precipitation.seasonality <- extract(r_bio15_current, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## solar radiation (kJ m-2 day-1)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_srad.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_srad/wc2.1_2.5m_srad_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    solar.rad <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## wind speed (m s-1)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_wind.zip
    # generate filename
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_wind/wc2.1_2.5m_wind_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    wind.speed <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## water vapor pressure (kPa)
    # http://biogeo.ucdavis.edu/data/worldclim/v2.1/base/wc2.1_2.5m_vapr.zip
    raster_filename <- paste(dir_env, "world_clim/wc2.1_2.5m_vapr/wc2.1_2.5m_vapr_", formatC(month, width=2, flag="0"), ".tif", sep="")
    # load raster file
    r <- raster(raster_filename)
    # get data
    water.vapor.pressure.set <- extract(r, d[sample_id, c("Long","Lat")], buffer=extract_buffer)[[1]]
    water.vapor.pressure <- mean(water.vapor.pressure.set[is.finite(water.vapor.pressure.set)])

    # compute relative humidity
    equilibrium.vapor.pressure <- 6.116441 * 10**(7.591386 * temp.ave.set / (temp.ave.set + 240.7263)) / 10
    humidity.set <- water.vapor.pressure.set / equilibrium.vapor.pressure
    humidity <- mean(humidity.set[is.finite(humidity.set)])

    ## Elevation
    elevation <- extract(r_elevation, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## Warming velocity
    # get current temperature
    bio1.current <- extract(r_bio1_current, d[sample_id, c("Long","Lat")], buffer=extract_buffer)
    # get past temperature
    bio1.past <- extract(r_bio1_past, d[sample_id, c("Long","Lat")], buffer=extract_buffer)
    # get spatial gradient of temperature
    bio1.slope <- extract(r_slope_bio1_current, d[sample_id, c("Long","Lat")], buffer=extract_buffer)
    # compute warming velocity
    warming.velocity <- (bio1.current[[1]] - bio1.past[[1]] * 0.1) / bio1.slope[[1]]
    # remove Inf and NA
    warming.velocity <- mean(warming.velocity[is.finite(warming.velocity)])

    ## Human footprint
    human.footprint <- extract(r_human_footprint, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## Population density
    population.density <- extract(r_population_density, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## Aridity index
    aridity <- extract(r_aridity, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## GDP per capita
    gdp.per.capita <- extract(r_GDP, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    ## Human development index
    human.dev.index <- extract(r_HDI, d[sample_id, c("Long","Lat")], buffer=extract_buffer, fun=mean)

    sum_data <- c(metadata, temp.min, temp.max, temp.ave, d.temp.range, temp.seasonality, precipitation, precipitation.seasonality, aridity, solar.rad, wind.speed, water.vapor.pressure, humidity, elevation, warming.velocity, human.footprint, gdp.per.capita, human.dev.index, population.density)
    # print
    cat(sum_data,"\n")
    # summarize the data
    environmental_data <- rbind(environmental_data, sum_data)
}

# generate header
header_names <- c("sample.id", names(d)[1:(first_idx_count - 1)])
param_names <- c("temp.min", "temp.max", "temp.ave", "d.temp.range", "temp.seasonality", "precipitation", "precipitation.seasonality", "aridity", "solar.rad", "wind.speed", "water.vapor.pressure", "humidity", "elevation", "warming.velocity", "human.footprint", "gdp.per.capita", "human.dev.index", "population.density")
header_names <- c(header_names, param_names)

# add header
environmental_data <- as.data.frame(environmental_data)
names(environmental_data) <- header_names

# output csv file
if(is.na(fixed_observational_period)){
    filename <- paste("data_latest/environmental_parameters_", dataset_type, "_period_first2peak_case_count_threshold_", case_count_threshold, "_radius_", extract_buffer, ".csv", sep="")
} else {
    filename <- paste("data_latest/environmental_parameters_", dataset_type, "_period_first_", fixed_observational_period, "_days_case_count_threshold_", case_count_threshold ,"_radius_", extract_buffer, ".csv", sep="")
}
write.csv(environmental_data, filename, row.names=F)


