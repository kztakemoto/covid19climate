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
library(ggplot2)
library(mixdist)
library(distcrete)
library(incidence)
library(earlyR)
library(R0)
source("estimate_power_law.R")

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
# random seed
set.seed(123)

# fix observational period (if NA, the period between the first date and peak date considered)
fixed_observational_period <- 15
# threshold of confirmed cases
case_count_threshold <- 30

# load data
d <-read.csv(filename, check.names = FALSE)
# check data table dimension
size <- dim(d)
# get date
date_reference <- as.Date(names(d)[first_idx_count:size[[2]]], tryFormats = c("%m/%d/%y"))

# computation for each sample
epidemic_data <- c()
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
        # conver to incidence format data
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

        first_date <- as.character(first_date)
        end_date <- as.character(end_date)

        # check whether #data points > 1
        if(length(which(count[first_idx:end_idx] > 0)) > 1){
            # compute #cases per day
            nb_cases_per_day <- mean(count[first_idx:end_idx])
            # compute power law exponent
            power.law.exponent <- estimate.power.law(count[first_idx:end_idx])
            # compute incidence
            flag <- try(incidence_result <- fit(i.1[1:observational_period]))
            if(class(flag)[[1]] == "try-error"){
                incidence_result <- fit(i.1)
            }
            # get growth rate
            growth_rate_incidence <- incidence_result$info$r
            # get doubling time
            if(growth_rate_incidence > 0){
                doubling_time <- incidence_result$info$doubling
            } else {
                doubling_time <- Inf
            }

            # early R0
            # based on Nishiura et al. https://www.ijidonline.com/article/S1201-9712(20)30119-3/fulltex
            param <- weibullpar(mu=4.8, sigma=2.3)
            si_nishiura <- distcrete("weibull", 1, shape=param$shape, scale=param$scale)
            flag <- try(earlyR0_nishiura <- get_R(i.1[1:observational_period], si=si_nishiura))
            if(class(flag)[[1]] == "try-error"){
                earlyR0_nishiura <- get_R(i.1, si=si_nishiura)
            }
            # based on Wang et al. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3551767
            param <- weibullpar(mu=7.4, sigma=5.2)
            si_wang <- distcrete("weibull", 1, shape=param$shape, scale=param$scale)
            flag <- try(earlyR0_wang <- get_R(i.1[1:observational_period], si=si_wang))
            if(class(flag)[[1]] == "try-error"){
                earlyR0_wang <- get_R(i.1, si=si_wang)
            }
            # based on Du et al. https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
            si_du <- distcrete("norm", 1, mean=4.0, sd=4.8)
            flag <- try(earlyR0_du <- get_R(i.1[1:observational_period], si=si_du))
            if(class(flag)[[1]] == "try-error"){
                earlyR0_du <- get_R(i.1, si=si_du)
            }

            # generation time distribution
            # based on Nishiura et al. https://www.ijidonline.com/article/S1201-9712(20)30119-3/fulltext
            mGT_nishiura <- generation.time("weibull", c(4.8, 2.3))
            # based on Wang et al. https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3551767
            mGT_wang <- generation.time("weibull", c(7.4, 5.2))
            # based on Du et al. https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article
            # NOTE: consider Weibull distribution here (see Truncated (>0) in Appendix Table 3) because normal distribution is unavailable in generation.time function
            mGT_du <- generation.time("weibull", c(5.6, 3.8))

            # estimate R from exponential growth rate
            R0_EG_nishiura <- est.R0.EG(count, mGT_nishiura, begin=first_idx, end=end_idx)
            R0_EG_wang <- est.R0.EG(count, mGT_wang, begin=first_idx, end=end_idx)
            R0_EG_du <- est.R0.EG(count, mGT_du, begin=first_idx, end=end_idx)

            # estimate R by maximum likelihood
            flag <- try(R0_ML_nishiura <- est.R0.ML(count, mGT_nishiura, begin=first_idx, end=end_idx, range=c(1e-6,100)))
            if(class(flag)[[1]] == "try-error"){
                R0_ML_nishiura$R <- NA
            }
            flag <- try(R0_ML_wang <- est.R0.ML(count, mGT_wang, begin=first_idx, end=end_idx, range=c(1e-6,100)))
            if(class(flag)[[1]] == "try-error"){
                R0_ML_wang$R <- NA
            }
            flag <- try(R0_ML_du <- est.R0.ML(count, mGT_du, begin=first_idx, end=end_idx, range=c(1e-6,100)))
            if(class(flag)[[1]] == "try-error"){
                R0_ML_wang$R <- NA
            }

            # combine the data
            sum_data <- c(metadata, nb_cases_per_day, power.law.exponent, growth_rate_incidence, doubling_time, earlyR0_nishiura$R_ml, earlyR0_wang$R_ml, earlyR0_du$R_ml, R0_EG_nishiura$r, R0_EG_nishiura$R, R0_EG_wang$R, R0_EG_du$R, R0_ML_nishiura$R, R0_ML_wang$R, R0_ML_du$R, sum(count[first_idx:end_idx]), first_date, end_date, observational_period)
        } else{
            sum_data <- c(metadata, 0, 0, -Inf, Inf, rep(0,3), -Inf, rep(0,6), sum(count[first_idx:end_idx]), first_date, end_date, observational_period)
        }
    } else {
        sum_data <- c(metadata, 0, 0, -Inf, Inf, rep(0,3), -Inf, rep(0,6), 0, NA, NA, NA)
    }
    # print
    cat(sum_data,"\n")
    # summarize the data
    epidemic_data <- rbind(epidemic_data, sum_data)
}

# generate header
header_names <- c("sample.id", names(d)[1:(first_idx_count - 1)])
param_names <- c("nb.cases.per.day","power.law.exponent","growth.rate.incidence.package","doubling.time","earlyR0.nishiura","earlyR0.wang","earlyR0.du","growth.rate.R0.package","R0.EG.nishiura","R0.EG.wang","R0.EG.du","R0.ML.nishiura","R0.ML.wang","R0.ML.du","total.nb.cases","first.date","end.date","observational.period")
header_names <- c(header_names, param_names)

# add header
epidemic_data <- as.data.frame(epidemic_data)
names(epidemic_data) <- header_names

# output csv file
if(is.na(fixed_observational_period)){
    filename <- paste("data_latest/epidemic_parameters_", dataset_type, "_period_first2peak_case_count_threshold_", case_count_threshold, ".csv",sep="")
} else {
    filename <- paste("data_latest/epidemic_parameters_", dataset_type, "_period_first_", fixed_observational_period, "_days_case_count_threshold_", case_count_threshold, ".csv", sep="")
}
write.csv(epidemic_data, filename, row.names=F)
