# load libraries
library(spdep)
library(spatialreg)
library(MuMIn)
options(na.action = "na.fail")

# first index for count data
first_idx_count <- 5

# directry containing data
dir_data <- "data_used_in_manuscript/"
#dir_data <- "data_latest/"

## load global data
# load epidemic parameter data
d_epi_global <- read.csv(paste(dir_data, "epidemic_parameters_global_period_first_15_days_case_count_threshold_30.csv", sep=""))
# load environmental parameter data
d_env_global <- read.csv(paste(dir_data, "environmental_parameters_global_period_first_15_days_case_count_threshold_30_radius_50000.csv", sep=""))

## load US data
# load epidemic parameter data
d_epi_US <- read.csv(paste(dir_data, "epidemic_parameters_US_period_first_15_days_case_count_threshold_30.csv", sep=""))
names(d_epi_US)[[11]] <- "Long"
names(d_epi_US)[[8]] <- "Province.State"
names(d_epi_US)[[9]] <- "Country.Region"
d_epi_US <- d_epi_US[,-c(2:7,12)]
# load environmental parameter data
d_env_US <- read.csv(paste(dir_data, "environmental_parameters_US_period_first_15_days_case_count_threshold_30_radius_10000.csv", sep=""))
d_env_US <- d_env_US[,-c(2:7,12)]

## combine global data and US data
d_global <- cbind(d_epi_global, d_env_global[first_idx_count:dim(d_env_global)[[2]]])
# remove US from global data
d_global <- d_global[-which(d_global$Country.Region == "US"),]
d_US <- cbind(d_epi_US, d_env_US[first_idx_count:dim(d_env_US)[[2]]])
d <- rbind(d_global, d_US)

## data preprocessing #################################
d$first.date <- as.Date(d$first.date)
# categorical variable for travel restrictions
d$ban <- ifelse(d$first.date < as.Date("2020-03-17"), 0, 1)
# remove short-term data
d <- d[d$observational.period >= 15, ]
# rescale temperature due to the quadratic relationship
d$temp.ave.peak <- sqrt((d$temp.ave - 7.8) ** 2)

epi_param <- c("nb.cases.per.day","power.law.exponent","growth.rate.incidence.package","doubling.time","earlyR0.nishiura","earlyR0.wang","earlyR0.du","growth.rate.R0.package","R0.EG.nishiura","R0.EG.wang","R0.EG.du","R0.ML.nishiura","R0.ML.wang","R0.ML.du","total.nb.cases")
env_param <- c("temp.min", "temp.max", "temp.ave", "temp.ave.peak", "d.temp.range", "temp.seasonality", "precipitation", "precipitation.seasonality", "aridity", "solar.rad", "wind.speed", "water.vapor.pressure", "humidity", "elevation", "warming.velocity", "gdp.per.capita", "human.dev.index", "human.footprint", "population.density", "ban", "pm25")

# simple pairwise correlation analysis
for(epi in epi_param){
    for(env in env_param){
        sum <- cor.test(d[[epi]], d[[env]], m="s")
        if(sum$p.value < 0.05){
            cat(epi, env, sum$estimate, sum$p.value, "\n")
        }
    }
}

d$population.density <- ifelse(d$population.density > 0, d$population.density, NA)
d$total.nb.cases <- ifelse(d$total.nb.cases > 0, d$total.nb.cases, NA)
d$growth.rate.incidence.package <- ifelse(is.finite(d$growth.rate.incidence.package), d$growth.rate.incidence.package, NA)

# exclude NA
d <- na.omit(d[c("sample.id","growth.rate.incidence.package", "total.nb.cases", "Long", "Lat", env_param)])

# check quadratic relationship between temperature and transmission rate
summary(nls(log(total.nb.cases)~b+c_temp*(temp.ave-temp_c)**2 + c_hum*(humidity - hum_c)**2, data=d, start=c(b=1,c_temp=1, c_hum=1, temp_c=1, hum_c=1)))

# Normalization
d_sub <- data.frame(
    d[c("Long", "Lat", "ban")],
    scale(d[c("temp.ave", "temp.ave.peak", "d.temp.range", "humidity", "temp.seasonality", "wind.speed", "elevation", "human.dev.index","solar.rad")]),
    scale(d[c("growth.rate.incidence.package")]),
    #scale(d[c("power.law.exponent")]),
    scale(sqrt(d[c("precipitation", "precipitation.seasonality")])),
    scale(log(d[c("aridity", "population.density", "gdp.per.capita")])),
    scale(abs(d[c("pm25")])**(1/3)),
    scale(log(d[c("total.nb.cases")])),
    #scale(log(d[c("nb.cases.per.day")])),
    #scale(log(d[c("doubling.time")])),
    scale(log(abs(d["warming.velocity"])))
)

cat(dim(d_sub),"\n")

## OLS regression analysis #######################################
res1 <- lm(growth.rate.incidence.package~temp.ave.peak+d.temp.range+temp.seasonality+wind.speed+precipitation+precipitation.seasonality+solar.rad+humidity+population.density+human.dev.index+warming.velocity+ban+pm25, data=d_sub)

res1_models <- dredge(res1,rank="AICc")
res1_best <- get.models(res1_models, subset = 1)[1]
res1_ave <-model.avg(get.models(res1_models, cumsum(weight) <= .95))

# OLS full model
summary(res1)
AICc(res1)
# OLS best model
summary(res1_best[[1]])
AICc(res1_best[[1]])
# OLS averaged model
summary(res1_ave)


## Spatial analysis based on SEVM modeling approach ##############
set.seed(123)
long <- d_sub$Long
lat <- d_sub$Lat
n_num <- length(d_sub[[1]])
longlat <- cbind(long+rnorm(n_num)/1000,lat+rnorm(n_num)/1000)
longlat_knn <- knearneigh(longlat,k=2,longlat=T)
longlat_nb <- knn2nb(longlat_knn)
longlat_listw <- nb2listw(longlat_nb)

# Moran's I test for the OLS full model
lm.morantest.exact(res1,longlat_listw)
# Moran's I test for the OLS best model
lm.morantest.exact(res1_best[[1]],longlat_listw)

# A spatial eigenvector mapping (SEVM) modelling approach removing spatial autocorrelation in the model residuals
lagcol <- SpatialFiltering(growth.rate.incidence.package~1, ~temp.ave.peak+d.temp.range+temp.seasonality+wind.speed+precipitation+precipitation.seasonality+solar.rad+humidity+population.density+human.dev.index+warming.velocity+ban+pm25-1, data=d_sub, style="W", nb=longlat_nb, ExactEV=T)

res1_corr <- lm(growth.rate.incidence.package~temp.ave.peak+d.temp.range+temp.seasonality+wind.speed+precipitation++precipitation.seasonality+solar.rad+humidity+population.density+human.dev.index+warming.velocity+ban+pm25+fitted(lagcol), data=d_sub)

res1_corr_models <- dredge(res1_corr,rank="AICc",fixed = "fitted(lagcol)")
res1_corr_best <- get.models(res1_corr_models, subset = 1)[1]
res1_corr_ave <- model.avg(get.models(res1_corr_models, cumsum(weight) <= .95))

# SEVM full model
summary(res1_corr)
AICc(res1_corr)
# SEVM best model
summary(res1_corr_best[[1]])
AICc(res1_corr_best[[1]])
# SEVM averaged model
summary(res1_corr_ave)

# Moran's I test for the SEVM full model
lm.morantest.exact(res1_corr,longlat_listw)
# Moran's I test for the SEVM best model
lm.morantest.exact(res1_corr_best[[1]],longlat_listw)
