# covid19climate
This repository contains the datasets and R codes for reproducing the results in [*Global COVID-19 transmission rate is influenced by precipitation seasonality and the speed of climate temperature warming*](https://www.medrxiv.org).

## Terms of use
MIT licensed.
Happy if you cite our study when using the data and codes:

Chiyomaru K & Takemoto K (2020) **Global COVID-19 transmission rate is influenced by precipitation seasonality and the speed of climate temperature warming.** *in preparation*.


## Datasets
Our datasets were generated based on [2019 Novel Coronavirus COVID-19 (2019-nCoV) Data Repository by Johns Hopkins CSSE (JHU CSSE COVID-19)](https://github.com/CSSEGISandData/COVID-19).

``data_used_in_manuscirpt`` directory contains the datasets of the epidemic parameters, climate parameters, etc., used in [our study](https://www.medrxiv.org).
The datasets were generated using the R codes (see "Codes" section).
We used the [JHU CSSE COVID-19 data](https://github.com/CSSEGISandData/COVID-19) for the period between January 22, 2020 - April 6, 2020.

``data_latest`` directory contains the datasets generated using the latest [JHU CSSE COVID-19 data](https://github.com/CSSEGISandData/COVID-19) (updated when requested).


## Codes

### Data analysis
``data_analysis.R`` analyses the data: ordinary least-squares regression and spatial analysis (spatial eigenvector mapping modeling).

```
Rscript data_analysis.R | tee output.txt
```

### Epidemic parameters
NOTE: need to clone [JHU CSSE COVID-19](https://github.com/CSSEGISandData/COVID-19).
```
git clone https://github.com/CSSEGISandData/COVID-19.git
```

``compute_epidemic_parameters.R`` computes the following parameters using the [time-series data](https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data/csse_covid_19_time_series) in JHU CSSE COVID-19.

```
Rscript compute_epidemic_parameters.R [dataset_type]
```
#### Argument
``dataset_type``: Dataset type, ``global`` or ``US``. 

#### List of the parameters
* number of confirmed cases per day
* total number of confirmed cases
* [power law exponent](https://www.medrxiv.org/content/10.1101/2020.03.31.20048827v1)
* [growth rate (incidence package)](https://www.repidemicsconsortium.org/incidence/)
* [doubling time](https://www.repidemicsconsortium.org/incidence/)
* [growth rate (R0 package)](https://bmcmedinformdecismak.biomedcentral.com/articles/10.1186/1472-6947-12-147)
* [R0 (R0 package)](https://bmcmedinformdecismak.biomedcentral.com/articles/10.1186/1472-6947-12-147)
  * R0 estimated from exponential growth rate
  * R0 estimated based on maximum likelihood
* [R0 (earlyR package)](https://www.repidemicsconsortium.org/earlyR/)

The parameters are computed from the data within a certain period (15 days since 30 cases have first recorded, in default).
For estimating R0, the serial interval distributions in [Nishiura et al. (2020)](https://www.ijidonline.com/article/S1201-9712(20)30119-3/fulltext), [Wang et al. (2020)](https://papers.ssrn.com/sol3/papers.cfm?abstract_id=3551767), and [Du et al. (2020)](https://wwwnc.cdc.gov/eid/article/26/6/20-0357_article) are used, respectively.
See also the comments in the codes for details.

### Climate parameters, etc.

``get_environmental_parameters.R`` obtains the following parameters from the databases based on the [JHU CSSE COVID-19 data](https://github.com/CSSEGISandData/COVID-19/tree/master/csse_covid_19_data).

```
Rscript get_environmental_parameters.R [dataset_type]
```
#### Argument
``dataset_type``: Dataset type, ``global`` or ``US``.

NOTE: need to manually download the raster data from the databases and to **modify** the code to specify the directories containing them.
See also the comments in the codes for details.

#### List of the parameters

* [average temperature](https://worldclim.org/data/worldclim21.html)*
* [maximum temperature](https://worldclim.org/data/worldclim21.html)*
* [minimum temperature](https://worldclim.org/data/worldclim21.html)*
* diurnal temperature range*: maximum temperature - minimum temperature
* [precipitation](https://worldclim.org/data/worldclim21.html)*
* [solar radiation](https://worldclim.org/data/worldclim21.html)*
* [wind speed](https://worldclim.org/data/worldclim21.html)*
* [water vapor pressure](https://worldclim.org/data/worldclim21.html)*
* relative humidity* (computed based on water vapor pressure and average temperature)
* [aridity index](https://figshare.com/articles/Global_Aridity_Index_and_Potential_Evapotranspiration_ET0_Climate_Database_v2/7504448/3)
* [temperature seasonality](https://worldclim.org/data/worldclim21.html)
* [precipitation seasonality](https://worldclim.org/data/worldclim21.html)
* [elevation](https://worldclim.org/data/worldclim21.html)
* [warming velocity](https://science.sciencemag.org/content/334/6056/660): computed based on [current annual mean temperature (AMT)](https://worldclim.org/data/worldclim21.html) and [past AMT (CCSM)](http://www.worldclim.com/past).
* [human footprint](https://sedac.ciesin.columbia.edu/data/set/wildareas-v2-human-footprint-geographic)
* [population density](https://sedac.ciesin.columbia.edu/data/set/gpw-v4-population-density-rev11)
* [GDP per capita](https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0)
* [Human development index](https://datadryad.org/stash/dataset/doi:10.5061/dryad.dk1j0)

*monthly parameters


## Contact
[Kazuhiro Takemoto](https://sites.google.com/site/kztakemoto/)
 * e-mail: takemoto@bio.kyutech.ac.jp
 * twitter: [@kztakemoto](https://twitter.com/kztakemoto)