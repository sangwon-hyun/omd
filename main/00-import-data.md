Importing data
================
Compiled at 2021-07-26 18:38:31 UTC

``` r
knitr::opts_chunk$set(fig.width=14, fig.height=8, echo=TRUE, eval=TRUE, cache=FALSE,
                      warning=FALSE, message=FALSE,
                      cache.lazy = FALSE)

## Load packages
library(tidyverse)
library(knitr)
library(here)
```

``` r
base = "00-import-data"
here::i_am("00-import-data.Rmd")
knitr::opts_chunk$set(fig.path = here::here("data", base, 'figures/'))
datadir = here::here("data", base)
## outputdir = here::here("data", base)
if(!dir.exists(datadir)) dir.create(datadir)
```

There are three datasets for this project, all placed
[/data/00-import-data]() [I’m a relative reference to a repository
file](../blob/master/LICENSE)

  - Climatology data, from
    [github](https://github.com/brorfred/ocean_clustering).

  - Non-climatology data, from [netcdf files in a web server
    directory](https://rsg.pml.ac.uk/shared_files/brj/proj/CBIOMES/Justin/).

  - Depth data, downloaded from CMAP (not done yet):

# Climatology data

The **climatology** chlorophyll dataset is a monthly dataset having
averaged 20 years worth of data for each month (January through
December).

There are multiple resolutions of this data; we mainly make use of the
1-degree resolution data.

``` r
## Read in data
filenames = c("tabulated_darwin_montly_clim_045_090_ver_0_2_6.csv",
              "tabulated_darwin_montly_clim_090_180_ver_0_2_6.csv", ## This will be used.
              "tabulated_darwin_montly_clim_180_360_ver_0_2_6.csv",
              "tabulated_darwin_montly_clim_360_720_ver_0_2_6.csv",

              "tabulated_geospatial_montly_clim_045_090_ver_0_2_5.csv",
              "tabulated_geospatial_montly_clim_090_180_ver_0_2_5.csv", ## This will be used.
              "tabulated_geospatial_montly_clim_180_360_ver_0_2.csv", ## The 2.5 version is missing.
              "tabulated_geospatial_montly_clim_360_720_ver_0_2_5.csv")
ddat = read.csv(file.path(datadir, filenames[2])) %>% as_tibble()
rdat = read.csv(file.path(datadir, filenames[2+4])) %>% as_tibble()
print(rdat)
```

    ## # A tibble: 194,400 x 19
    ##        X month   lat   lon   SST   Chl   PAR Kd490 euphotic_depth   mld  wind   EKE bathymetry Rrs412 Rrs443 Rrs490 Rrs510 Rrs555
    ##    <int> <int> <dbl> <dbl> <dbl> <dbl> <dbl> <dbl>          <dbl> <dbl> <dbl> <dbl>      <dbl>  <dbl>  <dbl>  <dbl>  <dbl>  <dbl>
    ##  1     0     1 -89.8 -180.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  2     1     1 -89.8 -178.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  3     2     1 -89.8 -176.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  4     3     1 -89.8 -174.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  5     4     1 -89.8 -172.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  6     5     1 -89.8 -170.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  7     6     1 -89.8 -168.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  8     7     1 -89.8 -166.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ##  9     8     1 -89.8 -164.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ## 10     9     1 -89.8 -162.    NA    NA    NA    NA             NA    NA    NA    NA          0     NA     NA     NA     NA     NA
    ## # … with 194,390 more rows, and 1 more variable: Rrs670 <dbl>

``` r
print(ddat)
```

    ## # A tibble: 130,692 x 19
    ##        X month   lat   lon   SST  SALT   Chl  Rrs412  Rrs443  Rrs490  Rrs510  Rrs555   Rrs670   MLD bathymetry  wind     TKE   PAR
    ##    <int> <int> <dbl> <dbl> <dbl> <dbl> <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>    <dbl> <dbl>      <dbl> <dbl>   <dbl> <dbl>
    ##  1   900     1 -79.8 -180. 0.726  34.0 0.495 0.00493 0.00532 0.00534 0.00453 0.00290 0.000342  19.3      172.   4.03 1.69e-4  532.
    ##  2   901     1 -79.8 -178. 0.671  34.0 0.486 0.00501 0.00539 0.00535 0.00451 0.00287 0.000337  19.8      140.   3.97 1.31e-4  522.
    ##  3   902     1 -79.8 -176. 0.643  33.9 0.480 0.00509 0.00545 0.00537 0.00451 0.00285 0.000333  19.5       75.0  3.96 1.32e-4  513.
    ##  4   903     1 -79.8 -174. 0.659  33.9 0.478 0.00515 0.00549 0.00537 0.00450 0.00283 0.000329  19.0       55    3.95 1.34e-4  508.
    ##  5   904     1 -79.8 -172. 0.702  33.8 0.478 0.00518 0.00551 0.00537 0.00448 0.00280 0.000322  18.3       65    3.95 1.32e-4  502.
    ##  6  1077     1 -79.8  174. 0.803  34.1 0.503 0.00484 0.00525 0.00531 0.00452 0.00290 0.000338  21.5      105.   3.85 1.50e-4  585.
    ##  7  1078     1 -79.8  176. 0.780  34.1 0.501 0.00485 0.00526 0.00531 0.00451 0.00290 0.000338  20.9      116.   3.94 1.67e-4  573.
    ##  8  1079     1 -79.8  178. 0.748  34.0 0.498 0.00488 0.00528 0.00532 0.00451 0.00289 0.000339  19.6      172.   3.98 1.70e-4  549.
    ##  9  1080     1 -77.8 -180. 0.317  34.0 0.461 0.00556 0.00585 0.00561 0.00465 0.00287 0.000327  21.1      635.   5.49 2.77e-4  529.
    ## 10  1081     1 -77.8 -178. 0.315  33.8 0.459 0.00563 0.00590 0.00559 0.00461 0.00283 0.000320  18.8      635.   5.51 2.65e-4  519.
    ## # … with 130,682 more rows, and 1 more variable: euphotic_depth <dbl>

# Non-climatology data

Non-climatology dataset was compiled by Bror Johnsson and shared as
netcdf files.

``` r
print('Hi world')
```

    ## [1] "Hi world"

# Depth data

The depth dataset has been downloaded from Simon’s CMAP, from the tables
"" and "". Here is the code to do compile this data.
