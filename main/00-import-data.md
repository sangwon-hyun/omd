Importing data
================
Compiled at 2022-07-28 22:19:48 UTC

``` r
knitr::opts_chunk$set(fig.width=14, fig.height=8, echo=TRUE, eval=TRUE, cache=FALSE,
                      warning=FALSE, message=FALSE,
                      cache.lazy = FALSE)

## Load packages
library(tidyverse, quietly = TRUE)
library(knitr, quietly = TRUE)
library(here, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidync, quietly = TRUE)
library(omd, quietly = TRUE)
```

    ## Warning: replacing previous import 'dplyr::union' by 'raster::union' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::select' by 'raster::select' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::intersect' by 'raster::intersect' when loading 'omd'

``` r
library(ggstar, quietly = TRUE)
sf::sf_use_s2(FALSE)
```

The R package `omd` to use is here
<https://github.com/sangwon-hyun/omd/>.

``` r
base = "00-import-data"
here::i_am("00-import-data.Rmd")
knitr::opts_chunk$set(fig.path = here::here("data", base, 'figures/'))
datadir = here::here("data", base)
if(!dir.exists(datadir)) dir.create(datadir)
figdir = here::here("figures")
source(here::here("00-helpers.R"))
```

There are three datasets for this project, all placed here
[/data/00-import-data](/data/00-import-data).

  - Climatology data, from
    [github](https://github.com/brorfred/ocean_clustering).

  - Non-climatology data, from [netcdf files in a web server
    directory](https://rsg.pml.ac.uk/shared_files/brj/proj/CBIOMES/Justin/).

  - Depth data, downloaded from CMAP. See
    [./depth-analysis-as-is](./depth-analysis-as-is).

## Climatology data

The **climatology** chlorophyll dataset is a monthly dataset having
averaged 20 years worth of data for each month (January through
December).

There are multiple resolutions of this data; we mainly make use of the
2-degree resolution data. Here is some sample code for loading this.

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

    ## # A tibble: 194,400 × 19
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

    ## # A tibble: 130,692 × 19
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

## Non-climatology data

Non-climatology dataset was compiled by Bror Johnsson and shared as
netcdf (`.nc`) files. This data was further cleaned (code shown further
below) and saved as RDS files containing data frames between 1998
through 2001. These RDS files are named like this:

`filename = paste0(type, "-nonclim-", year, "-", mo, ".RDS")`

where `type` is one of “darwin” or “real” (meaning remote sensing), and
`year` is one of {1998, 1999, 2000, .., 2006}, and `mo` is month number
1, …, 12. Here’s an example of one month’s data:

``` r
type = "real"
year = 2000
mo = 6
filename = paste0(type, "-nonclim-", year, "-", mo, ".RDS")
folder = paste0(type, "-nonclim")
one_month_dat = readRDS(file = file.path(datadir, folder, filename))
mo = unique(one_month_dat$mo)
year = unique(one_month_dat$year)
one_month_dat %>% 
  plot_dat(add_map = TRUE) +
  ggtitle(paste0(month.abb[mo], ",", year)) +
  theme(plot.title = element_text(size = rel(2), face = "bold"))
```

These data files (named
e.g. `data/darwin-nonclim/darwin-nonclim-2000-1.RDS`) were cleaned and
produced using this code (not run now):

``` r
## Latitude/longitude range 
lat = 19.8968
lon = -155.5828
boxsize = 40
lonrange = lon + c(-1,1) * boxsize
latrange = lat + c(-1,1) * boxsize

## Time indices, since this is too expensive to do in one swipe
source("nonclim-helpers.R")
for(dat_type in c("darwin", "real")){
  cat("data type is ", dat_type, fill = TRUE)
  if(dat_type == "darwin"){
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin_time.csv")
  }
  if(dat_type == "real"){
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_time.csv")
  }
  times = read.csv(timefile)
  maxt = nrow(times)
  x = 1:maxt
  n = 20
  indices = split(x, sort(x %% n))
  num_groups = length(indices)
  start.time = Sys.time()
  for(ii in 1:num_groups){
    ind = indices[[ii]]
    dat = get_nonclim_dat(type = dat_type, datadir = datadir,
                                 lonrange = lonrange,
                                 latrange = latrange,
                                 time_range = ind)
    filename = paste0(dat_type, "-nonclim-", ii, ".RDS")
    cat("saving to", filename, fill = TRUE)
    saveRDS(dat, file = file.path(datadir, filename))
  }
}


## Further process and save data.
for(dat_type in c("darwin", "real")){
  cat("data type is ", dat_type, fill = TRUE)
  if(dat_type == "darwin"){
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin_time.csv")
  }
  if(dat_type == "real"){
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_time.csv")
  }

  ## Read in individual files (equal chunks of dates)
  dat_list = list()
  num_groups = 20
  for(ii in 1:num_groups){
    filename = paste0(dat_type, "-nonclim-", ii, ".RDS")
    dat_list[[ii]] = readRDS(file = file.path(datadir, filename))
  }

  ## Combine and save
  alldat = data.table::rbindlist(dat_list)
  combined_filename = paste0(dat_type, "-nonclim-combined.RDS")
  saveRDS(alldat, file = file.path(datadir, combined_filename))

  ## Read in combined data
  alldat = readRDS(file = file.path(datadir, combined_filename))

  ## Split into individual months.
  dlist = alldat %>% group_by(mo, year) %>% group_split()

  ## Get 99th quantile of all chlorophyll values (to filter the extreme values)
  maxval = alldat %>% summarize(q=quantile(val, 0.9999)) %>% unlist()
  col_cuts = seq(from = 0, to = maxval, length=100)
  rm(alldat)

  ## Save each month-year into separate file.
  for(ii in 1:length(dlist)){

    ## Extract month information
    thismo = dlist[[ii]] %>% pull(mo) %>% unique()
    thisyear = dlist[[ii]] %>% pull(year) %>% unique()
    stopifnot((length(thismo) == 1) & (length(thisyear) == 1))

    ## Aggregate to unique month (at original resolution)
    one_month_dat = dlist[[ii]] %>%
      dplyr::filter(mo == thismo) %>%
      dplyr::filter(year == thisyear)

    ## Sanity check 1: Plotting the Chlorophyll map.
    plotfilename = paste0(dat_type, "-nonclim-", thisyear, "-", thismo, ".png")
    plot_dat(one_month_dat, add_map = TRUE) +
      ggtitle(paste0(month.abb[thismo], ",", thisyear)) +
      theme(plot.title = element_text(size = rel(2), face = "bold")) +
      ggsave(file.path(datadir, plotfilename), width=5, height=7, units="in")

    ## Sanity check 2: Plotting but with the same colors in all plots.
    plotfilename = paste0("samecol-", dat_type, "-nonclim-", thisyear, "-", thismo, ".png")
    plot_dat(one_month_dat, add_map = TRUE) +
      ggtitle(paste0(month.abb[thismo], ",", thisyear)) +
      scale_fill_gradientn(colours = c("blue", "red", "yellow"),
                           guide="colorbar", values = col_cuts) +
      theme(plot.title = element_text(size = rel(2), face = "bold")) +
      ggsave(file.path(datadir, plotfilename), width=5, height=7, units="in")

    filename = paste0(dat_type, "-nonclim-", thisyear, "-", thismo, ".RDS")
    saveRDS(one_month_dat, file = file.path(datadir, filename))
  }
}
```

This produces a Chlorophyll map that is larger, and shows two boxes.

<!-- Also, this should remove the need for the star in Figure 1. -->

``` r
## Form Smaller box
lat = 19.8968
lon = -155.5828
boxsize = 30
lonrange_small = lon + c(-1,1)*boxsize
latrange_small = lat + c(-.3,1)*boxsize

## Medium box
boxsize = 40
lonrange = lon + c(-1,1) * boxsize
latrange = lat + c(-1,1) * boxsize
latrange = pmin(latrange, 48.25) 
latrange = pmax(latrange, -20.5) 

## Medium box
boxsize = 40
lonrange_large = lon + c(-1,1) * boxsize
latrange_large = lat + c(-1,1) * boxsize

## Get January Darwin data
datadir = "/home/sangwonh/repos/omd/main/data/00-import-data"
ddat = read.csv(file.path(datadir, filenames[2])) %>% as_tibble()

## Successively narrow down the data
smallerdat = ddat %>% filter(month == 1) %>% dplyr::select(lon, lat, val = Chl)
smallerdat2 = smallerdat
smallerdat2$lon = smallerdat2$lon - 360
smallerdat = smallerdat %>% rbind(smallerdat2)
smallerdat = smallerdat %>% 
  dplyr::filter(lon < -50) %>%
  dplyr::filter(lon > -200) %>%
  dplyr::filter(lat < 85) %>%
  dplyr::filter(lat > -40)

p = smallerdat %>% 
  plot_dat(add_map = TRUE) +
  geom_rect(xmin = -180, xmax = lonrange[2],
            ymin = latrange[1], ymax = latrange[2],
            col = 'black', fill = NA, size = 1.3) + 
  geom_rect(xmin = -180, xmax = lonrange_small[2],
            ymin = latrange_small[1], ymax = latrange_small[2],
            col = 'white', fill = NA, size = 1.3, linetype = "dashed") +
  geom_rect(xmin = -180, xmax = lonrange_large[2],
            ymin = latrange_large[1], ymax = latrange_large[2],
            col = 'red', fill = NA, size = 1.3, linetype = "dotted") +
  scale_y_continuous(breaks = seq(from = -40,##ceiling(min(smallerdat$lat)),
                                  to = 80,   ##floor(max(smallerdat$lat)),
                                  by=20)) +
  geom_star(aes(x=lon, y=lat), data = data.frame(lat = 22.45, lon = -158), size=5,
            fill = 'yellow', col = 'black') 
plot(p)
```

![](/home/sangwonh/repos/omd/main/data/00-import-data/figures/larger-map-1.png)<!-- -->

``` r
## Plotting the map to a file.
if(FALSE){
  plotfilename = "orient-the-box.pdf"
  ggsave(filename = file.path(figdir, plotfilename), height = 5, width = 5)
}
```
