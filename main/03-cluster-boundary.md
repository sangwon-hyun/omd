Cluster boundary comparison using OMD
================
Compiled at 2022-07-29 22:26:28 UTC

``` r
knitr::opts_chunk$set(fig.width=14, fig.height=8, echo=TRUE, eval=TRUE, cache=TRUE,
                      warning=FALSE, message=FALSE)
                      ## cache.lazy = FALSE)

## Read in libraries
library(Matrix, quietly = TRUE)
library(transport, quietly = TRUE)
library(rworldmap, quietly = TRUE)
library(omd, quietly = TRUE)
library(maps, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(sf, quietly = TRUE)
library(rnaturalearth, quietly = TRUE)
library(rgeos, quietly = TRUE)
library(ggspatial, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(gridExtra, quietly = TRUE)
library(mvtnorm, quietly = TRUE)
library(omd, quietly = TRUE)
sf::sf_use_s2(FALSE)

## Helper function
go_down_lat_and_mark_changepoint <- function(onerow){
  onerow = as.numeric(onerow)
  newrow = rep(0, length(onerow))
  for(ii in 1:(length(onerow)-1)){
    if(onerow[ii+1] == 2)  break
  }
  newrow[ii] = 1
  return(newrow)
}
```

The R package `omd` to use is here
<https://github.com/sangwon-hyun/omd/>, commit ()\[\].

``` r
base = "03-cluster-boundary"
here::i_am("03-cluster-boundary.Rmd")
knitr::opts_chunk$set(fig.path = here::here("data", base, 'figures/'))
datadir = here::here("data", base)
datadir_orig = here::here("data", "00-import-data")
if(!dir.exists(datadir)) dir.create(datadir)
figdir = here::here("figures")
source(here::here("03-helpers.R"))
source(here::here("02-helpers.R"))
source(here::here("01-helpers.R"))
```

``` r
ddat = read.csv(file.path(datadir_orig, "tabulated_darwin_montly_clim_090_180_ver_0_2_6.csv")) %>% as_tibble()
rdat = read.csv(file.path(datadir_orig, "tabulated_geospatial_montly_clim_090_180_ver_0_2_5.csv")) %>% as_tibble()
lat = 19.8968
lon = -155.5828
boxsize = 30
lonrange = lon + c(-1,1)*boxsize
latrange = lat + c(-.3,1)*boxsize
restrictbox <- . %>% filter(lat > latrange[1],
                            lat < latrange[2],
                            lon > lonrange[1],
                            lon < lonrange[2])
```

**Supplemental Figure 1**

``` r
segment_pipe  = . %>% group_by(lon) %>%
    mutate(clustered_val = as.factor(kmeans(val, centers=c(0, 0.3))$cluster) ) %>%
    group_by(lon, clustered_val) %>% 
    mutate(meanval = mean(val)) %>%
    ungroup() %>% 
    group_by(lon) %>% 
    mutate(changepoint = go_down_lat_and_mark_changepoint(clustered_val)) %>% 
    arrange(lon)

plist = list()
for(mo in 1:12){

  ## Set up pipeline for data
  pipeline <- . %>%
    dplyr::filter(month == mo) %>%
    dplyr::select(lon, lat, val = Chl) %>%
    restrictbox() 

  ## Draw Darwin data boundary
  darwindat = ddat %>% pipeline()
  p = plot_dat(darwindat)
  border_darwin = darwindat %>% segment_pipe() %>% ungroup() %>%
    dplyr::select(lon, lat, val = changepoint) 
  p = p + geom_path(aes(x = lon, y = lat),
                    data = border_darwin %>% dplyr::filter(val == 1),
                    size = rel(2), col = rgb(1,1,0,0.8))##, lwd = 3)
  p = p + ylab("") + xlab("")
  p = p + theme(plot.margin=unit(-0.2 * c(1,2,1,2), "cm"))

  # Draw real data boundary
  realdat = rdat %>% pipeline()
  border_real = realdat %>% segment_pipe() %>% ungroup() %>%
    dplyr::select(lon, lat, val = changepoint)
  p = p + geom_path(aes(x = lon, y = lat),
                    data = border_real %>% dplyr::filter(val == 1),
                    size = rel(2), col = rgb(0,0,1,0.8))##, lwd = 3)
  p = p + ggtitle(month.name[mo])
  p = p + theme(plot.title=element_text(margin=margin(t=40,b=-20)))


  ## Obtain data
  twodat = inner_join(border_real,
                     border_darwin,
                     by = c("lat", "lon"),
                     suffix = c(".r", ".d")) %>% as_tibble()
  twodat$val.r = twodat$val.r + 1E-10
  twodat$val.d = twodat$val.d + 1E-10

  res = omd(M1_long = twodat %>% dplyr::select(lon, lat, val = val.d),
            M2_long = twodat %>% dplyr::select(lon, lat, val = val.r),
            p = 1) 
  all_rows = plot_omd_ggplot(res, return_arrows = TRUE)
  p = p + geom_segment(aes(x = fr_lon,
                           y = fr_lat,
                           xend = to_lon,
                           yend = to_lat,
                           size = mass),
                       col = rgb(1,1,0,0.8),
                       arrow = arrow(length = unit(0.05, "inches")),
                       data = all_rows) +
    scale_size_continuous(range = c(1E-3, 1)) 
  plist[[mo]] = p
}

do.call("grid.arrange", c(plist, ncol = 4, left="Latitude", bottom="Longitude")) 
```

![](03-cluster-boundary_files/figure-gfm/all-boundaries-1.png)<!-- -->

The first panel of **Figure 6** is here:

``` r
do.call("grid.arrange", c(plist[c(3,8)], ncol = 1, left="Latitude", bottom="Longitude")) 
```

![](03-cluster-boundary_files/figure-gfm/unnamed-chunk-1-1.png)<!-- -->

Now, the 24 x 24 distance matrix is calculated:

``` r
realdat_list = list()
darwindat_list = list()
for(mo in 1:12){

  pipeline <- . %>%
    dplyr::filter(month == mo) %>%
    dplyr::select(lon, lat, val = Chl) %>%
    restrictbox() 

  ## Obtain data
  realdat = rdat %>% pipeline()
  darwindat = ddat %>% pipeline()

  ## Create a boundary map
  realdat = realdat %>% segment_pipe() %>% ungroup() %>% dplyr::select(lon, lat, val=changepoint) 
  darwindat = darwindat %>% segment_pipe() %>% ungroup() %>% dplyr::select(lon, lat, val=changepoint) 
  darwinname = paste0("d", mo)
  realname = paste0("r", mo)
  twodat = merge(realdat, darwindat, by = c("lat", "lon"),
                 suffixes = c(".r", ".d")) %>% as_tibble() %>%
    rename(!!(realname) := val.r, !!(darwinname) := val.d) %>% 
    pivot_longer(cols = c(realname, darwinname), names_to = "dat_type", values_to = "val") 
  twodat$val = twodat$val + 1E-10

  realdat = twodat %>% filter(dat_type==realname) 
  darwindat = twodat %>% filter(dat_type==darwinname)

  realdat_list[[mo]] = realdat
  darwindat_list[[mo]] = darwindat
}
names(realdat_list) = paste0("r", 1:12)
names(darwindat_list) = paste0("d", 1:12)
alldat_list = c(realdat_list, darwindat_list)

distmat = matrix(NA, 24, 24)
for(ii in 1:24){
  for(jj in 1:24){
    if(ii < jj){

      ## Merge (which means getting rid of NAs)
      d1 = alldat_list[[ii]]
      d2 = alldat_list[[jj]]
      name1 = names(alldat_list)[ii]
      name2 = names(alldat_list)[jj]
      twodat = combine(d1, d2, name1, name2) ## from OMD

      ## Calculate OMD
      res = omd(M1_long = twodat %>% dplyr::filter(dat_type==name1),
                M2_long = twodat %>% dplyr::filter(dat_type==name2),
                geodesic = TRUE,
                p = 1)
      distmat[ii,jj] = distmat[jj,ii] = res$dist
    }
  }
}
rownames(distmat) = names(alldat_list)
colnames(distmat) = names(alldat_list)
saveRDS(distmat, file = file.path(datadir, "distmat-boundary.RDS"))
```

MDS plot **Figure 6**.

``` r
distmat = readRDS(file = file.path(datadir, "distmat-boundary.RDS"))
p = make_mds_plot(distmat, angle = -pi * 0.3) ## in 01-helpers.R
plot(p)
```

![](03-cluster-boundary_files/figure-gfm/mds-1.png)<!-- -->

MDS plot **Figure 6**.

``` r
p = plot_distmat_ggplot_clim(distmat)
plot(p)
```

![](03-cluster-boundary_files/figure-gfm/distmat-1.png)<!-- -->

``` r
plotfilename = "distmat-boundary.pdf"
```
