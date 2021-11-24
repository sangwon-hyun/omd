Lat/lon analysis of climatology Chlorophyll data
================
Compiled at 2021-11-24 23:02:51 UTC

``` r
knitr::opts_chunk$set(fig.width=14, fig.height=8, echo=TRUE, eval=TRUE, cache=FALSE,
                      warning=FALSE, message=FALSE,
                      cache.lazy = FALSE)

## Load packages
library(tidyverse, quietly = TRUE)
library(ggrepel, quietly = TRUE)
library(knitr, quietly = TRUE)
library(here, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(ggstar, quietly = TRUE)
library(omd, quietly = TRUE)
```

    ## Warning: replacing previous import 'dplyr::select' by 'raster::select' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::intersect' by 'raster::intersect' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::union' by 'raster::union' when loading 'omd'

``` r
source("01-helpers.R")
figdir = "~/Dropbox/Apps/Overleaf/OMD/figures" 
```

The R package `omd` to use is here
<https://github.com/sangwon-hyun/omd/>, commit ()\[\].

``` r
base = "01-climatology"
here::i_am("01-climatology.Rmd")
knitr::opts_chunk$set(fig.path = here::here("data", base, 'figures/'))
datadir = here::here("data", base)
datadir_orig = here::here("data", "00-import-data")
if(!dir.exists(datadir)) dir.create(datadir)
```

Here’s the box we’ll focus on.

``` r
lat = 19.8968
lon = -155.5828
boxsize = 40
lonrange = lon + c(-1,1) * boxsize
latrange = lat + c(-1,1) * boxsize
latrange = pmin(latrange, 48.25) 
latrange = pmax(latrange, -20.5) 
restrictbox <- . %>% filter(lat > latrange[1],
                            lat < latrange[2],
                            lon > lonrange[1],
                            lon < lonrange[2])
world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
```

## Load data

First, we load Darwin data (from
<https://github.com/brorfred/ocean_clustering/>) at 2-degree resolution.

``` r
ddat = read.csv(file.path(datadir_orig, "tabulated_darwin_montly_clim_090_180_ver_0_2_6.csv")) %>% as_tibble()
rdat = read.csv(file.path(datadir_orig, "tabulated_geospatial_montly_clim_090_180_ver_0_2_5.csv")) %>% as_tibble()
```

## Gather data

All twelve month’s worth of data is first gathered into two list
objects.

``` r
## Get 24 months' worth of data
realdat_list = lapply(1:12, function(mo){ 
  rdat %>% restrictbox() %>% dplyr::filter(month == mo) %>% 
    dplyr::select(month, lat, lon, val = Chl)
})
names(realdat_list) = paste0("r", 1:12)

darwindat_list = lapply(1:12, function(mo){
  ddat %>% restrictbox() %>% dplyr::filter(month == mo) %>% 
    dplyr::select(month, lat, lon, val = Chl)
})
names(darwindat_list) = paste0("d", 1:12)
alldat_list = c(realdat_list, darwindat_list)

## Removing landlocked points
remove_landlock <- . %>%
  dplyr::mutate(land = mark_land(lon, lat, overreact = TRUE, buffer_in_km = 250)) %>%
  dplyr::filter(land == FALSE)  %>%
  dplyr::select(-land)

## Making common set of lon, lats.
alldat_list = alldat_list %>%
  purrr::map(. %>%
             remove_landlock() %>%
             dplyr::select(lon, lat, val) %>% 
             complete_rectangle() %>%
             dplyr::arrange(lon, lat))

## Check that all lon, lats are exactly the same
all_lonlats = alldat_list %>% purrr::map(. %>% dplyr::select(lon, lat))
for(ii in 1:length(all_lonlats)){
  stopifnot(all(all_lonlats[[ii]] == all_lonlats[[1]]))
}
```

## Top row of Figure 2

Calculate OMD.

``` r
## Get the two datasets
name1 = "d4"
name2 = "r4"
name1_long = "Darwin data (April)"
name2_long = "Remote sensing (April)"
d1 = alldat_list[[name1]]
d2 = alldat_list[[name2]]
twodat = combine(d1, d2, name1_long, name2_long)
plot_dat(twodat, add_map = TRUE)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/calculate-omd-1.png)<!-- -->

``` r
## Calculate 1-Wasserstein
res = omd(M1_long = twodat %>% filter(dat_type == name1_long),
          M2_long = twodat %>% filter(dat_type == name2_long),
          p = 2,
          geodesic = TRUE)
```

Make plot of two maps from the two sources.

``` r
## 1. First map
limits = c(0, 0.00483)
colours = c("blue", "red", "yellow")
p1 = twodat %>% plot_dat(add_map = TRUE,
                         hide_legend = FALSE,
                         standardize = TRUE,
                         limits = limits,
                         colours = colours, legend_name = "Chlorophyll\n     PDF")
p1 = p1 + theme(legend.position="left") +
  theme(strip.text.x = element_text(size = rel(1.4))) 
plot(p1)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/clim-two-maps-1.png)<!-- -->

``` r
## ggsave(file.path(figdir, "climatology-maps.pdf"), width = 10, height=5)
```

Make plot of pixel-wise difference

``` r
## Second map
dat1 = twodat %>% filter(dat_type == name1_long)
dat2 = twodat %>% filter(dat_type == name2_long)
dat1$val = dat1$val / sum(dat1$val, na.rm = TRUE)
dat2$val = dat2$val / sum(dat2$val, na.rm = TRUE)
dat3 = dat2
dat3$val = dat1$val - dat2$val
dat3$dat_type = "Pixel-wise difference"
p2 = plot_dat(dat3, standardize = FALSE, add_map = TRUE) +
  theme(legend.position="left") +
  theme(strip.text.x = element_text(size = rel(1.4))) +
  theme(legend.title=element_text(size=rel(1.2))) + 
  scale_fill_gradientn(colours = c("blue", "white", "red"),
                       limits = c(-0.003, 0.003),
                       name = "Difference")
plot(p2)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/clim-pixelwise-diff-1.png)<!-- -->

``` r
## ggsave(file.path(figdir, "climatology-diff.pdf"), width = 6, height=6)
```

``` r
## 3. Transport plot
p = plot_omd_ggplot(res, add_map = TRUE, classify_quantile = TRUE)
p = p + ggtitle(label="Optimal transport pattern", subtitle="")
p = p + theme(plot.title = element_text(hjust = 0.5, size = rel(3.5)))
p = p + theme(strip.text = element_text(size = rel(2.5), colour = "black", face = "bold"))
p = p + theme(axis.text = element_text(size = rel(2.5)),
              axis.title = element_text(size = rel(2.5)))
p = p + theme(strip.text.x = element_text(size = rel(2.5)))
plot(p)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/clim-transport-plot-1.png)<!-- -->

``` r
## ggsave(file.path(figdir, "climatology-omd-mid.pdf"), width = 12, height=12)
```

## Bottom row of Figure 2

Calculate distance matrix.

``` r
distmat = matrix(NA, 24, 24)
for(ii in 1:24){
  printprogress(ii, 24)
  for(jj in 1:24){
    if(ii < jj){

      ## Calculate OMD
      res = omd(M1_long = alldat_list[[ii]],
                M2_long = alldat_list[[jj]],
                p = 2,
                geodesic = TRUE)
                ## costm = costm)
      plot_omd_ggplot(res)
      distmat[ii,jj] = distmat[jj,ii] = res$dist
    }
  }
} 
diag(distmat) = NA
colnames(distmat) = rownames(distmat) = names(alldat_list)
saveRDS(distmat, file = file.path(datadir, "geo-distmat-clim-p2.RDS"))
```

Plot distance matrix.

``` r
## Load distance matrix
datadir = "~/repos/omd/main/data/01-climatology"
distmat = readRDS(file = file.path(datadir, "geo-distmat-clim-p2-overreact.RDS"))
source("01-helpers.R")
p = plot_distmat_ggplot_clim(distmat)
p = p + theme(legend.position = "left") 
plot(p)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/distmat-plot-1.png)<!-- -->

``` r
plotfilename = "distmat-climatology-overreact.pdf"
ggsave(file.path(figdir, plotfilename), width = 6, height=6)
```

Plotting MDS plot.

``` r
## Load distance matrix
outputdir = "~/repos/omd/main/data/04-cluster-boundary"
distmat = readRDS(file = file.path(datadir, "geo-distmat-clim-p2-overreact.RDS"))
p = make_mds_plot(distmat)
plot(p)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/mds-plot-1.png)<!-- -->

``` r
## ggsave(file = file.path(figdir, "mds-climatology-overreact.pdf"), width = 6, height = 5) 
```

## Figure 1 (Toy example)

We’ll use a larger box, just for illustration.

``` r
lat = 19.8968
lon = -155.5828
boxsize = 60
lonrange = lon + c(-1,1)*boxsize
latrange = lat + c(-1,1)*boxsize
lonrange = lon + c(-1,1)*boxsize
latrange = lat + c(-1,1)*boxsize
restrictbox <- . %>% filter(lat > latrange[1],
                            lat < latrange[2],
                            lon > lonrange[1],
                            lon < lonrange[2])
```

Toy example map.

``` r
ddat = read.csv(file.path(datadir_orig, "tabulated_darwin_montly_clim_090_180_ver_0_2_6.csv")) %>% as_tibble()
rdat = read.csv(file.path(datadir_orig, "tabulated_geospatial_montly_clim_090_180_ver_0_2_5.csv")) %>% as_tibble()
mo = 1
darwindat = ddat %>% dplyr::filter(month == mo) %>%  dplyr::select(month, lat, lon, val = Chl)  
darwindat = darwindat %>% restrictbox()

p1 = plot_dat(darwindat %>% 
              add_blob(c(-150 + (10*0 ), -25),
                       rotate = 0),
              add_map = TRUE)

p1 = p1 + labs(title = "Map 1")
p1 = p1 + theme(plot.title = element_text(hjust = 0.5))
plot(p1)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/toy-map-1.png)<!-- -->

``` r
## ggsave(file = file.path(figdir, "toy-map-orig-without-star.png"), width = 4, height=5)

p1 = p1 + geom_star(aes(x=lon, y=lat), data = data.frame(lat = 22.45, lon = -158), size=5,
                    fill = 'yellow', col = 'black') 
plot(p1)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/toy-map-2.png)<!-- -->

``` r
## ggsave(file = file.path(figdir, "toy-map-orig.png"), width = 4, height=5)


p2 = plot_dat(darwindat %>%
              add_blob(c(-150 + (10*0 ), -25),
                       rotate = 0) %>%
              add_blob(c(-150 + (10*2), -25),
                       rotate = (5 * 2 / 20) * 90) %>% 
              add_blob(c(-150 + (10*4), -25),
                       rotate = (5 * 4 / 20) * 90),
              add_map = TRUE)
arrowdat1 = data.frame(x=c(-150, -150+40), y=c(-25, -25) + .2)
arrowdat2 = data.frame(x=c(-150, -150+20), y=c(-25, -25) - .2)
p2 = p2 + geom_line(aes(x=x, y=y), arrow = arrow(length = unit(0.05*2, "inches")), data = arrowdat1, lwd = rel(2), col = rgb(0,0,0,0.6))
p2 = p2 + geom_line(aes(x=x, y=y), arrow = arrow(length = unit(0.05*2, "inches")), data = arrowdat2, lwd = rel(2), col = rgb(0,0,0,0.6))
p2 = p2 + labs(title = "Map 2")
p2 = p2 + theme(plot.title = element_text(hjust = 0.5))
plot(p2)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/toy-map-3.png)<!-- -->

``` r
## ggsave(file = file.path(figdir, "toy-map.png"), width = 4, height=5)
```

Toy distances.

``` r
start.time = Sys.time()
dists = mclapply(1:20, function(ii){
  printprogress(ii, 20, start.time = start.time)

  dat1 = darwindat %>% add_blob(c(-150, -25))
  dat2 = darwindat %>% add_blob(c(-150 + 2*ii, -25), rotate = (ii / 20) * 90)

  val1 = dat1 %>% pull(val)
  val2 = dat2 %>% pull(val)

  costm = form_cost_matrix(dat1)
  res = omd(M1_long = dat1,
            M2_long = dat2,
            p = 1,
            ## costm = costm,
            geodesic = TRUE)

  diffval = val1 - val2
  twodists = c(omd = res$dist, rmse = sqrt(sum((diffval)^2)))

  return(twodists)
}, mc.cores = 3)

distmat = do.call(rbind, dists) %>% as_tibble()  %>%
  tibble::add_column(shift = 2*(1:20)) %>%
  mutate(omd = omd-min(omd), rmse = rmse-min(rmse)) %>% 
  mutate(omd = omd/max(omd), rmse = rmse/max(rmse)) %>% 
  tidyr::pivot_longer(cols = c('omd', 'rmse')) 

saveRDS(distmat, file = file.path(datadir, "distmat-toy.RDS"))
```

Plotting two distances in one plot, as two lines.

``` r
distmat = readRDS(file = file.path(datadir, "distmat-toy.RDS"))

distmat = distmat %>% mutate(name = str_replace(name, "omd", "Wasserstein distance"))
distmat = distmat %>% mutate(name = str_replace(name, "rmse", "RMSE"))

p = ggplot(distmat) +
  geom_line(aes(x = shift, y = value, col = name), lwd = 2) +
  theme_minimal() +
  ylab("Distance, scaled to (0, 1)") +
  xlim(c(0,40)) +
  guides(col = guide_legend(title="Distance measure")) +
  theme(legend.position = c(0.7, 0.25)) +
  xlab("Longitude shift of patch")
plot(p)
```

![](/home/sangwonh/repos/omd/main/data/01-climatology/figures/toy-dists-1.png)<!-- -->

``` r
p = p + labs(title = "Comparing Map1 and Map2")
p = p + theme(plot.title = element_text(hjust = 0.5))
## ggsave(file = file.path(figdir, "toy-dists.png"), width = 4, height =5)
```
