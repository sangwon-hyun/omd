Lat/lon analysis of all months from 1998-2006 (non-climatology)
Chlorophyll data
================
Compiled at 2022-07-29 01:10:50 UTC

``` r
knitr::opts_chunk$set(fig.width=14, fig.height=8, echo=TRUE, eval=TRUE, cache=TRUE,
                      warning=FALSE, message=FALSE)

## Load packages
library(ggrepel, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(knitr, quietly = TRUE)
library(here, quietly = TRUE)
library(ggplot2, quietly = TRUE)
library(tidyverse, quietly = TRUE)
library(scales, quietly = TRUE)
library(reshape2, quietly = TRUE)
library(tidync, quietly = TRUE)
library(omd, quietly = TRUE)
```

    ## Warning: replacing previous import 'dplyr::union' by 'raster::union' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::select' by 'raster::select' when loading 'omd'

    ## Warning: replacing previous import 'dplyr::intersect' by 'raster::intersect' when loading 'omd'

``` r
sf::sf_use_s2(FALSE)
```

The R package `omd` to use is here
<https://github.com/sangwon-hyun/omd/>, commit ()\[\].

``` r
base = "02-nonclimatology"
here::i_am("02-nonclimatology.Rmd")
knitr::opts_chunk$set(fig.path = here::here("data", base, 'figures/'))
datadir = here::here("data", base)
datadir_orig = here::here("data", "00-import-data")
if(!dir.exists(datadir)) dir.create(datadir)
figdur = here::here("figures")
source(here::here("01-helpers.R"))
source(here::here("02-helpers.R"))
```

Following up on the climatology data analysis in `01-climatology.Rmd`,
we analyze nonclimatology data (monthly data between 1999-2006).

# Load data

These files are required: - `yearmo_mat_large.RDS` -
`nonclim_omd_objects_large_distmat.RDS` -
`nonclim_rmse_objects_large_distmat.RDS`

Plot the distance matrix for EMD and RMSE (**Supplemental Figure 1,
2**).

``` r
datadir = "~/repos/omd/main/data/02-nonclimatology" ## Delete when done

## yearmo_mat = readRDS(file = file.path(datadir, "yearmo_mat.RDS")) 
yearmo_mat = readRDS(file = file.path(datadir, "yearmo_mat_large.RDS"))
yearmo_mat = rbind(yearmo_mat %>% tibble::add_column(dat_type = "darwin"),
                   yearmo_mat %>% tibble::add_column(dat_type = "real"))
  

three_distmats_by_type = list()
for(dist_type in c("omd", "rmse")){

  ## Read distance matrix
  filename = paste0("nonclim_", dist_type, "_objects_large_distmat.RDS")
  distmat = readRDS(file = file.path(datadir, filename))

  ## Reorder distance matrix by year and month
  reordered_yearmo_mat = yearmo_mat %>% arrange(dat_type, year, mo) 
  new_order = reordered_yearmo_mat %>% left_join(yearmo_mat %>% add_column(orig_row = 1:nrow(yearmo_mat)),
                                                 by = c("dat_type", "year", "mo")) %>% pull(orig_row)
  distmat_reordered = distmat[new_order, new_order]

  yearmo_mat %>% filter(dat_type == 'darwin') %>% arrange(mo, year) %>% head(1)
  yearmo_mat %>% filter(dat_type == 'darwin') %>% arrange(mo, year) %>% tail(1)
  yearmo_mat %>% head()

  ## Delete the columns and rows that have NA (these are the ones for which dates /dont/ overlap.
  delete_row = which(distmat_reordered %>% apply(1, function(a) all(is.na(a))))
  delete_col = which(distmat_reordered %>% apply(2, function(a) all(is.na(a))))
  distmat_reordered = distmat_reordered[-delete_row, -delete_col]
  drawmat_precise(distmat_reordered)## temporary


  ## Isolate our attention to a certain date range, between 3-98 and 12-06.
  ind_d_begin = which(rownames(distmat_reordered) == "d3-98")
  ind_d_end = which(rownames(distmat_reordered) == "d12-6")
  ind_r_begin = which(rownames(distmat_reordered) == "r3-98")
  ind_r_end = which(rownames(distmat_reordered) == "r12-6")
  inds = c(c(ind_d_begin:ind_d_end), c(ind_r_begin:ind_r_end))
  distmat_reordered = distmat_reordered[inds, inds]

  ## Take subsets to form within-source distance matrices
  three_distmats = list()
  inds = which(substr(rownames(distmat_reordered), 1,1)=="r")
  three_distmats[["r"]] = distmat_reordered[inds,inds]
  inds = which(substr(rownames(distmat_reordered), 1,1)=="d")
  three_distmats[["d"]] = distmat_reordered[inds,inds]

  ## Also store the entire thing
  three_distmats[["all"]] = distmat_reordered
  three_distmats_by_type[[dist_type]] = three_distmats
}


## Better combined plot
for(dist_type in c("omd", "rmse")){
  three_distmats = three_distmats_by_type[[dist_type]]
  for(matname in c("r", "d", "all")){
    if(matname == "all"){
      p = plot_distmat_ggplot_nonclim(three_distmats[["all"]])## * 100)
      if(dist_type == "rmse") p = p + scale_fill_gradientn(name = "RMSE", colors = hcl.colors(20, "RdYlGn"))
    } else {
      p = drawmat_precise_ggplot_custom(three_distmats[[matname]], colours = c("skyblue", "blue", "red", "yellow"),
                                        hcuts = 12 * (0:100),
                                        vcuts = 12 * (1:100) - 1.9)
    }
    p = p + coord_fixed()
    if(dist_type == "omd") mytitle = paste("Distance matrix: ", "Wasserstein distance")
    if(dist_type == "rmse") mytitle = paste("Distance matrix: ", "RMSE")
    p = p + labs(title = mytitle)
    plot(p)
    plotfilename = paste0("nonclim-distmat-", dist_type, "-", matname, "-large-common-dates-overreact.pdf")
    ## ggsave(file = file.path(figdir, plotfilename), width = 10, height = 10)
  }
}
```

![](02-nonclimatology_files/figure-gfm/distance-matrix-1.png)<!-- -->![](02-nonclimatology_files/figure-gfm/distance-matrix-2.png)<!-- -->![](02-nonclimatology_files/figure-gfm/distance-matrix-3.png)<!-- -->![](02-nonclimatology_files/figure-gfm/distance-matrix-4.png)<!-- -->![](02-nonclimatology_files/figure-gfm/distance-matrix-5.png)<!-- -->![](02-nonclimatology_files/figure-gfm/distance-matrix-6.png)<!-- -->

Calculate metric MDS (**Supplemental Figure 3, 4**).

``` r
set.seed(1001)
for(dist_type in c("omd", "rmse")){

  three_distmats = three_distmats_by_type[[dist_type]]
  for(matname in c("r", "d", "all")){

    ## Combined
    p = make_nonclim_mds_plot(three_distmats[[matname]], simple = FALSE)
    dt = make_nonclim_mds_plot(three_distmats[[matname]], simple = FALSE, return_dt = TRUE)
    p = p + ggtitle(paste0(dist_type, ": ", matname)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", matname, "-large-common-dates-overreact.png")
    ## ggsave(file = file.path(figdir, plotfilename), width = 8, height = 8)

    ## Combined, all years have their own grey line
    p0 = make_nonclim_mds_plot(three_distmats[[matname]], simple = FALSE, separate_grey = TRUE)
    p0 = p0 + ggtitle(paste0(dist_type, ": ", matname)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p0)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", matname, "-allconnect-large-common-dates-overreact.png")
    ## ggsave(file = file.path(figdir, plotfilename), width = 20, height = 20)

    ## Facet by month
    month_names <- lapply(1:12, function(ii) month.abb[ii])
    names(month_names) = 1:12
    month_labeller <- function(variable, value){ return(month_names[value])  }
    p2 = p + geom_text_repel(aes(x = x, y = y, label = yyshort, size = yy),max.overlaps = Inf,  data = dt %>% mutate(yyshort = substr(yy, 3,4))) +
      facet_wrap(~moname, labeller = month_labeller) +
      scale_radius(range = c(3,6)) +
      theme(strip.text.x = element_text(size = rel(3)))
    p2 = p2 + ggtitle(paste0(dist_type, ": ", matname)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p2)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", matname, "-large-by-month-common-dates-overreact.png")
    ## ggsave(file = file.path(figdir, plotfilename), width = 25, height = 20)

    ## Facet by year
    p = make_nonclim_mds_plot(three_distmats[[matname]], simple = FALSE, separate_grey = TRUE)
    p3 = p + facet_wrap(~year)## + ggsave("rmse-nonclim-big-by-month.png", width = 10, height = 8) 
    p3 = p3 + geom_text_repel(aes(x = x, y = y, label = mo, size = mo),
                              max_overlap = Inf,
                              data = dt %>% mutate(mo = as.numeric(mo)))## %>% mutate(yyshort = substr(yy, 3,4))) 
    p3 = p3 + scale_radius(range = c(3,6))
    p3 = p3 + theme(strip.text.x = element_text(size = rel(3)))
    p3 = p3 + ggtitle(paste0(dist_type, ": ", matname)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p3)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", matname, "-large-by-year-common-dates-overreact.png")
    ## ggsave(file = file.path(figdir, plotfilename), width = 25, height = 20)
  }
}
```

![](02-nonclimatology_files/figure-gfm/plot-mds-1.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-2.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-3.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-4.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-5.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-6.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-7.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-8.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-9.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-10.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-11.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-12.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-13.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-14.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-15.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-16.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-17.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-18.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-19.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-20.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-21.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-22.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-23.png)<!-- -->![](02-nonclimatology_files/figure-gfm/plot-mds-24.png)<!-- -->

# (OMD) vs (number of months apart)

**Figure 3** of the paper.

``` r
for(dist_type in c("omd", "rmse")){

  three_distmats = three_distmats_by_type[[dist_type]]

  plist = lapply(1:2, function(ii){

      matname = c("r", "d")[ii]
      distmat_small = three_distmats[[matname]] 
      mytitle = c("Between Remote-sensing",  "Between Darwin")[ii]

      print(paste0("Results from lm(); ", mytitle))
      p = trendplot_advanced(distmat_small, dist_type, mytitle) 
      if(dist_type == "omd") p = p + ylab("Wasserstein distance (in km)")
      if(dist_type == "rmse") p = p + ylab("RMSE")
      return(p)
  })

  g = gridExtra::grid.arrange(plist[[1]], plist[[2]], nrow=2)

  ## Save two ways
  plotfilename = paste0("nonclim-trend", "-", dist_type, "-common-dates-overreact-new.pdf")
  ## ggsave(file = file.path(figdir, plotfilename), width = 8, height = 8, g)
}
```

    ## [1] "Results from lm(); Between Remote-sensing"
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##    Min     1Q Median     3Q    Max 
    ## -7.927 -2.585 -0.483  2.115 14.862 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 18.315240   0.084848 215.860  < 2e-16 ***
    ## months_diff  0.027412   0.001904  14.396  < 2e-16 ***
    ## month_sep1  -2.137044   0.154704 -13.814  < 2e-16 ***
    ## month_sep2  -1.623129   0.109565 -14.814  < 2e-16 ***
    ## month_sep3  -0.171135   0.109609  -1.561    0.119    
    ## month_sep4   0.752885   0.109606   6.869 7.17e-12 ***
    ## month_sep5   1.114571   0.109604  10.169  < 2e-16 ***
    ## month_sep6   1.078375   0.109603   9.839  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 3.529 on 5557 degrees of freedom
    ## Multiple R-squared:  0.126,  Adjusted R-squared:  0.1249 
    ## F-statistic: 114.4 on 7 and 5557 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Results from lm(); Between Darwin"
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -18.9476  -4.8015  -0.0974   5.0952  20.4269 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  28.541133   0.170330 167.564  < 2e-16 ***
    ## months_diff   0.003200   0.003823   0.837    0.403    
    ## month_sep1  -12.893958   0.310565 -41.518  < 2e-16 ***
    ## month_sep2   -7.193511   0.219949 -32.705  < 2e-16 ***
    ## month_sep3   -1.196867   0.220037  -5.439 5.57e-08 ***
    ## month_sep4    2.986260   0.220031  13.572  < 2e-16 ***
    ## month_sep5    5.637017   0.220028  25.620  < 2e-16 ***
    ## month_sep6    6.267182   0.220026  28.484  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.083 on 5557 degrees of freedom
    ## Multiple R-squared:  0.4301, Adjusted R-squared:  0.4294 
    ## F-statistic: 599.1 on 7 and 5557 DF,  p-value: < 2.2e-16

![](02-nonclimatology_files/figure-gfm/trend-1.png)<!-- -->

    ## [1] "Results from lm(); Between Remote-sensing"
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.034065 -0.009422  0.000620  0.009330  0.044583 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.178e-01  3.053e-04 385.931  < 2e-16 ***
    ## months_diff  8.557e-05  6.851e-06  12.491  < 2e-16 ***
    ## month_sep1  -1.421e-02  5.566e-04 -25.532  < 2e-16 ***
    ## month_sep2  -1.090e-02  3.942e-04 -27.661  < 2e-16 ***
    ## month_sep3  -2.226e-03  3.944e-04  -5.645 1.73e-08 ***
    ## month_sep4   2.999e-03  3.943e-04   7.605 3.32e-14 ***
    ## month_sep5   6.271e-03  3.943e-04  15.903  < 2e-16 ***
    ## month_sep6   8.617e-03  3.943e-04  21.851  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.0127 on 5557 degrees of freedom
    ## Multiple R-squared:  0.2961, Adjusted R-squared:  0.2953 
    ## F-statistic:   334 on 7 and 5557 DF,  p-value: < 2.2e-16
    ## 
    ## [1] "Results from lm(); Between Darwin"
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.085371 -0.020658 -0.004272  0.021366  0.111235 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.561e-01  7.152e-04 218.239  < 2e-16 ***
    ## months_diff  4.946e-05  1.605e-05   3.081  0.00207 ** 
    ## month_sep1  -3.395e-02  1.304e-03 -26.035  < 2e-16 ***
    ## month_sep2  -2.054e-02  9.236e-04 -22.244  < 2e-16 ***
    ## month_sep3  -9.799e-04  9.240e-04  -1.060  0.28897    
    ## month_sep4   9.758e-03  9.240e-04  10.561  < 2e-16 ***
    ## month_sep5   1.461e-02  9.239e-04  15.812  < 2e-16 ***
    ## month_sep6   1.564e-02  9.239e-04  16.932  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02975 on 5557 degrees of freedom
    ## Multiple R-squared:  0.2325, Adjusted R-squared:  0.2316 
    ## F-statistic: 240.5 on 7 and 5557 DF,  p-value: < 2.2e-16

![](02-nonclimatology_files/figure-gfm/trend-2.png)<!-- -->

# Longer term trend, only on remote sensing data

**Supplemental Figure 5** of the paper, showing longer term trend in
longer time range (1998-2020), only on remote sensing data:

``` r
## Get the years and monnths
yearmo_mat = readRDS(file = file.path(datadir, "yearmo_mat_large.RDS"))
yearmo_mat = rbind(yearmo_mat %>% tibble::add_column(dat_type = "darwin"),
                   yearmo_mat %>% tibble::add_column(dat_type = "real"))

## Read distance matrix
for(dist_type in c("omd", "rmse")){
  for(dat_type in c("d", "r")){

    filename = paste0("nonclim_", dist_type, "_objects_large_distmat.RDS")
    distmat = readRDS(file = file.path(datadir, filename))
    
    ## Reorder distance matrix by year and month
    reordered_yearmo_mat = yearmo_mat %>% arrange(dat_type, year, mo) 
    new_order = reordered_yearmo_mat %>% left_join(yearmo_mat %>% add_column(orig_row = 1:nrow(yearmo_mat)),
                                                   by = c("dat_type", "year", "mo")) %>% pull(orig_row)
    distmat_reordered = distmat[new_order, new_order]
    
    ## Delete the columns and rows that have NA (these are the ones for which dates /dont/ overlap.
    delete_row = which(distmat_reordered %>% apply(1, function(a) all(is.na(a))))
    delete_col = which(distmat_reordered %>% apply(2, function(a) all(is.na(a))))
    distmat_reordered = distmat_reordered[-delete_row, -delete_col]
    
    ## Take subsets to form within-source distance matrices
    inds = which(substr(rownames(distmat_reordered), 1,1) == dat_type)
    distmat_small = distmat_reordered[inds,inds]
    
    ## Make plot
    ## p = trendplot(distmat_small, "OMD", "Between Remote Sensing (1998-2020)", size = rel(.5))
    if(dat_type == "r") longname = "Remote Sensing"
    if(dat_type == "d") longname = "Darwin"

    p = trendplot_advanced(distmat_small, toupper(dist_type),
                           paste0("Between ", longname," (1998-2020)"),
                           size = rel(.5), c(0,25))
    if(dist_type == "omd") p = p + ylab("Wasserstein distance (in km)")
    if(dist_type == "rmse") p = p + ylab("RMSE")
    plotfilename = paste0("nonclim-trend-", dist_type,"-", dat_type, "-longer-time-range-new.png")
    ## ggsave(file = file.path(figdir, plotfilename), width = 10, height = 4)
  }
} 
```

    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -19.1416  -4.7570  -0.1106   5.0266  20.5238 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  28.581622   0.165474 172.725  < 2e-16 ***
    ## months_diff   0.002302   0.003646   0.631    0.528    
    ## month_sep1  -12.860518   0.301809 -42.611  < 2e-16 ***
    ## month_sep2   -7.243774   0.213860 -33.872  < 2e-16 ***
    ## month_sep3   -1.359668   0.213852  -6.358  2.2e-10 ***
    ## month_sep4    2.825860   0.213847  13.214  < 2e-16 ***
    ## month_sep5    5.606575   0.213844  26.218  < 2e-16 ***
    ## month_sep6    6.425183   0.213843  30.046  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 7.015 on 5770 degrees of freedom
    ## Multiple R-squared:  0.4381, Adjusted R-squared:  0.4374 
    ## F-statistic: 642.7 on 7 and 5770 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -10.1613  -3.1972  -0.7847   2.5588  22.0304 
    ## 
    ## Coefficients:
    ##              Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept) 19.851718   0.039191 506.533  < 2e-16 ***
    ## months_diff  0.008707   0.000343  25.387  < 2e-16 ***
    ## month_sep1  -2.404494   0.070047 -34.327  < 2e-16 ***
    ## month_sep2  -1.816646   0.051268 -35.434  < 2e-16 ***
    ## month_sep3  -0.416410   0.051271  -8.122 4.74e-16 ***
    ## month_sep4   0.579980   0.051271  11.312  < 2e-16 ***
    ## month_sep5   1.174330   0.051271  22.904  < 2e-16 ***
    ## month_sep6   1.438006   0.051271  28.047  < 2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 4.275 on 37393 degrees of freedom
    ## Multiple R-squared:  0.1051, Adjusted R-squared:  0.1049 
    ## F-statistic: 627.2 on 7 and 37393 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.085083 -0.020409 -0.004057  0.020928  0.111003 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.564e-01  6.945e-04 225.186   <2e-16 ***
    ## months_diff  3.894e-05  1.530e-05   2.544    0.011 *  
    ## month_sep1  -3.377e-02  1.267e-03 -26.658   <2e-16 ***
    ## month_sep2  -2.049e-02  8.976e-04 -22.831   <2e-16 ***
    ## month_sep3  -9.833e-04  8.976e-04  -1.096    0.273    
    ## month_sep4   9.897e-03  8.975e-04  11.026   <2e-16 ***
    ## month_sep5   1.460e-02  8.975e-04  16.262   <2e-16 ***
    ## month_sep6   1.543e-02  8.975e-04  17.188   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.02945 on 5770 degrees of freedom
    ## Multiple R-squared:  0.2343, Adjusted R-squared:  0.2334 
    ## F-statistic: 252.3 on 7 and 5770 DF,  p-value: < 2.2e-16
    ## 
    ## 
    ## Call:
    ## lm(formula = sqrt(omd) ~ months_diff + month_sep, data = dists, 
    ##     contrasts = list(month_sep = contr.sum))
    ## 
    ## Residuals:
    ##       Min        1Q    Median        3Q       Max 
    ## -0.042601 -0.010216 -0.000920  0.008511  0.095154 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)    
    ## (Intercept)  1.258e-01  1.390e-04  905.57   <2e-16 ***
    ## months_diff  1.467e-05  1.216e-06   12.06   <2e-16 ***
    ## month_sep1  -1.258e-02  2.484e-04  -50.64   <2e-16 ***
    ## month_sep2  -9.782e-03  1.818e-04  -53.81   <2e-16 ***
    ## month_sep3  -2.503e-03  1.818e-04  -13.77   <2e-16 ***
    ## month_sep4   2.407e-03  1.818e-04   13.24   <2e-16 ***
    ## month_sep5   5.683e-03  1.818e-04   31.26   <2e-16 ***
    ## month_sep6   7.838e-03  1.818e-04   43.11   <2e-16 ***
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Residual standard error: 0.01516 on 37393 degrees of freedom
    ## Multiple R-squared:  0.1894, Adjusted R-squared:  0.1893 
    ## F-statistic:  1249 on 7 and 37393 DF,  p-value: < 2.2e-16

**Supplemental Figure 6** showing mds plot in longer time range.

``` r
par(bg="white")
for(dist_type in c("omd", "rmse")){
  for(dat_type in c("r", "d")){

    filename = paste0("nonclim_", dist_type, "_objects_large_distmat.RDS") 
    distmat = readRDS(file = file.path(datadir, filename))
    
    ## Reorder distance matrix by year and month
    reordered_yearmo_mat = yearmo_mat %>% arrange(dat_type, year, mo) 
    new_order = reordered_yearmo_mat %>% left_join(yearmo_mat %>% add_column(orig_row = 1:nrow(yearmo_mat)),
                                                   by = c("dat_type", "year", "mo")) %>% pull(orig_row)
    distmat_reordered = distmat[new_order, new_order]
    
    ## Delete the columns and rows that have NA (these are the ones for which dates /dont/ overlap.
    delete_row = which(distmat_reordered %>% apply(1, function(a) all(is.na(a))))
    delete_col = which(distmat_reordered %>% apply(2, function(a) all(is.na(a))))
    distmat_reordered = distmat_reordered[-delete_row, -delete_col]
    
    ## Take subsets to form within-source distance matrices
    inds = which(substr(rownames(distmat_reordered), 1,1) == dat_type)
    distmat_small = distmat_reordered[inds,inds]

    ## Combined
    p = make_nonclim_mds_plot( distmat_small, simple = FALSE)
    dt = make_nonclim_mds_plot(distmat_small, simple = FALSE, return_dt = TRUE)
    p = p + ggtitle(paste0(dist_type, ": ", dat_type)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", dat_type, "-large--overreact-longterm.pdf")
    ## ggsave(file = file.path(figdir, plotfilename), width = 8, height = 8)

    ## Combined, all years have their own grey line
    p0 = make_nonclim_mds_plot(distmat_small, simple = FALSE, separate_grey = TRUE)
    p0 = p0 + ggtitle(paste0(dist_type, ": ", dat_type)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p0)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", dat_type, "-allconnect-large-overreact-longterm.pdf")
    ## ggsave(file = file.path(figdir, plotfilename), width = 20, height = 20)

    ## Facet by month
    p2 = make_nonclim_mds_plot( distmat_small, simple = TRUE) +
      geom_text_repel(aes(x = x, y = y, label = yyshort, size = yy),
                      max.overlaps = Inf,
                      data = dt %>% mutate(yyshort = substr(yy, 3,4))) +
      facet_wrap(~moname)+
      scale_radius(range = c(3,6)) +
      theme(strip.text.x = element_text(size = rel(3)))
    p2 = p2 + ggtitle(paste0(dist_type, ": ", dat_type)) + theme(plot.title = element_text(size = rel(3), face = "bold"))
    plot(p2)
    plotfilename = paste0("nonclim-mds-", dist_type, "-", dat_type, "-large-by-month-overreact-longterm.pdf")
    ## ggsave(file = file.path(figdir, plotfilename), width = 15, height = 15)
  }
}
```

![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-1.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-2.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-3.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-4.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-5.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-6.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-7.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-8.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-9.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-10.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-11.png)<!-- -->![](02-nonclimatology_files/figure-gfm/long-term-trend-mds-12.png)<!-- -->

# Detailed analysis of distance matrix

**Figure 5 of the paper.** Comparing Jan 1998 Darwin data with other
months.

``` r
two_distmats = sapply(c("omd", "rmse"), function(dist_type){ 

  ## Read distance matrix
  filename = paste0("nonclim_", dist_type, "_objects_large_distmat.RDS")
  distmat = readRDS(file = file.path(datadir, filename))

  ## Reorder distance matrix by year and month
  reordered_yearmo_mat = yearmo_mat %>% arrange(dat_type, year, mo) 
  new_order = reordered_yearmo_mat %>% left_join(yearmo_mat %>% add_column(orig_row = 1:nrow(yearmo_mat)),
                                                 by = c("dat_type", "year", "mo")) %>% pull(orig_row)
  distmat_reordered = distmat[new_order, new_order]

  ## Delete the columns and rows that have NA (these are the ones for which dates /dont/ overlap.
  delete_row = which(distmat_reordered %>% apply(1, function(a) all(is.na(a))))
  delete_col = which(distmat_reordered %>% apply(2, function(a) all(is.na(a))))
  distmat_reordered = distmat_reordered[-delete_row, -delete_col]

  inds = which(substr(rownames(distmat_reordered), 1,1)=="d")
  return(distmat_reordered[inds,inds])

}, simplify = FALSE, USE.NAMES = TRUE)
  
oo = two_distmats[["omd"]]
rr = two_distmats[["rmse"]]

joint = full_join(reshape2::melt(oo[1,,drop=FALSE], na.rm=FALSE, as.is=TRUE) %>% as_tibble() %>% rename(omd = value),
                  reshape2::melt(rr[1,,drop=FALSE], na.rm=FALSE, as.is=TRUE) %>% as_tibble() %>% rename(rmse = value), by = c("Var1", "Var2"))

joint = joint %>% mutate(from_source = substr(Var1, 1,1),
                         from_mo = substr(Var1,2,3) %>% str_remove("-") %>% as.integer(),
                         from_yy = substr(Var1, 3, 5) %>% str_remove("-") %>% as.integer())

joint = joint %>% mutate(to_source = substr(Var2, 1,1),
                         to_mo = substr(Var2,2,3) %>% str_remove("-") %>% as.integer(),
                         to_yy = substr(Var2, 4, 6) %>% str_remove("-") %>% as.integer())

joint = joint %>% mutate(from_yyyy = case_when(from_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(from_yyyy = from_yyyy + from_yy) %>%
  add_column(day = 1) %>% 
  mutate(from_date = paste0(from_yyyy, "-", from_mo, "-", day)) %>%
  mutate(from_date = lubridate::ymd(from_date))

joint = joint %>% mutate(to_yyyy = case_when(to_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(to_yyyy = to_yyyy + to_yy) %>%
  mutate(to_date = paste0(to_yyyy, "-", to_mo, "-", day)) %>%
  mutate(to_date = lubridate::ymd(to_date))


joint$omd = joint$omd / max(joint$omd, na.rm=TRUE)
joint$rmse = joint$rmse / max(joint$rmse, na.rm=TRUE)

joint = joint %>% add_column(ind = 1:nrow(joint))
joint = joint %>% rename("Wasserstein distance \n(rescaled)" = omd, "RMSE (rescaled)" = rmse)
p = joint %>% pivot_longer(cols = c("Wasserstein distance \n(rescaled)", "RMSE (rescaled)")) %>%
  ggplot(aes(x = to_date, y = value, group = name, col = name)) +
  ggtitle("Comparison of Darwin data of January 1998, to other months") +
  geom_vline(xintercept = c(lubridate::ymd("1998-01-01"),
                            lubridate::ymd("2002-04-01"),
                            lubridate::ymd("2002-08-01")),
             col = rgb(0,0,0,0.5), lwd = rel(.8), lty="dotted") +
  geom_text(aes(x=lubridate::ymd("1998-01-20"), y = 1, label = "(I)"), col = "black", size = rel(5)) + 
  geom_text(aes(x=lubridate::ymd("2002-04-20"), y = 1, label = "(II)"), col = "black", size = rel(5)) + 
  geom_text(aes(x=lubridate::ymd("2002-08-20"), y = 1, label = "(III)"), col = "black", size = rel(5)) +
  geom_point() + geom_line() +
  xlab("Date of comparison") + 
  ylab("Rescaled distance")  +
  theme_bw() + 
  scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y")  +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) 
plot(p) 
```

![](02-nonclimatology_files/figure-gfm/one-slice-1.png)<!-- -->

``` r
## ggsave(file.path(figdir, "omd-darwin-slice-overreact.pdf"), width = 10, height = 3.5)
```

How to explain this? Compare January to April or August.

``` r
## Helper to remove coastline
remove_landlock <- . %>%
  dplyr::mutate(land = mark_land(lon, lat, overreact = TRUE)) %>%
  dplyr::filter(land == FALSE)  %>%
  dplyr::select(-land)

mae = c()
for(mo2 in c(4,8)){
  year2 = 2002
  type2 = "darwin"
  filename2 = paste0(type2, "-nonclim-", year2, "-", mo2, ".RDS")
  datadir2 = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data/darwin-nonclim"
  dat2 = readRDS(file = file.path(datadir2, filename2))

  ## First year-month
  year1 = 1998##yearmo_mat %>% .[ii1,] %>% pull(year)
  mo1   = 1##yearmo_mat %>% .[ii1,] %>% pull(mo)
  type1 = "darwin"##yearmo_mat %>% .[ii1,] %>% pull(dat_type)
  filename1 = paste0(type1, "-nonclim-", year1, "-", mo1, ".RDS")
  datadir1 = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data/darwin-nonclim"
  dat1 = readRDS(file = file.path(datadir1, filename1))
 
  ## Make combined data
  pipeline = . %>% dplyr::select(lat, lon, val) %>% dplyr::group_by(lat, lon) %>%
    dplyr::summarize(val = mean(val, na.rm = TRUE)) %>% ungroup() ## %>% coarsen(fac=4)  %>% 
  ## dplyr::filter(!is.na(val))
  dat1 = dat1 %>% pipeline()
  dat2 = dat2 %>% pipeline()
  name1 = "A"
  if(mo2 == 4){ name2 = "B";  name2_other = "(II)"}
  if(mo2 == 8){ name2 = "B";  name2_other = "(III)"}
  twodat = omd::combine(dat1, dat2, name1 = name1, name2 = name2)
  mycombine <- function(dlist, name1, name2){
    omd::combine(dlist[[1]], dlist[[2]],
                 name1 = name1,
                 name2 = name2)
  }
  twodat = twodat %>% group_by(dat_type) %>% group_split() %>%
    purrr::map(.%>% dplyr::select(-dat_type) %>% omd::coarsen(fac = 4)) %>%
    mycombine(name1, name2) %>% ##%>% restrictbox_further() 
    remove_landlock()

  level_key = c("(I)", name2_other)
  names(level_key) = c("A", name2)

  p1 = twodat %>% 
    mutate(dat_type = recode(dat_type, !!!level_key))  %>% 
    plot_dat(add_map = TRUE, limits = c(0, 0.0131), hide_legend = FALSE)
  p1 = p1 + scale_fill_gradientn(colours = c("blue", "red", "yellow"),
                                 name = "Chlorophyll \n     PDF",
                                 limits = c(0, 0.0131))
  p1 = p1 + theme(legend.position="bottom") +
    theme(strip.text.x = element_text(size = rel(1.4))) +
    theme(legend.key.width = unit(2, "cm"))

  ## Calculate optimal transport
  res = omd(M1_long = twodat %>% dplyr::filter(dat_type == name1),
            M2_long = twodat %>% dplyr::filter(dat_type == name2),
            p = 2)
  p2 = plot_omd_ggplot(res, add_map = TRUE, arrow_range = c(1E-3, 2), sample_arrows=TRUE)

  ## Getting the range for the pixel-wise difference
  val1 = twodat %>% filter(dat_type == name1) %>% arrange(lon, lat) %>% pull(val)
  val2 = twodat %>% filter(dat_type == name2) %>% arrange(lon, lat) %>% pull(val)
  val1 = val1/sum(val1)
  val2 = val2/sum(val2)
  df = abs(val1 - val2)
  print(range(df, na.rm=TRUE) %>% round(5))
  mymax = max(df, na.rm=TRUE) %>% ceiling()
  mymax = max(0.00936, 0.00310)
  ## print(mymax)
  p3 = twodat %>% pivot_wider(values_from = "val", names_from = "dat_type") %>%
    mutate(A = A/sum(A, na.rm = TRUE)) %>%
    mutate(B = B/sum(B, na.rm = TRUE)) %>%
    mutate(val = abs(A - B)) %>% dplyr::select(lon, lat, val) %>% 
    plot_dat(add_map = TRUE, standardize=FALSE, hide_legend = FALSE)
  p3 = p3 + scale_fill_gradientn(colours = c("white", "red", "darkred"),
                                 name = "",
                                 limits = c(0, mymax)) 
  p3 = p3 + ggtitle("Pixel-wise absolute difference")
  p3 = p3 + theme(strip.text.x = element_blank())
  p3 = p3 + theme(legend.position="bottom") +
    theme(strip.text.x = element_text(size = rel(1.4))) +
    theme(legend.key.width = unit(1, "cm"))
  p3 + theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 14))    # Center title position and size

  plot(p3)

  g = gridExtra::grid.arrange(p1, p3, p2, ncol=3, widths=c(0.6, 0.34,  0.3)) 
  ## ggsave(file.path(figdir, paste0("jan-to-", tolower(month.abb[mo2]), "-darwin-overreact-new.pdf")), g, width=15, height=5)
}
```

    ## [1] 0.00000 0.00935

![](02-nonclimatology_files/figure-gfm/one-slice-explanation-1.png)<!-- -->![](02-nonclimatology_files/figure-gfm/one-slice-explanation-2.png)<!-- -->

    ## [1] 0.00000 0.00309

![](02-nonclimatology_files/figure-gfm/one-slice-explanation-3.png)<!-- -->![](02-nonclimatology_files/figure-gfm/one-slice-explanation-4.png)<!-- -->

Instead of comparing 1998 to 2002, how about comparing 1998 to 2004.

``` r
## Helper to remove coastline
remove_landlock <- . %>%
  dplyr::mutate(land = mark_land(lon, lat, overreact = TRUE)) %>%
  dplyr::filter(land == FALSE)  %>%
  dplyr::select(-land)

mae = c()
for(mo2 in c(4,8)){
  year2 = 2004
  type2 = "darwin"
  filename2 = paste0(type2, "-nonclim-", year2, "-", mo2, ".RDS")
  datadir2 = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data/darwin-nonclim"
  dat2 = readRDS(file = file.path(datadir2, filename2))

  ## First year-month
  year1 = 1998##yearmo_mat %>% .[ii1,] %>% pull(year)
  mo1   = 1##yearmo_mat %>% .[ii1,] %>% pull(mo)
  type1 = "darwin"##yearmo_mat %>% .[ii1,] %>% pull(dat_type)
  filename1 = paste0(type1, "-nonclim-", year1, "-", mo1, ".RDS")
  datadir1 = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data/darwin-nonclim"
  dat1 = readRDS(file = file.path(datadir1, filename1))
 
  ## Make combined data
  pipeline = . %>% dplyr::select(lat, lon, val) %>% dplyr::group_by(lat, lon) %>%
    dplyr::summarize(val = mean(val, na.rm = TRUE)) %>% ungroup() ## %>% coarsen(fac=4)  %>% 
  ## dplyr::filter(!is.na(val))
  dat1 = dat1 %>% pipeline()
  dat2 = dat2 %>% pipeline()
  name1 = "A"
  if(mo2 == 4){ name2 = "B";  name2_other = "(II)"}
  if(mo2 == 8){ name2 = "B";  name2_other = "(III)"}
  twodat = omd::combine(dat1, dat2, name1 = name1, name2 = name2)
  mycombine <- function(dlist, name1, name2){
    omd::combine(dlist[[1]], dlist[[2]],
                 name1 = name1,
                 name2 = name2)
  }
  twodat = twodat %>% group_by(dat_type) %>% group_split() %>%
    purrr::map(.%>% dplyr::select(-dat_type) %>% omd::coarsen(fac = 4)) %>%
    mycombine(name1, name2) %>% ##%>% restrictbox_further() 
    remove_landlock()

  level_key = c("(I)", name2_other)
  names(level_key) = c("A", name2)

  p1 = twodat %>% 
    mutate(dat_type = recode(dat_type, !!!level_key))  %>% 
    plot_dat(add_map = TRUE, limits = c(0, 0.0131), hide_legend = FALSE)
  p1 = p1 + scale_fill_gradientn(colours = c("blue", "red", "yellow"),
                                 name = "Chlorophyll \n     PDF",
                                 limits = c(0, 0.0131))

  ## This is new!!
  p1 = p1 + theme(legend.position="bottom") +
    theme(strip.text.x = element_text(size = rel(1.4))) +
    theme(legend.key.width = unit(2, "cm"))

  ## Calculate optimal transport
  res = omd(M1_long = twodat %>% dplyr::filter(dat_type == name1),
            M2_long = twodat %>% dplyr::filter(dat_type == name2),
            p = 2)
  p2 = plot_omd_ggplot(res, add_map = TRUE, arrow_range = c(1E-3, 2), sample_arrows=TRUE)

  ## Getting the range for the pixel-wise difference
  val1 = twodat %>% filter(dat_type == name1) %>% arrange(lon, lat) %>% pull(val)
  val2 = twodat %>% filter(dat_type == name2) %>% arrange(lon, lat) %>% pull(val)
  val1 = val1/sum(val1)
  val2 = val2/sum(val2)
  df = abs(val1 - val2)
  print(range(df, na.rm=TRUE) %>% round(5))
  mymax = max(df, na.rm=TRUE) %>% ceiling()
  mymax = max(0.00936, 0.00310)
  print(mymax)
  p3 = twodat %>% pivot_wider(values_from = "val", names_from = "dat_type") %>%
    mutate(A = A/sum(A, na.rm = TRUE)) %>%
    mutate(B = B/sum(B, na.rm = TRUE)) %>%
    mutate(val = abs(A - B)) %>% dplyr::select(lon, lat, val) %>% 
    plot_dat(add_map = TRUE, standardize=FALSE, hide_legend = FALSE)
  p3 = p3 + scale_fill_gradientn(colours = c("white", "red", "darkred"),
                                 name = "",
                                 limits = c(0, mymax)) 
  p3 = p3 + ggtitle("Pixel-wise absolute difference")
  p3 = p3 + theme(strip.text.x = element_blank())

  ## this is new!!
  p3 = p3 + theme(legend.position="bottom") +
    theme(strip.text.x = element_text(size = rel(1.4))) +
    theme(legend.key.width = unit(1, "cm"))
  p3 + theme(plot.title = element_text(hjust = 0.5, vjust = 0.01, size = 14))    # Center title position and size

  plot(p3)

  g = gridExtra::grid.arrange(p1, p3, p2, ncol=3, widths=c(0.6, 0.34,  0.3)) 
  ## ggsave(file.path(figdir, paste0("2004-jan-to-", tolower(month.abb[mo2]), "-darwin-overreact-new.pdf")), g, width=15, height=5)
}
```

    ## [1] 0.00000 0.00757
    ## [1] 0.00936

![](02-nonclimatology_files/figure-gfm/alternate-one-slice-explanation-1.png)<!-- -->![](02-nonclimatology_files/figure-gfm/alternate-one-slice-explanation-2.png)<!-- -->

    ## [1] 0.00000 0.00297
    ## [1] 0.00936

![](02-nonclimatology_files/figure-gfm/alternate-one-slice-explanation-3.png)<!-- -->![](02-nonclimatology_files/figure-gfm/alternate-one-slice-explanation-4.png)<!-- -->

``` r
two_distmats = sapply(c("omd", "rmse"), function(dist_type){ 

  ## Read distance matrix
  filename = paste0("nonclim_", dist_type, "_objects_large_distmat.RDS")
  distmat = readRDS(file = file.path(datadir, filename))

  ## Reorder distance matrix by year and month
  reordered_yearmo_mat = yearmo_mat %>% arrange(dat_type, year, mo) 
  new_order = reordered_yearmo_mat %>% left_join(yearmo_mat %>% add_column(orig_row = 1:nrow(yearmo_mat)),
                                                 by = c("dat_type", "year", "mo")) %>% pull(orig_row)
  distmat_reordered = distmat[new_order, new_order]

  ## Delete the columns and rows that have NA (these are the ones for which dates /dont/ overlap.
  delete_row = which(distmat_reordered %>% apply(1, function(a) all(is.na(a))))
  delete_col = which(distmat_reordered %>% apply(2, function(a) all(is.na(a))))
  distmat_reordered = distmat_reordered[-delete_row, -delete_col]

  inds = which(substr(rownames(distmat_reordered), 1,1)=="d")
  return(distmat_reordered[inds,inds])

}, simplify = FALSE, USE.NAMES = TRUE)
  
oo = two_distmats[["omd"]]
rr = two_distmats[["rmse"]]

joint = full_join(reshape2::melt(oo[1,,drop=FALSE], na.rm=FALSE, as.is=TRUE) %>% as_tibble() %>% rename(omd = value),
                  reshape2::melt(rr[1,,drop=FALSE], na.rm=FALSE, as.is=TRUE) %>% as_tibble() %>% rename(rmse = value), by = c("Var1", "Var2"))

joint = joint %>% mutate(from_source = substr(Var1, 1,1),
                         from_mo = substr(Var1,2,3) %>% str_remove("-") %>% as.integer(),
                         from_yy = substr(Var1, 3, 5) %>% str_remove("-") %>% as.integer())

joint = joint %>% mutate(to_source = substr(Var2, 1,1),
                         to_mo = substr(Var2,2,3) %>% str_remove("-") %>% as.integer(),
                         to_yy = substr(Var2, 4, 6) %>% str_remove("-") %>% as.integer())

joint = joint %>% mutate(from_yyyy = case_when(from_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(from_yyyy = from_yyyy + from_yy) %>%
  add_column(day = 1) %>% 
  mutate(from_date = paste0(from_yyyy, "-", from_mo, "-", day)) %>%
  mutate(from_date = lubridate::ymd(from_date))

joint = joint %>% mutate(to_yyyy = case_when(to_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(to_yyyy = to_yyyy + to_yy) %>%
  mutate(to_date = paste0(to_yyyy, "-", to_mo, "-", day)) %>%
  mutate(to_date = lubridate::ymd(to_date))


joint$omd = joint$omd / max(joint$omd, na.rm=TRUE)
joint$rmse = joint$rmse / max(joint$rmse, na.rm=TRUE)

joint = joint %>% add_column(ind = 1:nrow(joint))
joint = joint %>% rename("Wasserstein distance \n(rescaled)" = omd, "RMSE (rescaled)" = rmse)
p = joint %>% pivot_longer(cols = c("Wasserstein distance \n(rescaled)", "RMSE (rescaled)")) %>%
  ggplot(aes(x = to_date, y = value, group = name, col = name)) +
  ggtitle("Comparison of Darwin data of January 1998, to other months") +
  geom_vline(xintercept = c(lubridate::ymd("1998-01-01"),
                            lubridate::ymd("2004-04-01"),
                            lubridate::ymd("2004-08-01")),
             col = rgb(0,0,0,0.5), lwd = rel(.8), lty="dotted") +
  geom_text(aes(x=lubridate::ymd("1998-01-20"), y = 1, label = "(I)"), col = "black", size = rel(5)) + 
  geom_text(aes(x=lubridate::ymd("2004-04-20"), y = 1, label = "(II)"), col = "black", size = rel(5)) + 
  geom_text(aes(x=lubridate::ymd("2004-08-20"), y = 1, label = "(III)"), col = "black", size = rel(5)) +
  geom_point() + geom_line() +
  xlab("Date of comparison") + 
  ylab("Rescaled distance")  +
  theme_bw() + 
  scale_x_date(date_breaks = "6 month", date_labels =  "%b %Y")  +
  theme(legend.title = element_blank()) +
  theme(axis.text.x=element_text(angle=90, hjust=1)) 
plot(p) 
```

![](02-nonclimatology_files/figure-gfm/alternative-one-slice-1.png)<!-- -->

``` r
## ggsave(file.path(figdir, "2004-omd-darwin-slice-overreact.pdf"), width = 10, height = 3.5)
```

What is causing this difference between the two distances? Compare
January to (1) April or (2) August.
