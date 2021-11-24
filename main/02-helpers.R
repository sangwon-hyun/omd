##' @param type Data type
##' @param datadir Directory of data
##' @param time_range Range of time (vector of two dates)
##' @param latrange Latitude range (vector of two numbers)
##' @param lonrange Longitude range (vector of two numbers)
##' @param verbose If TRUE, prints loudly
get_nonclim_dat <- function(type = c("real", "darwin"),
                            datadir = "/home/sangwonh/Dropbox/research/usc/ocean-provinces/data",
                            time_range = NULL,
                            latrange,
                            lonrange,
                            verbose = TRUE){


  ## Restrict data to a box
  restrictbox <- . %>% filter(lat >= latrange[1],
                              lat <= latrange[2],
                              lon >= lonrange[1],
                              lon <= lonrange[2])

  ## Setup
  library(ncdf4)
  type = match.arg(type)

  if(type=="real"){
    datfile = file.path(datadir, "Justin_pacific_box_halfdegree.nc")
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_time.csv")
  }
  if (type == "darwin"){
    datfile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin.nc")
    timefile = file.path(datadir, "Justin_pacific_box_halfdegree_Darwin_time.csv")
  }


  ## Read the data
  dat = tidync(datfile)
  nc_data <- nc_open(datfile)
  all_times = ncvar_get(nc_data, "time") %>% unique() %>% sort()
  if(is.null(time_range)){ times = all_times }
  if(!is.null(time_range)){ times = all_times[time_range] }

  ## Subset data into time frame of interest
  if(verbose){
    print("handling between")
    print(all_times[min(time_range)])
    print("and")
    print(all_times[max(time_range)])
  }
  dat = dat %>% hyper_filter(time = time %in% all_times[time_range])
  dat = dat %>% hyper_tibble()
  if(verbose){
    print("Size of data is")
    object.size(dat) %>% format("Mb") %>% print()
  }
  ## stopifnot((times %>% length()) ==
  ##           (read.csv(timefile) %>% .[,2] %>% length()))

  ## Get time strings from external file
  realtimes = read.csv(timefile) %>% .[,2] %>% .[time_range]
  conversion = tibble(time = all_times[time_range], realtime = realtimes)
  dat = dat %>% dplyr::filter(time %in% (conversion$time))
  ## object.size(dat) %>% format("Mb") %>% print()

  ## Add proper times
  dat = dat %>%
    dplyr::full_join(conversion, "time") %>%
    dplyr::filter(!is.na(realtime)) %>%
    dplyr::filter(!is.na(lat)) %>%
    dplyr::filter(!is.na(lon)) %>%
    dplyr::select(lat, lon, val = chl, time = realtime) %>%
    dplyr::mutate(time = lubridate::as_date(time)) %>%
    dplyr::mutate(day = lubridate::day(time),
                  mo = lubridate::month(time),
                  year = lubridate::year(time)) %>%
    ## remove_landlock() %>%
    dplyr::mutate(land = mark_land(lon, lat, overreact = TRUE)) %>%
    dplyr::filter(land == FALSE)  %>%
    dplyr::select(-land) %>%
    restrictbox()

  ## yy = 1998
  ## mm = 9
  ## this_month_dat = dat %>%
  ##   dplyr::filter(mo == mm) %>%
  ##   dplyr::filter(year == yy) %>%
  ##   group_by(lat, lon) %>%
  ##   summarize(val = mean(val, na.rm=TRUE))  %>%
  ##   ungroup() %>%
  ##   coarsen(4)
  ## this_month_dat %>% plot_dat(add_map = TRUE, hide_legend = FALSE)

  ## dplyr::mutate(val = pmin(val, 1))
  return(dat)
}





drawmat_precise_ggplot_custom <- function(distmat, colours = c("blue", "red"),
                                   limits = NULL,
                                   xlab = "", ylab = "", title = "",
                                   hcuts = NULL, vcuts = NULL){
  distmat = as.matrix(distmat)
  longData <- reshape2::melt(distmat[nrow(distmat):1,])##, varnames=c('Var1', 'Var2'))
  ## longData<-longData[longData$value!=0,]

  p = ggplot(longData, aes(x = Var2, y = Var1)) +
    geom_raster(aes(fill=value)) +
    ## scale_fill_gradient(low="grey90", high="red") +
    ## scale_fill_gradient(low="grey90", high="red") +
    scale_fill_gradientn(colours = colours, guide="colorbar", limits=limits) +
    labs(x = xlab, y = ylab, title = title) +
    theme_minimal() + theme(axis.text.x = element_text(size = rel(1), angle = 90, vjust = 0.3),
                            axis.text.y = element_text(size = rel(1)),
                            plot.title = element_text(size = rel(1.5)))

  ## Add cuts, if necessary.
  if(!is.null(vcuts)){
    p = p + geom_vline(xintercept = vcuts + 0.5, lwd = rel(1), col = rgb(0,0,0,0.5))
  }
  if(!is.null(hcuts)){
    p = p + geom_hline(yintercept = hcuts + 0.5, lwd = rel(1), col = rgb(0,0,0,0.5))
  }
  return(p)
}







make_nonclim_mds_plot <- function(distmat, simple = FALSE, separate_grey = FALSE, return_dt = FALSE){

  ## Prepare
  dmat <- as.dist(distmat, diag = TRUE)
  fit <- cmdscale(dmat, eig = TRUE, k = 2) # k is the number of dim
  x <- fit$points[,1]
  y <- fit$points[,2]
  nm = rownames(distmat)##[-(49:50)]
  nm_short = nm %>% sapply(., substr, 1,3) %>% stringr::str_split("-")  %>% map(. %>% .[1]) %>% unlist()
  dat_type = nm %>% sapply(.,substr, 0,1)
  moyy = nm %>% sapply(.,substr, 2,10) %>% stringr::str_split("-") %>% purrr::map(.%>%as.numeric()) %>%
    do.call(rbind,.) %>%
    as_tibble() %>%
    dplyr::select(mo = 1, yy = 2)
  mo = moyy %>% pull(mo)
  yy = moyy %>% pull(yy)
  yy = ifelse(yy>50, yy + 1900, yy + 2000)
  dt = data.frame(x = x, y = y,
             labels = nm_short,
             dat_type = dat_type,
             year = factor(yy, levels=unique(yy)),
             mo = mo)

  get_rotation_matrix <- function(aa){
    cbind(c(cos(aa), sin(aa)),
          c(-sin(aa), cos(aa)))
  }
  rot = get_rotation_matrix(-pi * 0.2)

  dt = data.frame(x = cbind(x,y) %*% rot %>% .[,1],
                  y = cbind(x,y) %*% rot %>% .[,2],
             labels = nm_short,
             dat_type = dat_type,
             dat_type_label = ifelse(dat_type=="r", "Remote Sensing", "Darwin model"),
             year = factor(yy, levels=unique(yy)),
             yy = yy,
             mo = mo,
             moname = month.abb[mo]) %>%
    mutate(mo = as.factor(mo),
           dat_type = as.factor(dat_type),
           dat_type_label = as.factor(dat_type_label),
           moname = factor(moname, levels = month.abb)) %>%
    as_tibble()

  ## dt = dt %>% mutate(group = interaction(dat_type, year)) %>% mutate(moname = mo)
  if(return_dt) return(dt)


  ## dt$labels = dt$labels %>% sapply(substr, 1,1)
  p = ggplot() +
    geom_point(aes(x = x, y = y, shape = dat_type_label, col = moname),  size = rel(3), data = dt) +
    xlab("Coordinate 1")+
    ylab("Coordinate 2") +
    scale_color_manual(values = RColorBrewer::brewer.pal(12, "Paired")) +
    theme_minimal() +
    guides(col = guide_legend(title = "Month of year"))  +
    guides(shape = guide_legend(title = "Data source"))  +
    coord_fixed()
    ## geom_text(aes(x=-5, y= -7, label = "Darwin model data", fontface = "plain"), cex=rel(4)) +
    ## geom_text(aes(x=0, y= 5, label = "Remote sensing data", fontface = "plain"), cex=rel(4))
    ## theme(legend.position = "hide")
  if(simple){
    return(p)
  }

  if(separate_grey){
    dt2 = dt %>% group_by(mo, dat_type_label)
    dt2 = dt2 %>% mutate(group = interaction(dat_type_label)) %>% mutate(moname = mo)
    dt2_extra = dt2 %>% dplyr::filter(mo == 1)

    p = p + geom_path(aes(x=x, y=y, group = group), col = rgb(0,0,0,0.2),
                  data = rbind(dt2, dt2_extra))

  } else {
    dt2 = dt %>% group_by(mo, dat_type_label) %>% summarise(x = mean(x), y = mean(y))
    dt2 = dt2 %>% mutate(group = interaction(dat_type_label)) %>% mutate(moname = mo)
    dt2_extra = dt2 %>% dplyr::filter(mo == 1)

    p = p + geom_path(aes(x=x, y=y, group = group), col = rgb(0,0,0,0.2),
                    data = rbind(dt2, dt2_extra))

    ## p + geom_point(aes(x=x, y=y, group = group, col = mo), data = rbind(dt2))##, dt2_extra)) ## In case we want to customize even more
  }

  return(p)
}





##' Make trend plot of (distances vs number of months apart)
##'
##' @param distmat_small Distance matrix
##' @param dist_type Use as y axis.
##' @param mytitle Title of plot
##'
##' @return ggplot object.
trendplot <- function(distmat_small, dist_type = "", mytitle = "", limits = c(NA, 25), size = rel(1)){
  tab = lapply(rownames(distmat_small), function(src){
    dat_type = src %>% substr(1,1)
    yearmo = src %>% substr(2,nchar(src)) %>% str_split("-") %>% unlist()
    year = yearmo[2] %>% as.numeric()
    mo = yearmo[1] %>% as.numeric()
    tibble(dat_type = dat_type, year = year, mo = mo)}) %>% bind_rows()

  nr0 = nrow(distmat_small)
  info_orig = lapply(1:nr0, function(jj){
    omds = distmat_small[,jj]
    time_apart = abs((1:nr0) - jj)
    tab %>% add_column(omd = omds) %>%
      add_column(x = time_apart) %>%
      add_column(src=rownames(distmat_small)[jj])
  }) %>% bind_rows()
  info = info_orig
  info$year = info$year %>% sapply(function(yr){if(yr>90) yr = yr + 1900 else yr = yr + 2000})
  info_mean = info %>% group_by(x) %>% summarise(omd = median(omd, na.rm=TRUE)) %>% ungroup()

  number_ticks <- function(n) {function(limits) pretty(limits, n)}

  p = info %>%
    ggplot() +
    geom_point(aes(x=x, y=omd, group = src), size = rel(1), col=rgb(0,0,0,0.1)) +
    geom_smooth(aes(x=x, y=omd), method="lm") +
    geom_line(aes(x=x, y=omd), info_mean, lwd = 2, col = rgb(1,0,0,0.5))  +
    theme_bw() +
    ylab(dist_type %>% toupper())  +
    xlab("Number of months apart")  +
    ggtitle(mytitle) +
    scale_x_continuous(breaks = seq(from=0,to=1000, by = 12))

  base_breaks <- function(n = 10){
    function(x) {
      axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
  }

  p = p +
    scale_y_continuous(trans = scales::log_trans(), breaks = base_breaks(), limits = limits) +
    theme(panel.grid.minor = element_blank())
  return(p)
}




##' Make trend plot of (distances vs number of months apart)
##'
##' @param distmat_small Distance matrix
##' @param dist_type Use as y axis.
##' @param mytitle Title of plot
##'
##' @return ggplot object.
trendplot_advanced <- function(distmat_small, ylab = "", mytitle = "", limits = c(NA, 25), size = rel(1)){


  ## Reformat into a long matrix
  dists = distmat_small %>% reshape2::melt() %>% as_tibble()
  dists = dists %>% mutate(from_source = substr(Var1, 1,1),
                           from_mo = substr(Var1, 2, 3) %>% str_remove("-") %>% as.integer(),
                           from_yy = substr(Var1, 4, 6) %>% str_remove("-") %>% as.integer())
  dists = dists %>% mutate(to_source = substr(Var2, 1,1),
                           to_mo = substr(Var2, 2, 3) %>% str_remove("-") %>% as.integer(),
                           to_yy = substr(Var2, 4, 6) %>% str_remove("-") %>% as.integer())
  dists = dists %>% mutate(from_yyyy = case_when(from_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(from_yyyy = from_yyyy + from_yy)
  dists = dists %>% mutate(to_yyyy = case_when(to_yy < 90 ~ 2000, TRUE ~ 1900)) %>% mutate(to_yyyy = to_yyyy + to_yy)
  dists = dists %>% dplyr::select(-Var1, -Var2, -from_yy, -to_yy)
  dists = dists %>% na.omit()

  dists = dists %>% add_column(day = 1) %>%
    mutate(from_date = paste0(from_yyyy, "-", from_mo, "-", day)) %>%
    mutate(from_date = lubridate::ymd(from_date)) %>%
    mutate(to_date = paste0(to_yyyy, "-", to_mo, "-", day)) %>%
    mutate(to_date = lubridate::ymd(to_date))

  dists = dists %>% filter(to_date > from_date)

  dists = dists %>% mutate(from_date = from_yyyy + (from_mo - 1) / 12,
                           to_date = to_yyyy + (to_mo - 1) / 12,
                           months_diff = round(12 * (to_date - from_date)),
                                        #         mod_diff = as.factor(months_diff %% 12),
                           month_sep = as.factor(pmin(12 - months_diff %% 12, months_diff %% 12)),
                           pair = paste(round(from_date, 3), round(to_date, 3))) %>%
    rename(omd = value)

  ## Fit a model
  fit <- lm(sqrt(omd) ~ months_diff + month_sep,
            contrasts = list(month_sep = contr.sum),
            data = dists)
  print(summary(fit))

  ## Make the predictions at each |i_j|
  pred <- dists %>%
    distinct(months_diff, month_sep) %>%
    mutate(as_tibble(predict(fit, data.frame(months_diff, month_sep),
                             interval = "prediction")))

  ## Linear prediction.
  b1 = fit %>% coef() %>% .["months_diff"]
  b0 = fit %>% coef() %>% .["(Intercept)"]
  modiff = dists$months_diff  %>% unique() %>% sort()
  linpred = (b0 + modiff * b1)^2
  df = tibble(linpred = linpred, x = modiff)

  ## Add distance.
  base_breaks <- function(n = 10){
    function(x) {
      axisTicks(log10(range(x, na.rm = TRUE)), log = TRUE, n = n)
    }
  }

  p = dists %>%
    ggplot(aes(x = months_diff, y = omd)) +
    geom_point(alpha = 0.1, size = size) +
    geom_line(data = pred, aes(x = months_diff, y = fit^2), color = "blue", lwd = 1) +
    ## geom_ribbon(data = pred, aes(x = months_diff, ymin = lwr^2, ymax = upr^2, y = fit^2),
    ##             alpha = 0.2)  +
    geom_line(aes(y=linpred, x=x), data = df, col = 'red', lwd = 1) +
    theme_bw() +
    ylab(ylab) +
    xlab("Number of months apart |ym(a)-ym(b)|")  +
    scale_x_continuous(breaks = seq(from=0,to=1000, by = 12)) +
    ## scale_y_continuous(trans = scales::log_trans(), breaks = base_breaks()) +
    scale_y_sqrt() +
    theme(panel.grid.minor = element_blank()) +
    ggtitle(mytitle)
    ## ylim(limits)
  return(p)
}
