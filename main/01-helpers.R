##' Helper for the hierarchical clustering.
labelCol <- function(x) {
  if (is.leaf(x)) {
    ## fetch label
    label <- attr(x, "label")
    ## set label color to red for A and B, to blue otherwise
    attr(x, "nodePar") <- list(lab.col=ifelse(sapply(label, substr, 1,4) == "real", "blue", "red"))
  }
  return(x)
}



## ##' Helper for making MDS plot.
## mds_plot <- function(distmat){
##   ## Todo: make the labels transparent, if possible.
##   fit <- cmdscale(distmat, eig=TRUE, k=2) # k is the number of dim

##   # plot solution
##   x <- fit$points[,1]
##   y <- fit$points[,2]
##   data.frame(x=x, y=y, labels=c(paste0(1:12),paste0(1:12)), dat_type=c(rep("Real", 12), rep("Darwin", 12))) %>%
##     ggplot() +
##     geom_label(aes(x=x, y=y, label=labels, col=dat_type), cex=5, fontface = "bold") +
##     xlab("Coordinate 1")+
##     ylab("Coordinate 2") +
##     ggtitle("Metric MDS") +
##     scale_color_manual(values=c("black","blue"))
## }



make_mds_plot <- function(distmat,
                          shortnames = c(paste0("real-", 1:12), paste0("darwin-", 1:12)),
                          longnames  = c(rep("Remote sensing", 12), rep("Darwin model", 12)),
                          angle = -pi * 0.2
                          ){

  ## Basic checks

  diag(distmat) = NA
  colnames(distmat) = rownames(distmat) = shortnames

  ## Turn into "dist" class object
  distmat <- as.dist(distmat, diag = TRUE)

  ## Calculating an MDS plot
  fit <- cmdscale(distmat, eig=TRUE, k=2) # k is the number of dim

  # Plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]

  ## Define rotation
  get_rotation_matrix = function(aa){
    cbind(c(cos(aa), sin(aa)),
        c(-sin(aa), cos(aa)))
  }
  rot = get_rotation_matrix(angle)

  ## Plot the points on a 2d plot
  dt = data.frame(x = x, y = y)
  dt = as.matrix(dt) %*% rot
  dt = data.frame(x= dt[,1], y = dt[,2],
                  labels = c(paste0(month.abb[1:12]),
                             paste0(month.abb[1:12])),
                  dat_type = longnames) %>%
    mutate(labels = factor(labels, levels = month.abb[1:12])) %>%
    arrange(labels)

  some_shortnames = shortnames[c(12, 1, 24, 13)]
  ## some_shortnames = c("real-12", "real-1", "darwin-12", "darwin-1"),])
  dt_extra = rbind(dt[some_shortnames,])

  p = dt %>% ggplot(aes(x = x, y = y, label = labels, col = dat_type)) +
    geom_point(cex=2) +
    geom_path(, alpha=.5, lty=2) +
    geom_path(aes(x = x, y = y, label = labels, col = dat_type), alpha=.5, lty=2, data = dt_extra) +
    geom_text_repel(cex = rel(5), ##fontface = "bold",
                     alpha = .8, show.legend = F) +
                     ## label.padding=.05) +
    labs(col = "Data source") +
    theme(legend.title= element_blank()) +
    theme_minimal() +
    xlab("Coordinate 1") +
    ylab("Coordinate 2") +
    ## ggtitle("Multidimensional scaling of climatology data (1998-2018)") +
    scale_color_manual(values=c("black","blue")) +
    ## scale_color_manual(values=c("black","blue")) +
    coord_fixed()  +
    theme(legend.text = element_text(size=rel(1), face=1))
  p = p + theme(legend.justification = c(1, 0), legend.position = c(1, 0))
  return(p)
}



plot_distmat_ggplot_clim <- function(distmat){

  ## Helper
  mold_distmat <- function(distmat){

    rownames(distmat) = rep(1:12, 2)
    colnames(distmat) = rep(1:12, 2)

    OMD_mat = list(remote_remote = distmat[1:12, 1:12],
                   remote_darwin = distmat[1:12, 13:24],
                   darwin_remote = distmat[13:24, 1:12],
                   darwin_darwin = distmat[13:24, 13:24])

    temp2 <- lapply(OMD_mat, function(x) {
      x <- x %>% as_tibble();
      x$from <- colnames(x)
      x <- tidyr::gather(x, key = 'to', value = 'OMD', -'from' )
      x
      })
    df_plot2 <- as_tibble()
    for (i in 1:length(temp2)){
      x <- temp2[[i]]
      temx <- strsplit(names(temp2)[i],'_')[[1]]
      x$type1 <- temx[1]
      x$type2 <- temx[2]
      df_plot2 <- bind_rows(list(df_plot2, x))
    }

    df_plot2$type1 <- factor(df_plot2$type1)
    df_plot2$type2 <- factor(df_plot2$type2)##, levels = c('obs','dar') )
    levels(df_plot2$type1) <- c('Darwin Model', 'Remote Sensing')
    levels(df_plot2$type2) <- c('Darwin Model', 'Remote Sensing')

    df_plot2$from = month.abb[as.numeric(df_plot2$from)]
    df_plot2$to   = month.abb[as.numeric(df_plot2$to)]

    df_plot2$from <- factor(df_plot2$from, levels = month.abb[1:12])
    df_plot2$to <-  factor(df_plot2$to, levels = month.abb[12:1])

    return(df_plot2)
  }

  ## Make data frame
  df_plot2 = mold_distmat(distmat)

  ## Plot data frame as heatmap with facets
  p = ggplot(df_plot2 , aes(from, to, fill= OMD)) +
    facet_grid(type2 ~ type1) +
    geom_tile() +
    scale_fill_gradientn(name = "Wasserstein \ndistance (in km)", colors = hcl.colors(20, "RdYlGn")) +
    theme_bw() + xlab("Month") + ylab('Month') +
    coord_fixed() +
    theme(axis.text.x = element_text(size = rel(1), angle = 90, vjust = 0.3)) +
    theme(legend.position="top")
  return(p)
}





plot_distmat_ggplot_nonclim <- function(distmat){

  ## distmat = distmat_reordered
  nr = nrow(distmat)
  moyr = substr(rownames(distmat), 2, 20)
  rownames(distmat) = colnames(distmat) = moyr

  ## Extra handling b/c the first two months are missing

  OMD_mat = list(remote_remote = distmat[(nr/2 + 1):nr, (nr/2 + 1):nr],
                 darwin_remote = distmat[1:(nr/2), (nr/2 + 1):nr],
                 remote_darwin = distmat[(nr/2 + 1):nr, 1:(nr/2)],
                 darwin_darwin = distmat[1:(nr/2), 1:(nr/2)])

  ## browser()
  ## distmat[1:(nr/2), (nr/2 + 1):nr] %>% rownames()
  ## distmat[1:(nr/2), (nr/2 + 1):nr] %>% colnames()

  temp2 <- lapply(OMD_mat, function(x) {
    x <- x %>% as_tibble();
    x$from <- colnames(x)
    x <- tidyr::gather(x, key = 'to', value = 'OMD', -'from' )
  })

  df_plot2 <- as_tibble()
  for (i in 1:length(temp2)){
    x <- temp2[[i]]
    temx <- strsplit(names(temp2)[i],'_')[[1]]
    x$type1 <- temx[1]
    x$type2 <- temx[2]
    df_plot2 <- bind_rows(list(df_plot2, x))
  }

  ## browser()
  df_plot2$type1 <- factor(df_plot2$type1)
  df_plot2$type2 <- factor(df_plot2$type2)##, levels = c('obs','dar') )
  levels(df_plot2$type1) <- c('Darwin Model', 'Remote Sensing')
  levels(df_plot2$type2) <- c('Darwin Model', 'Remote Sensing')

  levels_from = moyr[1:(nr/2)]
  levels_to = moyr[(nr/2):1]

  df_plot2$from <- factor(df_plot2$from, levels = levels_from)
  df_plot2$to <-  factor(df_plot2$to, levels = levels_to)

  ## temporary
  df_plot2$from = df_plot2$from %>% as.numeric()
  df_plot2$to = df_plot2$to %>% as.numeric()

  mo_from = levels_from %>% strsplit("-") %>% purrr::map(.%>%.[[1]] %>% as.numeric()) %>% unlist()
  mo_to = levels_to %>% strsplit("-") %>% purrr::map(.%>%.[[1]] %>% as.numeric()) %>% unlist()
  breaks = which(mo_from == 6)
  labels = paste0(1998:2006)
  breaks2 = which(mo_to == 6)
  labels2 = rev(labels)
  ## end of temporary

  ## Plot data frame as heatmap with facets
  p =
    ggplot(df_plot2 , aes(from, to, fill= OMD)) +
    facet_grid(type2 ~ type1) +
    geom_tile() +
    scale_fill_gradientn(name = "EMD (in km)", colors = hcl.colors(20, "RdYlGn")) +
    ## scale_fill_gradientn(name = "OMD", colors = c("skyblue", "blue", "red", "yellow")) +
    theme_bw() + xlab("Time") + ylab('Time') +
    coord_fixed() +
    theme(axis.text.x = element_text(size = rel(1), angle = 90, vjust = 0.3)) +
    theme(legend.position="left") +
    ## geom_vline(xintercept = breaks, col=rgb(0,0,0,0.5), lty='dotted')+
    ## geom_hline(yintercept = breaks2, col=rgb(0,0,0,0.5), lty = 'dotted')+
    geom_vline(xintercept = breaks+6.5, col=rgb(0,0,0,0.5))+
    geom_hline(yintercept = breaks2-6.5, col=rgb(0,0,0,0.5))+
    scale_y_continuous(breaks = breaks2,
                       labels = labels2,
                       expand = c(0, 0)) +
    scale_x_continuous(breaks = breaks,
                       labels = labels,
                       expand = c(0, 0))

  p = p + theme(panel.spacing = unit(0, "lines"))

  p = p + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                axis.ticks.x=element_blank(),
                axis.ticks.y=element_blank())
  p = p +
    theme(strip.text.x = element_text(size = rel(1))) +
    theme(strip.text.y = element_text(size = rel(1)))
  ## nr = nrow(three_distmats[["all"]])
  nrmax = ceiling(c(nr/2 / 12)) * 12
  p = p + theme(panel.spacing = unit(0, "lines"))
  return(p)
}
