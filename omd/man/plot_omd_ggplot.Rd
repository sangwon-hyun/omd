% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot.R
\name{plot_omd_ggplot}
\alias{plot_omd_ggplot}
\title{ggplot-based plotting of the optimal transport. Produces three figures: (1)
The "before" image, (2) the "after" image, and (3) the optimal transports.}
\usage{
plot_omd_ggplot(
  obj,
  plot_type = c("one", "four"),
  name_from = NULL,
  name_to = NULL,
  add_map = FALSE,
  classify_quantile = TRUE,
  sample_arrows = FALSE
)
}
\arguments{
\item{obj}{An omd object.}
}
\description{
ggplot-based plotting of the optimal transport. Produces three figures: (1)
The "before" image, (2) the "after" image, and (3) the optimal transports.
}
