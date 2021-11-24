Analysis of Chlorophyll data using Wasserstein distance
=============

This repository contains code for applying Wasserstein's distance for analyzing
ocean data, for this paper:
[https://arxiv.org/abs/2111.08736](https://arxiv.org/abs/2111.08736). The main
motivating application is that of Chlorophyll maps and distributions, over space
and time.

The name of the paper and this `R` package are "ocean mover's distance (omd)",
which is a play on "earth mover's distance".



## Installation

This R package can be installed using the following commands.

```{r}
library(devtools)
install_github("sangwon-hyun/omd", subdir = "omd")
```

## Usage

The code to conduct the experiments and produce the figures in the paper is in
[.\main](https://github.com/sangwon-hyun/omd/tree/master/main).
	
## Authors & Contributors

Sangwon Hyun, Aditya Mishra
