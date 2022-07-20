
# mSigHdp
  
The goal of mSigHdp is mutational signature discovery using 
hierarchical Dirichlet process (HDP) mixture models. mSigHdp
is only supported on Linux systems. Most users
will use the function RunHdpxParallel.

This package uses https://github.com/steverozen/hdpx for the
hierarchical Dirichlet process implementation.

Please also see our paper: 
Mo Liu, Yang Wu, Nanhai Jiang, Arnoud Boot, Steven G. Rozen,
mSigHdp: hierarchical Dirichlet process mixture modeling 
for mutational signature discovery, 
https://www.biorxiv.org/content/10.1101/2022.01.31.478587v1.


## Installation

### Latest stable version

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", build_vignettes = T)
```

### Get the development version

To use new features in the development version, you can install mSigHdp
from the master branch on [GitHub](https://github.com/), which may not
be stable:

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", ref = "master")
```

## Reference manual for development version

<https://github.com/steverozen/mSigHdp/blob/master/mSigHdp_2.0.1.0010.pdf>
