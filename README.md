
# mSigHdp

The goal of mSigHdp is mutational signature discovery using 
hierarchical Dirichlet process (HDP) mixture models. mSigHdp
is only supported on Linux systems. Most users
will use the function RunHdpxParallel.

This package uses https://github.com/steverozen/hdpx for the
hierarchical Dirichlet process implementation. 

## Installation

### Latest stable version

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", vignette = T)
```

### Get the development version

To use new features in the development version, you can install mSigHdp
from the master branch on [GitHub](https://github.com/), which may not
be stable:

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", ref = "master")
```

## Reference manual

<https://github.com/steverozen/mSigHdp/blob/v2.0.0-branch/mSigHdp_2.0.0.pdf>
