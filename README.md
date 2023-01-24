
# mSigHdp: Mutational Signature Discovery Using Hierarchical Dirichlet Processes
  
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
https://doi.org/10.1093/nargab/lqad005.


## Installation

### Singularity/Apptainer

`singularity pull library://rozen-lab/msighdp/msighdp:2.1.2`

A toy-example R script for using this container is available at 
https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/data-raw/container_scripts/test_mSigHdp.R.

### Latest stable version

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", ref = "v2.1.2-branch", build_vignettes = TRUE)
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

<https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/mSigHdp_2.1.2.pdf>
