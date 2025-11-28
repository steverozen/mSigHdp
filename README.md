
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

### Newest stable version

*Compatible with current Rcpp.* Corrects `Error in lower.to.upper.tri.inds(n) : 'n' must be >= 2`

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", ref = "v2.1.3-branch")
```

### Singularity/Apptainer

This is v2.1.2 and still has ``Error in lower.to.upper.tri.inds(n) : 'n' must be >= 2`  
on inputs with only a single signature.

`singularity pull library://rozen-lab/msighdp/msighdp:2.1.2`

For more details, see
https://github.com/steverozen/mSigHdp/blob/master/data-raw/container_scripts/mSigHdp_container_installation.pdf

A toy-example R script for using this container is available at 
https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/data-raw/container_scripts/test_mSigHdp.R.

### Previous stable version

*Incompatible with current Rcpp*

``` r
install.packages("remotes")
remotes::install_github(repo = "steverozen/mSigHdp", ref = "v2.1.2-branch")
```



## Reference manual

<https://github.com/steverozen/mSigHdp/blob/v2.1.2-branch/mSigHdp_2.1.2.pdf>
