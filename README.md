
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mSigHdp

This package extracts mutational signatures using the hdpx package This
branch of mSigHdp depends on a branch of hdpx that was forked from
Nicola Roberts’s hdp package, and that still implements the algorithm
from Section 4.2.2, “Extracting Consensus Signatures” from

Roberts, N. D. (2018). Patterns of somatic genome rearrangement in human
cancer. (PhD Thesis). Cambridge University, Cambridge, England, United
Kingdom. Retrieved from
<https://www.repository.cam.ac.uk/bitstream/handle/1810/275454/Roberts-2018-PhD.pdf>

We do not recommend using this version.

## Installation

*This version *must* use hdpx branch \`NR-version-plus-fixes*.

``` r
# Always install the version of hdpx and mSigHdp with the Nicola
# Roberts algorithm for combining raw clusters into signatures
hdpx.version <- "0.1.5.0099"
if (system.file(package = "hdpx") != "") {
  if (packageVersion("hdpx") != hdpx.version) {
    remove.packages("hdpx")
    remotes::install_github("steverozen/hdpx", ref = "NR-version-plus-fixes")
  }
} else {
  remotes::install_github("steverozen/hdpx", ref = "NR-version-plus-fixes")
}
message("hdpx version ", packageVersion("hdpx"))
stopifnot(packageVersion("hdpx") == hdpx.version)


mSigHdp.version <- "0.0.0.9019"
if (system.file(package = "mSigHdp") != "") {
  if (packageVersion("mSigHdp") != mSigHdp.version) {
    remove.packages("mSigHdp")
    remotes::install_github(repo = "steverozen/mSigHdp", 
                            ref = "for-NR-version-plus-fixes")
  }
} else {
  remotes::install_github(repo = "steverozen/mSigHdp", 
                          ref = "for-NR-version-plus-fixes")
}
```

## Reference manual

<https://github.com/steverozen/mSigHdp/blob/master/data-raw/mSigHdp_0.0.0.9015.pdf>
