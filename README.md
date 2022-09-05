
<!-- README.md is generated from README.Rmd. Please edit that file -->

# mSigHdp

The goal of mSigHdp is to do mutational signature extraction using the
hdpx package This branch of mSigHdp depends on a branch of hdpx that was
forked from Nicola Roberts’s hdp package, and that still implements the
algorithm from Section 4.2.2, “Extracting Consensus Signatures” from

Roberts, N. D. (2018). Patterns of somatic genome rearrangement in human
cancer. (PhD Thesis). Cambridge University, Cambridge, England, United
Kingdom. Retrieved from
<https://www.repository.cam.ac.uk/bitstream/handle/1810/275454/Roberts-2018-PhD.pdf>

We do not recommend using this version.

## Installation

``` r
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
