# [0.0.0.9005]- 2020-05-16
* Changed arg for `Runhdp4` and `RundAndEvalHdp4` from `input.catalog.file` to
  `input.catalog`; starting to retire `Runhpd5`, `RunAndEvalHdp5`, `RunhdpInternal5`.
* Added reference manual and README.
* Updated DESCRIPTION to include `biocViews` and `Imports` as field name.

# [0.0.0.9004]- 2020-05-15
* Still trying to make the code robust to errors in `mclapply` children
  run `hdp_posterior`.
