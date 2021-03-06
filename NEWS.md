# [1.0.2.003]
* MatchSigsAndRelabel was moved to ICAMSxtra

# [0.0.0.9030] - 2020-07-14
* No longer necessary to checkpoint clean.chlist.Rdata

# [0.0.0.9029] - 2020-07-13
* Added ability to specify parameters for gamma distribution prior on 
  Dirichlet concentration parameters for the top level Dirichlet process node.

# [0.0.0.9023] - 2020-06-27
* Added some checkpointing facilities in `MultipleSetupAndPosterior` and `RunHdpParallel`.

# [0.0.0.9015]- 2020-06-06
* Stable version with package source file available to download for installation. See [README.md](https://github.com/steverozen/mSigHdp/blob/master/README.md) for more details.

# [0.0.0.9005]- 2020-05-16
* Changed arg for `Runhdp4` and `RundAndEvalHdp4` from `input.catalog.file` to
  `input.catalog`; starting to retire `Runhpd5`, `RunAndEvalHdp5`, `RunhdpInternal5`.
* Added reference manual and README.
* Updated DESCRIPTION to include `biocViews` and `Imports` as field name.
* Fixed typos and updated WORDLIST for spelling check.
* Added test data for test-RunhdpInternal4-fast-2types.R.
* Fixed a bug in test-RunhdpInternal4-slow.R.

# [0.0.0.9004]- 2020-05-15
* Still trying to make the code robust to errors in `mclapply` children
  run `hdp_posterior`.
