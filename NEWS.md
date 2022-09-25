# [ 2.1.0.1]
 * Added function `show_downsample_curves`.

# [ 2.1.0 ]
 * Replaced `mSigAct` with `mSigTools` for writing and plotting exposure.
 * Updated required branch for hdpx.

# [ 2.0.1.0010 ]
 * Added downsampling capability.

# [ 2.0.1.0009 ]
 * 2022 06 24 Replaced hc.cutoff with merge.raw.cluster.args; passes check except for
   2 warnings unrelated to hc.cutoff, which I cannot figure out how to
   silence.

# [ 2.0.1 ]
 * Added pointers to the paper on bioRxiv

# [ 2.0.0 ] 
 * Removing "frozen nodes" (still in branch "last-with-NR-frozen-nodes")

# [ 1.2.7 ]
 * Dependency on new version of hdpx

# [ 1.2.6 ]
 * Added more tests
 * Major documentation cleanup
 * Some function renaming
 * Made some functions internal

# [ 1.2.1 ]
 * Removed "moderate" classification of confidence in aggregated clusters / signatures.
 * Incremented version of hdpx required.
 * Removed 'multi.types' from CombineChainsAndExtractSigs
 * Add test for AnalyzeAndPlotretval
 * Update some documentations

# [ 1.2.0 ]
 * Many changes

# [ 1.1.4.9001]
* Improved error handling and reporting from mclapply child processes.

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
