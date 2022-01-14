
test_that("ChainBurnin", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/2.type.ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.GibbsSamplingAfterBurnin.Rdata",
       envir = reg)
  load("RunhdpInternal.testdata/test.ExtractAfterBurnin.Rdata",
       envir = reg)

  hdp.state.1 <- SetupAndActivate(input.catalog[1:10 , 1:15],
                                  seedNumber          = (44 + 2e6),
                                  K.guess             = 5,
                                  multi.types         = TRUE,
                                  verbose             = TRUE,
                                  gamma.alpha         = 1,
                                  gamma.beta          = 20)

  hdp.state.2 <- SetupAndActivate(input.catalog[1:10 , 1:15],
                                  seedNumber          = (44 + 3e6),
                                  K.guess             = 5,
                                  multi.types         = TRUE,
                                  verbose             = TRUE,
                                  gamma.alpha         = 1,
                                  gamma.beta          = 20)

  ##Generate burnin checkpoint for testing
  burnin.output.1 <- ChainBurnin(hdp.state          = hdp.state.1,
                                 seedNumber          = (44 + 2e6),
                                 burnin              = 100,
                                 cpiter              = 3,
                                 burnin.verbosity    = 0,
                                 burnin.multiplier   = 1,
                                 checkpoint          = T)

  burnin.output.2 <- ChainBurnin(hdp.state          = hdp.state.2,
                                 seedNumber          = (44 + 3e6),
                                 burnin              = 100,
                                 cpiter              = 3,
                                 burnin.verbosity    = 0,
                                 burnin.multiplier   = 1,
                                 checkpoint          = T)

  # The names of the burnin checkpoint files are based on
  # the seedNumber arguments to the two calls to ChainBurnin, above.
  sample.chain.1 <- GibbsSamplingAfterBurnin(burnin.output  = "mSigHdp.burnin.checkpoint.2000044.Rdata",
                                             post.n         = 10,
                                             post.space     = 5)
  sample.chain.2 <- GibbsSamplingAfterBurnin(burnin.output  = "mSigHdp.burnin.checkpoint.3000044.Rdata",
                                             post.n         = 10,
                                             post.space     = 5)

  retvalx <-
    CombineChainsAndExtractSigs(clean.chlist = c(list(sample.chain.1),
                                                 list(sample.chain.2)),
                                input.catalog = input.catalog[1:10 , 1:15],
                                verbose = TRUE)

  if (FALSE) { # Set to TRUE to regenerate the baseline data
    save(sample.chain.1,sample.chain.2,
         file = "RunhdpInternal.testdata/test.GibbsSamplingAfterBurnin.Rdata")
    save(retvalx,
         file = "RunhdpInternal.testdata/test.ExtractAfterBurnin.Rdata")
  }

  expect_equal(sample.chain.1, reg$sample.chain.1)
  expect_equal(sample.chain.2, reg$sample.chain.2)
  expect_equal(retvalx, reg$retvalx )

  unlink( "mSigHdp.burnin.checkpoint.2000044.Rdata")
  unlink( "mSigHdp.burnin.checkpoint.3000044.Rdata")

})
