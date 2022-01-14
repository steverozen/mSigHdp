
test_that("ChainBurnin", {

  load("misc.test.data/1000044.activated.hdp.state.Rdata")

  load("RunhdpInternal.testdata/test.ChainBurnin.Rdata")
  reg <- new.env()
  load("RunhdpInternal.testdata/test.GibbsSamplingAfterBurnin.Rdata",
       envir = reg)

  ##Generate burnin checkpoint for testing
  burnin.output.1 <- ChainBurnin(hdp.state          = hdp.state,
                                 seedNumber          = (44 + 3e6),
                                 burnin              = 100,
                                 cpiter              = 3,
                                 burnin.verbosity    = 0,
                                 burnin.multiplier   = 1,
                                 checkpoint          = T)

  burnin.output.2 <- ChainBurnin(hdp.state          = hdp.state,
                                 seedNumber          = (44 + 2e6),
                                 burnin              = 100,
                                 cpiter              = 3,
                                 burnin.verbosity    = 0,
                                 burnin.multiplier   = 1,
                                 checkpoint          = T)

  sample.chain.1 <- GibbsSamplingAfterBurnin(burnin.output     = "mSigHdp.burnin.checkpoint.3000044.Rdata",
                                             post.n         = 10,
                                             post.space     = 5)
  sample.chain.2 <- GibbsSamplingAfterBurnin(burnin.output     = "mSigHdp.burnin.checkpoint.2000044.Rdata",
                                             post.n         = 10,
                                             post.space     = 5)



  # In we need to regenerate the baseline data
  save(sample.chain.1,sample.chain.2, file = "RunhdpInternal.testdata/test.GibbsSamplingAfterBurnin.Rdata")

  expect_equal(sample.chain.1, reg$sample.chain.1)
  expect_equal(sample.chain.2, reg$sample.chain.2)

  unlink( "mSigHdp.burnin.checkpoint.2000044.Rdata")
  unlink( "mSigHdp.burnin.checkpoint.3000044.Rdata")

})
