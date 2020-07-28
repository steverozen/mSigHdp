
test_that("ChainBurnin", {

  load("1000044.activated.hdp.state.Rdata")
  reg <- new.env()
  load("RunhdpInternal.testdata/test.ChainBurnin.Rdata",
       envir = reg)

  retvalx <- ChainBurnin(hdp.state = hdp.state,
                         seedNumber          = (44 + 3e6),
                         burnin              = 100,
                         cpiter              = 3,
                         burnin.verbosity    = 0,
                         burnin.multiplier = 2,
                         burnin.checkpoint = T)

  #save(retvalx, file = "RunhdpInternal.testdata/test.ChainBurnin.Rdata")

  expect_equal(retvalx, reg$retvalx)

})
