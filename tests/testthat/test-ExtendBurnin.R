
test_that("ExtendBurnin", {

  input <- new.env()
  load("RunhdpInternal.testdata/mSigHdp.burnin.checkpoint-input.Rdata",
       envir = input)

  reg <- new.env()
  load("RunhdpInternal.testdata/test.ExtendBurnin.Rdata",
       envir = reg)

  retvalx <- ExtendBurnin(previous.burnin.output = input$retvalx,
                           burnin                = 99,
                           cpiter                = 3)

  # If we need to regenerate the test data:
  if (FALSE) {
    save(retvalx,
         file = "RunhdpInternal.testdata/mSigHdp.burnin.checkpoint-input.Rdata")
  }
  expect_equal(retvalx, reg$retvalx)

  retvalx <- ExtendBurnin(
    previous.burnin.output =
      "RunhdpInternal.testdata/mSigHdp.burnin.checkpoint-input.Rdata",
    burnin                = 99,
    cpiter                = 3)
  expect_equal(retvalx, reg$retvalx)

})
