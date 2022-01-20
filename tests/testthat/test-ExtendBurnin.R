
test_that("ExtendBurnin", {

  input <- new.env()
  load("RunhdpInternal.testdata/ExtendBurnin.Rdata",
       envir = input)

  reg <- new.env()
  load("RunhdpInternal.testdata/test.ExtendBurnin.Rdata",
       envir = reg)

  retvalx  <- ExtendBurnin(hdplist    = input$retvalx$hdplist,
                           seedNumber  = 1,
                           burnin      = 100,
                           cpiter      = 3,
                           verbosity   = 0)

  # If we need to regenerate the test data:
  # save(retvalx, file = "RunhdpInternal.testdata/test.ExtendBurnin.Rdata")

  expect_equal(retvalx, reg$retvalx)

})
