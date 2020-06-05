
test_that("PrepInit-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.PrepInit.Rdata",
       envir = reg)


  retvalx <- PrepInit(input.catalog = input.catalog[1:10 , 1:15],
                      multi.types = F,
                      K.guess = 5,
                      verbose = T)

  #save(retvalx, file = "RunhdpInternal.testdata/test.PrepInit.Rdata")

  expect_equal(retvalx$signature, reg$retvalx$signature)
  expect_equal(retvalx$exposure, reg$retvalx$exposure)

})
