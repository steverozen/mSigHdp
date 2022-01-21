
test_that("ParallelGibbsSample-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  reg <- new.env()
  load("RunhdpInternal.testdata/test.MultipleSetupAndPosterior.Rdata",
       envir = reg)

  retvalx <- ParallelGibbsSample( input.catalog = input.catalog[1:10 , 1:15],
                                  seedNumber          = 44,
                                  K.guess             = 5,
                                  multi.types         = FALSE,
                                  verbose             = TRUE,
                                  burnin              = 50,
                                  post.n              = 5,
                                  post.space          = 50,
                                  post.cpiter         = 3,
                                  post.verbosity      = 0,
                                  CPU.cores           = 2,
                                  num.child.process   = 2,
                                  checkpoint          = FALSE
  )


  #save(retvalx, file = "RunhdpInternal.testdata/test.MultipleSetupAndPosterior.Rdata")

  expect_equal(retvalx, reg$retvalx)
})
