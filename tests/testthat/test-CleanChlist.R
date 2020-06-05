test_that("SetupAndPosterior-fast", {

  reg <- new.env()
  load("RunhdpInternal.testdata/test.CleanChlist.Rdata",
       envir = reg)

  seedNumber    = 44
  num.child.process =  2
  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")
  input.catalog = input.catalog[1:10,1:15]
  seedNumber          = 44
  K.guess = 5
  multi.types         = FALSE
  verbose             = TRUE
  post.burnin         = 50
  post.n              = 50
  post.space          = 5
  post.cpiter         = 1
  post.verbosity      = 0
  CPU.cores           = 2
  num.child.process   = 2

  run.setup.and.posterior <- function(seedNumber) {

    if (verbose) message("Runing run.setup.and.posterior on ", seedNumber)
    sample.chain <-SetupAndPosterior(
      input.catalog,
      seedNumber     = seedNumber,
      K.guess        = K.guess,
      multi.types    = multi.types,
      verbose        = verbose,
      post.burnin    = post.burnin,
      post.n         = post.n,
      post.space     = post.space,
      post.cpiter    = post.cpiter,
      post.verbosity = post.verbosity)
    return(sample.chain)
  }

  chlist <- parallel::mclapply(
    # Must choose a different seed for each of the chains
    X = (seedNumber + 1:num.child.process * 10^6) ,
    FUN = run.setup.and.posterior,
    mc.cores = CPU.cores)

  retvalx <- CleanChlist(chlist)

  # save(retvalx, file = "RunhdpInternal.testdata/test.CleanChlist.Rdata")

  expect_equal(retvalx,reg$retvalx)
})
