
test_that("RunhdpInternal4-fast", {

  input.catalog <-
    ICAMS::ReadCatalog("SBS96.ground.truth/ground.truth.syn.catalog.csv")

  input.catalog <- input.catalog[1:10, 1:6]
  colnames(input.catalog) <-c(paste0("a::", 1:3), paste0("b::", 4:6))

  reg <- new.env()
  load("RunhdpInternal.testdata/test.RunhdpInternal4-2types.Rdata",
       envir = reg)

  retvalx <- RunhdpInternal4(
    input.catalog = input.catalog,
    CPU.cores     = 1,
    seedNumber    = 44,
    K.guess       = 5,
    multi.types   = TRUE, # <--- !
    verbose       = TRUE,
    post.burnin   = 50,
    num.posterior = 1
  )

  # save(retvalx, file = "RunhdpInternal.testdata/test.RunhdpInternal4-2types.Rdata")

  expect_equal(retvalx$signature, reg$retvalx$signature)
  expect_equal(retvalx$exposure,  reg$retvalx$exposure)
})
