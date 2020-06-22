test_that("GenerateAverageCluster", {

  input <- new.env()
  load("RunhdpInternal.testdata/test.CleanChlist.Rdata",
       envir = input)
  reg <- new.env()
  load("RunhdpInternal.testdata/test.GenerateAverageCluster.Rdata",
       envir = reg)


  retvalx <- GenerateAverageCluster(input$retvalx)

  #save(retvalx, file = "RunhdpInternal.testdata/test.GenerateAverageCluster.Rdata")

  expect_equal(retvalx,reg$retvalx)


})
