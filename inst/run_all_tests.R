devtools::load_all()
Sys.setenv(MSIGHDP_LONG = "Y")
devtools::test()
cat("finished")

