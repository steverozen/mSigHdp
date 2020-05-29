
# input.catalog <- ICAMS::ReadCatalog("cat.csv")
load("input.catalog.Rdata")

colnames(input.catalog)[1:7] <- paste0("a::", colnames(input.catalog)[1:7])
colnames(input.catalog)[8:ncol(input.catalog)] <- paste0("b::", colnames(input.catalog)[8:ncol(input.catalog)])

input.catalog[3, ] <- 0

input.catalog[5, 1:7] <- 0
input.catalog[6, 8:15] <- 0
input.catalog[10, 8:ncol(input.catalog)] <- 3

retvalx <- hdpx::TestScaffold1(
  input.catalog = input.catalog,
  CPU.cores     = 1,
  seedNumber    = 1000,
  K.guess       = 5,
  multi.types   = TRUE,
  verbose       = TRUE,
  post.burnin   = 50,
  post.n        = 50,
  post.space    = 25,
  post.cpiter   = 3,
  num.posterior = 1
)
# nice valgrind --leak-check=yes ~/debug.hdp/R-4.0.0/bin/exec/R --vanilla < s5.R >s5.out &>s5.err &
