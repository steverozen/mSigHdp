# Example script for running Runhdp{4,5}

library(mSigHdp)

prog.number <- 4
if (prog.number == 4) {
  my.fn    <- Runhdp4
} else if (prog.number == 5) {
  my.fn <- Runhdp5
}

root.dir <- "."
seed     <- 7744
num.jobs <- 10
burnin   <- 1000 # SUPER LOW FOR TESTING FOR BUGS!
K.guess  <- 25
n        <- 200 # SUPER LOW FOR TESTING FOR BUGS!

out.dir  <-
  file.path(root.dir,
            paste("WES.id.out.dir",
                  "seed", seed,
                  "burnin", burnin,
                  "K", K.guess,
                  "n", n,
                  "prog", prog.number, sep = "."))
dir.create(out.dir, recursive = TRUE)
of <- file(file.path(out.dir, "out.txt"), open = "wt")
sink(file = of, type = "output")
ef <- file(file.path(out.dir, "err.txt"), open = "wt")
sink(file = ef, type = "message")
my.fn

cat1 <- ICAMS::ReadCatalog("WES_TCGA.indels.csv")
cat2 <- ICAMS::ReadCatalog("catKumar+Murphy.csv")
in.cat <- cbind(cat1, cat2)

system.time(
  retval <- my.fn(
    input.catalog      = in.cat,
    out.dir            = out.dir,
    CPU.cores          = num.jobs,
    seedNumber         = seed,
    K.guess            = K.guess,
    post.burnin        = burnin,
    post.n             = n,
    post.space         = 10,
    multi.types        = FALSE,
    overwrite          = TRUE,
    num.posterior      = num.jobs,
    plot.extracted.sig = TRUE))

# nice R --vanilla < WES.id.R > out.WES.id.K25.b10000n2000.txt &> err.WES.idK25.b10000n2000.txt &

