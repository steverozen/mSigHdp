# Example script for running RunAndEvalHdp5 or RunAndEvalHdp4

library(mSigHdp)

prog.number <- 5
if (prog.number == 4) {
  my.fn    <- RunAndEvalHdp4
} else if (prog.number == 5) {
  my.fn <- RunAndEvalHdp5
}

root.dir <- "."
seed     <- 7744
num.jobs <- 10
burnin   <- 1000
K.guess  <- 25

out.dir  <-
  file.path(root.dir,
            paste("out.dir",
                  "seed", seed,
                  "burnin", burnin,
                  "K", K.guess,
                  "prog", prog.number, sep = "."))
dir.create(out.dir, recursive = TRUE)
of <- file(file.path(out.dir, "out.txt"), open = "wt")
sink(file = of, type = "output")
ef <- file(file.path(out.dir, "err.txt"), open = "wt")
sink(file = ef, type = "message")
my.fn

retval <- my.fn(
  input.catalog.file         = file.path(root.dir,
                                         "ground.truth.syn.catalog.csv"),
  ground.truth.exposure.file = file.path(root.dir,
                                         "ground.truth.syn.exposures.csv"),
  ground.truth.sig.file      = file.path(root.dir,
                                         "ground.truth.syn.sigs.csv"),
  out.dir                    = out.dir,
  CPU.cores          = num.jobs,
  seedNumber         = seed,
  K.guess            = K.guess,
  post.burnin        = burnin,
  post.n             = 750,
  post.space         = 10,
  multi.types        = FALSE,
  overwrite          = TRUE,
  num.posterior      = num.jobs)

# nice R --vanilla < x5K25.R > out.x5K25.txt &> err.x5K25.txt &

