# Toy example for using RunHdpxParallel

# To run this at the command line

# R --vanilla < RunHdpxParallel.example.R  > out.txt 2> err.txt &

# See also https://github.com/steverozen/mSigHdp/raw/master/data-raw/Hierarchical.dirichlet.process.for.mutational.signatures.docx


library(mSigHdp)
library(ICAMS)
sessionInfo()

out.dir <- "./output.from.RunHdpxParallel.example"
dir.create(out.dir)
setwd(out.dir)
# So the checkpoint files will go into out.dir

retval <- RunHdpxParallel (

  input.catalog      = mSigHdp::test.spectra,
  # Can also be a file name, but the file has to be in ICAMS format.

  ground.truth.sig   = NULL,
  # If comparing to previously estimated signatures or signatures in synthetic
  # data, can use mSigHdp::test.ground.truth.sig for testing.

  ground.truth.exp   = NULL,
  # If comparing to previously estimated exposures or exposures in synthetic
  # data, can use mSigHdp::test.ground.truth.exposure for testing.

  out.dir            = ".",
  num.child.process  = 2, # We recommend >=20 for real data
  CPU.cores          = 2,
  seedNumber         = 123,
  K.guess            = 10,
  burnin.checkpoint  = TRUE,
  burnin             = 1000,

  burnin.multiplier  = 2,
  # We recommend >= 10,000 burn-in iterations in total for real data. This toy
  # example uses 1000 x 2 = 2000 iterations in total.

  post.n             = 20, # 200 is recommended for real data

  post.space         = 10, # 100 is recommended for real data

  multi.types        = TRUE,
  overwrite          = TRUE,
  gamma.alpha = 1,
  gamma.beta  = 20,
  high.confidence.prop = 0.9,
  moderate.confidence.prop = 0.9)

save(retval, file = "RunHdpxParallel.retval.Rdata")
