# Toy example for using RunHdpxParallel
library(mSigHdp)
library(ICAMS)
sessionInfo()

out.dir <- "./output.from.RunHdpxParallel.example2"
dir.create(out.dir)
setwd(out.dir)
# This way the checkpoint files will go into out.dir

retval <- RunHdpxParallel (

  input.catalog      = mSigHdp::test.spectra, # Can also be a file name, but the file has to be in ICAMS format.

  ground.truth.sig   = NULL, # If comparing to signatures in synthetic data, can use mSigHdp::test.ground.truth.sig for testing.

  ground.truth.exp   = NULL, # If comparing to previously estimated exposures or exposures in synthetic data.
                             # Can use mSigHdp::test.ground.truth.exposure for testing.

  out.dir            = ".",

  num.child.process  = 2, # We recommend >=20 for real data

  CPU.cores          = 2,
  seedNumber         = 123,
  K.guess            = 10,

  burnin.checkpoint  = T, # Periodically checkpoint the burnin chain

  post.burnin        = 1000,

  burnin.multiplier  = 2, # We recommend >= 10,000 burn-in iterations in total for real data.
                          # This toy example uses 1000 x 2 = 2000 iterations in total.

  post.n             = 20, # 200 is recommended for real data

  post.space         = 10, # 100 is recommended for real data

  multi.types        = TRUE,
  overwrite          = TRUE,
  gamma.alpha = 1,
  gamma.beta  = 20,
  cos.merge = 0.90,
  confident.prop = 0.9,
  noise.prop = 0.5)

# To run this at the command line
# R --vanilla < RunHdpxParallel.example.R  > out.txt 2> err.txt &
