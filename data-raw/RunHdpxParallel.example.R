# Toy example for using RunHdpxParallel
#

library(mSigHdp)

out.dir <- "./output.from.RunHdpxParallel.example2"
dir.create(out.dir)
setwd(out.dir)
# This way the checkpoint files will go into out.dir

retval <- RunHdpxParallel (

  input.catalog      = mSigHdp::test.spectra[,1:3], # Can also be a file name, but the file has to be in ICAMS format.

  ground.truth.sig   = NULL, # If comparing to signatures in synthetic data, can use this "mSigHdp::test.ground.truth.sig" for testing.

  ground.truth.exp   = NULL, # If comparing to previously estimated exposures or exposures in synthetic data.
                             # Can use mSigHdp::test.ground.truth.exposure for testing.

  out.dir            = ".",

  num.child.process  = 2, #number of independent MCMC chain initiated. 20 and more is recommended. 2 is for a quick test

  CPU.cores          = 2, #equal as num.child.process

  seedNumber         = 123, # Random number seed

  K.guess            = 10, # Set to twice your guess of how many signatures will be extracted; does not need to be precise.

  post.burnin        = 1000, # Number of burn-in iterations

  burnin.checkpoint  = T, # Periodically checkpoint the burnin chain

  burnin.multiplier  = 2, # after every post.burnin iterations, a burnin.checkpoint will be saved. This settings means 1000x2 = 2000 burn-in iterations.
                          # This setting depends on the data size. We recommended a total of 10,000 burn-in iterations for real data

  post.n             = 20, #number of posterior samples collected for each chain in Gibbs sampling

  post.space         = 10, #interval between collection of two posterior samples

  multi.types        = TRUE, #if TRUE, hdp sets the structure of grandparent node for whole dataset and parent node for each tumor type.

  ##The following parameters
  overwrite          = TRUE,
  gamma.alpha = 1,
  gamma.beta  = 20,
  cos.merge = 0.90,
  confident.prop = 0.9,
  noise.prop = 0.5)

# To run this at the command line
# R --vanilla < RunHdpxParallel.example.R  > hdp.test.run.out 2> hdp.test.run.err &
