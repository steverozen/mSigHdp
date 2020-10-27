# Toy example for using RunHdpxParallel
#

library(mSigHdp)

out.dir <- "./output.from.RunHdpxParallel.example4"
dir.create(out.dir)
setwd(out.dir)
# This way the checpoint files will go into out.dir

seed <- 1234

retval <- RunHdpxParallel (

  input.catalog      = mSigHdp::test.spectra, # Can also be a file name, but the file has to be in ICAMS format.

  ground.truth.sig   = NULL, # UPDATE this "~/mSigHdp/tests/ground.truth.syn.sigs.csv", #for evaluation purpose,can be NULL

  ground.truth.exp   = NULL, # If comparing to previously estimated exposures or exposures in synthetic data.
                             # Can use mSigHdp/tests/test.SP.Syn.Bladder-TCC.ground.truth.exposure.csv for testing.

  out.dir            = ".",

  num.child.process  = 2, #number of independent MCMC chain initiated. 20 and more is recommended. 2 is for a quick test

  CPU.cores          = 2, #equal as num.child.process

  seedNumber         = seed, #seed of random generator

  K.guess            = 10, # Set to twice your guess of how many signatures will be extracted; does not need to be precise.

  post.burnin        = 1000, # For real data recommend 10,000

  burnin.checkpoint  = T, # Periodically checkpoint the burnin chain

  burnin.multiplier  = 2, # after every post.burnin iterations, a burnin.checkpoint will be saved. This settings means 1000x2 = 2000 burn-in iterations. This setting depends on the data size.

  post.n             = 20, #number of posterior samples collected for each chain in Gibb's sampling

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
