# Toy example for using RunHdpxParallel
#

library(mSigHdp)

seed <- 1234

retval <- RunHdpxParallel (

  input.catalog = "~/mSigHdp/tests/test.SP.Syn.Bladder-TCC.catalog.csv", #Path to your data. Please transform your data to the same format. And sample names as ‘TumorType::sample” (e.g. “Bladder-TCC::SP1234”)

  ground.truth.sig  = "~/mSigHdp/tests/ground.truth.syn.sigs.csv", #for evaluation purpose,can be NULL

  ground.truth.exp = "~/mSigHdp/tests/test.SP.Syn.Bladder-TCC.ground.truth.exposure.csv",#for evaluation purpose,can be NULL

  out.dir            = out.dir, #the directory to store results

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
# R --vanilla < hdp.test.run.R  > hdp.test.run.out 2> hdp.test.run.err &
