# Toy example for using RunHdpxParallel

# To run this at the command line

# R --vanilla < RunHdpxParallel.example.R  > out.txt 2> err.txt &

# See also https://github.com/steverozen/mSigHdp/raw/master/data-raw/Hierarchical.dirichlet.process.for.mutational.signatures.docx


# library(mSigHdp)
library(ICAMS)
library(PCAWG7)

out.dir <- "./output.from.RunHdpxParallel.example"
dir.create(out.dir)

retval <- RunHdpxParallel (

  input.catalog            = PCAWG7::spectra$PCAWG$SBS96[ , 1:25],
  # Can also be a file name, but the file has to be in ICAMS format.

  out.dir                  = out.dir,
  num.child.process        = 2, # We recommend >=20 for real data
  CPU.cores                = 2,
  seedNumber               = 123,
  K.guess                  = 10,
  burnin.checkpoint        = FALSE,
  burnin                   = 500, # Very very short

  burnin.multiplier        = 2,
  # We recommend >= 10,000 burn-in iterations in total for real data. This toy
  # example uses 1000 x 2 = 2000 iterations in total.

  post.n                   = 20, # 200 is recommended for real data

  post.space               = 10, # 100 is recommended for real data

  multi.types              = TRUE,
  overwrite                = TRUE,
  gamma.alpha              = 1,
  gamma.beta               = 20,
  high.confidence.prop     = 0.9,
  moderate.confidence.prop = 0.9,
  checkpoint.chlist        = FALSE,
  checkpoint.1.chain       = FALSE,
  posterior.checkpoint     = FALSE
  )

# Comment out for profiling
save(retval, file = file.path(out.dir, "RunHdpxParallel.retval.Rdata"))
