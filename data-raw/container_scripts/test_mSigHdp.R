# This is designed to run in a singularity/apptainer container.
# Do not call it directly from the command line.
#
# Suggested call for simple example
# singularity exec <my_R_container> Rscript --vanilla test_mSigHdp.R 123
#
# But normally something like:
# nice singularity exec <my_R_container> Rscript --vanilla test_mSigHdp.R 123 >& log.txt &

# We need to tell R to look for packages first inside the container
# but also to be able to look outside the container (the "host"
# system). This works only if the host OS is the same as
# the container OS, or if the libraries in the host system do not
# have compiled code. In the last case, you would have to
# convert the container to a sandbox, add the necessary libraries,
# and then convert back to a .sif.
.libPaths(c("/usr/local/lib/R/site-library", .libPaths()))

# args have the arguments passed to Rscript from bash script
args <- commandArgs(trailingOnly = TRUE)
if (length(args) == 0) {
  stop("Please provide the random seed as an argument")
}
seed <- as.numeric(args[1])

library(ICAMS)
library(mSigHdp)

# All inputs and outputs are in the current working directory

input_file <- "test_catalog.csv"
input_catalog <- ICAMS::ReadCatalog(file = input_file)

output_home <- paste0("seed_", seed)

if (!dir.exists(output_home)) {
  if (!dir.create(path = output_home, recursive = TRUE)) {
    stop("unable to create ", output_home)
  }
}

checkpoint_dir <- file.path(output_home, "checkpoint")
if (!dir.exists(checkpoint_dir)) {
  if (!dir.create(path = checkpoint_dir, recursive = TRUE)) {
    stop("unable to create ", checkpoint_dir)
  }
}

# Set working directory to checkpoint_dir directory so that the checkpoint Rdata
# files will be saved there
setwd(checkpoint_dir)

K_guess <- 26

message("Start running mSigHdp")
message("Using seed ", seed)
message("Using K_guess ", K_guess)
message("mSigHdp version: ")
print(packageVersion("mSigHdp"))

# The parameters used here are for testing purpose only to save time.
# For recommended parameters to use, please refer to this paper
# Mo Liu, Yang Wu, Nanhai Jiang, Arnoud Boot, Steven G. Rozen,
# mSigHdp: hierarchical Dirichlet process mixture modeling
# for mutational signature discovery,
# https://doi.org/10.1093/nargab/lqad005

mSigHdp::RunHdpxParallel(
  input.catalog = input_catalog,
  seedNumber = seed,
  K.guess = K_guess,
  out.dir = output_home,
  multi.types = FALSE,
  burnin = 1000,
  burnin.multiplier = 2,
  post.n = 5,
  post.space = 10,
  num.child.process = 2,
  CPU.cores = 2,
  high.confidence.prop = 0.9,
  gamma.alpha = 1,
  gamma.beta = 20,
  checkpoint = TRUE,
  verbose = FALSE,
  downsample_threshold = 3000
)

message("Finished running mSigHdp")
