#!/usr/bin/env Rscript

use_container <- TRUE

if (use_container) {
  # If using singularity container to run this script, need to tell R to look
  # for packages in the path inside the container. Otherwise, R will use the
  # library path in the host system, which can cause package dependency errors
  .libPaths("/usr/local/lib/R/site-library")

  # It is recommended to use bash script to pass arguments to this R script
  # when running mSigHdp using singularity container

  # args have the arguments passed to Rscript from bash script
  args <- commandArgs(trailingOnly = TRUE)
  project_dir <- args[1]
  seed <- as.numeric(args[2])
  K_guess <- as.numeric(args[3])
  output_dir <- args[4]
} else {
  project_dir <- "/home/e0012078/container/msighdp"
  seed <- 519
  K_guess <- 26
  output_dir <- "/home/e0012078/output/msighdp"
}

library(ICAMS)
library(mSigHdp)

input_file <- file.path(project_dir, "test_catalog.csv")
input_catalog <- ICAMS::ReadCatalog(file = input_file)

output_home <- file.path(output_dir, paste0("seed_", seed))

if (!dir.exists(output_home)) {
  dir.create(path = output_home, recursive = TRUE)
}

checkpoint_dir <- file.path(output_home, "checkpoint")
if (!dir.exists(checkpoint_dir)) {
  dir.create(path = checkpoint_dir, recursive = TRUE)
}
# Set working directory to checkpoint_dir directory so that the checkpoint Rdata
# files will be saved there
setwd(checkpoint_dir)

message("Start running mSigHdp")
message("Using seed ", seed)
message("Using K_guess ", K_guess)
message("mSigHdp version: ")
print(packageVersion("mSigHdp"))

time_used <- system.time(expr = {
  # The parameters used here are for testing purpose only to save time.
  # For recommended parameters to use, please refer to this paper
  # Mo Liu, Yang Wu, Nanhai Jiang, Arnoud Boot, Steven G. Rozen,
  # mSigHdp: hierarchical Dirichlet process mixture modeling
  # for mutational signature discovery,
  # https://www.biorxiv.org/content/10.1101/2022.01.31.478587v1

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
})

saveRDS(time_used, file = file.path(output_home, "time_used.Rds"))
message("Time used: ")
print(time_used)
message("Finished running mSigHdp")
