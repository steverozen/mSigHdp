# Template for script to Gibbs sampling after burnin

# Run this script from within a directory that will also
# contain the output


load(paste0(last.burnin, ".burnin.output.Rdata"))

# Loaded variable is burnin.output
# this should be a well burned-in MCMC chain that available for just Gibbs sampling

num.child.process <- 20
CPU.cores         <- 20
n.post.sample     <- 10
n.space           <- 10
seedNumber        <- 44

run.hdp.posterior.sample <- function(seed){
  one_sample_chains <-
    hdpx::hdp_posterior_sample(burnin.output    = burnin.output,
                               n                = n.post.sample,
                               space            = n.space,
                               cpiter           = 3,
                               verbosity        = 0,
                               seed             = seed)
  return(one_sample_chains)
}

multiple_sample_chains <- parallel::mclapply(
  # Must choose a different seed for each of the chains
  X = (seedNumber + 1:num.child.process * 10^6) ,
  FUN = run.hdp.posterior.sample,
  mc.cores = CPU.cores)

clean.chlist <- CleanChlist(multiple_sample_chains)

save(clean.chlist, file = "posterior.sample.chains.output.Rdata")


##the clean.chlist can be used as input of CombinePosteriorChains, then AnalyzeandPlotretval
