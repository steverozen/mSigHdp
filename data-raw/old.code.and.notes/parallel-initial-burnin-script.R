# Template for script to start a new burnin

# Run this script from within a directory that will also
# contain the output


library(mSigHdp)
library(hdpx)


if (file.exists("parameters.csv")) {
  params <- read.csv("parameters.csv")
  cpiter <- params$cpiter
  burnin <- params$burnin
} else {
  cpiter <- 3
  burnin <- 15000
  write.csv(data.frame(cpiter = cpiter, burnin = burnin),
            "parameters.csv")
}

setwd("~/panc.test/panc.syn.binom.size.100/") ##the dir to get catalog and output
seedNumber    <- 20616
no.chains     <- 1 ##number of independent initial burn-ins
K.guess       <- 400
input.catalog <- ICAMS::ReadCatalog("ground.truth.syn.catalog.csv") # Update this



run.activateandburn <- function(seed){
  message(paste0("Start ",burnin, " burn-in iterations with seed ",seed))
  burnin.output <- ActivateAndBurnin(input.catalog    = input.catalog,
                                     seedNumber       = seed,
                                     K.guess          = K.guess, # Change this
                                     multi.types      = TRUE,
                                     verbose          = TRUE,
                                     burnin           = burnin,
                                     cpiter           = cpiter,
                                     burnin.verbosity = 0,
                                     gamma.alpha = 1,
                                     gamma.beta  = 1)

  burnin.output$number <- 1
  
  burnin.output$all.lik <- burnin.output$likelihood ##this will be used in continue burnin
  
  
  pdf(file = paste0(seed,".1.burnin.lik.pdf"), paper = "a4")
  plot(burnin.output$likelihood, pch=16, cex=0.5) ##change to a bigger dot
  dev.off()
  
  save(burnin.output, file = paste0(seed,".1.burnin.alpha1beta1output.k400.Rdata"))
}

if(no.chains ==1){
  run.activateandburn(seedNumber)
}else{
  parallel::mclapply(
    # Must choose a different seed for each of the chains
    X = (seedNumber + 1:no.chains * 10^6) ,
    FUN = run.activateandburn,
    mc.cores = no.chains)
}




