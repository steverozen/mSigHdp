# Template for script to start a new burnin

# Run this script from within a directory that will also
# contain the output

# ** Your code here to get the input catalog

if (file.exists("parameters.csv")) {
  params <- read.csv("parameters.csv")
  cpiter <- params$cpiter
  burnin <- params$burnin
} else {
  cpiter <- 3
  burnin <- 100
  write.csv(data.frame(cpiter = cpiter, burnin = burnin),
            "parameters.csv")
}

input.catalog <- ICAMS::ReadCatalog("ground.truth.syn.catalog.csv") # Update this
input.catalog <- input.catalog[1:10, 1:5] # FOR TESTING ONLY!

burnin.output <-
  ActivateAndBurnin(input.catalog    = input.catalog,
                    seedNumber       = 1,
                    K.guess          = 20, # Change this
                    multi.types      = FALSE,
                    verbose          = TRUE,
                    burnin           = burnin,
                    cpiter           = cpiter,
                    burnin.verbosity = 0,
                    gamma.alpha      = 1,
                    gamma.beta       = 1)

burnin.output$number <- 1
burnin.output$all.lik <- burnin.output$likelihood

save(burnin.output, file = "1.burnin.output.Rdata")

pdf(file = "1.burnin.lik.pdf", paper = "a4")
plot(burnin.output$likelihood, pch=16, cex=0.5) ##change to a bigger dot
dev.off()
