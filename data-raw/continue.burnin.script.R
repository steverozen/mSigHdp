# Template for script to continue a new burnin

# Run this script from within a directory that will also
# contain the output

load(paste0(last.burnin, ".burnin.output.Rdata"))
# Loaded variable is burnin.output
last.burnin <- burnin.output$number
this.burnin <- last.burnin + 1

if (file.exists("parameters.csv")) {
  params <- read.csv("parameters.csv")
  cpiter <- params$cpiter
  burnin <- params$burnin
} else {
  stop("Need parameters in parameters.csv")
}

burnin.output2 <-
  hdpx::hdp_burnin(hdp       = hdpx:::as.hdpState(burnin.output$hdplist),
                   burnin    = burnin,
                   cpiter    = cpiter,
                   verbosity = 0)

burnin.output2$number <- this.burnin
burnin.output2$all.lik <-
  c(burnin.output$all.lik, burnin.output2$likelihood)
burnin.output <- burnin.output2

save(burnin.output, file = paste0(this.burnin, ".burnin.output.Rdata"))

pdf(file = paste0(this.burnin, ".burnin.lik.pdf"), paper = "a4")
plot(burnin.output$likelihood, pch = ".")
plot(burnin.output$all.lik, pch = ".")
dev.off()
