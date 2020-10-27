#' Diagnostic plot for a hdpSampleMulti object. This function is
#' compatible with the return object from Liu's extract_components_from_clusters
#' @param retval output from \code{\link{CombineChainsAndExtractSigs}}
#'
#' @param input.catalog ground truth catalog
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#' @importFrom grDevices dev.off pdf
#'
#' @importFrom graphics par
#'
#' @export



ComponentDiagnosticPlotting <- function(retval,
                                   input.catalog,
                                   out.dir,
                                   verbose){
  multi <- retval$extracted.retval[["multi.chains"]] # class hdpSampleMulti
  chains <- hdpx::chains(multi)      # list of hdpSampleChain

  if (verbose) message("Writing HDP diagnostics")

  pdf(file = file.path(out.dir,"diagnostics.likelihood.pdf"))
  par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
  lapply(chains, hdpx::plot_lik, bty = "L")
  dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.numcluster.pdf"))
  # This is the number of raw clusters sampled along each chain
  par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
  lapply(chains, hdpx::plot_numcluster, bty = "L")
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.data.assigned.pdf"))
  # This is the number of mutations assigned as a function of
  # the number of raw clusters
  par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
  lapply(chains, hdpx::plot_data_assigned, bty = "L")
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.hdp.signature.exposure.each.sample.pdf"))
  myCol <- grDevices::rainbow(ncol(retval$signature), alpha = 1)
  graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))

  mSigHdp::PlotSamplesHighSigExp(retval           = retval,
                                  hdpsample        = multi,
                                  input.catalog    = input.catalog)

  grDevices::dev.off()

  #because extract_component_from_clusters extract informations across chains
  #so the chain information was not recorded as NR's hdp_extract_components
  #therefore, we use extract_ccc_cdc_from_hdp to seek for the information of
  #raw clusters that highly similar as components (most likely these clusters
  #contribute to the components during extract_component_from_clusters).

  ccc_0 <- lapply(chains, function(ch){
    lapply(hdpx::clust_categ_counts(ch), function(x){
      ans <- cbind(x)
      return(ans[, -ncol(ans)])
    })
  })

  cdc_0 <- lapply(chains, function(ch){
    lapply(hdpx::clust_dp_counts(ch), function(x){
      ans <- cbind(x)
      return(ans[, -ncol(ans)])
    })
  })
  sigmatchretval <- apply(retval$signature,2,function(x){
    hdpx::extract_ccc_cdc_from_hdp(x,
                                   ccc_0 = ccc_0,
                                   cdc_0 = cdc_0,
                                   cos.merge = 0.90)})

  grDevices::pdf(file = file.path(out.dir,"diagnostics.signatures.pdf"))
  graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
  # This plots the component (signature) profiles with
  # 95% credibility intervals
  hdpx::plot_component_with_credint(retval =sigmatchretval)
  grDevices::dev.off()



  grDevices::pdf(file.path(out.dir,"component.distribution.in.posterior.samples.pdf"))
  graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
  hdpx::plot_component_posterior_samples(components = retval$signature,
                                         retval = sigmatchretval)
  grDevices::dev.off()


}



