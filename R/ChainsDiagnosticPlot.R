#' Diagnostic plot for hdp multi sample chains (output from CombinePosteriorChains)
#' @param retval output from CombinePosteriorChains.A list with the following elements:\describe{
#'    \item{signature}{The extracted signature profiles as a matrix;
#'             rows are mutation types, columns are
#'             samples (e.g. tumors).}
#'    \item{exposure}{The inferred exposures as a matrix of mutation counts;
#'            rows are signatures, columns are samples (e.g. tumors).}

#'
#'    \item{multi.chains}{A \code{\link[hdpx]{hdpSampleMulti-class}} object.
#'     This object has the method \code{\link[hdpx]{chains}} which returns
#'     a list of \code{\link[hdpx]{hdpSampleChain-class}} objects. Each of these
#'     sample chains objects has a method \code{\link[hdpx]{final_hdpState}}
#'     (actually the methods seems to be just \code{hdp})
#'     that returns the \code{hdpState} from which it was generated.}
#' }
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#'
#'
#' @export



ChainsDiagnosticPlot <- function(retval,
                                 out.dir,
                                 verbose){
  multi <- retval[["multi.chains"]] # class hdpSampleMulti
  chains <- hdpx::chains(multi)      # list of hdpSampleChain

  if (verbose) message("Writing HDP diagnostics")
  par(mfrow=c(2,2), mar=c(4, 4, 2, 1))
  pdf(file = file.path(out.dir,"diagnostics.likelihood.pdf"))
  lapply(chains, hdpx::plot_lik, bty = "L")
  dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.numcluster.pdf"))
  # This is the number of raw clusters sampled along each chain
  lapply(chains, hdpx::plot_numcluster, bty = "L")
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.data.assigned.pdf"))
  # This is the number of mutations assigned as a function of
  # the number of raw clusters
  lapply(chains, hdpx::plot_data_assigned, bty = "L")
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.comp.size.pdf"))
  graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
  # Components were already extracted, so this call will work
  hdpx::plot_comp_size(multi, bty="L")
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.signatures.pdf"))
  graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
  # This plots the component (signature) profiles with
  # 95% credibility intervals
  hdpx::plot_comp_distn(multi)
  grDevices::dev.off()

  # TODO, need argument dpindices and col_comp;
  # Need to return the hdp object (perhaps) from RunhdpInternal
  # to get the required values.
  if (FALSE) { # Not finished
    num.dpindices <- length(chains[[1]]@hdp@ppindex)
    hdpx::plot_dp_comp_exposure(
      multi, dpindices = 3:num.dpindices,
      col_comp = myCol[1:ncol(retval$signature)],
      dpnames = colnames(retval$exposure))
  }
}

