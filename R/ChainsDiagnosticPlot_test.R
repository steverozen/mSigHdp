#' Diagnostic plot for a hdpSampleMulti object
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
#' @param input.catalog ground truth catalog
#'
#' @inheritParams AnalyzeAndPlotretval
#'
#' @importFrom grDevices dev.off pdf
#'
#' @importFrom graphics par
#'
#' @export



ChainsDiagnosticPlotMo <- function(retval,
                                   input.catalog,
                                   out.dir,
                                   verbose){
  multi <- retval[["multi.chains"]] # class hdpSampleMulti
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

  grDevices::pdf(file = file.path(out.dir,"diagnostics.comp.size.pdf"))
  graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))
  hdpx::mo_plot_comp_size(retval = retval$diagnostic.retval)
  grDevices::dev.off()


  grDevices::pdf(file = file.path(out.dir,"diagnostics.signatures.pdf"))
  graphics::par(mfrow=c(8, 1), mar = c(1, 1, 1, 1))
  # This plots the component (signature) profiles with
  # 95% credibility intervals
  hdpx::mo_plot_comp_distn_with_credint(retval = retval$diagnostic.retval)
  grDevices::dev.off()

  grDevices::pdf(file = file.path(out.dir,"diagnostics.hdp.signature.exposure.each.sample.pdf"))
  myCol <- grDevices::rainbow(ncol(retval$signature), alpha = 1)
  graphics::par(mfrow=c(1,1), mar=c(5, 4, 4, 2))

  hdpx::mo_plot_sig_exposure_for_dp(retval           = retval$diagnostic.retval,
                                    hdpsample        = multi,
                                    input.catalog    = input.catalog)

  grDevices::dev.off()

  #grDevices::pdf(file = file.path(out.dir,"tsne.sig.vs.tumortype.pdf"))
  #graphics::par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
  #hdpx::plot_tsne_sigs_tumortype(exposure = retval$exposure)
  #grDevices::dev.off()

  #grDevices::pdf(file = file.path(out.dir,"pca.sig.vs.tumortype.pdf"))
  #graphics::par(mfrow=c(1,1), mar=c(4, 4, 2, 1))
  #hdpx::plot_pca_sigs_tumortype(exposure = retval$exposure)
  #grDevices::dev.off()
}

