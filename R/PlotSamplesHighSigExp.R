#' Plot hdp signature exposure in each sample. This function returns the plot of
#' top 5 samples with the highest exposure to a signature.
#' Each spectrum's title is in the form of: SampleName(Proportion of Signature Assginment)
#' This function is here because it is specific for signature extraction application.
#'
#' @param hdpsample  A \code{\link[hdpx]{hdpSampleChain-class}} or \code{\link[hdpx]{hdpSampleMulti-class}} object including output
#'  from \code{\link[hdpx]{extract_components_from_clusters}}
#' @param retval An object return from \code{\link{extract_ccc_cdc_from_hdp}}
#' @param input.catalog Input spectra catalog as a matrix or
#' in \code{\link[ICAMS]{ICAMS}} format.
#' @param col_comp Colours of each component, from 0 to the max number. If NULL,                          default colors will be used
#' @param incl_numdata_plot Logical - should an upper barplot indicating the number of
#'  data items per DP be included? (Default TRUE)
#' @param ylab_numdata Vertical axis label for numdata plot
#' @param ylab_exp Vertical exis label for exposure plot
#' @param leg.title Legend title
#' @param cex.names Expansion factor for bar labels (dpnames) in exposure plot
#' @param cex.axis Expansion factor for vertical-axis annotation
#' @param mar See ?par
#' @param oma See ?par
#' @importFrom graphics barplot
#' @importFrom reshape2 melt
#' @export
PlotSamplesHighSigExp <- function(retval, hdpsample, input.catalog,
                                  col_comp = NULL,
                                  incl_numdata_plot=F,
                                  ylab_numdata="Number of data items",
                                  ylab_exp="Component exposure",
                                  leg.title="Component", cex.names=0.6,
                                  cex.axis=0.7, mar=c(1, 4, 2, 0.5),
                                  oma=c(1.5, 1.5, 1, 1)){

  # input checks

  exposures <- retval$exposureProbs

  signature <- retval$signature

  ndp <- hdpsample@chains[[1]]@hdp@numdp
  ncomp <- ncol(signature)

  dpnames <- colnames(input.catalog)

  dpindices <- (ndp-length(dpnames)+1):ndp

  if(is.null(col_comp)){
    col_comp <- grDevices::rainbow(ncol(signature), alpha = 1)
  }


  if (!is.numeric(dpindices) | any(dpindices %% 1 != 0) |
      any(dpindices < 1) | any(dpindices > ndp)) {
    stop(paste("dpindices must be integers between 1 and", ndp))
  }


  if(class(incl_numdata_plot) != "logical") {
    stop("incl_numdata_plot must be TRUE or FALSE")
  }

  # save pre-existing par conditions, and reset on exit
  par_old <- par(no.readonly=TRUE)
  on.exit(par(par_old), add=TRUE)

  # Number of data items per DP
  if (class(hdpsample) == "hdpSampleChain") {
    dps <- hdpx::dp(hdpx::final_hdpState(hdpsample))[dpindices]
  } else if (class(hdpsample) == "hdpSampleMulti") {
    dps <- hdpx::dp(hdpx::final_hdpState(hdpx::chains(hdpsample)[[1]]))[dpindices]

  }

  numdata <- sapply(dps, function(x) x@numdata)
  dp_order <- order(numdata, decreasing=TRUE)

  # if incl_numdata_plot TRUE, throw error if one DP has no data associated
  if (incl_numdata_plot & any(numdata == 0)) {
    stop("Can't have incl_numdata_plot TRUE if
         one or more dpindices have no associated data item/s")
  }


  # which components to include in this plot
  inc <- which(rowSums(exposures, na.rm=T)>0)

  num_leg_col <- floor(sqrt(length(inc)))

  ##These are not informative. Removed before submit the manuscript
  if(F){
    if (incl_numdata_plot){
      par(mfrow=c(2, 1), mar=mar, oma=oma, cex.axis=cex.axis, las=2)

      graphics::barplot(numdata[dp_order], main="Mutations of each tumor", col="gray", space=0, border=NA,
                        names.arg="", ylab=ylab_numdata,
                        legend.text=names(inc),
                        args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                                         ncol=num_leg_col))

      graphics::barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc], border=NA,
                        ylim=c(0, 1), names.arg=dpnames[dp_order], ylab=ylab_exp,
                        cex.names=cex.names)
    } else {

      par(cex.axis=cex.axis, las=2)
      # don't understand why legend.text needs rev() here and not in above case,
      # but seems to work?
      graphics::barplot(as.matrix(exposures[inc, dp_order, drop=FALSE]), space=0, col=col_comp[inc],
                        border=NA, ylim=c(0, 1),
                        xlim=c(0, length(dpindices) + num_leg_col + 1),
                        names.arg=dpnames[dp_order],
                        ylab=ylab_exp, cex.names=cex.names,
                        legend.text=rev(names(inc)),
                        args.legend=list(fill=col_comp[inc], bty="n", title=leg.title,
                                         ncol=num_leg_col))
    }

  }

  data.exposures <- t(numdata*t(exposures))
  colnames(data.exposures) <- colnames(input.catalog)
  row.names(data.exposures) <- colnames(retval$signature)
  numdata.cutoff <- 0.5*nrow(input.catalog)##exclude low mutated samples
  data.exposures <- data.exposures[,which(colSums(input.catalog)>numdata.cutoff),drop=F]##exclude extremely low samples
  exposures <- exposures[,which(colSums(input.catalog)>numdata.cutoff),drop=F]
  input.catalog <- input.catalog[,which(colSums(input.catalog)>numdata.cutoff),drop=F]
  x <- barplot(rowSums(data.exposures), las=2,cex.names = 0.8) # Do not plot any axes


  Signature <- Sample <- Exposure <- Tumor <- NULL

  df <- reshape2::melt(data.exposures)
  colnames(df) <- c("Signature","Sample","Exposure")

  if(any(grepl("::",colnames(data.exposures)))){

    df$Tumor <- apply(df,1,function(x){
      x["Tumor"] <- unlist(strsplit(x["Sample"],"::"))[1]
    })
  }else{
    df$Tumor <- "No specific tumor type"
  }

  ##for each signature, plotting the top 5 tumors
  for(i in 1:nrow(data.exposures)){
    dp_order_sig <- order(exposures[i,], decreasing=TRUE)
    this.par <- par(mfrow=c(2, 1), mar=mar, oma=oma, cex.axis=cex.axis, las=1)

    this.col <- rep("grey",length(col_comp))
    this.col[i] <- col_comp[i]


    sig <- row.names(data.exposures)[i]

    plot(exposures[i,]~data.exposures[i,],xlab="Exposure",ylab="Proportion of exposure",main=sig,pch=16)

    #sig.df <- df[df$Signature==sig,]

    on.exit(par(this.par))

    old.par <- par(mfrow = c(6, 1), mar = c(2, 2, 2, 2), oma = c(2, 2, 2, 2))
    on.exit(par(old.par))

    ICAMS::PlotCatalog(ICAMS::as.catalog(signature[,i,drop=FALSE],infer.rownames = T,catalog.type = "counts.signature"))

    max.sample <- min(5, ncol(input.catalog))

    this.catalog <- input.catalog[,dp_order_sig[1:max.sample], drop=FALSE]
    this.prop <- exposures[i,dp_order_sig[1:max.sample]]
    this.prop <- round(this.prop,3)
    colnames(this.catalog) <- paste0(colnames(this.catalog),"(",this.prop,")")
    for (j in 1:ncol(this.catalog)) {
      ICAMS::PlotCatalog(ICAMS::as.catalog(this.catalog[ , j, drop = FALSE]))
    }
  }
}
