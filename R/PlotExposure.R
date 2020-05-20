#' Plot exposures in multple plots each with a manageable number of samples.
#'
#' @param exp Exposures as a numerical matrix (or data.frame) with
#'    signatures in rows and samples in columns. Rownames are taken
#'    as the signature names and column names are taken as the
#'    sample IDs. If you want \code{exp} sorted from largest to smallest
#'    use \code{\link{SortExp }}. Do not use column names that start
#'    with multipe underscores. The exposures will often be mutation
#'    counts, but could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proprotions rather than counts.
#'
#' @param num.per.line Number of samples to show in each plot.
#'
#' @param ... Other arguments pased to \code{\link{PlotExposure}}. If \code{ylab}
#'    is not included, it defaults to a value depending on \code{plot.proportion}.
#'    If \code{col} is not supplied the function tries to do something
#'    reasonable.
#'
#' @export
#'
PlotExposureByRange <- function(
  exp,
  num.per.line    = 30,
  plot.proportion = FALSE,
  ...
) {

  args <- list(...)
  ylab <- args$ylab
  if (is.null(ylab)) {
    ylab <- ifelse(plot.proportion,
                   "Proportion of mutations",
                   "Number of mutations")
  }

  n.sample <- ncol(exp)
  num.ranges <- n.sample %/% num.per.line
  size.of.last.range <- n.sample %% num.per.line
  if (size.of.last.range > 0) {
    padding.len <- num.per.line - size.of.last.range
    padding <- matrix(0,nrow = nrow(exp), ncol = padding.len)
    # The column names starting with lots of underscore
    # will not be plotted in the final output.
    colnames(padding) <- paste("_____", 1:ncol(padding), sep = "_")
    exp <- cbind(exp, padding)
    starts <- 0:num.ranges * num.per.line + 1
  } else {
    starts <- 0:(num.ranges - 1) *num.per.line + 1
  }
  ends   <- starts + num.per.line - 1

  plot.legend <- TRUE
  for (i in 1:length(starts)) {
    PlotExposure(exp[ , starts[i]:ends[i]],
                 plot.proportion = plot.proportion,
                 plot.legend    = plot.legend,
                 ...)
    plot.legend <- FALSE
  }
}


#' Plot a single exposure plot
#'
#' @param exp Exposures as a numerical matrix (or data.frame) with
#'    signatures in rows and samples in columns. Rownames are taken
#'    as the signature names and column names are taken as the
#'    sample IDs. If you want \code{exp} sorted from largest to smallest
#'    use \code{\link{SortExp }}. Do not use column names that start
#'    with multipe underscores. The exposures will often be mutation
#'    counts, but could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proprotions rather than counts.
#'
#' @param plot.legend If \code{TRUE} plot a legend.
#'
#' @param ... Parameters passed to \code{\link[graphics]{barplot}}.
#'
#' @export

PlotExposure <-
  function(s.weights, # This is actually the exposure "counts"
           plot.proportion = FALSE,
           plot.legend     = TRUE,
           ...
  ) {

    # note - might be reals > 1, not necessary colSum==1
    s.weights <- as.matrix(s.weights) # in case it is a data frame
    num.sigs  <- nrow(s.weights)
    num.samples <- ncol(s.weights)

    col <- list(...)$col
    if (is.null(col)) {
      if (num.sigs <= 8) {
        col <- # c('skyblue', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'red')
          c('red', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'skyblue')

      } else {
        # lots of signatures; use shaded lines to differentiate
        col <- grDevices::rainbow(num.sigs, alpha = 1)
      }
    }
    if (num.sigs <= 12) {
      p.dense = -1 # will repeat as needed, -1 = solid
      p.angle = 0  # ditto
    } else {
      # lots of signatures; use shaded lines to differentiate
      p.dense = c(-1,35,35,50,50) # will repeat as needed, -1 = solid
      p.angle = c(0,0,135,45,135)  # ditto
    }
    # for legend, we need to put in reverse order. make sure this matches
    # order used in barplot
    num.repeats = ceiling(num.sigs/length(p.dense))
    p.dense.rev = rev(rep(p.dense,num.repeats)[1:num.sigs])
    p.angle.rev = rev(rep(p.angle,num.repeats)[1:num.sigs])

    l.cex = if (num.sigs > 15) .5 else 1 # char expansion for legend (was 0.7)
    direction = 2 # 1=always horizontal, 2=perpendicular to axis

    if (plot.proportion) {
      # Matrix divided by vector goes col-wise, not row-wise, so transpose twice
     plot.what <- t(t(s.weights)/colSums(s.weights))
    } else {
      plot.what <- s.weights
    }

    # ignore column names; we'll plot them separately to make them fit
    bp = barplot(plot.what,
                 las      = 1,
                 yaxt     = 's',
                 xaxt     = 'n', # Do not plot the X axis
                 density  = p.dense,
                 angle    = p.angle,
                 border   = F, # ifelse(num.samples>200,NA,1),
                 cex.main = 1.2,
                 col      = col,
                 ...) # removed xlab, main, ylim, ylab

    # get max y values for plot region, put legend at top right
    dims = par('usr') # c(x.min, x.max, y.min, y.max)
    y.max = dims[4]

    if (plot.legend) {
      # less space between rows (y.intersp), and between box & label (x.intersp)
      # reverse the order, so sig 1 is at bottom (to match bargraph)
      legend.x <- ncol(s.weights) * .7
      legend.y <- y.max * 0.8
      legend(x=legend.x, y=legend.y,
             rev(row.names(s.weights)),
             density=p.dense.rev, angle=p.angle.rev,
             bg=NA, xpd=NA,
             fill=col[num.sigs:1],
             x.intersp=.4, y.intersp=.8,
             bty='n', cex=l.cex * 0.9)
      text(x=legend.x, y = legend.y, "Signature", adj=-0.09)
    }

    # now add sample names, rotated to hopefully fit
    # don't even try to show all if there are too many
    if (num.samples <= 200) {
      if (length(bp)<50) size.adj = .75
      else if (length(bp)<80) size.adj = .65
      else if (length(bp)<100) size.adj = 0.4 # .5
      else if (length(bp)<120) size.adj = .4
      else if (length(bp)<150) size.adj = .3
      else size.adj = .3
      cnames <- colnames(s.weights)
      cnames <- sub("_____.*", "", cnames)
      mtext(cnames, side=1, at=bp, las=direction, cex=size.adj)
    }
   }

#' Sort columns of an exposure matrix from largest to smaller (or vice versa).
#'
#' @param exposures The exposures to sort; columns are samples.
#'
#' @param decreasing If \code{TRUE} sort from largest to smallest.
#'
#' @export
#'
SortExp <- function(exposures, decreasing = TRUE) {
  retval <- exposures[   , order(colSums(exposures), decreasing = decreasing)]
  return(retval)
}

# This one is trivial -- not needed
PDFExposureByRange <- function(path,   # Out file path
                               exp,
                               num.per.line,
                               plot.proportion,
                               col=NULL,
                               main=NULL,
                               xlab=NULL,
                               ylab=NULL
) {
  pdf(path, width = 8.2677, height = 11.6929, # for A4
      onefile = TRUE, useDingbats = FALSE)
  PlotExposureByRange(
    exp             = exp,
    num.per.line    = num.per.line,
    plot.proportion = plot.proportion,
    col = col,  main = main,
    xlab = xlab,
    ylab = ylab)
  dev.off()
}

