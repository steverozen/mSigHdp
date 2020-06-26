#' Plot exposures in multiple plots each with a manageable number of samples.
#'
#' @param exposures Exposures as a numerical matrix (or data.frame) with
#'    signatures in rows and samples in columns. Rownames are taken
#'    as the signature names and column names are taken as the
#'    sample IDs. If you want \code{exposures} sorted from largest to smallest
#'    use \code{\link{SortExp}}. Do not use column names that start
#'    with multiple underscores. The exposures will often be mutation
#'    counts, but could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param num.per.line Number of samples to show in each plot.
#'
#' @param ... Other arguments passed to \code{\link{PlotExposure}}. If \code{ylab}
#'    is not included, it defaults to a value depending on \code{plot.proportion}.
#'    If \code{col} is not supplied the function tries to do something
#'    reasonable.
#'
#' @export
#'
PlotExposureByRange <- function(
  exposures,
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

  n.sample <- ncol(exposures)
  num.ranges <- n.sample %/% num.per.line
  size.of.last.range <- n.sample %% num.per.line
  if (size.of.last.range > 0) {
    padding.len <- num.per.line - size.of.last.range
    padding <- matrix(0,nrow = nrow(exposures), ncol = padding.len)
    # The column names starting with lots of underscore
    # will not be plotted in the final output.
    colnames(padding) <- paste("_____", 1:ncol(padding), sep = "_")
    exposures <- cbind(exposures, padding)
    starts <- 0:num.ranges * num.per.line + 1
  } else {
    starts <- 0:(num.ranges - 1) *num.per.line + 1
  }
  ends   <- starts + num.per.line - 1

  plot.legend <- TRUE
  for (i in 1:length(starts)) {
    PlotExposure(exposures[ , starts[i]:ends[i]],
                 plot.proportion = plot.proportion,
                 plot.legend    = plot.legend,
                 ...)
    plot.legend <- FALSE
  }
}


#' Plot a single exposure plot
#'
#' @param exposures Exposures as a numerical matrix (or data.frame) with
#'    signatures in rows and samples in columns. Rownames are taken
#'    as the signature names and column names are taken as the
#'    sample IDs. If you want \code{exp} sorted from largest to smallest
#'    use \code{\link{SortExp}}. Do not use column names that start
#'    with multiple underscores. The exposures will often be mutation
#'    counts, but could also be e.g. mutations per megabase.
#'
#' @param plot.proportion Plot exposure proportions rather than counts.
#'
#' @param plot.legend If \code{TRUE} plot a legend.
#'
#' @param ... Parameters passed to \code{\link[graphics]{barplot}}.
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom graphics barplot legend mtext par text
#'
#' @export

PlotExposure <-
  function(exposures, # This is actually the exposure "counts"
           plot.proportion = FALSE,
           plot.legend     = TRUE,
           ...
  ) {

    # note - might be reals > 1, not necessary colSum==1
    exposures <- as.matrix(exposures) # in case it is a data frame
    num.sigs  <- nrow(exposures)
    num.samples <- ncol(exposures)

    three.dots <- list(...)
    if (is.null(three.dots$col)) {
      if (num.sigs <= 8) {
        three.dots$col <- # c('skyblue', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'red')
          c('red', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'skyblue')
      } else {
        # lots of signatures; use shaded lines to differentiate
        three.dots$col <- grDevices::rainbow(num.sigs, alpha = 1)
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
     plot.what <- t(t(exposures)/colSums(exposures))
    } else {
      plot.what <- exposures
    }

    # ignore column names; we'll plot them separately to make them fit
    bp = do.call(
      barplot,
      args = c(list(height = plot.what,
                    las      = 1,
                    yaxt     = 's',
                    xaxt     = 'n', # Do not plot the X axis
                    density  = p.dense,
                    angle    = p.angle,
                    border   = F, # ifelse(num.samples>200,NA,1),
                    cex.main = 1.2),
               three.dots))

    # get max y values for plot region, put legend at top right
    dims = par('usr') # c(x.min, x.max, y.min, y.max)
    y.max = dims[4]

    if (plot.legend) {
      # less space between rows (y.intersp), and between box & label (x.intersp)
      # reverse the order, so sig 1 is at bottom (to match the stacked bar graph)
      legend.x <- ncol(exposures) * .7   # Nanhai, we could pass in legend.x and legend.y as optional arguments from the caller
      legend.y <- y.max * 0.8
      legend(x         = legend.x,
             y         = legend.y,
             legend    = rev(row.names(exposures)),
             density   = p.dense.rev,
             angle     = p.angle.rev,
             bg        = NA,
             xpd       = NA,
             fill      = three.dots$col[num.sigs:1],
             x.intersp = .4,
             y.intersp = .8,
             bty       = 'n',
             cex       = l.cex * 0.9)
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
      cnames <- colnames(exposures)
      cnames <- sub("_____.*", "", cnames)
      mtext(cnames, side=1, at=bp, las=direction, cex=size.adj)
    }

    return(bp)

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
