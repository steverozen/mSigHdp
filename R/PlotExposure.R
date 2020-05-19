### Plot exposures broken up by ranges
PlotExposureByRange <- function(spectrum,
                             sigs,
                             exp,
                             ranges, # A list of vectors
                             col=NULL,
                             main=NULL,
                             xlab=NULL,
                             ylab=NULL,
                             mfrow=c(3,1) # 3 exposure graphs per page by default
) {

  par(
    mfrow=mfrow,
    # mar  =c(1.5,1,2,2), # space between plot and figure, for axis labels etc
    oma  =rep(3, 4) # outer margin between figure and device.
  )
  if (is.null(ylab)) ylab <- 'Number of mutations'
  first=T
  for (range in ranges) {
    if (first) {
      PlotExposure(exp[ ,range], signatures=sigs,
                     input.genomes = spectrum[ , range],
                     plot.proprtion = F,
                     col=col,
                     main=main,
                     ylab=ylab, xlab=xlab)
      first=F
    } else {
      PlotExposure(exp[ ,range], signatures=sigs,
                     input.genomes = spectrum[ , range],
                     plot.proprtion = F, plot.legend=F,
                     main=main,
                     col=col,
                     ylab=ylab, xlab=xlab)
    }
  }
}


PDFExposureByRange <- function(path,   # Out file path
                            spectrum,
                            sigs,
                            exp,
                            ranges, # A list of vectors
                            col=NULL,
                            main=NULL,
                            xlab=NULL,
                            ylab=NULL
) {
  pdf(path, width=8.2677, height=11.6929, # for A4
      onefile=T, useDingbats=F)
  PlotExposureByRange(spectrum, sigs, exp, ranges, col, main, xlab, ylab)
  dev.off()
}


PlotExposure <-
  function(s.weights, # This is actually the exposure "counts"
           # (or floats approximating the exposure counts)
           signames=NULL,
           scale.num=NULL,
           signatures=NULL,
           input.genomes=NULL,
           plot.proprtion=T,
           plot.legend=T,
           ylim=NULL,
           main=NULL,
           ylab=NULL,
           xlab=NULL,
           col=NULL
  ) {

    # note - might be reals > 1, not necessary colSum==1
    s.weights <- as.matrix(s.weights) # in case it is a data frame
    signature_proportions <- t(s.weights)
    num.sigs = dim(s.weights)[1]
    num.samples = dim(s.weights)[2]

    if (is.null(col)) {
      if (num.sigs <= 8) {
        col = # c('skyblue', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'red')
          c('red', 'black', 'grey', 'yellow', 'blue', 'brown', 'green4', 'skyblue')

      } else {
        # lots of signatures; use shaded lines to differentiate
        col = rainbow(num.sigs)
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

    # add names (if not already set as row.names in the original input frame)
    # for sorting. (needs a "Sample" column in the signature_proportions frame)
    # if (length(colnames(s.weights))==0) colnames(s.weights) = signature_proportions$Sample
    if (is.null(scale.num)) {
      ylabel = '# mutations'
    } else {# show as rate instead of count
      if (scale.num < 3000 || scale.num > 3300) {
        warning('assuming "scale.num" should be divisor for human genome. using 3000')
        scale.num = 3000
      }
      s.weights = s.weights/scale.num
      ylabel = '# mutations/Mbase'
    }
    l.cex = if (num.sigs > 15) .5 else 1 # char expansion for legend (was 0.7)
    direction = 2 # 1=always horizontal, 2=perpendicular to axis

    # if we weights file and counts file have samples in different order
    if (!all(colnames(s.weights) == colnames(input.genomes))) {
      input.genomes = input.genomes[,colnames(s.weights)]
      warning('weights file and counts file are ordered differently; re-ordering counts.')
    }

    # ignore column names; we'll plot them separately to make them fit
    bp = barplot(s.weights,
                 ylim=ylim,
                 las=1,
                 col=col,
                 ylab=ylab,
                 yaxt='s',
                 xaxt='n',
                 xlab=xlab,
                 density=p.dense, angle=p.angle,
                 border=ifelse(num.samples>200,NA,1),
                 main=main, cex.main=1.2)

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
      text(x=legend.x, y = legend.y, "Mutational signature", adj=-0.09)
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
      mtext(colnames(s.weights), side=1, at=bp, las=direction, cex=size.adj)
    }

    if (plot.proprtion) {
      # add proportion panel; eg col should sum() to 1. matrix divided by
      # a vector goes col-wise, not row-wise, so transpose twice :(
      bp = barplot(t(t(s.weights)/colSums(s.weights)), las=1, col=col, ylab='Proportion',
                   density=p.dense, angle=p.angle.rev,
                   main='', axisnames=F, border=NA)
    }

  }
