library(ggplot2)
library(lattice)
library(grid)
setwd('../output')
################################################################################
################################################################################
################################################################################
load('../data/Figure5B_DNase_PRE_NNNCNN_log_fold_change_mk2.Rdata')






x = DNase_PRE_NNNCNN_log_fold_change

panel.violin.hack <-
  function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio),
            horizontal = TRUE, alpha = plot.polygon$alpha, border =  
              plot.polygon$border,
            lty = plot.polygon$lty, lwd = plot.polygon$lwd, col = plot.polygon 
            $col,
            varwidth = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
            window = NULL, width = NULL, n = 50, from = NULL, to = NULL,
            cut = NULL, na.rm = TRUE, ...){
    if (all(is.na(x) | is.na(y)))
      return()
    x <- as.numeric(x)
    y <- as.numeric(y)
    plot.polygon <- trellis.par.get("plot.polygon")
    darg <- list()
    darg$bw <- bw
    darg$adjust <- adjust
    darg$kernel <- kernel
    darg$window <- window
    darg$width <- width
    darg$n <- n
    darg$from <- from
    darg$to <- to
    darg$cut <- cut
    darg$na.rm <- na.rm
    my.density <- function(x) {
      ans <- try(do.call("density", c(list(x = x), darg)),
                 silent = TRUE)
      if (inherits(ans, "try-error"))
        list(x = rep(x[1], 3), y = c(0, 1, 0))
      else ans
    }
    numeric.list <- if (horizontal)
      split(x, factor(y))
    else split(y, factor(x))
    levels.fos <- as.numeric(names(numeric.list))
    d.list <- lapply(numeric.list, my.density)
    dx.list <- lapply(d.list, "[[", "x")
    dy.list <- lapply(d.list, "[[", "y")
    max.d <- sapply(dy.list, max)
    if (varwidth)
      max.d[] <- max(max.d)
    xscale <- current.panel.limits()$xlim
    yscale <- current.panel.limits()$ylim
    height <- box.width
    if (horizontal) {
      for (i in seq_along(levels.fos)) {
        if (is.finite(max.d[i])) {
          pushViewport(viewport(y = unit(levels.fos[i],
                                         "native"), height = unit(height, "native"),
                                yscale = c(max.d[i] * c(-1, 1)), xscale = xscale))
          grid.polygon(x = c(dx.list[[i]], rev(dx.list[[i]])),
                       y = c(dy.list[[i]], -rev(dy.list[[i]])),  
                       default.units = "native",
                       # this is the point at which the index is added
                       gp = gpar(fill = col[i], col = border, lty = lty,
                                 lwd = lwd, alpha = alpha))
          popViewport()
        }
      }
    }
    else {
      for (i in seq_along(levels.fos)) {
        if (is.finite(max.d[i])) {
          pushViewport(viewport(x = unit(levels.fos[i],
                                         "native"), width = unit(height, "native"),
                                xscale = c(max.d[i] * c(-1, 1)), yscale = yscale))
          grid.polygon(y = c(dx.list[[i]], rev(dx.list[[i]])),
                       x = c(dy.list[[i]], -rev(dy.list[[i]])),  
                       default.units = "native",
                       # this is the point at which the index is added
                       gp = gpar(fill = col[i], col = border, lty = lty,
                                 lwd = lwd, alpha = alpha))
          popViewport()
        }
      }
    }
    invisible()
  }


x[,3] = as.character(x[,3])
x[which(x[,3] == 'Rules Ensemble'),3] = 'Rule Ensemble'

x[,3] <- as.factor(x[,3])
x[,3] = factor(x[,3], levels = c('Unscaled', 'seqOutBias', 'Rule Ensemble'))

#
length(which(x$Difference < 0.8 & x$Difference > -0.80 & x$Treatment == 'seqOutBias')) / length(which(x$Treatment == 'seqOutBias'))
length(which(x$Difference < 2.1 & x$Difference > -2.1)) / length(x$Treatment)



#Figure out which TFs are most improved by rule ensemble based on variance:
RE_var = NULL
sOB_var = NULL
for (i in 1:length(unique(x$Factor))) {
  RE_var[i] = var(x$Difference[which(x$Factor == unique(x$Factor)[i] & x$Treatment == 'Rule Ensemble')])
  sOB_var[i] = var(x$Difference[which(x$Factor == unique(x$Factor)[i] & x$Treatment == 'seqOutBias')])
}

var_compare = cbind(unique(x$Factor), RE_var, sOB_var)
var_compare = as.data.frame(var_compare)
var_compare[,2:3] = sapply(var_compare[,2:3],as.numeric)
var_compare[,4] = var_compare[,3] - var_compare[,2]
#These are the top 9 improved TFs
var_compare$V1[sort(var_compare$V4, decreasing = TRUE, index.return = TRUE)[[2]]][1:9]

#Subset x to only include top 9 and make supplemental for other 9
supp_x = x[-c(which(x$Factor %in% var_compare$V1[sort(var_compare$V4, decreasing = TRUE, index.return = TRUE)[[2]]][1:9])),]
x = x[which(x$Factor %in% var_compare$V1[sort(var_compare$V4, decreasing = TRUE, index.return = TRUE)[[2]]][1:9]),]


log2plot = bwplot(Difference ~ Treatment | Factor , data = x,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-2.10, 2.10),
                  ylab =expression("log"[2]~frac("Unbiased Signal","Model Signal Output")),
                  horizontal =FALSE,  col= 'black',
                  aspect = 2,
                  par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                                    par.ylab.text=list(cex=1.2,font=1),
                                    par.main.text=list(cex=1.2, font=1),
                                    plot.symbol = list(col='black', lwd=0, pch =19, cex = 0.4)),
                  par.strip.text = list(cex = 1.2),
                  strip = function(..., which.panel, bg) {
                    bg.col = c("white")
                    strip.default(..., which.panel = which.panel,
                                  bg = rep(bg.col, length = which.panel)[which.panel])
                  },
                  panel = function(..., box.ratio, col) {
                    panel.abline(h = 0, col = 'grey45', lty = 2)
                    panel.violin.hack(..., col = c("#00ffaf", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })

pdf('Figure5B_DNase_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()



log2plot = bwplot(Difference ~ Treatment | Factor , data = supp_x,
                  between=list(y=0.5, x = 0.5),
                  scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1.2,font=1),
                  #xlab = '',
                  ylim = c(-2.10, 2.10),
                  ylab =expression("log"[2]~frac("Unbiased Signal","Model Signal Output")),
                  horizontal =FALSE,  col= 'black',
                  aspect = 2,
                  par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                                    par.ylab.text=list(cex=1.2,font=1),
                                    par.main.text=list(cex=1.2, font=1),
                                    plot.symbol = list(col='black', lwd=0, pch =19, cex = 0.4)),
                  par.strip.text = list(cex = 1.2),
                  strip = function(..., which.panel, bg) {
                    bg.col = c("white")
                    strip.default(..., which.panel = which.panel,
                                  bg = rep(bg.col, length = which.panel)[which.panel])
                  },
                  panel = function(..., box.ratio, col) {
                    panel.abline(h = 0, col = 'grey45', lty = 2)
                    panel.violin.hack(..., col = c("#00ffaf", "#FF0000", "#0000FF"),
                                      varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                    panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                    panel.bwplot(..., pch = '|', do.out = FALSE)
                    
                  })

pdf('Supplemental_Figure5B_DNase_log2_comparison.pdf', useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939300", lwd=2),
                box.rectangle = list(col = '#93939300', lwd=1.6),
                plot.symbol = list(col='#93939300', lwd=1.6, pch ='.'))

################################################################################
print(update(log2plot, layout=c(9,1)))
dev.off()
