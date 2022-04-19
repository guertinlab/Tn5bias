library(ggplot2)
library(lattice)
library(grid)

load('Figure5D_Tn5_PRE_log2_fold.Rdata')

Maskvarplot <- function(){
  ggplot(log_fold, aes(x=factor(Factor), y=Difference, fill = Treatment, color = Treatment)) +
    xlab('') + ylab(expression(log[2]*('Fold Change'))) +
    geom_violin() +
    scale_fill_manual(values=c('#808080', '#FF0000', '#0000FF')) +
    scale_color_manual(values=c('#808080', '#FF0000', '#0000FF')) +
    #scale_x_discrete(position = "top") +
    geom_hline(yintercept=0, linetype="dashed", color = "black") +
    #geom_point(size = 1, shape = 20) +
    theme_classic() + 
    theme(axis.text.x = element_text(angle = 45, vjust = 0.75, hjust=0.75, size = 16, face = 'bold', color = 'black'), 
          axis.title.y=element_text(size=16, color = 'black', face = 'bold'),
          #axis.title.x = element_text(vjust = -0.5),
          axis.text.y = element_text(size = 12, color = 'black', face = 'bold'),
          axis.line = element_line(colour = 'black', size = 1),
          legend.title = element_text(size=0),
          legend.text = element_text(size=16, face = 'bold')) 
}

pdf("Figure5D_PRE_unscaled_test_log2.pdf", width=12, height=7)
Maskvarplot()
dev.off()




bwplot(log_fold[,1] ~ log_fold[,3] | log_fold[,2], aspect = 2.5)

panel.violin.hack <-
function (x, y, box.ratio = 1, box.width = box.ratio/(1 + box.ratio),
     horizontal = TRUE, alpha = plot.polygon$alpha, border =  
plot.polygon$border,
     lty = plot.polygon$lty, lwd = plot.polygon$lwd, col = plot.polygon 
$col,
     varwidth = FALSE, bw = NULL, adjust = NULL, kernel = NULL,
     window = NULL, width = NULL, n = 50, from = NULL, to = NULL,
     cut = NULL, na.rm = TRUE, ...)
{
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



factor="DNase"

x = log_fold
x[,3] = gsub("Tn5 NNNCNN", "seqOutBias", x[,3])
x[,3] = gsub("Tn5 Rules Ensemble", "Rule Ensemble", x[,3])

x[,3] <- as.factor(x[,3])


for (i in c('Rule Ensemble', 'seqOutBias', 'Unscaled' )) {
  x[,3] = relevel(x[,3],ref=i)
}

pdf(paste0(factor,'_bias_width.pdf'), useDingbats = FALSE, width=10.83, height=6)

trellis.par.set(box.umbrella = list(lty = 1, col="#93939380", lwd=2),
                box.rectangle = list(col = '#93939380', lwd=1.6),
                plot.symbol = list(col='#93939380', lwd=1.6, pch ='.'))

#levels(log_fold$Treatment) = c("Unscaled", "Tn5 NNNCNN","Tn5 Rules Ensemble" )

print(bwplot(Difference ~ Treatment | Factor , data = x,
             between=list(y=0.5, x = 0.5),
             scales=list(x=list(draw=TRUE),rot = 45, alternating=c(1,1,1,1),cex=1,font=1),
            #xlab = '',
             ylim = c(-4.5, 3.5),
         #   main = "Pause Index (PI) Ratio",
             ylab =expression("log"[2]~"(deviation from random cleavage)"),
             horizontal =FALSE,  col= 'black',
             aspect = 2,
             par.settings=list(par.xlab.text=list(cex=1.2,font=1),
                               par.ylab.text=list(cex=1.2,font=1),
                               par.main.text=list(cex=1.2, font=1),
                               plot.symbol = list(col='black', lwd=1.6, pch =19, cex = 0.25)),
             strip = function(..., which.panel, bg) {
               bg.col = c("grey90")
               strip.default(..., which.panel = which.panel,
                             bg = rep(bg.col, length = which.panel)[which.panel])
             },
             panel = function(..., box.ratio, col) {
                 panel.abline(h = 0, col = 'grey45', lty = 2)
                 panel.violin.hack(..., col = c("#fff4c7", "#ff5fd2", "#50ffe6"),
                              varwidth = FALSE, box.ratio = box.ratio, outer = FALSE)
                 panel.stripplot(..., col='#54545380', do.out=FALSE, jitter.data=TRUE, amount = 0.2, pch = 16)
                 panel.bwplot(..., pch = '|', do.out = FALSE)
           
       }))
dev.off()
