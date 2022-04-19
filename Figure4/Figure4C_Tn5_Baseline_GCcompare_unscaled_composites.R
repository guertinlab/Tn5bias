load('Figure4C_Tn5_GCcon_baseline.Rdata')

pdf("Figure4C_Tn5_unscaled_BaselineVsGCcon.pdf", width=6, height=6)
plot(V3 ~ V1, data = Tn5_GCcon_baseline, cex = 0, xlab="Baseline Signal", ylab = "Motif GC%", main = "Tn5 baseline vs. motif GC content",
     cex.lab=1.45, cex.axis=1.45)
text(V3 ~ V1, labels=V4,data=Tn5_GCcon_baseline, cex=0.75, font=2)
dev.off()
