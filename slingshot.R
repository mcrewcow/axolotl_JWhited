library(slingshot)
dimred <- axolotl_intact@reductions$umap@cell.embeddings
clustering <- axolotl_intact$RNA_snn_res.1
counts <- as.matrix(axolotl_intact@assays$RNA@counts[axolotl_intact@assays$RNA@var.features, ])



BiocManager::install("DelayedMatrixStats")
set.seed(1)
lineages <- getLineages(data = dimred, clusterLabels = clustering)

lineages

pto<-getLineages(dimred, clustering)
sds<-as.SlingshotDataSet(pto) 
plot(dimred, col = clustering,asp=1) 
lines(sds,type='l',lwd=3)


par(mfrow = c(1, 2))
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
for (i in levels(clustering)) {
    text(mean(dimred[clustering == i, 1]), mean(dimred[clustering == i, 2]), labels = i, font = 2)
}
plot(dimred[, 1:2], col = pal[clustering], cex = 0.5, pch = 16)
lines(sds,type='l',lwd=3)

curves <- getCurvesptoapprox_points = 300, thresh = 0.01, stretch = 0.8, allow.breaks = FALSE, shrink = 0.99)
curves

sds<-as.SlingshotDataSet(curves)

plot(dimred, col = pal[clustering], asp = 1, pch = 16)
lines(sds, lwd = 3, col = "black")
