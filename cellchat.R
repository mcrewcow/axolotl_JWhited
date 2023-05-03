axolotl <- readRDS("G:/DPayz_sysAct_seurat_20230327.rds")
axolottest <- axolotl
library(CellChat)
library(Seurat)
axolotl_contra <- subset(axolottest, subset = sample == c('contra'))
axolotl_intact <- subset(axolottest, subset = sample == c('intact'))
interaction_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/interaction_input_CellChatDB.csv', row.names = 1)
complex_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/complex_input_CellChatDB.csv', row.names = 1)
cofactor_input <- read.csv(file = 'C://Users/Emil/10X/cellchat/cofactor_input_CellChatDB.csv', row.names = 1)
geneInfo <- read.csv(file = 'C://Users/Emil/10X/cellchat/geneInfo_input_CellChatDB.csv', row.names = 1)
CellChatDB <- list()
CellChatDB$interaction <- interaction_input
CellChatDB$complex <- complex_input
CellChatDB$cofactor <- cofactor_input
CellChatDB$geneInfo <- geneInfo
CellChatDB.mouse <- CellChatDB
features = c('IGF-IR-IGF1R-AMEX60DD003508','HGF.L-HGF-AMEX60DD006229','MET-AMEX60DD006320','CIB84-004984-IGF1-AMEX60DD008414','EGF-AMEX60DD003122',
'EGF-AMEX60DD044302',
'EGF-AMEX60DD051379',
'ERBB1-EGFR-AMEX60DD022047',
'FGFR1-AMEX60DD002586',
'FGFRL1-AMEX60DD046097',
'FGFR2-AMEX60DD053280',
'FGFR-2-FGFR2-AMEX60DD053280',
'KGFR-FGFR2-AMEX60DD053280',
'FGFR3-AMEX60DD046037',
'FGFR4-AMEX60DD027855',
'FGF1-AMEX60DD028820',
'FGF-2-FGF2-AMEX60DD044865',
'FGF3-AMEX60DD005046',
'FGF4-AMEX60DD005045',
'FGF4-FGF6-AMEX60DD005045',
'FGF6-AMEX60DD006587',
'FGF7-AMEX60DD003767',
'FGF8-AMEX60DD052891',
'FGF9-AMEX60DD049338',
'FGF10-AMEX60DD042028',
'FGF11-AMEX60DDU001008401',
'FGF11-AMEX60DDU001008413',
'FGF11-AMEX60DDU001008413',
'FGF12-AMEX60DD003211',
'FGF12-AMEX60DD003212',
'FGF13-AMEX60DD037407',
'FGF14-AMEX60DD048162',
'DBR06-SOUSAS5610051-FGF14-AMEX60DD048159',
'BN2614-LOCUS1-FGF14-AMEX60DD048160',
'FGF16-AMEX60DD037078',
'FGF17-AMEX60DD002799',
'FGF20-AMEX60DD045469',
'FGF18-AMEX60DD046672',
'FGF21-AMEX60DD017409',
'FGF23-AMEX60DD006585',
'INS-AMEX60DD004535',
'INS-AMEX60DD004533',
'INSR-AMEX60DD014048',
'INSR-AMEX60DDU001028323',
'NGF-AMEX60DD008715',
'NGFR-AMEX60DD009964',
'NGFR-AMEX60DD026029',
'TGFBR1-AMEX60DD038999',
'TGFBR2-AMEX60DD022702',
'TGFBR2-AMEX60DD056061',
'TGFBR3-AMEX60DD019130',
'TGFBR3L-TGFBR3-AMEX60DD031769',
'TGFB1I1-AMEX60DD027949',
'TGFB1I1-PXN-AMEX60DD000525',
'TGFB2-AMEX60DD036126',
'TGFB2-AMEX60DD025328',
'TGFB3-AMEX60DDU001011686',
'TGFB3-AMEX60DDU001011687',
'TGFB3-AMEX60DD011162',
'TGFB3-AMEX60DD011161',
'NPY-AMEX60DD009927',
'NPY-AMEX60DD022071',
'PPY-AMEX60DD009929',
'TAC1-AMEX60DD022344',
'NPY1R-AMEX60DD001252',
'NPY1R-AMEX60DD045196',
'NPY2R-AMEX60DD030157',
'NPY2R-AMEX60DD045112',
'NPY5R.L-NPY5R-AMEX60DD045202',
'GPR83.L-AMEX60DD036978',
'GPR83.2.L-GPR83-AMEX60DD049565',
'GPR83.2-GPR83-AMEX60DD049565',
'TACR1-AMEX60DD002928',
'TACR1-AMEX60DD002932',
'TACR1-AMEX60DD002935',
'WNT11-AMEX60DD049747',
'WNT9A-AMEX60DD020504',
'WNT10B-AMEX60DD029981',
'FZD1-AMEX60DD022211',
'FZD10-AMEX60DD000016',
'FZD2-AMEX60DD009878',
'FZD3-AMEX60DD032782',
'LLAP-16351-FZD4-AMEX60DD049634',
'FZD5-AMEX60DD055073',
'FZD6-AMEX60DD040248',
'FZD7-AMEX60DD054983',
'FZD8-AMEX60DD021741',
'FZD9-AMEX60DD054233',
'FZD2-AMEX60DD009877',
'FZD4-AMEX60DD049635',
'LRP5-AMEX60DD004352',
'LRP6-AMEX60DD006973',
'WNT11-AMEX60DD036922')


cellchat <- createCellChat(object = axolotl_contra, group.by = "labels")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
identifyOverExpressedGenes <- function(object, data.use = NULL, group.by = NULL, idents.use = NULL, invert = FALSE, group.dataset = NULL, pos.dataset = NULL, features.name = "features",  only.pos = TRUE, features = NULL, return.object = TRUE,
thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05) {
if (!is.list(object@var.features)) {
stop("Please update your CellChat object via `updateCellChat()`")
}
if (is.null(data.use)) {
X <- object@data.signaling
if (nrow(X) < 3) {stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")}
} else {
X <- data.use
}
if (is.null(features)) {
features.use <- row.names(X)
} else {
features.use <- intersect(features, row.names(X))
}
data.use <- X[features.use,]
data.use <- as.matrix(data.use)
if (is.null(group.by)) {
labels <- object@idents
if (!is.factor(labels)) {
message("Use the joint cell labels from the merged CellChat object")
labels <- object@idents$joint
}
} else {
labels <- object@meta[[group.by]]
}
if (!is.factor(labels)) {
labels <- factor(labels)
}
level.use <- levels(labels)[levels(labels) %in% unique(labels)]
if (!is.null(idents.use)) {
if (invert) {
level.use <- level.use[!(level.use %in% idents.use)]
} else {
level.use <- level.use[level.use %in% idents.use]
}
}
numCluster <- length(level.use)
if (!is.null(group.dataset)) {
labels.dataset <- as.character(object@meta[[group.dataset]])
if (!(pos.dataset %in% unique(labels.dataset))) {
cat("Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
stop()
}
}
my.sapply <- ifelse(
test = future::nbrOfWorkers() == 1,
yes = pbapply::pbsapply,
no = future.apply::future_sapply
)
mean.fxn <- function(x) {
return(log(x = mean(x = expm1(x = x)) + 1))
}
labels <- as.character(labels)
genes.de <- vector("list", length = numCluster)
for (i in 1:numCluster) {
features <- features.use
if (is.null(group.dataset)) {
cell.use1 <- which(labels == level.use[i])
cell.use2 <- base::setdiff(1:length(labels), cell.use1)
} else {
cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
}
# feature selection (based on percentages)
thresh.min <- 0
pct.1 <- round(
x = rowSums(data.use[features, cell.use1, drop = FALSE] > thresh.min) /
length(x = cell.use1),
digits = 3
)
pct.2 <- round(
x = rowSums(data.use[features, cell.use2, drop = FALSE] > thresh.min) /
length(x = cell.use2),
digits = 3
)
data.alpha <- cbind(pct.1, pct.2)
colnames(x = data.alpha) <- c("pct.1", "pct.2")
alpha.min <- apply(X = data.alpha, MARGIN = 1, FUN = max)
names(x = alpha.min) <- rownames(x = data.alpha)
features <- names(x = which(x = alpha.min > thresh.pc))
if (length(x = features) == 0) {
#stop("No features pass thresh.pc threshold")
next
}
# feature selection (based on average difference)
data.1 <- apply(X = data.use[features, cell.use1, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
data.2 <- apply(X = data.use[features, cell.use2, drop = FALSE],MARGIN = 1,FUN = mean.fxn)
FC <- (data.1 - data.2)
if (only.pos) {
features.diff <- names(which(FC > thresh.fc))
} else {
features.diff <- names(which(abs(FC) > thresh.fc))
}
features <- intersect(x = features, y = features.diff)
if (length(x = features) == 0) {
#  stop("No features pass thresh.fc threshold")
next
}
data1 <- data.use[features, cell.use1, drop = FALSE]
data2 <- data.use[features, cell.use2, drop = FALSE]
pvalues <- unlist(
x = my.sapply(
X = 1:nrow(x = data1),
FUN = function(x) {
# return(wilcox.test(data1[x, ], data2[x, ], alternative = "greater")$p.value)
return(wilcox.test(data1[x, ], data2[x, ])$p.value)
}
)
)
pval.adj = stats::p.adjust(
p = pvalues,
method = "bonferroni",
n = nrow(X)
)
genes.de[[i]] <- data.frame(clusters = level.use[i], features = as.character(rownames(data1)), pvalues = pvalues, logFC = FC[features], data.alpha[features,, drop = F],pvalues.adj = pval.adj, stringsAsFactors = FALSE)
}
markers.all <- data.frame()
for (i in 1:numCluster) {
gde <- genes.de[[i]]
if (!is.null(gde)) {
gde <- gde[order(gde$pvalues, -gde$logFC), ]
gde <- subset(gde, subset = pvalues < thresh.p)
if (nrow(gde) > 0) {
markers.all <- rbind(markers.all, gde)
}
}
}
if (only.pos & nrow(markers.all) > 0) {
markers.all <- subset(markers.all, subset = logFC > 0)
}
if (!is.null(group.dataset)) {
markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
markers.all$datasets <- factor(markers.all$datasets, levels = levels(factor(object@meta[[group.dataset]])))
markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
}
markers.all$features <- as.character(markers.all$features)
features.sig <- markers.all$features
object@var.features[[features.name]] <- features.sig
features.name <- paste0(features.name, ".info")
object@var.features[[features.name]] <- markers.all
if (return.object) {
return(object)
} else {
return(markers.all)
}
}
cellchat <- identifyOverExpressedGenes(cellchat, features = features, thresh.p = 0.05)
cellchat@var.features
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.mouse)
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0, raw.use = FALSE, population.size = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
netVisual_aggregate(cellchat, signaling = c('HGF'), layout = "chord")
saveRDS(cellchat, file = "C://Users/Emil/10X/cellchat_total_contra")
netVisual_aggregate(cellchat, signaling = c('Wnt10b'), layout = "chord")
