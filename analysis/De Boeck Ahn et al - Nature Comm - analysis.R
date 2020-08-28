#De Boeck Ahn et al - Nature Communications 2020 - R 3.6.3
#working directory contains folder with CellRanger aggr output of 6 libraries (3 pcDNA and 3 IL-33 xenografts)

#Load libraries
library(Seurat)
library(dplyr)
library(ggplot2)
library(viridis)

#Load data
fdata<-Read10X("filtered_feature_bc_matrix_U87FULL/")
U87<-CreateSeuratObject(counts=fdata, min.cells=3, min.features=200, project="CON", names.field = 2, names.delim = "-")

#Add metadata for conditions and sequencing batches
samplename = U87@meta.data$orig.ident
batchid = rep("pcDNA-1",length(samplename))
batchid[samplename %in% c("2")] = "pcDNA-2"
batchid[samplename %in% c("3")] = "pcDNA-3"
batchid[samplename %in% c("4")] = "IL-33-1"
batchid[samplename %in% c("5")] = "IL-33-2"
batchid[samplename %in% c("6")] = "IL-33-3"
names(batchid) = rownames(U87@meta.data)
U87 <- AddMetaData(
  object = U87,
  metadata = batchid,
  col.name = "batchid")
xeno = rep("CON",length(samplename))
xeno[samplename %in% c("2")] = "CON"
xeno[samplename %in% c("3")] = "CON"
xeno[samplename %in% c("4")] = "IL-33"
xeno[samplename %in% c("5")] = "IL-33"
xeno[samplename %in% c("6")] = "IL-33"
names(xeno) = rownames(U87@meta.data)
U87 <- AddMetaData(
  object = U87,
  metadata = xeno,
  col.name = "xeno")

#Run Seurat (adapted from:https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html)
U87[["percent.mt"]] <- PercentageFeatureSet(U87, pattern = "^mt-")
RPS.genes <- grep(pattern = "^Rps", x = rownames(x = U87@assays$RNA), value = TRUE)
RPL.genes <- grep(pattern = "^Rpl", x = rownames(x = U87@assays$RNA), value = TRUE)
MRPS.genes <- grep(pattern = "^Mrps", x = rownames(x = U87@assays$RNA), value = TRUE)
MRPL.genes <- grep(pattern = "^Mrpl", x = rownames(x = U87@assays$RNA), value = TRUE)
VlnPlot(U87, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
U87 <- subset(U87, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
U87 <- NormalizeData(U87, normalization.method = "LogNormalize", scale.factor = 10000)
U87 <- FindVariableFeatures(U87, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(U87)
U87 <- ScaleData(U87, features = all.genes)
U87 <- RunPCA(U87, features = VariableFeatures(object = U87))
hU87<-U87
U87 <- U87 %>%
  RunUMAP(reduction = "pca", dims = 1:15) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution = 0.3)
UMAPPlot(U87, pt.size=0.5, label=TRUE, label.size=8)+NoLegend()+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
UMAPPlot(U87, pt.size=0.5, group.by="xeno", cols=c('#00BFC4', '#F8766D'))+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
U87markers <- FindAllMarkers(U87, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

#Evaluate batch effects with Seurat integrate (adapted from:https://satijalab.org/seurat/v3.1/immune_alignment.html)
sU87<-U87
sU87 <- SplitObject(sU87, split.by = "xeno")
sU87 <- lapply(X = sU87, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})
U87.anchors <- FindIntegrationAnchors(object.list = sU87, dims = 1:20)
U87.combined <- IntegrateData(anchorset = U87.anchors, dims = 1:20)
DefaultAssay(U87.combined) <- "integrated"
U87.combined <- ScaleData(U87.combined, verbose = FALSE)
U87.combined <- RunPCA(U87.combined, npcs = 30, verbose = FALSE)
U87.combined <- RunUMAP(U87.combined, reduction = "pca", dims = 1:20)
U87.combined <- FindNeighbors(U87.combined, reduction = "pca", dims = 1:20)
U87.combined <- FindClusters(U87.combined, resolution = 0.5)
UMAPPlot(U87.combined, pt.size=0.5, label=TRUE, label.size=8)+NoLegend()+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
UMAPPlot(U87.combined, pt.size=0.5, group.by="xeno", cols=c('#00BFC4', '#F8766D'))+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')

#Evaluate batch effects with Harmony (adapted from:http://htmlpreview.github.io/?https://github.com/immunogenomics/harmony/blob/master/docs/SeuratV3.html)
library(Harmony)
p1 <- DimPlot(object = hU87, reduction = "pca", pt.size = .1, group.by = "xeno")
p2 <- VlnPlot(object = hU87, features = "PC_1", group.by = "xeno", pt.size = .1)
plot_grid(p1,p2)
U87 <- U87 %>%
  RunHarmony("xeno", plot_convergence = TRUE)
harmony_embeddings <- Embeddings(U87, 'harmony')
harmony_embeddings[1:5, 1:5]
p1 <- DimPlot(object = hU87, reduction = "harmony", pt.size = .1, group.by = "xeno")
p2 <- VlnPlot(object = hU87, features = "harmony_1", group.by = "xeno", pt.size = .1)
plot_grid(p1,p2)
ElbowPlot(hU87, reduction='harmony')
hU87 <- hU87 %>%
  RunUMAP(reduction = "harmony", dims = 1:15) %>%
  FindNeighbors(reduction = "harmony", dims = 1:15) %>%
  FindClusters(resolution = 0.375)

UMAPPlot(hU87, pt.size=0.5, label=TRUE, label.size=8)+NoLegend()+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
UMAPPlot(hU87, pt.size=0.5, group.by="xeno", cols=c('#00BFC4', '#F8766D'))+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')

#Remove the cluster that is CD45 (PTPRC) negative
U87<-subset(U87, idents=10, invert=TRUE)
U87 <- FindVariableFeatures(U87, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(U87)
U87 <- ScaleData(U87, features = all.genes)
U87 <- RunPCA(U87, features = VariableFeatures(object = U87))
U87 <- U87 %>%
  RunUMAP(reduction = "pca", dims = 1:15) %>%
  FindNeighbors(reduction = "pca", dims = 1:15) %>%
  FindClusters(resolution = 0.375)
UMAPPlot(U87, pt.size=0.5, label=TRUE, label.size=8)+NoLegend()+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
UMAPPlot(U87, pt.size=0.5, group.by="xeno", cols=c('#00BFC4', '#F8766D'))+theme(axis.line = element_line(color = "black", size = 1, linetype = "solid"), axis.text.x = element_blank(),axis.text.y = element_blank(), axis.ticks = element_blank()) + xlab('UMAP 1')+ylab('UMAP 2')
U87markers <- FindAllMarkers(U87, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
FeaturePlot(U87, features = 'Gpr34', pt.size=0.5)+ scale_color_viridis(option = 'D') + NoAxes() + theme(plot.title = element_text(size= 36, face='italic', hjust=0.55))

#convert to loom format; export to scanpy (workflow adapted from: https://scanpy.readthedocs.io/en/stable/tutorials.html)
library(loomR)
U87@graphs <- list()
U87loom<-as.loom(u87, file='U87.loom')

#SCORPIUS trajectory analysis of monocytes, differentiating monocytes, and glioma-associated BMDM (adapted from:https://github.com/rcannood/SCORPIUS/blob/master/vignettes/seurat.md)
library(SCORPIUS)
library(plyr)
macro<-subset(U87, idents=c(3,7,11))
srt<-subset(macro, subset=xeno=='IL-33')
srt$seurat_clusters<-factor(srt$seurat_clusters)
expression <- t(as.matrix(srt@assays$RNA@data))
group_name <- srt@meta.data$seurat_clusters
group_name<-revalue(group_name, c("3"="BMDM2", "7"="Monocytes", "11"="BMDM3"))
space <- reduce_dimensionality(expression, "spearman", ndim = 3)
draw_trajectory_plot(space, progression_group = group_name, contour = TRUE)
traj <- infer_trajectory(space)
draw_trajectory_plot(
  space, 
  progression_group = group_name,
  path = traj$path,
  contour = TRUE
)
gimp <- gene_importances(expression, traj$time, num_permutations = 0, num_threads = 6)
gene_sel <- gimp[1:50,]
expr_sel <- expression[,gene_sel$gene]
draw_trajectory_heatmap(expr_sel, traj$time, group_name, show_labels_row = TRUE)
modules <- extract_modules(scale_quantile(expr_sel), traj$time, verbose = FALSE)
draw_trajectory_heatmap(expr_sel, traj$time, group_name, modules, show_labels_row = TRUE, color=viridis(5), fontsize=14)

#Assign identitites
fullid <- c("Microglia1", "BMDM1", "Microglia2", "BMDM2", "Microglia3", "Microglia4", "Microglia5", "Monocytes", "Microglia6", "NK Cells","Microglia7", "BMDM3", "BMDM4", "Granulocytes")
names(fullid) <- levels(U87)
U87 <- RenameIdents(U87, fullid)

#subset microglia clusters
micro<-subset(U87, idents=c(0,2,4,5,6,8,10))
newid <- c("Microglia1", "Microglia2", "Microglia3", "Microglia4", "Microglia5", "Microglia6", "Microglia7")
names(newid) <- levels(micro)
micro <- RenameIdents(micro, newid)
microgliamarkers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10micro <- microgliamarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(micro, features = top10micro$gene)+ theme(text = element_text(size = 10))

#subset macrophage clusters
macro<-subset(U87, idents=c(1,3,11,12))
newid <- c("BMDM1", "BMDM2", "BMDM3", "BMDM4")
names(newid) <- levels(macro)
macro <- RenameIdents(macro, newid)
macromarkers <- FindAllMarkers(macro, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top10macro <- macromarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(macro, features = top10macro$gene)+ theme(text = element_text(size = 10))

#Gene ontology analysis (adapted from:https://github.com/mahmoudibrahim/genesorteR/wiki/From-Cluster-to-Pathway-Enrichment-in-Large-scRNA-Seq-Data)
library(genesorteR)
library(msigdbr)
library(clusterProfiler)
sg = genesorteR::sortGenes(U87@assays$RNA@data, Idents(U87))
msig = msigdbr::msigdbr(species = "Mus musculus", category = "C5")
msig = data.frame(msig$gs_name[which(msig$gs_subcat == "CP:CC")], msig$gene_symbol[which(msig$gs_subcat == "CP:CC")])
topMarkerGenes = genesorteR::plotTopMarkerHeat(sg, top_n = 200, outs = TRUE, plotheat = FALSE)
cp = clusterProfiler::compareCluster(geneCluster = topMarkerGenes, fun = "enricher", TERM2GENE = msig)
write.csv(cp,'goanalysis.csv')
clusterProfiler::dotplot(cp, showCategory = 4, font.size = 8)