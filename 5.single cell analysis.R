pacman::p_load('dplyr','data.table','fst','plyr','tidyr','Hmisc','tableone','zoo','mgcv','splines',BiocManager,stringr,dplyr,progress,survival,openxlsx,plotRCS,SingleR,Matrix,Seurat)

rm(list = ls())
dir <- 'D:\5.single_cell\\IBD'
setwd(dir)

folder_list <- list.dirs(path = dir, full.names = TRUE, recursive = FALSE)
seurat_list <- list()
for (folder in folder_list) {
  print(folder)
  data <- Read10X(data.dir = folder)
  seurat_obj <- CreateSeuratObject(counts = data, project = basename(folder), min.cells = 3, min.features = 200)
  seurat_list[[basename(folder)]] <- seurat_obj
}
View(seurat_list)
seurat_merge <- merge(seurat_list[[1]], y = seurat_list[c(2:18)], add.cell.ids = names(seurat_list))
seurat <- seurat_merge
allgene <- rownames(seurat)

seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

seurat <- subset(seurat, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 20 & nCount_RNA < 30000)
VlnPlot(seurat, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
seurat <- NormalizeData(seurat, normalization.method = "LogNormalize", scale.factor = 10000)

seurat<- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = 2000)
top10 <- head(VariableFeatures(seurat), 10)
plot1 <- VariableFeaturePlot(seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot1
plot2

all.genes <- rownames(seurat)
seurat <- ScaleData(seurat, features = all.genes)
GetAssayData(seurat,slot="counts",assay="RNA") 
GetAssayData(seurat,slot="data",assay="RNA")
GetAssayData(seurat,slot="scale.data",assay="RNA")

seurat <- RunPCA(seurat, features = VariableFeatures(object = seurat))
plot1 <- DimPlot(seurat, reduction = "pca", group.by="orig.ident")
plot2 <- ElbowPlot(seurat, ndims=20, reduction="pca") 
plot1+plot2

VizDimLoadings(seurat, dims = 1:2, reduction = "pca")
DimPlot(seurat, reduction = "pca") + NoLegend()
DimHeatmap(seurat, dims = 1:10, cells = 500, balanced = TRUE)
ElbowPlot(seurat,ndims = 40)

pc.num <- 1:15
seurat <- FindNeighbors(seurat, dims = pc.num)
seurat <- FindClusters(seurat, resolution = 0.5)
table(seurat@meta.data$seurat_clusters)
metadata <- seurat@meta.data
cell_cluster <- data.frame(cell_ID=rownames(metadata), cluster_ID=metadata$seurat_clusters)

seurat = RunTSNE(seurat, dims = pc.num)
plot1 = DimPlot(seurat, reduction = "tsne") 
plot1

seurat <- RunUMAP(seurat, dims = pc.num)
plot2 = DimPlot(seurat, reduction = "umap") 
plot2

b <- FetchData(seurat, vars = c("CCL20",	'IL18','IL12B','TNFRSF9', 'CCL7','IL1RL1','IL19','CCL13','CD72','EPHB4','SEPT8','IFNG','CD6','LY9','IL10RA','IGLC2','MXRA8','ITGAV','CXCL9','NOV','RSPO3','NOS3'))
colSums(b)

seurat.merg <- JoinLayers(seurat)
head(seurat.merg@meta.data)
table(seurat.merg$orig.ident)
sum(table(Idents(seurat.merg)))
data <- GetAssayData(seurat.merg,slot='counts',assay='RNA')
data[1:5,1:5]
seurat.merg@meta.data[["sample"]] <-seurat.merg@meta.data[["orig.ident"]]
clinical <- c("GSM6614348_HC-1" = "HC",
              "GSM6614349_HC-2" = "HC",
              "GSM6614350_HC-3" = "HC",
              "GSM6614351_HC-4" = "HC",
              "GSM6614352_HC-5"="HC",
              "GSM6614353_HC-6"="HC",
              "GSM6614354_UC-1" = "IBD",
              "GSM6614355_UC-2" = "IBD",
              "GSM6614356_UC-3" = "IBD",
              "GSM6614357_UC-4" = "IBD",
              "GSM6614358_UC-5"="IBD",
              "GSM6614359_UC-6"="IBD",
              "GSM6614360_CD-1" = "IBD",
              "GSM6614361_CD-2" = "IBD",
              "GSM6614362_CD-3" = "IBD",
              "GSM6614363_CD-4" = "IBD",
              "GSM6614364_CD-5"="IBD",
              "GSM6614365_CD-6"="IBD")
seurat.merg[['clinical']] <- unname(clinical[seurat.merg@meta.data$sample])
Idents(seurat.merg) <- "clinical"
table(seurat.merg$clinical)

a <- FetchData(seurat.merg, vars = c("CCL20",	'IL18','IL12B','TNFRSF9', 'CCL7','IL1RL1','IL19','CCL13','CD72','EPHB4','SEPT8','IFNG','CD6','LY9','IL10RA','IGLC2','MXRA8','ITGAV','CXCL9','NOV','RSPO3','NOS3'))
colSums(is.na(a))
colSums(a)

VlnPlot(seurat.merg, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeaturePlot(subset(seurat.merg,clinical=="HC"), reduction = "tsne",
            features = c("CD6"))
FeaturePlot(subset(seurat.merg,clinical=="IBD"), reduction = "tsne",
            features = c("CD6"))

library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()
hpca.se

testdata <- GetAssayData(seurat.merg, slot="data")
testdata[1:5,1:5]
clusters <- seurat.merg@meta.data$seurat_clusters

pred.seurat.merg <- SingleR(test =  testdata, 
                            ref = hpca.se, 
                            labels = hpca.se$label.fine,  
                            method = "clusters",
                            clusters = clusters, 
                            assay.type.test = "logcounts", 
                            assay.type.ref = "logcounts")
pred.seurat.merg$labels <- sub(":.*", "", pred.seurat.merg$labels)

plotScoreHeatmap(pred.seurat.merg,clusters=pred.seurat.merg@rownames, fontsize.row = 9,show_colnames = T)

celltype = data.frame(ClusterID=rownames(pred.seurat.merg), 
                      celltype=pred.seurat.merg$labels,
                      , stringsAsFactors = F)
celltype
table(celltype$celltype)

seurat.merg@meta.data$celltype = "NA"
table(seurat.merg@meta.data$celltype)
for(i in 1:nrow(celltype)){
  seurat.merg@meta.data[which(seurat.merg@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]
}


umap_celltype<-DimPlot(seurat.merg,
                       group.by="celltype",
                       reduction="tsne",
                       label = "T", 
                       pt.size = 0.2,
                       label.size = 5)
umap_celltype<-DimPlot(subset(seurat.merg,clinical=="IBD"),
                       group.by="celltype",
                       reduction="tsne",
                       label = "T", 
                       pt.size = 0.2,
                       label.size = 5)
umap_celltype<-DimPlot(subset(seurat.merg,clinical=="HC"),
                       group.by="celltype",
                       reduction="tsne",
                       label = "T", 
                       pt.size = 0.2,
                       label.size = 5)

Idents(seurat.merg) <- "seurat_clusters"

umap_celltype_num <- DimPlot(
  seurat.merg,
  reduction = "tsne",
  label = TRUE,         
  pt.size = 0.2,
  label.size = 5
)
genes <- c('CCL20', 'IL18', 'IL12B', 'TNFRSF9', 'CCL7', 'IL1RL1', 'IL19', 'CCL13', 'CD72', 'EPHB4',
           'SEPT8' ,'IFNG' ,'CD6' ,'LY9', 'IL10RA', 'IGLC2', 'MXRA8', 'ITGAV', 'CXCL9', 'NOV',
           'RSPO3' ,'NOS3')

list1 <- list()
list2 <- list()
list3 <- list()

for (gene in genes) {
  print(gene)
  gene_data <- GetAssayData(seurat.merg, assay = "RNA", slot = "data")[gene, ]
  seurat.merg$gene_expression <- gene_data
  filtered_data_IBD <- subset(seurat.merg, clinical == "IBD" & gene_expression > 0)
  p1 <- FeaturePlot(
    object = filtered_data_IBD, 
    features = gene,
    reduction = "tsne",
    cols = c("#FF7078", "#CC1C2F"),
    pt.size = 0.1
  )+ xlim(-50, 50) + ylim(-50, 50) 
  list1[[gene]] <- p1  

  if (gene %in% c("IL12B", "IL19")) {
    p2 <- p1  
  } else {
    filtered_data_HC <- subset(seurat.merg, clinical == "HC" & gene_expression > 0)
    p2 <- FeaturePlot(
      object = filtered_data_HC, 
      features = gene,
      reduction = "tsne",
      cols = c("#79ABE2", "#0072B4"),
      pt.size = 0.1
    )+ xlim(-50, 50) + ylim(-50, 50) 
  }
  list2[[gene]] <- p2  
  filtered_data_zero <- subset(seurat.merg, gene_expression == 0)
  p3 <- FeaturePlot(
    object = filtered_data_zero, 
    features = gene,
    reduction = "tsne",
    cols = c("#ADB6B6FF", "#ADB6B6FF"),
    pt.size = 0.1
  )+ xlim(-50, 50) + ylim(-50, 50) 
  list3[[gene]] <- p3  
}
wrap_plots(list1, ncol = 4)
wrap_plots(list2, ncol = 4)
wrap_plots(list3, ncol = 4)

seurat.merg$celltype <- factor(seurat.merg$celltype, levels = rev(sort(unique(seurat.merg$celltype))))
head(seurat.merg)
Idents(seurat.merg) <-'celltype'
p3 <- DotPlot(seurat.merg, split.by = "clinical",cols = c("#0072B4", "#CC1C2F"), 
              scale.by = "size",
              features = genes) + 
  RotatedAxis()


p4 <- DotPlot(subset(seurat.merg,clinical=="IBD"), 
              cols = c("grey","#CC1C2F"),
              scale.by = "size",
              features = genes) +
  RotatedAxis()


p5 <- DotPlot(subset(seurat.merg,clinical=="HC"), 
              cols = c("grey","#0072B4"),
              scale.by = "size",
              features = genes) +
  RotatedAxis()

seurat.merg_subset <- seurat.merg
expression_data <- GetAssayData(seurat.merg, assay = "RNA", slot = "data")
filtered_data <- expression_data[genes, ]
seurat.merg_subset <-  SetAssayData(seurat.merg, assay = "RNA", slot = "data", new.data = filtered_data)
table(seurat.merg_subset$celltype)
table(seurat.merg_subset$clinical)
head(seurat.merg_subset)
seurat.merg_subset.markers <- FindAllMarkers( seurat.merg_subset, only.pos = TRUE, min.pct = 0.1 , logfc.threshold = 0.25)
seurat.merg_subset.markers <- seurat.merg_subset.markers[seurat.merg_subset.markers$p_val_adj<0.05,]
seurat.merg_subset.markers$cluster <- factor(seurat.merg_subset.markers$cluster, levels = rev(sort(unique(seurat.merg_subset.markers$cluster))))
genes_order <-  c('CCL20', 'IL18', 'TNFRSF9', 'CCL7','CCL13', 'CD72', 'EPHB4',
                    'SEPT8' ,'IFNG' ,'CD6' ,'LY9', 'IL10RA', 'IGLC2', 'MXRA8', 'ITGAV', 'CXCL9', 'NOV',
                    'RSPO3' ,'NOS3')
seurat.merg_subset.markers$gene <- factor(seurat.merg_subset.markers$gene, levels = genes_order)

p6 <- ggplot(seurat.merg_subset.markers, aes(interaction(group,gene ), avg_log2FC,  fill = cluster)) + 
  geom_bar(stat = "identity", position = "stack", na.rm = FALSE) + 
  scale_x_discrete(drop = FALSE) +  
  ylab("Average log2FC") +
  xlab("Gene") +
  theme_classic() +
  theme(
    legend.title = element_blank(),
    legend.key.width = unit(1.2, "cm"),
    legend.text = element_text(size = 16, face = "bold"),
    
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 16, face = "bold"),
    
    strip.text = element_text(size = 16, face = "bold"),
    strip.background = element_blank(),
    
    panel.grid.major.y = element_line(size = 0.5, linetype = "dashed", colour = "grey"),
    panel.spacing = unit(1.2, "lines")
  )
p6

library(dplyr)

mymatrix2 <- GetAssayData(seurat.merg, assay = "RNA", slot = "data")
mymatrix2[1:20,1:2]
mymatrix2 <- mymatrix2[rownames(mymatrix2) %in% genes,]

library(SeuratData)
library(Seurat)
library(ggSCvis)
head(seurat.merg)
mymatrix2 <- t(mymatrix2)%>% as.data.frame()
mymatrix2$celltype <- seurat.merg$celltype  

mymatrix2$sample <- seurat.merg$sample  

mymatrix2$clinical <- seurat.merg$clinical 
head(mymatrix2)

results_list <- list()

for (i in unique(mymatrix2$celltype)) {
  print(paste("Processing cell type:", i))

  sub_data <- subset(mymatrix2, celltype == i)

  hc_group <- sub_data[sub_data$clinical == "HC", genes, drop = FALSE]
  ibd_group <- sub_data[sub_data$clinical == "IBD", genes, drop = FALSE]

  logFC_results <- c()
  pval_results <- c()
  i_type <- c()
  hc_means <- c()  
  ibd_means <- c() 

  for (gene in genes) {
    if (gene %in% colnames(hc_group) && gene %in% colnames(ibd_group)) {
      hc_values <- hc_group[[gene]]
      ibd_values <- ibd_group[[gene]]

      test_result <- wilcox.test(hc_values, ibd_values, exact = FALSE) 

      logFC <- log2(mean(ibd_values + 0.001) / mean(hc_values + 0.001))

      hc_mean <- mean(hc_values, na.rm = TRUE)
      ibd_mean <- mean(ibd_values, na.rm = TRUE)

      logFC_results <- c(logFC_results, logFC)
      pval_results <- c(pval_results, test_result$p.value)
      i_type <- c(i_type, i)
      hc_means <- c(hc_means, hc_mean)  
      ibd_means <- c(ibd_means, ibd_mean) 
    }
  }

  result_df <- data.frame(
    Gene = genes,
    type = i_type,
    logFC = logFC_results,
    p.value = pval_results,
    HC_Mean = hc_means,
    IBD_Mean = ibd_means
  )

  results_list[[i]] <- result_df
}

final_results <- do.call(rbind, results_list)

head(final_results)
cell_gene_pairs <- data.frame(
  CellType = c("Tissue_stem_cells", "Smooth_muscle_cells", "Macrophage", "T_cell",
               "Smooth_muscle_cells", "T_cell", "NK_cell", "Macrophage", "DC", 
               "Macrophage", "Endothelial_cells", "NK_cell", "B_cell", "T_cell", 
               "Macrophage", "DC", "Epithelial_cells", "Macrophage", 
               "Endothelial_cells", "Smooth_muscle_cells", "Tissue_stem_cells",
               "Macrophage", "B_cell", "DC", "Tissue_stem_cells", "Smooth_muscle_cells", 
               "Endothelial_cells", "Endothelial_cells", "Smooth_muscle_cells",
               "Tissue_stem_cells", "Smooth_muscle_cells", "Endothelial_cells", "NK_cell"),
  Gene = c("CCL13", "CCL13", "CCL20", "CCL20", "CCL7", "CD6", "CD6", "CD72", "CD72","CXCL9",
           "EPHB4", "IFNG", "IGLC2", "IL10RA", "IL10RA", "IL10RA", "IL18", 
           "IL18", "ITGAV", "ITGAV", "ITGAV", "ITGAV", "LY9", "LY9", "MXRA8", 
           "MXRA8", "NOS3", "NOV", "RSPO3", "SEPT8", "SEPT8", "SEPT8", "TNFRSF9")
)
filtered_results <- final_results %>%
  inner_join(cell_gene_pairs, by = c("Gene" = "Gene", "type" = "CellType"))

filtered_results$adjusted_pvalue <- p.adjust(filtered_results$p.value, method = "fdr")
library(fdrtool)
fdr <- fdrtool(filtered_results$p.value, statistic="pvalue")
filtered_results$qval <- fdr$qval

dim(seurat.merg)

table(seurat.merg$orig.ident)
Idents(seurat.merg) <-'celltype'
table(Idents(seurat.merg))
head(seurat.merg)
VlnPlot(seurat.merg,features = 'CCL13',split.by = 'clinical',idents = 'Tissue_stem_cells')
VlnPlot(seurat.merg,features = 'CCL13',split.by = 'clinical',idents = 'Smooth_muscle_cells')

VlnPlot(seurat.merg,features = 'CD6',split.by = 'clinical',idents = 'T_cell')

cell_gene_pairs <- data.frame(
  CellType = c("Tissue_stem_cells", "Smooth_muscle_cells", "Macrophage", "T_cell",
               "Smooth_muscle_cells", "T_cell", "NK_cell", "Macrophage", "DC", 
               "Macrophage", "Endothelial_cells", "NK_cell", "B_cell", "T_cell", 
               "Macrophage", "DC", "Epithelial_cells", "Macrophage", 
               "Endothelial_cells", "Smooth_muscle_cells", "Tissue_stem_cells",
               "Macrophage", "B_cell", "DC", "Tissue_stem_cells", "Smooth_muscle_cells", 
               "Endothelial_cells", "Endothelial_cells", "Smooth_muscle_cells",
               "Tissue_stem_cells", "Smooth_muscle_cells", "Endothelial_cells", "NK_cell"),
  Gene = c("CCL13", "CCL13",
           "CCL20", "CCL20", 
           "CCL7", 
           "CD6", "CD6", 
           "CD72", "CD72",
           "CXCL9",
           "EPHB4", "IFNG", "IGLC2", "IL10RA", "IL10RA", "IL10RA", "IL18",
           "IL18", "ITGAV", "ITGAV", "ITGAV", "ITGAV", "LY9", "LY9", "MXRA8", 
           "MXRA8", "NOS3", "NOV", "RSPO3", "SEPT8", "SEPT8", "SEPT8", "TNFRSF9")
)
cell_gene_pairs$Gene <- factor(cell_gene_pairs$Gene,levels = c('CCL20', 'IL18',  'TNFRSF9', 'CCL7', 'CCL13', 'CD72', 'EPHB4',
                                                               'SEPT8' ,'IFNG' ,'CD6' ,'LY9', 'IL10RA', 'IGLC2', 'MXRA8', 'ITGAV', 'CXCL9', 'NOV',
                                                               'RSPO3' ,'NOS3'))

gene <- cell_gene_pairs$Gene[1]

plots <- lapply(levels(cell_gene_pairs$Gene), function(gene) {

  current_cells <- cell_gene_pairs[cell_gene_pairs$Gene == gene, "CellType"]

  plot_data <- FetchData(seurat.merg, vars = c(gene, 'clinical', 'celltype')) %>%
    filter(celltype %in% current_cells) %>%
    group_by(clinical, celltype) %>%
    dplyr::summarize(mean_value = mean(get(gene), na.rm = TRUE), 
                     sd_value = sd(get(gene), na.rm = TRUE), 
                     .groups = 'drop')  

  print(paste("Processing gene:", gene))
  print(head(plot_data)) 

  p <- VlnPlot(seurat.merg, features = gene, split.by = 'clinical',
               cols = c("#0072B4", "#CC1C2F"), alpha = 1,idents = current_cells) +
    ggtitle(gene) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          axis.title.x = element_blank(),
          legend.position = 'left') 

  return(p)
})
