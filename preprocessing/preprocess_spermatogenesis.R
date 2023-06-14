##load mouse data
rm(list = ls())
library(Seurat)

cnts_mouse <- Matrix::readMM("../datasets/raw/Spermatogenesis/mouse/raw_counts.mtx")

gene_mouse <- read.table("../datasets/raw/Spermatogenesis/mouse/genes.tsv", sep='\t', header=TRUE)
idx <- as.vector(tapply(seq_along(gene_mouse$Symbol), gene_mouse$Symbol, function(x){x[1]})[unique(gene_mouse$Symbol)])
gene_mouse <- gene_mouse[idx, ]
cnts_mouse <- cnts_mouse[idx, ]
rownames(cnts_mouse) <- as.vector(gene_mouse$Symbol)

meta_mouse <- read.table("../datasets/raw/Spermatogenesis/mouse/cell_metadata.txt")
idx <- as.vector(tapply(seq_along(meta_mouse$Barcode), meta_mouse$Barcode, function(x){x[1]})[unique(meta_mouse$Barcode)])
meta_mouse <- meta_mouse[idx, ]
cnts_mouse <- cnts_mouse[, idx]
rownames(meta_mouse) <- as.vector(meta_mouse$Barcode)
colnames(cnts_mouse) <- as.vector(meta_mouse$Barcode)
meta_mouse <- meta_mouse[colnames(cnts_mouse), ]

obj_mouse <- CreateSeuratObject(cnts_mouse)
obj_mouse@meta.data$AnnotatedClusters <- as.vector(meta_mouse$AnnotatedClusters)
obj_mouse@meta.data$BroadClusters <- as.vector(meta_mouse$BroadClusters)
obj_mouse@meta.data$library_id <- as.vector(meta_mouse$Sample)
obj_mouse@meta.data$sub_cell_type <- as.vector(meta_mouse$AnnotatedClusters)

save(obj_mouse, file="../datasets/raw/Spermatogenesis/obj_mouse.RData")



##mouse annotation
rm(list = ls())
library(Seurat)
set.seed(1234)
load("../datasets/raw/Spermatogenesis/obj_mouse.RData")

cnts <- obj_mouse[["RNA"]]@counts

germs <- c("eP1","eP2","mP","lP1","lP2","D","MI","MII","S1","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11","Spermatogonia")

idx <- which((as.vector(obj_mouse$library_id) %in% c("P30")) & 
             (as.vector(obj_mouse$AnnotatedClusters) %in% germs))

cnts <- cnts[, idx]

obj_mouse_v2 <- CreateSeuratObject(cnts)
obj_mouse_v2 <- NormalizeData(obj_mouse_v2)
obj_mouse_v2 <- FindVariableFeatures(obj_mouse_v2, nfeatures = 2000, selection.method = "vst", verbose = FALSE)
obj_mouse_v2 <- ScaleData(obj_mouse_v2, features = rownames(obj_mouse_v2), verbose = FALSE)
hvg_mouse_v2 <- obj_mouse_v2@assays$RNA@var.features
obj_mouse_v2 <- RunPCA(obj_mouse_v2, features = hvg_mouse_v2, npcs = 50, verbose = FALSE)
obj_mouse_v2 <- RunUMAP(obj_mouse_v2, reduction = "pca", dims = 1:50, umap.method = "umap-learn", metric = "correlation", verbose = FALSE)
obj_mouse_v2 <- FindNeighbors(obj_mouse_v2, dims = 1:50)
obj_mouse_v2 <- FindClusters(obj_mouse_v2, resolution = 2.5)
DimPlot(obj_mouse_v2, reduction = "umap", label=T)
DotPlot(obj_mouse_v2, features = c("Sycp1","Uchl1","Crabp1","Stra8","Sohlh2","Dazl","Scml2","Rpa2","Rad51",#SPG
                                   "Hormad1","Piwil1","Pttg1","Insl6","Spag6","Tbpl1", #Scytes
                                   "Tssk1","Acrv1","Spaca1","Tsga8", #STids
                                   "Prm1","Prm2","Tnp1","Tnp2")) + RotatedAxis() 

barcode <- colnames(obj_mouse_v2)
cell_type_new <- as.vector(Idents(obj_mouse_v2))
# cell_type_new[which(cell_type_new %in% c("4","7"))] <- "Spermatogonia"
# cell_type_new[which(cell_type_new %in% c("18","13","20","19","15","14","3","5","10"))] <- "Spermatocytes"
# cell_type_new[which(cell_type_new %in% c("12","17","11","0","1","8","2","16"))] <- "Round spermatids"
# cell_type_new[which(cell_type_new %in% c("6","9"))] <- "Elongating spermatids"
cell_type_new[which(cell_type_new %in% c("4","7"))] <- "Spermatogonia"
cell_type_new[which(cell_type_new %in% c("18","13","20","19","15","14","3","5","10"))] <- "Spermatocyte"
cell_type_new[which(cell_type_new %in% c("12","17","11","0","1","8","2","16"))] <- "RoundSpermatid"
cell_type_new[which(cell_type_new %in% c("6","9"))] <- "Elongating"

newanno <- data.frame(barcode = barcode, cell_type_new = cell_type_new, louvain = as.vector(Idents(obj_mouse_v2)))
write.csv(newanno, file = "../datasets/raw/Spermatogenesis/mouseP30-anno.csv")

load("../datasets/raw/Spermatogenesis/obj_mouse.RData")
idx <- which((as.vector(obj_mouse$library_id) %in% c("P30")) & 
             (as.vector(obj_mouse$BroadClusters) == "Germ") &
             (as.vector(obj_mouse$AnnotatedClusters) != "Outliers"))
cnts <- obj_mouse[["RNA"]]@counts
cnts <- cnts[, idx]
rownames(newanno) <- as.vector(newanno$barcode)
newanno <- newanno[colnames(cnts), ]
celltype <- as.vector(newanno$cell_type_new)
rm(obj_mouse)
obj_mouse <- CreateSeuratObject(cnts)
obj_mouse@meta.data$celltype <- celltype
save(obj_mouse, file="../datasets/raw/Spermatogenesis/obj_mouse_v2.RData")



##load human data
rm(list = ls())
library(Seurat)

cnts_human <- read.table("../datasets/raw/Spermatogenesis/human/GSE142585_MergedHumanTestis4_DGE.txt")
anno_human <- read.table("../datasets/raw/Spermatogenesis/human/GSE142585_MergedHumanTestis4_PerCellAttributes.txt")
anno_human <- anno_human[colnames(cnts_human), ]

idx <- ((!(as.vector(anno_human$CellType) %in% c("Macrophage","m-Pericyte","Myoid","f-Pericyte","Endothelial","ImmLeydig","Tcell"))) &
        (as.vector(anno_human$orig.ident) %in% c("Human1.1","Human1.2","Human1.3","Human1.4","Human1.5")))

obj_human <- CreateSeuratObject(cnts_human[, idx])
save(obj_human, file="../datasets/raw/Spermatogenesis/obj_human.RData")



##load macaque data
rm(list = ls())
library(Seurat)

cnts_macaque <- read.table("../datasets/raw/Spermatogenesis/macaque/GSE142585_MergedMonkeyTestis5_DGE.txt")
anno_macaque <- read.table("../datasets/raw/Spermatogenesis/macaque/GSE142585_MergedMonkeyTestis5_PerCellAttributes.txt")
anno_macaque <- anno_macaque[colnames(cnts_macaque), ]

idx <- ((!(as.vector(anno_macaque$CellType) %in% c("Macrophage","m-Pericyte","Myoid","f-Pericyte","Endothelial","ImmLeydig","Tcell"))) &
        (as.vector(anno_macaque$orig.ident) == "Monkey2"))

obj_macaque <- CreateSeuratObject(cnts_macaque[, idx])
save(obj_macaque, file="../datasets/raw/Spermatogenesis/obj_macaque.RData")



##find shared genes
rm(list = ls())
library(Seurat)
set.seed(1234)
load("../datasets/raw/Spermatogenesis/obj_human.RData")
load("../datasets/raw/Spermatogenesis/obj_macaque.RData")
load("../datasets/raw/Spermatogenesis/obj_mouse_v2.RData")

#1-1-1
df_macaque <- read.table("orthologues_human_macaque.txt", fill=TRUE, header=TRUE, sep="\t")
df_mouse <- read.table("orthologues_human_mouse.txt", fill=TRUE, header=TRUE, sep="\t")
df_macaque <- df_macaque[which(df_macaque$Macaque.homology.type == "ortholog_one2one"), ]
df_mouse <- df_mouse[which(df_mouse$Mouse.homology.type == "ortholog_one2one"), ]

cnt_mouse <- obj_mouse[["RNA"]]@counts
cnt_macaque <- obj_macaque[["RNA"]]@counts

macaque_set <- intersect(df_macaque$Macaque.gene.name, rownames(cnt_macaque))
mouse_set <- intersect(df_mouse$Mouse.gene.name, rownames(cnt_mouse))

df_macaque <- df_macaque[which(df_macaque$Macaque.gene.name %in% macaque_set), ]
df_mouse <- df_mouse[which(df_mouse$Mouse.gene.name %in% mouse_set), ]

cnt_macaque <- cnt_macaque[which(rownames(cnt_macaque) %in% macaque_set), ]
cnt_mouse <- cnt_mouse[which(rownames(cnt_mouse) %in% mouse_set), ]

#transfer gene names
genes_dict_mouse2human <- as.vector(df_mouse$Gene.name)
names(genes_dict_mouse2human) <- as.vector(df_mouse$Mouse.gene.name)

genes_dict_macaque2human <- as.vector(df_macaque$Gene.name)
names(genes_dict_macaque2human) <- as.vector(df_macaque$Macaque.gene.name)

rownames(cnt_mouse) <- as.vector(genes_dict_mouse2human[rownames(cnt_mouse)])
rownames(cnt_macaque) <- as.vector(genes_dict_macaque2human[rownames(cnt_macaque)])
cnt_human <- obj_human[["RNA"]]@counts

#shared genes
shared_genes <- intersect(rownames(cnt_human), rownames(cnt_mouse))
shared_genes <- intersect(shared_genes, rownames(cnt_macaque))

#reform obj_mouse
cell_type <- as.vector(obj_mouse$celltype)
rm(obj_mouse)
obj_mouse <- CreateSeuratObject(cnt_mouse[shared_genes, ])
obj_mouse@meta.data$cell_type <- cell_type

#reform obj_human
obj_human <- CreateSeuratObject(cnt_human[shared_genes, ])

#reform obj_macaque
obj_macaque <- CreateSeuratObject(cnt_macaque[shared_genes, ])

save(obj_human, file="../datasets/raw/Spermatogenesis/obj_human_sharedgenes.RData")
save(obj_macaque, file="../datasets/raw/Spermatogenesis/obj_macaque_sharedgenes.RData")
save(obj_mouse, file="../datasets/raw/Spermatogenesis/obj_mouse_sharedgenes.RData")


