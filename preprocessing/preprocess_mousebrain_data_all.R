rm(list = ls())
library(Seurat)
library(ggplot2)
library(DropSeq.util)



##CB
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- rownames(dge)
cnt <- dge


##FC
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##PC
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##ENT
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60EntoPeduncular.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##GP
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60GlobusPallidus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##Hippo
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##STR
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Striatum.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##SN
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60SubstantiaNigra.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##TH
dge.path <- "../datasets/raw/MouseBrain/Drop-seq/GSE116470_F_GRCm38.81.P60Thalamus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])



sce <- cnt
save(sce, file="../datasets/raw/MouseBrain/Drop-seq_all_withNAanno.RData")



load("../datasets/raw/MouseBrain/Drop-seq_meta_all.RData")
sce <- sce[, as.vector(meta_all$cellid)]
print(dim(sce))
save(sce, file="../datasets/raw/MouseBrain/Drop-seq_all.RData")




