library(dplyr)
library(Seurat)
# Load the PBMC dataset
pbmc.data <- Read10X(data.dir = "./data/pbmc3k/filtered_gene_bc_matrices/hg19/")

# Initialize the Seurat object with the raw (non-normalized data).
pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
pbmc

pbmc[["percent.mt"]] <- PercentageFeatureSet(object = pbmc, pattern = "^MT-")
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

pbmc <- subset(x = pbmc, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & percent.mt < 5)
VlnPlot(object = pbmc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)




all.genes <- rownames(x = pbmc)
pbmc <- ScaleData(object = pbmc, features = all.genes)
pbmc <- RunPCA(object = pbmc, features = VariableFeatures(object = pbmc))
DimPlot(object = pbmc, reduction = "pca")





saveRDS(pbmc, file = "../output/pbmc_tutorial.rds")
cluster1.markers <- FindMarkers(object = pbmc, ident.1 = 1, min.pct = 0.25)
head(x = cluster1.markers, n = 5)



if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("monocle")

library(Biobase)
library(knitr)
library(reshape2)
library(ggplot2)
library(HSMMSingleCell)
library(monocle)
data(HSMM_expr_matrix) ## RPKM 矩阵,271个细胞，47192个基因
data(HSMM_gene_annotation)
data(HSMM_sample_sheet)
HSMM_expr_matrix[1:10,1:5] 
head(HSMM_gene_annotation)

pd <- new("AnnotatedDataFrame", data = HSMM_sample_sheet)
fd <- new("AnnotatedDataFrame", data = HSMM_gene_annotation)

# First create a CellDataSet from the relative expression levels

## 这里仅仅是针对rpkm表达矩阵的读取
HSMM <- newCellDataSet(as.matrix(HSMM_expr_matrix),   
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.1,
                       expressionFamily=tobit(Lower=0.1))

# Next, use it to estimate RNA counts
rpc_matrix <- relative2abs(HSMM)
rpc_matrix[1:10,1:5] 
## rpkm格式的表达值需要转换成reads counts之后才可以进行下游分析

# Now, make a new CellDataSet using the RNA counts
HSMM <- newCellDataSet(as(as.matrix(rpc_matrix), "sparseMatrix"),
                       phenoData = pd, 
                       featureData = fd,
                       lowerDetectionLimit=0.5,
                       expressionFamily=negbinomial.size())


HSMM_myo <- HSMM[,pData(HSMM)$CellType == "Myoblast"]   
HSMM_myo <- estimateDispersions(HSMM_myo)

if(F){
  diff_test_res <- differentialGeneTest(HSMM_myo[expressed_genes,],
                                        fullModelFormulaStr="~Media")
  ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
}
HSMM_myo <- setOrderingFilter(HSMM_myo, ordering_genes)
plot_ordering_genes(HSMM_myo)
