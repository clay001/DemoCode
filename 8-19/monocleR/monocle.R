devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(dplyr)
library(Seurat)

# Load the data
expression_matrix <- readRDS('cao_l2_expression.rds')
cell_metadata <- readRDS('cao_l2_colData.rds')
gene_annotation <- readRDS('cao_l2_rowData.rds')

# Make the CDS object
cds <- new_cell_data_set(expression_matrix,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)

# pre-process
cds <- preprocess_cds(cds, num_dim = 100)

# plot pc weight
plot_pc_variance_explained(cds)

# dimension reduction
cds <- reduce_dimension(cds)
plot_cells(cds)
plot_cells(cds, color_cells_by="cao_cell_type")
plot_cells(cds, genes=c("cpna-2", "egl-21", "ram-2", "inos-1"))

# tsne
cds <- reduce_dimension(cds, reduction_method="tSNE")
plot_cells(cds, reduction_method="tSNE", color_cells_by="cao_cell_type")
