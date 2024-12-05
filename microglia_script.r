library(Seurat)

# plotting and data science packages
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)


# co-expression network analysis packages:
library(WGCNA)
library(hdWGCNA)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("/users/PAS2694/eblige2/allCellTypes_garnettCleaned_finalObj_294114nuclei.rds") 

seurat_obj <- UpdateSeuratObject(object=seurat_obj)

p <- DimPlot(seurat_obj,group.by='seurat_clusters',  label=TRUE) +
   umap_theme() + ggtitle('all_cells') + NoLegend()

p

seurat_obj@meta.data$cellType

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction", # the gene selection approach
  fraction = 0.05, # fraction of cells that a gene needs to be expressed in order to be included 5%
  wgcna_name = "brasse_etal" # the name of the hdWGCNA experiment
)

# construct metacells  in each group
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("Sample_ID",'cellType'), # specify the columns in seurat_obj@meta.data to group by
  reduction = 'pca', # select the dimensionality reduction to perform KNN on
  max_shared = 10, 
   k = 50 ,# maximum number of shared cells between two metacells
  ident.group = 'cellType'
  )

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)

seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Micro"),# the name of the group of interest in the group.by column
  group.by='cellType',# the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'SCT', 
  slot = 'data' # using normalized data
)

# Test different soft powers:
seurat_obj <- TestSoftPowers(
  seurat_obj,
  networkType = 'signed' # you can also use "unsigned" or "signed hybrid"
)

# plot the results:
plot_list <- PlotSoftPowers(seurat_obj)

# assemble with patchwork
wrap_plots(plot_list, ncol=2)

seurat_obj <- ConstructNetwork(
  seurat_obj,
  tom_name = 'microglia' # name of the topoligical overlap matrix written to disk
)

PlotDendrogram(seurat_obj, main='Microglia hdWGCNA Dendrogram')

# need to run ScaleData first or else harmony throws an error:
seurat_obj <- ScaleData(seurat_obj, features=VariableFeatures(seurat_obj))

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleEigengenes(
 seurat_obj,
 group.by.vars="Sample_ID"
)

# compute eigengene-based connectivity (kME):
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = 'cellType', group_name = 'Micro'
)


# get the module assignment table:
modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

# show the first 6 columns:
head(modules[,1:6])

# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 30)

head(hub_df)

saveRDS(seurat_obj, file='hdWGCNA_object.rds')

seurat_obj <- readRDS("hdWGCNA_object.rds")

# Remanimg modules from color
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "micro" # the base name for the new modules
)

# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)

hub_df 

write.csv(hub_df , "hub_genes.csv")

# Replace "cell_type_column" with the actual column name that has the cell type labels
seurat_obj <- SetIdent(seurat_obj, value = "cellType")

# Visualize the gene expression across clusters/cells
FeaturePlot(seurat_obj, features = "PSEN1", label =TRUE)




# Get the gene expression
gene_expression <- FetchData(seurat_obj, vars = "PSEN1")

# Get module eigengenes (replace with the correct location in your object)
module_eigengenes <- seurat_obj@reductions$pca@cell.embeddings

# Correlate the gene expression with each module eigengene
correlations <- cor(gene_expression, module_eigengenes)

# View correlations
print(correlations)


write.csv(correlations , "ME_cor.csv")

# Fetch metadata for cell types (assuming you have a cell type column)
cell_types <- seurat_obj@meta.data$cellType

# Compare gene expression across cell types
gene_expression_by_cell_type <- aggregate(FetchData(seurat_obj, vars = "PSEN1"), by = list(cell_types), FUN = mean)

# View gene expression by cell type
print(gene_expression_by_cell_type)



gene_module <- seurat_obj@meta.data[rownames(seurat_obj@meta.data) == "PSEN1", "module_assignments"]

# View the module assignment for the specific gene
print(gene_module)


# Violin plot for a specific gene across clusters
VlnPlot(seurat_obj, features = "PSEN1", group.by = "cellType")







