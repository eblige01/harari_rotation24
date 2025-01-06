## The purpose of this scripp is to construct coexpression networks using data from brasse et al.
# Loading libraries
library(Seurat)
library(tidyverse)
library(cowplot)
library(patchwork)
library(devtools)
library(WGCNA)
library(hdWGCNA)
library(CEMiTool)

# using the cowplot theme for ggplot
theme_set(theme_cowplot())

# set random seed for reproducibility
set.seed(12345)

# optionally enable multithreading
enableWGCNAThreads(nThreads = 8)

seurat_obj <- readRDS("/users/PAS2694/eblige2/allCellTypes_garnettCleaned_finalObj_294114nuclei.rds") 

seurat_obj <- UpdateSeuratObject(object=seurat_obj)
# Has to be done BEFORE setupforwdwgcna

DefaultAssay(seurat_obj) <- "RNA"

seurat_obj <- subset(seurat_obj, subset = cellType == "Micro")

raw_counts <- GetAssayData(seurat_obj, slot = "counts", assay = "RNA")

# Save as CSV
write.csv(as.data.frame(raw_counts), "raw_counts.csv")


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
  reduction = 'pca',
    k =36 ,# select the dimensionality reduction to perform KNN on
  max_shared = 10, # maximum number of shared cells between two metacells
  ident.group = 'cellType'
  )

# normalize metacell expression matrix:
seurat_obj <- NormalizeMetacells(seurat_obj)


seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = c("Micro"),# the name of the group of interest in the group.by column
  group.by='cellType',# the metadata column containing the cell type info. This same column should have also been used in MetacellsByGroups
  assay = 'RNA', 
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
  tom_name = 'microstate3_Presympt' # name of the topoligical overlap matrix written to disk
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



# get hub genes
hub_df <- GetHubGenes(seurat_obj, n_hubs = 30)



# get hMEs from seurat object
MEs <- GetMEs(seurat_obj, harmonized=TRUE)
modules <- GetModules(seurat_obj)
mods <- levels(modules$module); mods <- mods[mods != 'grey']

# add hMEs to Seurat meta-data:
seurat_obj@meta.data <- cbind(seurat_obj@meta.data, MEs)

# Remanimg modules from color
seurat_obj <- ResetModuleNames(
  seurat_obj,
  new_name = "micro_state3_Presympt" # the base name for the new modules
)

saveRDS(seurat_obj, file="hdWGCNA_object_microState4.rds")



# make a featureplot of hMEs for each module
plot_list <- ModuleFeaturePlot(
  seurat_obj,
  features='hMEs', # plot the hMEs
  order=TRUE # order so the points with highest hMEs are on top
)

# stitch together with patchwork
wrap_plots(plot_list, ncol=6)


write.csv(hub_df , "hub_genes_micro3_Presympt.csv")

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


write.csv(correlations , "ME_cor_micro3_Presympt.csv")

# Fetch metadata for cell types (assuming you have a cell type column)
cell_types <- seurat_obj@meta.data$cellType

# Compare gene expression across cell types
gene_expression_by_cell_type <- aggregate(FetchData(seurat_obj, vars = "PSEN1"), by = list(cell_types), FUN = mean)

# View gene expression by cell type
print(gene_expression_by_cell_type)



gene_module <- seurat_obj@meta.data[rownames(seurat_obj@meta.data) == "PSEN1", "module_assignments"]

modules <- GetModules(seurat_obj) %>% subset(module != 'grey')

write.csv(modules , "modules_micro3_Presympt.csv")



