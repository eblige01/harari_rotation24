# The purpose of this script is to compare and contrast module genes to psen1 substrates found in Guner and Lichtenthaler 
library(biomaRt)
library(dplyr)
library(readr)
library(readxl)
library(venn)
library(ggplot2)
library(sf)
library(ggVennDiagram)


# Creating mouse and human protien list
humanlist <- c("acp2","calsyntenin-1","calsyntenin-2","calsyntenin-3","aplp1","aplp2","app","axl","bcma","betacellulin","betaglycan","car","cd147","cd200","cd43","cd44","cd99","cx3cl1","cxcl16","desmoglein-2","dner","dr6","dystroglycan","e-cadherin","epcam","epha2","epha5","epha7","ephb3","ephb4","ephb6","ephrin-b2","erbb4","f11r","fgfr4","glg1","hla-a2","ifnar2","igf-1r","il11r","il-1r1","il-1r2","il6r","ir","irela","ilrelb","jagged2","kcne2","kcne1","klotho","lar","ldlr","lrp1","lrp1b","lrp6","mer","met","muc1","musk","n-cadherin","neogenin","neurexin-3-b","NLRR3","p75NTR","pianp","plxdc2","pmel17","pdpn","polycystin-1","polyductin","protocadherin-12","ptk7","rage","robo1","sez6","sez6l","sez6l2","sirpa","sorcs1","sorla","sortilin","synedecan-1","synedecan-2","synedecan-3","tie1","tmeff2","tnfr1","trem2","trka","tyro3","vasorin","vegfr3")

mouselist <- c("apoER2","CACHD1","CADM1","CSF1R","Delta1","dscam","dscaml1","ephb2","ephrin-b2","fgfr3","ghr","l1","nectin1a","nectin3","ng2","notch1","notch2","notch3","notch4","nprc","nradd","pianp","prima","pvrl2","sez6","tkrB","trop2","tyrp1","tyrp2","vgscb1","vgscb2","vgscb3","vgscb4","vldlr")

mouse_mart <- useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
     

# Get UniProt IDs for these proteins
uniprot_mouse <- getBM(attributes = c("uniprotswissprot", "mgi_symbol"),
                      filters = "external_gene_name",
                      values = mouselist,
                      mart = mouse_mart)

mouse_ids <- uniprot_mouse$uniprotswissprot
# converting from mouse to human genes
uniprot_mouse <- getBM(attributes = c("uniprotswissprot", "mgi_symbol"),
                      filters = "uniprotswissprot",
                      values = mouse_ids,
                      mart = mouse_mart)

mouse_genes <- uniprot_mouse$mgi_symbol
     
mapping <- read.delim("HMD_HumanPhenotype.rpt", header = FALSE)
     
colnames(mapping) <- c("Human_Gene","hum_num", "Mouse_Gene", "Mouse_MGI_ID", "HPO_Term","last")
     
result <- mapping %>%
  filter(Mouse_Gene %in% mouse_genes) %>%
  select(Mouse_Gene, Human_Gene)
     
converted_mouse <- result$Human_Gene
     

ensembl <- useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")

     

# Fetches hgnc_symbol and unirpot ID
results_hum <- getBM(
  attributes = c("hgnc_symbol","uniprotswissprot"), # Attributes: UniProt ID and Gene Symbol
  filters = "hgnc_symbol",                      # Filter by UniProt IDs
  values = humanlist,                                 # Input your list of UniProt IDs
  mart = ensembl                                      # Use the Ensembl mart connection
)
     
# saving results to a list
human_genes <- results_hum$hgnc_symbol
     

human_genes <-c(human_genes,mouse_genes)
     
     

ad_genes <- read_excel("AD_genes.xlsx", col_names = FALSE)

ad_genes <- ad_genes$"...1"




# Initialize a list to store module lists from each file
module_lists <- c()
csv_filepaths <- c( "oligo_modules/oligoState1.csv", "oligo_modules/oligoState2.csv","oligo_modules/oligoState3.csv","oligo_modules/oligoState4.csv","oligo_modules/oligoState5.csv","oligo_modules/oligoState6.csv","oligo_modules/oligoState7 (1).csv","oligo_modules/oligoState8 (1).csv")
for (filepath in csv_filepaths) {
# Loop through each CSV file

  # Read the CSV file
  data <- read.csv(filepath)
  
  psen1_module <- data$module[data$gene_name == "PSEN1"]
  data <- subset(data, module == psen1_module)  
  # Add the module list to the overall list
  gene_sets <- list(
  "AD Genes" = ad_genes,
  "Module Genes" = data$gene_name,
  "PSEN1 substrates" = human_genes
)
# Create Venn Diagram with custom colors
venn <- ggVennDiagram(
  gene_sets,
  label = "both" # Show both counts and genes
) +
  # Adjust fill colors for better contrast
  scale_fill_gradient(low = "grey", high = "lightblue") +
  # Adjust outline colors for clarity
  scale_color_manual(values = c("black", "black", "black")) +
  # Tweak theme for visibility
  theme(
    legend.position = "none", # Remove legend
    panel.background = element_rect(fill = "white"), # White background
    plot.background = element_rect(fill = "white")  # White plot background
  )
     

psen_module_intersection <- intersect(data$gene_name,human_genes)
all_intersection <- intersect(psen_module_intersection,ad_genes) 
print(psen_module_intersection)
print( all_intersection)
}


