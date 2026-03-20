# Title: Perform functional enrichment using GSEA based on GO ontology
# Author: Jingxun Chen
# Date: code compiled on 20241209
# Related publication: Emma K. Costa, and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Figure 1F, 1G, Extended Data Fig. 2, Extended Data Fig. 3

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory
setwd("/Users/jingxun/Dropbox/BrunetLab/Data/AtlasData/GSEA/")

# Load packages
library("DOSE")
library("clusterProfiler") 
library("org.Hs.eg.db")

# To convert NCBI ids to human entrez ids. This is needed to run the package. There are ways to adapt it for nfur only, but for now I do everything based on human orthologs
hSymbols = read.table("Input/NCBI-Human-orthologs.txt", head = T, sep = "\t")

# Select one tissue at a time to run
# t = 'Bone'
# t = 'Brain'
# t = 'Fat'
# t = 'Gonad'
# t = 'Gut'
# t = 'Heart'
# t = 'Kidney'
# t = 'Liver'
# t = 'Muscle'
# t = 'Skin'
# t = 'SpinalCord'
# t = 'Spleen'
t = 'Eye'

# Select one sex at a time to run
# sex = 'female'
sex = 'male'

# ------------------------------------------------------------------
# Select input
# ------------------------------------------------------------------
# Select a proper input for GSEA by removing the '#' sign. 
# Add a '#' sign for input not analyzed currently. 
# Run each input one at a time.

data <- read.csv(paste0("Input/List_by_tissue/241208_allCorrresults", t ,"_", sex,"Only.csv"))

# ------------------------------------------------------------------
# Generate GSEA input, which is a ranked list of genes.
# ------------------------------------------------------------------
# calculate ranks based on (-log10(p-value of Spearman's correlation) x Spearman's correlation)

data$mlog10QvalxFC <- (-log10(data$p_spear))*(data$cor_spear)
data <- subset(data, select = c(gene, mlog10QvalxFC))
colnames(data) <- c('Gene', 'mlog10QvalxFC')
head(data)


# ------------------------------------------------------------------------------------------------
# Get human ortholog symbols based on the BLAST results file using org.Hs.eg.db package
# ------------------------------------------------------------------------------------------------

# Some Ids will fail to map and will be ignored
dataH = merge(hSymbols, data, by.x = "ncbi", by.y = "Gene") 
entrezIds = bitr(as.character(dataH[,2]), fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db") # Get entrez ids for annotation
### Note: "Warning message: 4.13% of input gene IDs are fail to map" ###
dataHE = merge(dataH, entrezIds, by.x = "human", by.y = "SYMBOL") # Get human symbols
head(dataHE)

# There can be duplicate values because of paralogs. I take the average of those for quantitative score
dataHE$mlog10QvalxFC <- as.numeric(dataHE$mlog10QvalxFC)
unique = aggregate(dataHE[,3], list(dataHE$human), mean)
dataHEU = merge(unique, entrezIds, by.x = "Group.1", by.y = "SYMBOL")
colnames(dataHEU) = c("human", "mlog10QvalxFC", "entrez")
head(dataHEU)

geneList = dataHEU[,2]  # gene list for GO 
names(geneList) = as.character(dataHEU[,1]) # with entrez ids as names
  
# *** Sort the gene list based on quantitative score in decreasing order. This is critical for GSEA  
geneList = sort(geneList, decreasing = TRUE)
  
head(geneList)
tail(geneList)
  


# ------------------------------------------------------------------------------------------------
# Do different enrichment analyses
# ------------------------------------------------------------------------------------------------

# --------------------- Gene Ontology -------------------
ego3 <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              keyType      = 'SYMBOL',
              ont          = c("ALL"),
              pvalueCutoff = 1)

head(ego3)

# Select proper outfile name based on input
write.table(ego3, paste0("Output/241209_allCorrresults", t ,"_", sex,"Only_GOGSEA.csv"), sep = ",", quote = T, row.names = F)


# ------------------------------------------------------------------
# Generate the killifish & human gene name conversion file  
# ------------------------------------------------------------------
# Select the proper output name based on input
write.table(dataHE, paste0("Output/241209_allCorrresults", t ,"_", sex,"Only_GOGSEA_HumanName.csv"), sep = ",", quote = T, row.names = F)


# ------------------------------------------------------------------
# Clear list to run the script again  
# ------------------------------------------------------------------

rm(list=ls()) 

sessionInfo() 