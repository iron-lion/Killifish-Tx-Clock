# Title: Plot selected GO term genes as heat maps
# Author: Jingxun Chen
# Date: code compiled on 20241215
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025



# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/GSEA/")
set.seed(1234)

# load package
library(pheatmap)
library(dplyr)


# ---------- Select tissue to run ----------
# Select a proper input for GSEA by removing the '#' sign. 
# Add a '#' sign for input not analyzed currently. 
# Run one tissue at a time.

# t = 'Brain'
# t = 'Liver'
# t = 'Muscle'
t = 'Gut'

# Select a GO term to plot. Pick one termID at a time (keep the other ones commented out):
# --------- Brain section ---------
# termID = "GO:0045087"
# termName = 'InnateImmuneResponse'
# 
# termID = "GO:0051306"
# termName = 'Mitotic sister chromatid separation'


# --------- Liver section ---------
# termID = "GO:0051306"
# termName = 'Mitotic sister chromatid separation'
# 
# termID = 'GO:0002253'
# termName = 'Activation of immune response'
#

# --------- Muscle section ---------
# termID = "GO:0016236"
# termName = 'Macroautophagy'
#   
# termID = 'GO:0001503'
# termName = 'Ossification'
# 

# --------- Gut section ---------
# termID = "GO:0006111"
# termName = 'Regulation of gluconeogenesis'
#    
termID = "GO:0005759"
termName = "Mitochondrial matrix"
#  

# ---------------------------------------------------------------------------
# Load data 
# ---------------------------------------------------------------------------
# Load Spearman's rank correlation for male and female separately
data_m <- read.csv(paste0("Input/List_by_tissue/241208_allCorrresults", t ,"_maleOnly.csv"))
data_f <- read.csv(paste0("Input/List_by_tissue/241208_allCorrresults", t ,"_femaleOnly.csv"))

# Pick male or female GSEA results to plot 
# Male results:
dataHE_m <- read.csv(paste0("Output/241209_allCorrresults", t ,"_maleOnly_GOGSEA_HumanName.csv"), header = T)
GOterm_m <- read.csv(paste0("Output/241209_allCorrresults", t ,"_maleOnly_GOGSEA.csv"), header = T)

# Female results:
dataHE_f <- read.csv(paste0("Output/241209_allCorrresults", t ,"_femaleOnly_GOGSEA_HumanName.csv"), header = T)
GOterm_f <- read.csv(paste0("Output/241209_allCorrresults", t ,"_femaleOnly_GOGSEA.csv"), header = T)


# ------------------------------------------------------------------
# Clean up data
# ------------------------------------------------------------------
# Keep only relevant columns
data_m_1 <- subset(data_m, select = c(gene, cor_spear, padj_spear))
data_f_1 <- subset(data_f, select = c(gene, cor_spear, padj_spear))

# Rename columns to keep track of male and female terms
colnames(data_m_1) <- c('gene', 'cor_spear_m', 'padj_spear_m')
colnames(data_f_1) <- c('gene', 'cor_spear_f', 'padj_spear_f')

# Rename the killifish gene column to "gene", matching the correlation file
colnames(dataHE_m) <- c('human_m', "gene", 'mlog10QvalxFC_m', 'ENTREZID_m')
colnames(dataHE_f) <- c('human_f', "gene", 'mlog10QvalxFC_f', 'ENTREZID_f')

# Join the two input files by the killifish gene name
all_m <- left_join(dataHE_m, data_m_1, by= "gene")
all_f <- left_join(dataHE_f, data_f_1, by= "gene")

# Check for gene duplicates
dup_m <- sum(duplicated(all_m$gene)) # 4 genes
dup_f <- sum(duplicated(all_f$gene)) # 4 genes

dup.m.df <- all_m[duplicated(all_m$gene), ] # Identify duplicate rows, same genes as in females
dup.f.df <- all_f[duplicated(all_f$gene), ] # Identify duplicate rows, same genes as in males

all_m_unique <- all_m[!duplicated(all_m$gene), ] # Drop rows with duplicated values in column 'gene' (keep first occurrence)
all_f_unique <- all_f[!duplicated(all_f$gene), ] # Drop rows with duplicated values in column 'gene' (keep first occurrence)

# Join male and female files
all <- inner_join(all_f_unique, all_m_unique, by = 'gene')

# Remove unnecessary columns
all_1 <- subset(all, select = -c(human_m, ENTREZID_m, ENTREZID_f, mlog10QvalxFC_m, mlog10QvalxFC_f))
colnames(all_1) [1] <- 'human'

# ------------------------------------------------------------------
# Find the genes for each GO term selected for plotting
# ------------------------------------------------------------------
# Find the enriched genes associated with the GO term for both males and females

# Female
GOterm_f_go <- GOterm_f[GOterm_f$ID == termID, ] # subset to include the row of the termID
GOterm_f_go_1 <- subset(GOterm_f_go, select = c(ID, core_enrichment)) # remove unnecessary terms
GeneGO_f <- strsplit(GOterm_f_go_1$core_enrichment, split = "/") # split the enriched genes associated with each GO terms
GeneGO_f_df <- as.data.frame(GeneGO_f[[1]]) # make this list as a dataframe
colnames(GeneGO_f_df) <- 'human'

# Male
GOterm_m_go <- GOterm_m[GOterm_m$ID == termID, ] # subset to include the row of the termID
GOterm_m_go_1 <- subset(GOterm_m_go, select = c(ID, core_enrichment)) # remove unnecessary terms
GeneGO_m <- strsplit(GOterm_m_go_1$core_enrichment, split = "/") # split the enriched genes associated with each GO terms
GeneGO_m_df <- as.data.frame(GeneGO_m[[1]]) # make this list as a dataframe
colnames(GeneGO_m_df) <- 'human'

# Combine the male and female data frames into one
GeneGO_fm_df <- rbind(GeneGO_m_df, GeneGO_f_df)  
GeneGO_fm_df$tracker <- 1 # add a counter for plotting later

# Drop rows with duplicated values in column 'human' (keep first occurrence). There are duplicated terms because some genes are common between males and females.
GeneGO_fm_df_u <- GeneGO_fm_df[!duplicated(GeneGO_fm_df$human), ] 


# ------------------------------------------------------------------
# Make heatmap for Spearman's rank correlation
# ------------------------------------------------------------------
# Join the data frames to add all the counts to the gene list, and then remove rows with tracker = NA
allspear_GeneGO <- full_join(all_1, GeneGO_fm_df_u, by="human")
allspear_GeneGO <- filter(allspear_GeneGO, !is.na(allspear_GeneGO$tracker))

# Get the data frame ready for plotting
allspear_GeneGO$GeneHuman <- paste0(allspear_GeneGO$gene, ' / ', allspear_GeneGO$human) # make killifiish/human names
rownames(allspear_GeneGO) <- allspear_GeneGO$GeneHuman # set killifiish/human names as row names
allspear_GeneGO_plot = subset(allspear_GeneGO, select = -c(human, tracker, gene, GeneHuman, padj_spear_m, padj_spear_f))
allspear_GeneGO_plot <- unique(allspear_GeneGO_plot) # plot only unique genes



# Plot heatmap for the selected genes
pdf(paste0('Plots/Heatmap/241209_allCorrresults', t ,'_', termName, '_GOGSEA_heatmap_241228.pdf'), width = 8, height = 12)
pheatmap(allspear_GeneGO_plot,
         cluster_rows = TRUE,
         color = colorRampPalette(c("blue", "white", "red"))(50), 
         breaks = seq(-1, 1, length.out = 51),
         fontsize = 9,         # Adjust overall text size
         fontsize_row = 3,     # Adjust row label text size
         fontsize_col = 12,       # Adjust column label text size# Blue-to-red gradient
         main = paste0(t ,' ', termName, ' ', termID, ' GOGSEA heatmap')
)
dev.off()


rm(list = ls())

# ------------------------------------------------------------------
sessionInfo() 

