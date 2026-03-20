# Title: Plot Hypergeometric GO analysis results
# Author: Jingxun Chen
# Date: code compiled on 20250120
# Related publication: Emma K. Costa and Jingxun Chen, Nature Aging, 2025

# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxun/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/HyperGO/")
set.seed(1234)

# Load packages
library('ggplot2')
library('dplyr')


# ------------------------------------------------------------------
# Load data, clean up
# ------------------------------------------------------------------
# Upregulated GO terms (age-correlated genes shared by 5 or more tissues)
up <- read.csv('Output/GOBP_Atlas_POScorrgene_AtLeast5Tissues.csv', header = TRUE) # upregulated

# Downregulated GO terms (age-correlated genes shared by 5 or more tissues)
down <- read.csv('Output/GOBP_Atlas_NEGcorrgene_AtLeast5Tissues.csv', header = TRUE) # downregulated

# Extract only the enrichment and padj terms
up <- up[c('GOBPID', 'Term', 'enrichment', 'padj')]
down <- down[c('GOBPID', 'Term', 'enrichment', 'padj')]

# Make the GO term names as row names
rownames(up) <- up$GOBPID 
rownames(down) <- down$GOBPID


# ------------------------------------------------------------------
# Subset specific GO terms, define plotting order
# ------------------------------------------------------------------
# Subset specific GO terms
up_go <- up[c('GO:0050776', 'GO:0050778', 'GO:0045087', 'GO:0019221'), ]
down_go <- down[c('GO:0043062','GO:0030198', 'GO:0071230', 'GO:0007584'), ]

# order the dataframe based on decreasing padj values
up_go <- up_go[order(up_go$padj), ]
down_go <- down_go[order(down_go$padj), ]

# Change the level of 'Term' based on the order I defined here.
up_go$Term <- factor(up_go$Term, levels=rev(unique(up_go$Term)))
up_go['Condition']<- 'upregulated'

down_go$Term <- factor(down_go$Term, levels=rev(unique(down_go$Term)))
down_go['Condition']<- 'downregulated'

# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------
pdf('Plots/GOBP_Atlas_POScorrgene_AtLeast5Tissues_bubble_250120.pdf', width = 8, height = 10)
ggplot(data=up_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "red", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for upregulated genes shared by 5 or more tissues', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 

pdf('Plots/GOBP_Atlas_NEGcorrgene_AtLeast5Tissues_bubble_250120.pdf', width = 8, height = 10)
ggplot(data=down_go, aes(x=Condition, y=Term)) + 
  geom_point(aes(color=enrichment, size = -log10(padj))) + 
  scale_colour_gradient2(high = "dark blue", space = "Lab", na.value = "grey50") +
  labs(title = 'GO terms for downregulated genes shared by 5 or more tissues', x='Conditions', y = '',size='-log(FDR)', color='Enrichment') 
dev.off() 



sessionInfo() 
