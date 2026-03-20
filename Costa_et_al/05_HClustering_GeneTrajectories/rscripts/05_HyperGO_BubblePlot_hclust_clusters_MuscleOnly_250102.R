# Title: Plot Hypergeometric GO analysis results
# Author: Jingxun Chen
# Date: code compiled on 20220805
# Related publication: Andrew McKay, Emma K. Costa, and Jingxun Chen, eLife, 2022


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
indir = '/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/05_HClustering_GeneTrajectories/Output/'

# Set wd to the current directory and seed (for reproducible results)
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/05_HClustering_GeneTrajectories/')
set.seed(1234)

# Load packages
library('ggplot2')
library('dplyr')
library(stringr)


# ------------------------------------------------------------------
# Muscle
# ------------------------------------------------------------------
# Muscle clustering data - GOBP
tissue = 'Muscle'
Cluster_GOBP <- list()
for(clust in 1:10){
  Cluster_GOBP_temp <- read.delim(paste0("Output/GOBP_250101_Atlas_",tissue,"hClust",clust,".txt"), row.names = 1)
  Cluster_GOBP_temp$GOBPID <- rownames(Cluster_GOBP_temp)
  Cluster_GOBP_temp$cluster <- clust
  Cluster_GOBP[[clust]] <- Cluster_GOBP_temp
}

GOBP_allclust = plyr::ldply(Cluster_GOBP, rbind) #note there are duplicate GO terms, so I cannot simply rename the rows

# Extract only the enrichment and padj terms
GOBP_allclust <- GOBP_allclust[c('GOBPID', 'Term', 'enrichment', 'padj', 'cluster')]

# Remove GO terms with padj > 0.05
GOBP_allclust_sig <- GOBP_allclust %>% filter(padj < 0.05)

# ------------------------------------------------------------------
# Muscle clustering data - Gene GOBP
Cluster_Gene_GOBP <- list()
for(clust in 1:10){
  Cluster_Gene_GOBP_temp <- read.delim(paste0("Output/GeneGOBP_240101_Atlas_",tissue,"hClust",clust,".txt"), header = F)
  colnames(Cluster_Gene_GOBP_temp) = 'Gene'
  Cluster_Gene_GOBP_temp[c('GOBPID', 'Genes')] <- str_split_fixed(Cluster_Gene_GOBP_temp$Gene, ' ', 2)
  Cluster_Gene_GOBP_temp$Gene = NULL
  Cluster_Gene_GOBP_temp$cluster <- clust
  Cluster_Gene_GOBP[[clust]] <- Cluster_Gene_GOBP_temp
}

Gene_GOBP_allclust = plyr::ldply(Cluster_Gene_GOBP, rbind) #note there are duplicate GO terms, so I cannot simply rename the rows
# ------------------------------------------------------------------
clusters <- read.csv(paste0('Output/Atlas_hclust_loessfits_clusterdata_',tissue,'_confinterval_241231.csv'), row.names =1)
table(clusters$Cluster)
# ------------------------------------------------------------------
# Interesting terms by eye, by cluster

#Muscle
GOBP_allclust_select <- GOBP_allclust_sig %>% filter(cluster == '7')
IDs.bycluster <- c('GO:0044786','GO:0000278','GO:0044770','GO:0006260',#cell cycle
                   'GO:0006457', 'GO:0006487', 'GO:0061448', 'GO:0006941', 'GO:0043403')

# ------------------------------------------------------------------
# Subset specific GO terms, define plotting order
# ------------------------------------------------------------------
#rm(GOBP_allclust_go)
IDs.top = IDs.bycluster

GOBP_allclust_go <- GOBP_allclust_sig %>% filter(GOBPID %in% IDs.top)

GOBP_allclust_go <- GOBP_allclust_go[order(GOBP_allclust_go$GOBPID), ]
GOBP_allclust_go$cluster <- as.numeric(GOBP_allclust_go$cluster)

order.df <- data.frame(IDs.bycluster)
colnames(order.df) = c('GOBPID')
rownames(order.df) = order.df$GOBPID

order.df$hallmark = c(rep(1,4), 2,3,4,5,6)

GOBP_allclust_go$hallmark <- NA

for(i in 1:nrow(GOBP_allclust_go)){
  row.temp = GOBP_allclust_go[i,]
  ID = row.temp$GOBPID
  h = order.df[ID,]$hallmark
  GOBP_allclust_go[i,]$hallmark = h
}


# order the dataframe based on hallmark
GOBP_allclust_go <- GOBP_allclust_go[order(GOBP_allclust_go$hallmark), ]
GOBP_allclust_go$cluster <- as.numeric(GOBP_allclust_go$cluster)

GOBP_allclust_go <- GOBP_allclust_go

GOBP_allclust_go$hallmark <- factor(GOBP_allclust_go$hallmark, levels = sort(unique(GOBP_allclust_go$hallmark )))
GOBP_allclust_go$Term <-factor(GOBP_allclust_go$Term, levels = unique(GOBP_allclust_go$Term))


# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------
library(viridis)
library(hrbrthemes)

#GOBP_allclust_go <- GOBP_allclust_go %>% filter(cluster == '7')

### ExtData Fig 5D ###
pdf('Output/Plots/Muscle_hClust_GO_bubble_250102.pdf', width = 4, height = 3)
ggplot(data=GOBP_allclust_go, aes(x=cluster, y=Term, size = -log10(padj), fill = enrichment)) + 
  geom_point(shape=21, color="black") + 
  scale_fill_viridis_c(limits = c(0, 15), oob = scales::squish) +
  theme_classic() +
  scale_y_discrete(limits = rev) +
  labs(x='cluster', y = '',size='-log10(FDR)', color='Enrichment') + 
  theme(legend.position = 'right')
dev.off() 
 