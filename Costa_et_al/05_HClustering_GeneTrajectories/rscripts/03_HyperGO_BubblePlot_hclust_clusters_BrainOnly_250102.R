# Title: Plot Hypergeometric GO analysis results
# Author: Jingxun Chen
# Date: code compiled on
# Related publication: 


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
# Load data, clean up
# ------------------------------------------------------------------
# Brain clustering data - GOBP
tissue = 'Brain'
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
# Brain clustering data - Gene GOBP
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
IDs.bycluster <- c('GO:0007049', 'GO:0030182', 'GO:0048699','GO:0048858', 'GO:0007409','GO:0007623','GO:0019953','GO:0006695','GO:0046165', #cluster1
                   'GO:0007399', 'GO:0099177', 'GO:0097485', 'GO:0007219', 'GO:0048813',#cluster2
                   'GO:0006397', 'GO:0006357', 'GO:0000398','GO:0048024', 'GO:0031123', #cluster3
                   'GO:0034440', 'GO:0019395',#cluster4
                   'GO:0006955','GO:0042119', 'GO:0002274', 'GO:0002250', 'GO:0001816', 'GO:0006915', #cluster5
                   'GO:0032373', 'GO:0032374', 'GO:0010875','GO:0032373', #cluster6
                   'GO:0006099','GO:0072524','GO:0009117', 'GO:0009259', #cluster7
                   'GO:0022900','GO:0007005','GO:0006119','GO:0009060', #cluster8
                   'GO:0043547','GO:0022008', #cluster9
                   'GO:0030509','GO:0016226'#cluster10 **not significant
                   )

# ------------------------------------------------------------------
# Subset specific GO terms, define plotting order
# ------------------------------------------------------------------
# Subset specific GO terms
# IDs.top <- c()
# for(clust in unique(GOBP_allclust_sig$cluster)){
#   subset <- GOBP_allclust_sig %>% filter(cluster == clust)
#   subset.ordered <- subset[,order('cluster')] 
#   top <- subset.ordered[1]
#   IDs.top <- c(IDs.top, top)
# }

IDs.top = IDs.bycluster

GOBP_allclust_go <- GOBP_allclust_sig %>% filter(GOBPID %in% IDs.top)

GOBP_allclust_go <- GOBP_allclust_go[order(GOBP_allclust_go$GOBPID), ]
GOBP_allclust_go$cluster <- as.numeric(GOBP_allclust_go$cluster)

GOBP_allclust_go$Genes <- NA

# add the genes associated with the go term - this is actually the full set size for the GO term, will need to overlap with cluster genes
for(row in 1:nrow(GOBP_allclust_go)){
  temp <- GOBP_allclust_go[row,]
  lookup <- Gene_GOBP_allclust %>% filter(GOBPID == temp$GOBPID & cluster == temp$cluster)
  GOBP_allclust_go[row,]$Genes <- lookup$Genes
}

# Get the actual genes assigning to each cluster
clust.id <- unique(clusters$cluster)
GOBP_allclust_go$geneset <- NA

for(row in 1:nrow(GOBP_allclust_go)){
  temp <- GOBP_allclust_go[row,]
  clust <- temp$cluster
  select.cluster <- clusters %>% filter(Cluster == clust)
  cluster.genes <- rownames(select.cluster)
  full.genelist <- strsplit(temp$Genes, " ")
  GO.genes <- intersect(cluster.genes, full.genelist[[1]])
  GOBP_allclust_go[row,]$geneset <- list(GO.genes) 
}

GOBP_allclust_go$Genes <- NULL

# Convert the LOC's
nfur.ortho <- read.csv(file = 'Supporting_files/nfur-ncbi_orthologs_20170302_pps.csv')
nfur.ortho <- as.data.frame(nfur.ortho)
nfur.ortho$NCBI.Definition <- str_replace(nfur.ortho$NCBI.Definition,"\\s[\\(][A-Za-z0-9]*[\\)]" ,"")
names(nfur.ortho)[1] <- "Killifish"

GOBP_allclust_go <- as.data.frame(GOBP_allclust_go)

# Initialize the 'Human' column as a character vector 
GOBP_allclust_go$Human <- character(nrow(GOBP_allclust_go))
GOBP_allclust_go$Mouse <- character(nrow(GOBP_allclust_go))
GOBP_allclust_go$Zebrafish <- character(nrow(GOBP_allclust_go))

for(row in 1:nrow(GOBP_allclust_go)){
  temp <- GOBP_allclust_go[row,]
  genelist <- temp$geneset[[1]]  # Extract the list of genes for this row
  human_genes <- c()  # Temporary character vector to store human gene names for this row
  for(gene in genelist){
    if(gene %in% nfur.ortho$Killifish){
      conv <- nfur.ortho[nfur.ortho$Killifish == gene,]$Human
      if(length(conv) > 0) {
        human_genes <- c(human_genes, conv)  # Concatenate converted gene names
      }
    }
  }
  # Collapse the human_genes vector into a single, comma-separated string if needed
  GOBP_allclust_go$Human[row] <- paste(human_genes, collapse = ", ")
}

for(row in 1:nrow(GOBP_allclust_go)){
  temp <- GOBP_allclust_go[row,]
  genelist <- temp$geneset[[1]]  # Extract the list of genes for this row
  ms_genes <- c()  # Temporary character vector to store mouse gene names for this row
  for(gene in genelist){
    if(gene %in% nfur.ortho$Killifish){
      conv <- nfur.ortho[nfur.ortho$Killifish == gene,]$Mouse
      if(length(conv) > 0) {
        ms_genes <- c(ms_genes, conv)  # Concatenate converted gene names
      }
    }
  }
  # Collapse the ms_genes vector into a single, comma-separated string if needed
  GOBP_allclust_go$Mouse[row] <- paste(ms_genes, collapse = ", ")
}


for(row in 1:nrow(GOBP_allclust_go)){
  temp <- GOBP_allclust_go[row,]
  genelist <- temp$geneset[[1]]  # Extract the list of genes for this row
  dr_genes <- c()  # Temporary character vector to store zebrafish gene names for this row
  for(gene in genelist){
    if(gene %in% nfur.ortho$Killifish){
      conv <- nfur.ortho[nfur.ortho$Killifish == gene,]$Zebrafish
      if(length(conv) > 0) {
        dr_genes <- c(dr_genes, conv)  # Concatenate converted gene names
      }
    }
  }
  # Collapse the dr_genes vector into a single, comma-separated string if needed
  GOBP_allclust_go$Zebrafish[row] <- paste(dr_genes, collapse = ", ")
}


# Convert list columns to character columns by collapsing each list to a single string
GOBP_allclust_go[] <- lapply(GOBP_allclust_go, function(col) {
  if (is.list(col)) {
    sapply(col, function(x) paste(x, collapse = ", "))  # Collapse lists to comma-separated strings
  } else {
    col  # Leave non-list columns unchanged
  }
})

# save table
write.table(GOBP_allclust_go, file = 'Output/250102_Brain_topGOterms-per-hClust_genes.csv', sep = ",", row.names = FALSE)


# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------
library(viridis)
library(hrbrthemes)

pdf('Output/Plots/Brain_hClust_GO_bubble_250102.pdf', width = 8, height = 3)
  ggplot(data=GOBP_allclust_go, aes(x=cluster, y=Term, size = -log10(padj), fill = enrichment)) + 
    geom_point(shape=21, color="black") + 
    scale_fill_viridis_c(limits = c(0, 15), oob = scales::squish) +
    theme_classic() +
    labs(x='cluster', y = '',size='-log10(FDR)', color='Enrichment') + 
    theme(legend.position = 'right')
dev.off() 


# ------------------------------------------------------------------
# Subset specific GO terms, define plotting order
# ------------------------------------------------------------------
rm(GOBP_allclust_go)
IDs.byhallmarks <- IDs.bycluster <- c('GO:0007049', 'GO:0030182', 'GO:0048699','GO:0048858', 'GO:0007409','GO:0007623','GO:0019953','GO:0006695','GO:0046165', #cluster1
                                      'GO:0007399', 'GO:0099177', 'GO:0097485', 'GO:0007219', 'GO:0048813',#cluster2
                                      'GO:0006397', 'GO:0006357', 'GO:0000398','GO:0048024', 'GO:0031123', #cluster3
                                      'GO:0034440', 'GO:0019395',#cluster4
                                      'GO:0006955','GO:0042119', 'GO:0002274', 'GO:0002250', 'GO:0001816', 'GO:0006915', #cluster5
                                      'GO:0032373', 'GO:0032374', 'GO:0010875', #cluster6
                                      'GO:0006099','GO:0072524','GO:0009117', 'GO:0009259', #cluster7
                                      'GO:0022900','GO:0007005','GO:0006119','GO:0009060', #cluster8
                                      'GO:0043547','GO:0022008', #cluster9
                                      'GO:0030509','GO:0016226'#cluster10 **not significant
)

IDs.byhallmarks <- unique(IDs.byhallmarks)
GOBP_allclust_go <- GOBP_allclust_sig %>% filter(GOBPID %in% IDs.byhallmarks) #cluster 10 terms nonsig
clust10_nonsig <- GOBP_allclust %>% filter(GOBPID %in% c('GO:0030509','GO:0016226'))
clust10_nonsig <- clust10_nonsig %>% filter(cluster =='10')
GOBP_allclust_go <- rbind(GOBP_allclust_go, clust10_nonsig)

order.df <- data.frame(IDs.byhallmarks)
colnames(order.df) = c('GOBPID')
rownames(order.df) = order.df$GOBPID

order.df$hallmark = c(rep(1,9),
                      rep(2,5),
                      rep(3,5),
                      rep(4,2),
                      rep(5,6),
                      rep(6,3),
                      rep(7,4),
                      rep(8,4),
                      rep(9,2),
                      rep(10,2))

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

pdf('Output/Plots/Brain_hClust_selectclusterterms_GO_bubble_250102.pdf', width = 8, height = 5)
  ggplot(data=GOBP_allclust_go, aes(x=cluster, y=Term, size = -log10(padj), fill = enrichment)) + 
    geom_point(shape=21, color="black") + 
    scale_fill_viridis_c(limits = c(0, 20), oob = scales::squish) +
    theme_classic() +
    scale_y_discrete(limits = rev) +
    #scale_x_discrete(limits = c(1,10)) +
    labs(x='cluster', y = '',size='-log10(FDR)', color='Enrichment') + 
    theme(legend.position = 'right')
dev.off() 


sessionInfo() 
