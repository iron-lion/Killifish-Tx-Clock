# Hierarchical clustering  

library(tidyverse)
#vignette - https://tavareshugo.github.io/data-carpentry-rnaseq/04b_rnaseq_clustering.html#Clustering_basics
# https://scienceparkstudygroup.github.io/rna-seq-lesson/08-cluster-analysis/index.html

# load data
# Load necessary libraries
library(DESeq2)
library(dplyr)
library(ggplot2)
library(stats)

#Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#93cca8','#93cca8','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#ab5673', '#f1a8a4','#ef9ac2' ) 

setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/05_HClustering_GeneTrajectories/')

#------------------------------------------------------------------
# Load Data
# ------------------------------------------------------------------
#First, we load our Deseq2 object as it contains all the counts and metadata & the corr res
object.indir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))
load(paste0(object.indir,'dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')) #sexes combined corr results

indir = '/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/Output/Tables/'
num_Expressed <- read.csv(paste0(indir, '240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv'), row.names = 1)

outdir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/05_HClustering_GeneTrajectories/Output/"

#load merged tpm matrix
tpm <- read.csv('path/to/TPM_Atlas_allbatches_merged_v3.csv', row.names = 1)

#------------------------------------------------------------------
# Define genes expressed in each tissue, defined by TPM cutoff
# ------------------------------------------------------------------
tissues = names(BulkSeq_Agecorrelation_results_list)

ExpressedGenes <- list()
for(tissue in tissues){
  ExpressedGenes_geneset <- read.csv(paste0(outdir,"240712_expressedgeneset-corranalysisinput_TPMcutoff_universe_",tissue,".csv"), row.names = 1)
  ExpressedGenes[[tissue]] <- ExpressedGenes_geneset$ExpressedGenes
}

Expr.geneset = plyr::ldply(ExpressedGenes, rbind)
rownames(Expr.geneset) = Expr.geneset$.id
Expr.geneset$.id = NULL
Expr.geneset = t(Expr.geneset)

x <- Reduce(intersect, ExpressedGenes)
y <- data.frame(GeneID = x) #10,847 genes expressed in all tissues, from the gene set of genes expressed in each tissue

# save the genes expressed in all tissues
write.csv(y,file = paste0(outdir,'240905_geneexpressedinalltissues.csv'))

#------------------------------------------------------------------
# Tissue-by-tissue gene expression trajectory clustering for 10,847 genes expr in all tissues - zscaling ncounts w/o median expression per age group
# ------------------------------------------------------------------
clust.outdir.v1 = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Clustering_Hierarch/output/2024_12_31_outs/'

#gene.set = 10,847 genes expressed in all tissues
gene.set <- read.csv(file = paste0(outdir,'240905_geneexpressedinalltissues.csv'), row.names = 1)

tissue.list = names(dds_TPS_list)[2:14]

### adjusted code
tissue.list2 = c('Liver', 'Kidney', 'Muscle', 'Brain', 'Gut')
#tissue = 'Liver' #debugging

for(tissue in tissue.list2){
  #dir.create(file.path(paste0(clust.outdir.v1, tissue)), showWarnings = T)
  
  dds <- dds_TPS_list[[tissue]]
  
  sampleTable <- as.data.frame(dds@colData)
  countTable <- counts(dds, normalized = TRUE) #normalize by tissue
  
  #Ensure samples are ordered by age
  sampleTable_filtered <- sampleTable %>% filter(!age_days %in% c(102,103)) #remove the middle time point
  sampleTable.ordered <- sampleTable_filtered %>% arrange(age_days)
  ordered_samples <- sampleTable.ordered$lib
  
  # Norm count subset and ordering
  norm_counts <- countTable[, ordered_samples]
  norm_counts <- norm_counts[gene.set$GeneID,] #only keep genes expressed in all tissues
  
  # Z-scale the normalized counts
  z_scaled_counts <- t(scale(t(norm_counts)))
  
  # Estimate gene trajectories using LOESS regression
  ages <- as.numeric(as.character(unique(sampleTable.ordered$age_days)))
  loess_fits <- data.frame(matrix(nrow = nrow(z_scaled_counts), ncol = length(ages)))
  rownames(loess_fits) <- rownames(z_scaled_counts)
  colnames(loess_fits) <- ages
  
  for (gene in rownames(z_scaled_counts)) {
    gene_data <- data.frame(Age = sampleTable.ordered$age_days, Expression = z_scaled_counts[gene, ])
    gene_data$Age <- as.numeric(as.character(gene_data$Age))
    loess_fit <- loess(Expression ~ Age, data = gene_data)
    loess_fits[gene, ] <- predict(loess_fit, ages)
  }
  
  # Compute a distance matrix between gene trajectories using Euclidean distance
  dist_matrix <- dist(loess_fits, method = "euclidean")
  
  # Perform hierarchical clustering using the complete method
  hc <- hclust(dist_matrix, method = "complete")
  
  # Cut the dendrogram into clusters
  num_clusters <- 10 # Adjust this number as needed
  clusters <- cutree(hc, k = num_clusters)
  
  # Add cluster information to the data frame
  loess_fits$Cluster <- as.factor(clusters)
  
  write.csv(as.data.frame(loess_fits), file = paste0(clust.outdir.v1, tissue, "/Atlas_hclust_loessfits_clusterdata_", tissue, "_confinterval_241231.csv"))
  
  # Visualization of gene expression dynamics for each cluster
  for (cluster in levels(loess_fits$Cluster)) {
    genes_in_cluster <- rownames(loess_fits[loess_fits$Cluster == cluster, ])
    cat("Cluster", cluster, ": Genes =", paste(genes_in_cluster, collapse = ", "), "\n")
    
    cluster_data <- loess_fits[loess_fits$Cluster == cluster, -ncol(loess_fits)]
    cluster_data_long <- tidyr::pivot_longer(as.data.frame(t(cluster_data)), 
                                             cols = everything(), 
                                             names_to = "Gene", 
                                             values_to = "Expression")
    cluster_data_long$Age <- rep(ages, each = nrow(cluster_data))
    
    # Calculate mean and SEM for the cluster
    summary_data <- cluster_data_long %>%
      group_by(Age) %>%
      summarise(
        MeanExpression = mean(Expression, na.rm = TRUE),
        SEM = sd(Expression, na.rm = TRUE) / sqrt(n())
      )
    summary_data <- summary_data %>%
      mutate(
        LowerCI = MeanExpression - SEM,
        UpperCI = MeanExpression + SEM
      )
    
    # Plot with average trajectory and confidence intervals
    pdf(file = paste0(clust.outdir.v1, tissue, "/Atlas_hclust_cluster", cluster, "_", tissue, "_confinterval_241231.pdf"), width = 2, height = 2)
    print(
      ggplot() +
        # Plot individual gene trajectories with transparency
        geom_line(data = cluster_data_long, aes(x = Age, y = Expression, group = Gene, color = Gene), alpha = 0.3) +
        # Plot average trajectory
        geom_line(data = summary_data, aes(x = Age, y = MeanExpression), color = "black", linewidth = 1) +
        # Add confidence intervals
        geom_ribbon(data = summary_data, aes(x = Age, ymin = LowerCI, ymax = UpperCI), fill = "blue", alpha = 0.2) +
        # Labels and theme
        labs(
          title = paste("Gene Expression Dynamics for Cluster", cluster),
          x = "Age",
          y = "Expression (Z-scaled)"
        ) +
        theme_minimal() +
        theme(legend.position = "none")
    )
    dev.off()
  }
}



