##########################################################################################################
###
###              RNAseq analysis  - Patricia Moran Losada
###                   adapted by Emma Costa 
###
##########################################################################################################

# ------------------------------------------------------------------
# Please run 000_merge_counts_TPM.R to get merged TPM matrices
# ------------------------------------------------------------------

setwd("/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/02_QC/")

library(reshape) #0.8.9
library(ggpubr) #0.6.0
library(Hmisc)
library(WGCNA) 
library(ggplot2)
library(reshape)
library(ggrepel)
library(DESeq2)
library(vsn)
library(pheatmap)
library(edgeR)
library(biomaRt)
library(magrittr) 
library(dplyr) 
library(DEGreport)
library(ggpubr)
library (plyr)
library(stringr)
library(nlme)
library(gprofiler2)
library(gridExtra)


tpm <- read.csv('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/02_QC/Output/TPM_Atlas_allbatches_merged_v3.csv', stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1)
metadata <- read.csv(file="/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2/Input/ExperimentDesign_allbatches_combined_v7.csv", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1 )

col_order <- rownames(metadata)
tpm <- tpm[,col_order]

table(rownames(metadata)==colnames(tpm)) #check order columns/rows in tpm and metadata file

################################################################
## Remove genes with low TPM expression values
################################################################

# filter genes with TPM 0.5 in 15% of samples
pres = apply(tpm>=0.5,1,sum) 
Expressed = (pres >= 0.15*ncol(tpm))
tpm_expressed = tpm[Expressed,]

gene_names <- as.data.frame(rownames(tpm_expressed))
colnames(gene_names)[1] <- "gene"


################################################################
#Remove outlier samples based on gene expression connectivity
################################################################
all_outliers = character()
par(mfrow=c(1,1), mar=c(5,4,2,2))
metadata$ConnectivityZscore = NA

idx = match(colnames(tpm_expressed), rownames(metadata))
normadj <- (0.5+0.5*bicor(tpm_expressed))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku)) ## Z-score the connectivity
metadata$ConnectivityZscore[idx] = z.ku


#Bone, Brain, Eye, Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 

tissues <- sort(unique(metadata$tissue))

for(i in 1:length(tissues)){
  t = tissues[i] #tissue name
  metadata$Colour[metadata$tissue == t] = colTPS[i]
}

outliers = (z.ku < -2)
table(outliers)
all_outliers = c(all_outliers, colnames(tpm_expressed)[outliers])

#label the outliers
for(j in 1:length(all_outliers)){
  o = all_outliers[j]
  metadata$outlier_status[metadata$lib == o] = "Outlier"
}

to_change <- which(metadata$outlier_status == 'Outlier')
metadata$outlier_label = NA
for(k in to_change){
  metadata[k,]$outlier_label = metadata[k,]$lib
}

plot(1:length(z.ku),z.ku,col=metadata$Colour,pch=19,main=paste("Outliers"),
     ps=30, ylab="Connectivity Z score",xlab="", ylim=c(-11,11))
abline(h=-2, lty=2)
text(1:length(z.ku),z.ku, labels = metadata$outlier_label)

plot.new()


write.csv(metadata, file="Output/ExperimentDesign_allbatches_combined_with_Connectivity.csv")


##ggplot
# metadata <- metadata[order(metadata$tissue),] #order metadata by tissue
# metadata$RNA_ID <- factor(metadata$RNA_ID,levels=metadata$RNA_ID)
# data=metadata[metadata$ConnectivityZscore <= -2 | metadata$ConnectivityZscore >= 2,]
# 
# g1 = ggplot(metadata, aes(x=RNA_ID, y=ConnectivityZscore, color=tissue)) + geom_point(size = 3) + 
#   geom_hline(yintercept = -2, lty=2) +  
#   geom_hline(yintercept = 2, lty=2) +
#   theme(axis.text.x = element_blank(),axis.ticks.x = element_blank()) +
#   xlab("Samples") + ggtitle("Outlier Sample ") +ylim(-3,2.5)+
#   theme(panel.background = element_rect(fill = NA))+
#   theme(panel.border = element_rect(color = "black", fill = NA, size = 1),
#         strip.background = element_rect(color = "black", size = 1))+
#   ylab("Network Connectivity (Z score)")  +geom_text_repel(data=data,
#                                                            aes(x=data$RNA_ID, y=data$ConnectivityZscore,label = RNA_ID),box.padding = unit(0.5, "lines"))
# 
# g1
