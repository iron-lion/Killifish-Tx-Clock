
##Here, we perform will correlate expression for each gene in each tissue/region with age and quantify the number of significant age-correlated genes

library(gprofiler2)
library(ComplexHeatmap)
library(UpSetR)
library(hrbrthemes)
library(dplyr)
library("RColorBrewer")
library("circlize")


setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/')

#Sex: male, female 
colsex <-c('#068ec9','#ba1e2d')

#Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#7961a2','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#ab5673', '#f1a8a4','#ef9ac2' ) 

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))

colage <- c('#395ea4','#395ea4','#395ea4','#3789c0','#3789c0','#3789c0',  '#dfeaf0','#dfeaf0', '#f3d4ac','#f3d4ac', '#f7b77b','#f7b77b','#f7b77b', '#a92b46','#a92b46')



#First, we load our Deseq2 object as it contains all the counts and metadata
object.indir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))


dds_CA1_list <- dds_TPS_list # there was an instance where the list of deseq objects was called dds_CA1_list. 
#------------------------------------------------------------------
# Age-corr genes from TPM filtering method
# ------------------------------------------------------------------
#load merged tpm matrix
tpm <- read.csv('path/to/TPM_Atlas_allbatches_merged_v3.csv', row.names = 1)

#Designate with sex you'll be using for this analysis, options = 'Both', 'Male', 'Female'
sex <- 'Both'

#We'll setup a results list in which we'll store the output of the following analysis
BulkSeq_Agecorrelation_results_list <- list()

#We will rotate over each tissue/region (except the 'All') object
for (tissue in everything_but(names(dds_CA1_list), 'All')) {
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]
  
  #We use the DESEQ2-normalized data, which we do by getting the counts them from the DESeq2 object
  count_table <- counts(dds_tissue, normalized=T)
  count_table <- as.data.frame(count_table)
  
  #To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
  #we consider as 'expressed' - in this case these will be genes within the tissue that have TPM > 0.5 for 80% of the samples
  
  #find sample names for tissue
  coldata_tissue <- dds_tissue@colData

  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples 
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(coldata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue)) 
  print(paste0(tissue," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]

  ExpresedGenes <- rownames(expressed)
  #We define the numeric values we want to correlate against - in this case, the age
  Score_to_compare <- as.numeric(as.character(as.data.frame(colData(dds_tissue))[,'age_days']))
  
  #Perform gene-vs-age score correlataion
  #To do that we set up a dataframe that can hold the results
  #We want to record both the spearman and pearson correlation coefficient, as well as the results from cor.test to see that we have statistical
  #confidence in our correlation coefficents
  correlation_frame <- data_frame(gene='A', cor_spear=1, p_spear=1, cor_pears=1, p_pears=1)
  
  #We set up a counter so that we can have a progressbar for our analysis
  i <- 1
  pb = txtProgressBar(min = 1, max = length(ExpresedGenes), initial = 1, style = 3)
  #We will now iterate over very expresed gene and perform the following steps
  for (gene in ExpresedGenes) {
    setTxtProgressBar(pb,i)
    #We set up a temporary dataframe that holds age and normalized counts 
    cor_frame <- data.frame(expr=t(count_table[gene,]), concentration= Score_to_compare)
    colnames(cor_frame) <- c('expr', 'concentration')
    #We wil only retain samples where the count is not 0
    cor_frame <- cor_frame %>% filter(expr > 0)
    #After filtering for non-expressing samples, we will only continue the analysis if we have at least 20 samples letf - otherwise we end up
    #with very spurious correlations
    if (nrow(cor_frame) >= 20) {
      #Now we run the correlation tests for spearman and pearson
      cor_vector_spearman <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'spearman')
      cor_vector_pearson <- cor.test(cor_frame$expr, cor_frame$concentration, method = 'pearson')
      #and we now appending to the correlation_frame a new row that contains all the coeffiients and pvalues
      correlation_frame[i, 'gene'] <- gene
      correlation_frame[i, 'cellsExpr'] <- nrow(cor_frame)
      correlation_frame[i, 'cor_spear'] <- cor_vector_spearman$estimate
      correlation_frame[i, 'p_spear'] <- cor_vector_spearman$p.value
      correlation_frame[i, 'cor_pears'] <- cor_vector_pearson$estimate
      correlation_frame[i, 'p_pears'] <- cor_vector_pearson$p.value
      
      #if the gene had fewer than 20 samples expressing it, we'll add a dummy row
    } else {
      correlation_frame[i, 'gene'] <- gene
      correlation_frame[i, 'cellsExpr'] <- nrow(cor_frame)
      correlation_frame[i, 'cor_spear'] <- NA
      correlation_frame[i, 'p_spear'] <- NA
      correlation_frame[i, 'cor_pears'] <- NA
      correlation_frame[i, 'p_pears'] <- NA
      
    }
    #We add 1 to the conter for the progress bar
    i <- i + 1
  }
  #Cleanup results and run the multiple testing correction
  correlation_frame <- na.omit(correlation_frame)
  correlation_frame$padj_spear <- p.adjust(correlation_frame$p_spear, method = 'BH')
  correlation_frame$padj_pear <- p.adjust(correlation_frame$p_pears, method = 'BH')
  #Relabel
  tissue_correlation_fullTable <- correlation_frame
  rm(correlation_frame)
  #Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibhitng an absolute correlation value of 0.5)
  print(paste('TopCorGenes: ', nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
  nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05))
  #Store data the results in our correlation results list
  BulkSeq_Agecorrelation_results_list[[tissue]][['resall']] <- tissue_correlation_fullTable
  #Filter and store data for high correlates in our correlation results list
  BulkSeq_Agecorrelation_results_list[[tissue]][['topCor']]  <- (tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))
  
  
  
}


####

#We will save the results list for later use
save(BulkSeq_Agecorrelation_results_list, file='Robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')


#Load if you are starting from here
#load('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')


####
#####save the geneset of ExpressedGenes 
outdir = '/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/Output/'

#tissue = 'Gonad' #debugging
Expressed_Gene_list <- data.frame(id = 'A')
for(tissue in names(BulkSeq_Agecorrelation_results_list)){
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]
  #We use the DESEQ2-normalized data, which we do by getting the counts them from the DESeq2 object
  count_table <- counts(dds_tissue, normalized=T)
  count_table <- as.data.frame(count_table)
  #find sample names for tissue
  coldata_tissue <- dds_tissue@colData
  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples 
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(coldata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue)) 
  print(paste0(tissue," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]
  ExpresedGenes <- rownames(expressed)
  expr.gene.df <- data.frame(id = ExpresedGenes)
  Expressed_Gene_list <- rbind(Expressed_Gene_list,expr.gene.df)
  #expr.gene.df <- data.frame(ExpressedGenes = ExpresedGenes)
  #write.csv(expr.gene.df, paste0(outdir,"240714_expressedgeneset-corranalysisinput_TPMcutoff_universe_",tissue,".csv"))
}

Expressed_Gene_list <- Expressed_Gene_list[-1,]
Expressed_Gene_list <- as.data.frame(Expressed_Gene_list)
Expressed_Gene_list <- Expressed_Gene_list %>% distinct()
colnames(Expressed_Gene_list) <- c('id')

write.csv(Expressed_Gene_list, paste0(outdir, "Tables/240714_allexpressedgenes_TPMcutoff_universe.csv"))

#------------------------------------------------------------------
# TPM filtering / Further analysis of age-correlated genes
# ------------------------------------------------------------------
outdir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/Output/"
load('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/Output/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')


#---------------- Bar plot of total # age corr genes ----------------
#Now that we've run the correlation analysis we want to plot how many genes correlated significantly with age in each tissue

#We will set a range of pvalue cutoffs, to see if the results remain stable 
padj_range <-c(0.001, 0.01, 0.05, 0.1)
#We have to set a correlation coefficient cutoff (absolute cutoff of 0.5)
cor_cutoff <- 0.5

#We will have to establish a 'master data frame' with some dummy data. We will iteratively add data to this data frame
plot_df_crossTissue <- data.frame(padj_cutoff=0.05, Type='A', Comparison='A', Numb_of_genes=1, tissue='A')

results_list_CA1_bothSex <- BulkSeq_Agecorrelation_results_list
#We will iterate over each of the tissues, get the correlation results to build the bargraph plots
for (tissue in names(results_list_CA1_bothSex)) {
  
  print(paste('Processing tissue:', tissue, sep = ' '))
  #Extract the relevant resultstable
  results_list_tissue <- results_list_CA1_bothSex[[tissue]]
  #We'll have to set up a list onto which we'll attach the extracted number of DEGs to build the bargraph plots
  barplot_list <- list()
  #We'll repeat this step for each p value cutoff
  for (padj_cutoff in c(0.001,0.01,0.05,0.1)) {
    barplot_frame <- data.frame()
    comparison <- 'correlation'
    #We extract the significantly postiively and negatively correlated genes seperately and store their total number in the bargraph dataframe
    results_list_tissue_temp <- as.data.frame(results_list_tissue$resall)
    down_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear < -1*cor_cutoff))
    up_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear > 1*cor_cutoff))
    barplot_frame["down",comparison] <- data.frame(down_regulated)
    barplot_frame["up",comparison] <- data.frame(up_regulated)
    #After having the respective information assembled, we'sll store the outcome
    barplot_list[[as.character(padj_cutoff)]] <- barplot_frame
    
  }
  #Having assemled the relevant information across a range of pvalue cutoffs, we can now aggregate the results
  barplot_list <- bind_rows(barplot_list,.id = "id")
  #We'll use simple labels of to define the direction. Since we strictly followed the design of getting first down/up-regulated genes, we'll just apply the rep function
  barplot_list$Type <- rep(c('down','up'), length(padj_range))
  #Time to 'melt' the dataframe into the long format
  plot_df <- melt(barplot_list, id.vars = c('id','Type'))
  #Relabel columns
  colnames(plot_df) <- c("padj_cutoff", 'Type',"Comparison", "Numb_of_genes")
  
  plot_df$tissue <- tissue
  #Having gathered and prepared all the relevant information from this tissue/region, we're now appending that to the dataframe created previously
  plot_df_crossTissue <- rbind(plot_df_crossTissue, plot_df)
  
  
}

plot_df_crossTissue <- plot_df_crossTissue[-1,] #Remove the first row we needed to set up the dataframe

#we'll also have to set the order of the factors in the 'Type' column that contains informatino about  up-/down-regualtion
plot_df_crossTissue$Type <- factor(plot_df_crossTissue$Type, levels = c("up","down"), ordered = T)

#Now order the tissues according to their total number of age-correlated genes
plot_df_crossTissue$tissue <- factor(plot_df_crossTissue$tissue, levels = rev(as.character((plot_df_crossTissue  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
                                     ordered = T)
#For the resulting plot, we'll focus on the genes passing the padj cutoff of 0.05 
plot_df_crossTissue <- plot_df_crossTissue %>% filter(padj_cutoff == 0.05)

# #####Barplot of total # age correlated genes
# pdf(paste0(outdir,"Plots/Atlas_corrgene_barplot_alltissue_TPMcutoff_240714.pdf"), width = 15, height = 10)
# ggplot(plot_df_crossTissue, aes(x=tissue, y=Numb_of_genes, fill=Type), color=NA) + geom_bar(stat = "identity") +theme_classic()
# dev.off()

### Fig 1D ###
#####Barplot of corr genes as % expressed genes
num_Expressed <- read.csv(paste0(outdir,"Tables/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv")) #made in correlation_checkconfounds.R
adjust_barplot_df <- (merge(plot_df_crossTissue, num_Expressed))
adjust_barplot_df$Numb_of_genes_perc <- (adjust_barplot_df$Numb_of_genes / adjust_barplot_df$total.expressed) * 100

tissue_levels <- c('Muscle', 'Eye', 'Skin', 'SpinalCord', 'Fat', 'Brain', 'Heart', 'Spleen', 'Kidney', 'Liver', 'Gut', 'Gonad', 'Bone')
adjust_barplot_df$tissue <- factor(adjust_barplot_df$tissue, levels = tissue_levels)

pdf(paste0(outdir,"Plots/Atlas_ALLcorrgenes_percofexpressed_barplot_alltissue_Gonadcombo_TPMcutoff_240822.pdf"), width = 15, height = 10)
ggplot(adjust_barplot_df, aes(x=tissue, y=Numb_of_genes_perc, fill = Type), color=NA) + geom_bar(stat = "identity") +theme_classic() + 
  ggtitle('All Correlated Genes (Abs(Spearman R > 0.5))') + 
  ylab('Percent of Genes Expressed in Tissue') +
  ylim(c(0,20))
dev.off()

# The following are helpful controls, but aren't in the figures, scroll down for more code
# library(psych)
# #Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
# #colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#93cca8','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#ab5673', '#f1a8a4','#ef9ac2' ) 
# 
# 
# #Muscle, Eye, Skin, SpinalCord, Fat, Brain, Heart, Spleen, Kidney, Liver, Gut,  Bone, Ovary
# colTPS <- c('#f4c489','#b8b8c0','#ab5673', '#f1a8a4','#eee09b','#00a550', '#f0932e','#ef9ac2','#fcd328', '#6cc0ee', '#010101', '#93cca8','#d8413f') 
# 
# #v1 - tissues colored, up vs down is diff shape
# ggplot(adjust_barplot_df, aes(x=Numb_of_genes, y=total.expressed, color = tissue)) + 
#   geom_point(size = 5, aes(shape = Type)) +
#   scale_color_manual(values = colTPS) +
#   ggtitle('Age-associated genes vs total num genes expressed in tissue') +
#   scale_y_continuous(breaks = seq(0,22000, by = 2000), limits = c(0,22000)) +
#   scale_x_continuous(limits = c(0,1500)) +
#   theme_classic()
# 
# #v2 - tissues shape, up vs down is diff color
# ggplot(adjust_barplot_df, aes(x=Numb_of_genes, y=total.expressed, color = Type)) + 
#   geom_point(size = 5, aes(shape = tissue)) +
#   ggtitle('Age-associated genes vs total num genes expressed in tissue') +
#   scale_y_continuous(breaks = seq(0,22000, by = 2000), limits = c(0,22000)) +
#   scale_x_continuous(limits = c(0,1500)) +
#   theme_classic()
# 
# adjust_barplot_df.up <- adjust_barplot_df %>% filter(Type == 'up') 
# adjust_barplot_df.up$tissue <- factor(adjust_barplot_df.up$tissue, levels = rev(as.character((plot_df_crossTissue  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
#                                       ordered = T)
# 
# adjust_barplot_df.down <- adjust_barplot_df %>% filter(Type == 'down')
# adjust_barplot_df.down$tissue <- factor(adjust_barplot_df.down$tissue, levels = rev(as.character((plot_df_crossTissue  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
#                                         ordered = T)
# #Eye, Muscle, Skin, SpinalCord, Brain, Fat, Heart, Kidney, Spleen, Testis, Gut, Liver, Bone, Ovary
# colTPS <- c('#b8b8c0','#f4c489','#ab5673', '#f1a8a4','#00a550', '#eee09b', '#f0932e','#fcd328',  '#ef9ac2', '#93cca8', '#010101', '#6cc0ee', '#d8413f', '#93cca8') 
# 
# pdf(paste0(outdir,"Plots/Alltissue_POScorrgeneNumber_vs_TotalGenesExpressed_TPMcutoff_240822.pdf"), width = 15, height = 10)
# ggplot(adjust_barplot_df.up, aes(x=Numb_of_genes, y=total.expressed, color = tissue)) + 
#   geom_point(size = 7) +
#   scale_color_manual(values = colTPS) +
#   ggtitle('Positively associated genes vs total num genes expressed in tissue') +
#   scale_y_continuous(breaks = seq(0,22000, by = 2000), limits = c(0,22000)) +
#   scale_x_continuous(limits = c(0,1500)) +
#   theme_classic()
# dev.off()
# p <-corr.test(adjust_barplot_df.up$Numb_of_genes,  adjust_barplot_df.up$total.expressed) # Pearson
# print(p, short = F)
# 
# pdf(paste0(outdir,"Plots/Alltissue_NEGcorrgeneNumber_vs_TotalGenesExpressed_TPMcutoff_240822.pdf"), width = 15, height = 10)
# ggplot(adjust_barplot_df.down, aes(x=Numb_of_genes, y=total.expressed, color = tissue)) + 
#   geom_point(size = 7) +
#   scale_color_manual(values = colTPS) +
#   ggtitle('Negatively associated genes vs total num genes expressed in tissue') +
#   scale_y_continuous(breaks = seq(0,22000, by = 2000), limits = c(0,22000)) +
#   scale_x_continuous(limits = c(0,1500)) +
#   theme_classic()
# dev.off()
# corr.test(adjust_barplot_df.down$Numb_of_genes,  adjust_barplot_df.up$total.expressed) # Pearson

#---------------- Upset plot of all corr (>0.5) and anticorr (<-0.5) genes in each tissue ----------------

####Make sure to run UpSetR_helperfunctions_EC.R first, there are some essential functions there

upset.input.pos <- list()
upset.input.neg <- list()

#tissue = everything_but(names(dds_CA1_list), 'All')[1] #debugging
for (tissue in names(BulkSeq_Agecorrelation_results_list)) {
  print(tissue)
  table.topcorr <- BulkSeq_Agecorrelation_results_list[[tissue]]$topCor 
  table.topcorr.pos <- table.topcorr %>% filter(cor_spear > 0)
  tissue.pos.geneset <- table.topcorr.pos$gene
  upset.input.pos[[tissue]] = tissue.pos.geneset 
  table.topcorr.neg <- table.topcorr %>% filter(cor_spear < 0)
  tissue.neg.geneset <- table.topcorr.neg$gene
  upset.input.neg[[tissue]] = tissue.neg.geneset
}

##save the positively correlated and negatively correlated lists
pos.geneset = plyr::ldply(upset.input.pos, rbind)
rownames(pos.geneset) = pos.geneset$.id
pos.geneset$.id = NULL
pos.geneset = t(pos.geneset)

neg.geneset = plyr::ldply(upset.input.neg, rbind)
rownames(neg.geneset) = neg.geneset$.id
neg.geneset$.id = NULL
neg.geneset = t(neg.geneset)

# save your files
write.csv(pos.geneset, paste0(outdir,"Tables/240714_TPMcutoff_alltissues_positivelycorrelatedgenes.csv"))
write.csv(neg.geneset, paste0(outdir,"Tables/240714_TPMcutoff_alltissues_negativelycorrelatedgenes.csv"))

#----------------Correlated genes ----------------
upset(fromList(upset.input.pos), order.by = "freq",nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")

#pdf(paste0(outdir,"Plots/Atlas_POScorrgene_upsetplot_alltissue_TPMcutoff_240714.pdf"), width = 25, height = 5)
upset(fromList(upset.input.pos), order.by = c("degree", "freq"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
#dev.off()

intersect.pos.all <- fromList2(upset.input.pos)
intersect.pos.all <- intersect.pos.all %>% mutate(tissue.overlap = rowSums(.)) 
write.csv(intersect.pos.all, file = paste0(outdir,'Tables/Atlas_POScorrgene_intersections_alltissue_TPMcutoff_240714.csv'))

# subset genes that are shared among 6 or more tissues
intersect.pos.all.high <- intersect.pos.all %>% filter(tissue.overlap >= 6)
write.csv(intersect.pos.all.high, file = paste0(outdir,'Tables/Atlas_POScorrgene_intersections_sharedby6ormore_alltissue_TPMcutoff_240714.csv'))

intersect.pos.all.high$tissue.overlap <- NULL
#pdf(paste0(outdir,"Plots/Atlas_POScorrgene_upsetplot_sharedby6ormore_alltissue_TPMcutoff_240714.pdf"), width = 15, height = 10)
upset(intersect.pos.all.high, order.by = c("degree"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
#dev.off()

#----------------Anticorrelated genes ----------------
upset(fromList(upset.input.neg), order.by = "freq",nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")

#pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_upsetplot_alltissue_TPMcutoff_240714.pdf"), width = 25, height = 5)
upset(fromList(upset.input.neg), order.by = c("degree", "freq"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
#dev.off()

intersect.neg.all <- fromList2(upset.input.neg)
intersect.neg.all <- intersect.neg.all %>% mutate(tissue.overlap = rowSums(.))
#write.csv(intersect.neg.all, file = paste0(outdir,'Tables/Atlas_NEGcorrgene_intersections_alltissue_TPMcutoff_240714.csv'))

intersect.neg.all.high <- intersect.neg.all %>% filter(tissue.overlap >= 6)
#write.csv(intersect.neg.all.high, file = paste0(outdir,'Tables/Atlas_NEGcorrgene_intersections_sharedby6ormore_alltissue_TPMcutoff_240714.csv'))

intersect.neg.all.high$tissue.overlap <- NULL
#pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_upsetplot_sharedby6ormore_alltissue_TPMcutoff_240714.pdf"), width = 15, height = 10)
upset(intersect.neg.all.high, order.by = c("degree"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
#dev.off()


#---------------- Small Heatmap of correlation values for overlapping top correlated and anti-correlated genes in each tissue - Spearman Correlation ----------------

### Fig 2A ###
#intersect.pos.all.high
corr.heatmap.pos <- intersect.pos.all.high
corr.heatmap.pos$gene <- rownames(corr.heatmap.pos)

for(i in 1:length(rownames(corr.heatmap.pos))){
  gene <- rownames(corr.heatmap.pos)[i]
  for(tissue in colnames(corr.heatmap.pos)[1:14]){
    tissue.corres <- as.data.frame(BulkSeq_Agecorrelation_results_list[[tissue]]$resall)
    rownames(tissue.corres) <- tissue.corres$gene
    gene.corval <- tissue.corres[gene,]$cor_spear
    corr.heatmap.pos[gene,][[tissue]] <- gene.corval
  }
}

corr.heatmap.pos$gene <- NULL
corr.heatmap.pos$tissue.overlap <- NULL
#corr.heatmap.pos$Ovary <- NULL

col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
corr.heatmap.pos <- as.matrix(corr.heatmap.pos)

pdf(paste0(outdir,"Plots/Atlas_POScorrgene_heatmapSpear_sharedby6ormore_alltissue_TPMcutoff_corlabels_240903.pdf"), width = 3, height = 3)
Heatmap(corr.heatmap.pos, col = col_fun, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5), cell_fun = function(j, i, x, y, w, h, col) {
  if(!is.na(corr.heatmap.pos[i,j])){if(corr.heatmap.pos[i, j] > 0.5) {grid.text(round(corr.heatmap.pos[i, j], digits = 2), x, y, gp = gpar(fontsize = 4))}}  
})
dev.off()

#intersect.neg.all.high
corr.heatmap.neg <- intersect.neg.all.high
corr.heatmap.neg$gene <- rownames(corr.heatmap.neg)

for(i in 1:length(rownames(corr.heatmap.neg))){
  gene <- rownames(corr.heatmap.neg)[i]
  for(tissue in colnames(corr.heatmap.neg)[1:14]){
    tissue.corres <- as.data.frame(BulkSeq_Agecorrelation_results_list[[tissue]]$resall)
    rownames(tissue.corres) <- tissue.corres$gene
    gene.corval <- tissue.corres[gene,]$cor_spear
    corr.heatmap.neg[gene,][[tissue]] <- gene.corval
  }
}

corr.heatmap.neg$gene <- NULL
corr.heatmap.neg$tissue.overlap <- NULL
#corr.heatmap.neg$Ovary <- NULL

corr.heatmap.neg <- as.matrix(corr.heatmap.neg)

pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_heatmapSpear_sharedby6ormore_alltissue_TPMcutoff_corlabels_240903.pdf"), width = 3, height = 3)
Heatmap(corr.heatmap.neg, col = col_fun, show_heatmap_legend = F, row_names_gp = gpar(fontsize = 5), column_names_gp = gpar(fontsize = 5), cell_fun = function(j, i, x, y, w, h, col) {
  if(!is.na(corr.heatmap.neg[i,j])){if(corr.heatmap.neg[i, j] < -0.5) {grid.text(round(corr.heatmap.neg[i, j], digits = 2), x, y, gp = gpar(fontsize = 4))}}  
})
dev.off()



