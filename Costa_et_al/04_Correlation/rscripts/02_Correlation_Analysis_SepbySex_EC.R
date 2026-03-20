

##Here, we perform will correlate expression for each gene in each tissue/region with age and quantify the number of significant age-correlated genes
library(gprofiler2)
library(ComplexHeatmap)
library(UpSetR)
library(hrbrthemes)
library(dplyr)
library("RColorBrewer")


setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/')


#Sex: male, female 
colsex <-c('#068ec9','#ba1e2d')

#Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#7961a2','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#ab5673', '#f1a8a4','#ef9ac2' ) 


#First, we load our Deseq2 object as it contains all the counts and metadata
object.indir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

dds_CA1_list <- dds_TPS_list
#------------------------------------------------------------------
# Age-corr genes from TPM filtering method
# ------------------------------------------------------------------
#load merged tpm matrix
tpm <- read.csv('path/to/TPM_Atlas_allbatches_merged_v3.csv', row.names = 1)

#We'll setup a results list in which we'll store the output of the following analysis
BulkSeq_Agecorrelation_results_list <- list()

#We will rotate over each tissue/region (except the 'All') object
#tissue = 'Fat' #debugging
for (tissue in everything_but(names(dds_CA1_list), c('All', 'Gonad', 'Eye', 'Bone'))) {
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]

  #To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
  #we consider as 'expressed' - in this case these will be genes within the tissue that have TPM > 0.5 for 80% of the samples
  #find sample names for tissue
  coldata_tissue <- dds_tissue@colData
  coldata_tissue <- as.data.frame(coldata_tissue)
  
  #We use the DESEQ2-normalized data, which we do by getting the counts them from the DESeq2 object
  count_table <- counts(dds_tissue, normalized=T)
  count_table <- as.data.frame(count_table)
  
  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples 
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(coldata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue)) 
  print(paste0(tissue," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]

  ExpresedGenes <- rownames(expressed)

  for(s in c("M", "F")){
    #subset dds and ncount table for the tissue by sex
    coldata_tissue_sex <- coldata_tissue %>% filter(sex == s)
    count_table_sex <- count_table[,rownames(coldata_tissue_sex)]
  
    #We define the numeric values we want to correlate against - in this case, the age
    Score_to_compare <- as.numeric(as.character(coldata_tissue_sex[,'age_days']))
    
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
      cor_frame <- data.frame(expr=t(count_table_sex[gene,]), concentration= Score_to_compare)
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
    #Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibiting an absolute correlation value of 0.5)
    print(paste('TopCorGenes: ', nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
    nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05))
    #Store data the results in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['resall']] <- tissue_correlation_fullTable
    #Filter and store data for high correlates in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['topCor']]  <- (tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))
  }
}

#Gonad, Eye, Bone alone because missing not enough samples for 20 sample cutoff
for (tissue in c('Gonad', 'Eye', 'Bone')) {
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]
  
  #To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
  #we consider as 'expressed' - in this case these will be genes within the tissue that have TPM > 0.5 for 80% of the samples
  #find sample names for tissue
  coldata_tissue <- dds_tissue@colData
  coldata_tissue <- as.data.frame(coldata_tissue)
  
  #We use the DESEQ2-normalized data, which we do by getting the counts them from the DESeq2 object
  count_table <- counts(dds_tissue, normalized=T)
  count_table <- as.data.frame(count_table)
  
  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples 
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(coldata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue)) 
  print(paste0(tissue," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]
  
  ExpresedGenes <- rownames(expressed)
  
  for(s in c("M", "F")){
    #subset dds and ncount table for the tissue by sex
    coldata_tissue_sex <- coldata_tissue %>% filter(sex == s)
    count_table_sex <- count_table[,rownames(coldata_tissue_sex)]
    
    #We define the numeric values we want to correlate against - in this case, the age
    Score_to_compare <- as.numeric(as.character(coldata_tissue_sex[,'age_days']))
    
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
      cor_frame <- data.frame(expr=t(count_table_sex[gene,]), concentration= Score_to_compare)
      colnames(cor_frame) <- c('expr', 'concentration')
      #We wil only retain samples where the count is not 0
      cor_frame <- cor_frame %>% filter(expr > 0)
      #After filtering for non-expressing samples, we will only continue the analysis if we have at least 20 samples letf - otherwise we end up
      #with very spurious correlations
      if (nrow(cor_frame) >= 10) {
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
    #Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibiting an absolute correlation value of 0.5)
    print(paste('TopCorGenes: ', nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
    nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05))
    #Store data the results in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['resall']] <- tissue_correlation_fullTable
    #Filter and store data for high correlates in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['topCor']]  <- (tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))
  }
}

for (tissue in c('Bone')) {
  print(tissue)
  #We extract the respective Deseq2 object
  dds_tissue <- dds_CA1_list[[tissue]]
  
  #To make sure that we are not correlating flat noise with flat noise, we will focus on genes that
  #we consider as 'expressed' - in this case these will be genes within the tissue that have TPM > 0.5 for 80% of the samples
  #find sample names for tissue
  coldata_tissue <- dds_tissue@colData
  coldata_tissue <- as.data.frame(coldata_tissue)
  
  #We use the DESEQ2-normalized data, which we do by getting the counts them from the DESeq2 object
  count_table <- counts(dds_tissue, normalized=T)
  count_table <- as.data.frame(count_table)
  
  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples 
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(coldata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue)) 
  print(paste0(tissue," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]
  
  ExpresedGenes <- rownames(expressed)
  
  for(s in c("M", "F")){
    #subset dds and ncount table for the tissue by sex
    coldata_tissue_sex <- coldata_tissue %>% filter(sex == s)
    count_table_sex <- count_table[,rownames(coldata_tissue_sex)]
    
    #We define the numeric values we want to correlate against - in this case, the age
    Score_to_compare <- as.numeric(as.character(coldata_tissue_sex[,'age_days']))
    
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
      cor_frame <- data.frame(expr=t(count_table_sex[gene,]), concentration= Score_to_compare)
      colnames(cor_frame) <- c('expr', 'concentration')
      #We wil only retain samples where the count is not 0
      cor_frame <- cor_frame %>% filter(expr > 0)
      #After filtering for non-expressing samples, we will only continue the analysis if we have at least 20 samples letf - otherwise we end up
      #with very spurious correlations
      if (nrow(cor_frame) >= 10) {
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
    #Inspect the results and select the top correlating genes (defined by passing the statistical testing and exhibiting an absolute correlation value of 0.5)
    print(paste('TopCorGenes: ', nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))))
    nrow(tissue_correlation_fullTable %>% filter(padj_spear < 0.05))
    #Store data the results in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['resall']] <- tissue_correlation_fullTable
    #Filter and store data for high correlates in our correlation results list
    BulkSeq_Agecorrelation_results_list[[tissue]][[s]][['topCor']]  <- (tissue_correlation_fullTable %>% filter(padj_spear < 0.05)  %>% filter(abs(cor_spear) > .5))
  }
}

####
#We will save the results list for later use
save(BulkSeq_Agecorrelation_results_list, file='Robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_bysex_TPMcutoff_240714.bin')
####

# NOTE:
# geneset of ExpressedGenes will be the same as the sex combined tissue-level analysis

#------------------------------------------------------------------
# TPM filtering / Further analysis of age-correlated genes
# ------------------------------------------------------------------
outdir = "Output/"

load(paste0(object.indir, 'dds_BulkSeq_Aging_CorrelationResults_allsamples_bysex_TPMcutoff_240714.bin'))

# Skip the following section commented out - does make into figs but useful, scroll down further
#---------------- Bar plot of total # age corr genes ----------------
# #Now that we've run the correlation analysis we want to plot how many genes correlated significantly with age in each tissue
# 
# #We will set a range of pvalue cutoffs, to see if the results remain stable 
# padj_range <-c(0.001, 0.01, 0.05, 0.1)
# #We have to set a correlation coefficient cutoff (absolute cutoff of 0.5)
# cor_cutoff <- 0.5
# 
# 
# #We will have to establish a 'master data frame' with some dummy data. We will iteratively add data to this data frame
# plot_df_crossTissue_M <- data.frame(padj_cutoff=0.05, Type='A', Comparison='A', Numb_of_genes=1, tissue='A')
# plot_df_crossTissue_F <- data.frame(padj_cutoff=0.05, Type='A', Comparison='A', Numb_of_genes=1, tissue='A')
# 
# 
# #We will iterate over each of the tissues, get the correlation results to build the bargraph plots
# for (tissue in names(BulkSeq_Agecorrelation_results_list)) {
#   
#   print(paste('Processing tissue:', tissue, sep = ' '))
#   for(s in c('M', 'F')){
#     #Extract the relevant resultstable
#     results_list_tissue <- BulkSeq_Agecorrelation_results_list[[tissue]][[s]]
#     barplot_list <- list()
#     #We'll repeat this step for each p value cutoff
#     for (padj_cutoff in c(0.001,0.01,0.05,0.1)) {
#       barplot_frame <- data.frame()
#       comparison <- 'correlation'
#       #We extract the significantly postiively and negatively correlated genes seperately and store their total number in the bargraph dataframe
#       results_list_tissue_temp <- as.data.frame(results_list_tissue$resall)
#       down_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear < -1*cor_cutoff))
#       up_regulated <- nrow(results_list_tissue_temp %>% filter(padj_spear < padj_cutoff) %>% filter(cor_spear > 1*cor_cutoff))
#       barplot_frame["down",comparison] <- data.frame(down_regulated)
#       barplot_frame["up",comparison] <- data.frame(up_regulated)
#       #After having the respective information assembled, we'sll store the outcome
#       barplot_list[[as.character(padj_cutoff)]] <- barplot_frame
#       
#     }
#     #Having assemled the relevant information across a range of pvalue cutoffs, we can now aggregate the results
#     barplot_list <- bind_rows(barplot_list,.id = "id")
#     #We'll use simple labels of to define the direction. Since we strictly followed the design of getting first down/up-regulated genes, we'll just apply the rep function
#     barplot_list$Type <- rep(c('down','up'), length(padj_range))
#     #Time to 'melt' the dataframe into the long format
#     plot_df <- melt(barplot_list, id.vars = c('id','Type'))
#     #Relabel columns
#     colnames(plot_df) <- c("padj_cutoff", 'Type',"Comparison", "Numb_of_genes")
#     
#     plot_df$tissue <- tissue
#     if(s == 'M'){
#       #Having gathered and prepared all the relevant information from this tissue/region, we're now appending that to the dataframe created previously
#       plot_df_crossTissue_M <- rbind(plot_df_crossTissue_M, plot_df)
#     }
#     else{
#       #Having gathered and prepared all the relevant information from this tissue/region, we're now appending that to the dataframe created previously
#       plot_df_crossTissue_F <- rbind(plot_df_crossTissue_F, plot_df)
#     }
#   }
#   
# }
# 
# plot_df_crossTissue_M <- plot_df_crossTissue_M[-1,] #Remove the first row we needed to set up the dataframe
# plot_df_crossTissue_F <- plot_df_crossTissue_F[-1,] #Remove the first row we needed to set up the dataframe
# 
# #we'll also have to set the order of the factors in the 'Type' column that contains informatino about  up-/down-regualtion
# plot_df_crossTissue_M$Type <- factor(plot_df_crossTissue_M$Type, levels = c("up","down"), ordered = T)
# plot_df_crossTissue_F$Type <- factor(plot_df_crossTissue_F$Type, levels = c("up","down"), ordered = T)
# 
# #Now order the tissues according to their total number of age-correlated genes
# plot_df_crossTissue_M$tissue <- factor(plot_df_crossTissue_M$tissue, levels = rev(as.character((plot_df_crossTissue_M  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
#                                      ordered = T)
# plot_df_crossTissue_F$tissue <- factor(plot_df_crossTissue_F$tissue, levels = rev(as.character((plot_df_crossTissue_F  %>% dplyr::filter(padj_cutoff == 0.05) %>% dplyr::group_by(tissue, padj_cutoff)  %>% dplyr::summarise(total=sum(Numb_of_genes)) %>% dplyr::arrange(total))$tissue)),
#                                      ordered = T)
# 
# #For the resulting plot, we'll focus on the genes passing the padj cutoff of 0.05 
# plot_df_crossTissue_M <- plot_df_crossTissue_M %>% filter(padj_cutoff == 0.05)
# plot_df_crossTissue_F <- plot_df_crossTissue_F %>% filter(padj_cutoff == 0.05)
# 
# 
# #####Barplot of total # age correlated genes
# pdf(paste0(outdir,"Plots/Atlas_corrgene_barplot_alltissue_MaleOnly_TPMcutoff_240712.pdf"), width = 15, height = 10)
# ggplot(plot_df_crossTissue_M, aes(x=tissue, y=Numb_of_genes, fill=Type), color=NA) + geom_bar(stat = "identity") +theme_classic()
# dev.off()
# 
# pdf(paste0(outdir,"Plots/Atlas_corrgene_barplot_alltissue_FemaleOnly_TPMcutoff_240712.pdf"), width = 15, height = 10)
# ggplot(plot_df_crossTissue_F, aes(x=tissue, y=Numb_of_genes, fill=Type), color=NA) + geom_bar(stat = "identity") +theme_classic()
# dev.off()
# 
# 
# #####Barplot of corr genes as % expressed genes
# num_Expressed <- read.csv("/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Correlation/Output/Tables/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv")
# 
# adjust_barplot_df_M <- (merge(plot_df_crossTissue_M, num_Expressed))
# adjust_barplot_df_M$Numb_of_genes_perc <- (adjust_barplot_df_M$Numb_of_genes / adjust_barplot_df_M$total.expressed) * 100
# adjust_barplot_df_M <- adjust_barplot_df_M %>% filter(!tissue == 'Eye')
# 
# tissue_levels <- c('Skin','Muscle','SpinalCord', 'Heart', 'Kidney', 'Fat', 'Brain', 'Spleen', 'Gonad', 'Liver', 'Gut', 'Bone')
# adjust_barplot_df_M$tissue <- factor(adjust_barplot_df_M$tissue, levels = tissue_levels)
# pdf(paste0(outdir,"Plots/Atlas_ALLcorrgenes_percofexpressed_barplot_alltissue_MaleOnly_TPMcutoff_240714.pdf"), width = 15, height = 10)
# ggplot(adjust_barplot_df_M, aes(x=tissue, y=Numb_of_genes_perc, fill = Type), color=NA) + geom_bar(stat = "identity") +theme_classic() + 
#   ggtitle('All Correlated Genes (Spearman R > 0.5)') + 
#   ylab('Percent of Expressed Genes') +
#   ylim(c(0,20))
# dev.off()
# 
# adjust_barplot_df_F <- (merge(plot_df_crossTissue_F, num_Expressed))
# adjust_barplot_df_F$Numb_of_genes_perc <- (adjust_barplot_df_F$Numb_of_genes / adjust_barplot_df_F$total.expressed) * 100
# adjust_barplot_df_F <- adjust_barplot_df_F %>% filter(!tissue == 'Eye')
# 
# tissue_levels <- c('Muscle','Fat','SpinalCord', 'Skin','Brain', 'Kidney', 'Liver', 'Gonad', 'Spleen','Gut', 'Heart', 'Bone')
# adjust_barplot_df_F$tissue <- factor(adjust_barplot_df_F$tissue, levels = tissue_levels)
# pdf(paste0(outdir,"Plots/Atlas_ALLcorrgenes_percofexpressed_barplot_alltissue_FemaleOnly_TPMcutoff_240714.pdf"), width = 15, height = 10)
# ggplot(adjust_barplot_df_F, aes(x=tissue, y=Numb_of_genes_perc, fill = Type), color=NA) + geom_bar(stat = "identity") +theme_classic() + 
#   ggtitle('All Correlated Genes (Spearman R > 0.5)') + 
#   ylab('Percent of Expressed Genes') +
#   ylim(c(0,20))
# dev.off()


#---------------- Upset plot of all corr (>0.5) and anticorr (<-0.5) genes in each tissue ----------------
sex = 'M' # to be done for males and females
#Make sure to run UpSetR_helperfunctions_EC.R first
upset.input.pos <- list()
upset.input.neg <- list()

#tissue = everything_but(names(dds_CA1_list), 'All')[1] #debugging
for (tissue in names(BulkSeq_Agecorrelation_results_list)) {
  print(tissue)
  table.topcorr <- BulkSeq_Agecorrelation_results_list[[tissue]][[sex]]$topCor 
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

write.csv(pos.geneset, paste0(outdir,"Tables/240714_TPMcutoff_alltissues_",sex,"Only_positivelycorrelatedgenes.csv"))
write.csv(neg.geneset, paste0(outdir,"Tables/240714_TPMcutoff_alltissues_",sex,"Only_negativelycorrelatedgenes.csv"))


##correlated
upset(fromList(upset.input.pos), order.by = "freq",nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")

pdf(paste0(outdir,"Plots/Atlas_POScorrgene_upsetplot_alltissue",sex,"Only_TPMcutoff_240714.pdf"), width = 25, height = 5)
upset(fromList(upset.input.pos), order.by = c("degree", "freq"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
dev.off()

intersect.pos.all <- fromList2(upset.input.pos)
intersect.pos.all <- intersect.pos.all %>% mutate(tissue.overlap = rowSums(.)) 
write.csv(intersect.pos.all, file = paste0(outdir,'Tables/Atlas_POScorrgene_intersections_alltissue_',sex,'Only_TPMcutoff_240714.csv'))

intersect.pos.all.high <- intersect.pos.all %>% filter(tissue.overlap >= 5)
write.csv(intersect.pos.all.high, file = paste0(outdir,'Tables/Atlas_POScorrgene_intersections_sharedby5ormore_alltissue',sex,'Only_TPMcutoff_240712.csv'))

###output for GO
write.csv(data.frame('id' = rownames(intersect.pos.all.high)), file = paste0(outdir,'Tables/POScorrgene_5plusoverlap_alltissues_',sex,'Only_TPMcutoff_240714.csv'))

intersect.pos.all.high$tissue.overlap <- NULL
pdf(paste0(outdir,"Plots/Atlas_POScorrgene_upsetplot_sharedby5ormore_alltissue",sex,"_TPMcutoff_240714.pdf"), width = 15, height = 10)
upset(intersect.pos.all.high, order.by = c("degree"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
dev.off()


##anticorrelated
upset(fromList(upset.input.neg), order.by = "freq",nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")

pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_upsetplot_alltissue_",sex,"TPMcutoff_240714.pdf"), width = 25, height = 5)
upset(fromList(upset.input.neg), order.by = c("degree", "freq"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
dev.off()

intersect.neg.all <- fromList2(upset.input.neg)
intersect.neg.all <- intersect.neg.all %>% mutate(tissue.overlap = rowSums(.))
write.csv(intersect.neg.all, file = paste0(outdir,'Tables/Atlas_NEGcorrgene_intersections_alltissue',sex,'Only_TPMcutoff_240714.csv'))

intersect.neg.all.high <- intersect.neg.all %>% filter(tissue.overlap >= 5)
write.csv(intersect.neg.all.high, file = paste0(outdir,'Tables/Atlas_NEGcorrgene_intersections_sharedby5ormore_alltissue',sex,'_TPMcutoff_240712.csv'))

###output for GO
write.csv(data.frame('id' = rownames(intersect.neg.all.high)), file = paste0(outdir,'Tables/NEGcorrgene_5plusoverlap_alltissues_',sex,'Only_TPMcutoff_240714.csv'))


intersect.neg.all.high$tissue.overlap <- NULL
pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_upsetplot_sharedby5ormore_alltissue",sex,"Only_TPMcutoff_240714.pdf"), width = 15, height = 10)
upset(intersect.neg.all.high, order.by = c("degree"),nintersects = NA,nsets = 13, number.angles = 30, mainbar.y.label = "# of Genes Intersecting", sets.x.label = "# of Age-Correlated (>0.5) Genes per Tissue")
dev.off()



#---------------- Small Heatmap of correlation values for overlapping top correlated and anti-correlated genes in each tissue - Spearman Correlation ----------------
# #intersect.pos.all.high
# corr.heatmap.pos <- intersect.pos.all.high
# corr.heatmap.pos$gene <- rownames(corr.heatmap.pos)
# corr.heatmap.pos$tissue.overlap <- NULL
# 
# for(i in 1:length(rownames(corr.heatmap.pos))){
#   gene <- rownames(corr.heatmap.pos)[i]
#   for(tissue in colnames(corr.heatmap.pos)[1:13]){
#     tissue.corres <- as.data.frame(BulkSeq_Agecorrelation_results_list[[tissue]][[sex]]$resall)
#     rownames(tissue.corres) <- tissue.corres$gene
#     gene.corval <- tissue.corres[gene,]$cor_spear
#     corr.heatmap.pos[gene,][[tissue]] <- gene.corval
#   }
# }
# 
# corr.heatmap.pos$gene <- NULL
# 
# 
# library("RColorBrewer")
# library("circlize")
# col_fun <- colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
# 
# pdf(paste0(outdir,"Plots/Atlas_POScorrgene_heatmapSpear_sharedby5ormore_alltissue",sex,"_TPMcutoff_240714.pdf"), width = 7, height = 7)
# Heatmap(corr.heatmap.pos, col = col_fun)
# dev.off()
# 
# #intersect.neg.all.high
# corr.heatmap.neg <- intersect.neg.all.high
# corr.heatmap.neg$gene <- rownames(corr.heatmap.neg)
# corr.heatmap.neg$tissue.overlap <- NULL
# 
# for(i in 1:length(rownames(corr.heatmap.neg))){
#   gene <- rownames(corr.heatmap.neg)[i]
#   for(tissue in colnames(corr.heatmap.neg)[1:13]){
#     tissue.corres <- as.data.frame(BulkSeq_Agecorrelation_results_list[[tissue]][[sex]]$resall)
#     rownames(tissue.corres) <- tissue.corres$gene
#     gene.corval <- tissue.corres[gene,]$cor_spear
#     corr.heatmap.neg[gene,][[tissue]] <- gene.corval
#   }
# }
# 
# corr.heatmap.neg$gene <- NULL
# 
# pdf(paste0(outdir,"Plots/Atlas_NEGcorrgene_heatmapSpear_sharedby5ormore_alltissue",sex,"_TPMcutoff_240714.pdf"), width = 7, height = 7)
# Heatmap(corr.heatmap.neg, col = col_fun)
# dev.off()

#------------------------------------------------------------------
# Upset Plot - Kidney GO enrichment input
# ------------------------------------------------------------------
### Fig 4C: input for GO enrichment ###

# Get the genes for input into GO enrichment 
sex <- "M" # Males
kidney.topcorr <- BulkSeq_Agecorrelation_results_list$Kidney[[sex]]$topCor
kidney.topcorr.up <- kidney.topcorr %>% filter(cor_spear > 0.5)
kidney.topcorr.down <- kidney.topcorr %>% filter(cor_spear < -0.5)

###output for GO
write.csv(data.frame('id' = kidney.topcorr.up$gene), file = paste0(outdir,'Tables/POScorrgene_Kidney_sex',sex,'_Only_TPMcutoff_240715.csv'))
write.csv(data.frame('id' = kidney.topcorr.down$gene), file = paste0(outdir,'Tables/NEGcorrgene_Kidney_sex',sex,'_Only_TPMcutoff_240715.csv'))

sex <- "F" # Females
kidney.topcorr <- BulkSeq_Agecorrelation_results_list$Kidney[[sex]]$topCor
kidney.topcorr.up <- kidney.topcorr %>% filter(cor_spear > 0.5)
kidney.topcorr.down <- kidney.topcorr %>% filter(cor_spear < -0.5)

###output for GO
write.csv(data.frame('id' = kidney.topcorr.up$gene), file = paste0(outdir,'Tables/POScorrgene_Kidney_sex',sex,'_Only_TPMcutoff_240715.csv'))
write.csv(data.frame('id' = kidney.topcorr.down$gene), file = paste0(outdir,'Tables/NEGcorrgene_Kidney_sex',sex,'_Only_TPMcutoff_240715.csv'))

###universe
#For GO enrichment - need to get the universe of all genes
dds_tissue <- dds_CA1_list$Kidney

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

write.csv(data.frame('id' = rownames(expressed)), file = paste0(outdir,'Tables/GOenrich_universe_corrgene_Kidney_bothSexes_TPMcutoff_240715.csv'))

#------------------------------------------------------------------
# Bubble plot of Kidney Up and Down genes
# ------------------------------------------------------------------

### Fig 4B ###
genes.to.plot <- c('IRF4', 'RAG1', 'RAG2', 'PAX5', 'TCF7', 'CD247', 'CD79B', 'LYZ', 'CSF1R', 'MARCO', 'TPO', 'EPX', 'CD3E', 'CD74')
topcorr.kidney.M$sex <- 'M'
topcorr.kidney.F$sex <- 'F'
toplot <- rbind(topcorr.kidney.F, topcorr.kidney.M)
plotdata <- toplot %>% filter(Human %in% genes.to.plot)
ggplot(data=plotdata, aes(x=sex, y=gene)) + 
  geom_point(aes(color=cor_spear, size = -log10(padj_spear))) + 
  scale_colour_gradientn(colours=c('blue','white','red'), 
                         limits=c(-1,1)) +
  labs(title = 'Age-Association of Immune Genes in Kidney', x='Sex', y = '',size='-log10(padj)', color='Spearman_R') 


# #------------------------------------------------------------------
# # Volcano plot of Kidney Up and Down genes
# # ------------------------------------------------------------------
# library(EnhancedVolcano)
# topcorr.kidney.M <- as.data.frame(BulkSeq_Agecorrelation_results_list$Kidney$M$resall)
# #topcorr.kidney.M$abs.cor_spear <- abs(topcorr.kidney.M$cor_spear)
# #topcorr.kidney.M <- topcorr.kidney.M %>% arrange(desc(abs.cor_spear)) %>% top_n(30, abs.cor_spear)
# 
# topcorr.kidney.F <- as.data.frame(BulkSeq_Agecorrelation_results_list$Kidney$F$resall)
# #topcorr.kidney.F$abs.cor_spear <- abs(topcorr.kidney.F$cor_spear)
# #topcorr.kidney.F <- topcorr.kidney.F %>% arrange(desc(abs.cor_spear)) %>% top_n(30, abs.cor_spear)
# 
# rownames(topcorr.kidney.M) <- topcorr.kidney.M$gene
# rownames(topcorr.kidney.F) <- topcorr.kidney.F$gene
# 
# #####
# EnhancedVolcano(topcorr.kidney.M,
#                 lab = rownames(topcorr.kidney.M),
#                 x = 'cor_spear',
#                 y = 'padj_spear',
#                 selectLab = c('irf4', 'LOC107383908'),
#                 FCcutoff = 0.5,
#                 pCutoff = 10e-2) +
#   scale_x_continuous(breaks=seq(-1,1, 0.25))
# 
# EnhancedVolcano(topcorr.kidney.F,
#                 lab = rownames(topcorr.kidney.F),
#                 x = 'cor_spear',
#                 y = 'padj_spear',
#                 selectLab = c('rag1', 'rag2'),
#                 FCcutoff = 0.5,
#                 pCutoff = 10e-2) +
#   scale_x_continuous(breaks=seq(-1,1, 0.25))
# 
# #####
# library(stringr)
# # Human genes for volcano plot
# nfur.ortho <- read.csv(file = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Rscripts/Ortholog_Conversion/nfur-ncbi_orthologs_20170302_pps.csv')
# nfur.ortho <- as.data.frame(nfur.ortho)
# nfur.ortho$NCBI.Definition <- str_replace(nfur.ortho$NCBI.Definition,"\\s[\\(][A-Za-z0-9]*[\\)]" ,"")
# names(nfur.ortho)[1] <- "Killifish"
# 
# topcorr.kidney.M <- as.data.frame(topcorr.kidney.M)
# topcorr.kidney.M$Gene <- rownames(topcorr.kidney.M)
# topcorr.kidney.M$Human <- NA
# topcorr.kidney.M$Mouse <- NA
# topcorr.kidney.M$Zebrafish <- NA
# 
# 
# for(row in 1:nrow(topcorr.kidney.M)){
#   g <- topcorr.kidney.M[row,]$Gene
#   print(g)
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Human
#     topcorr.kidney.M[row,]$Human<- ifelse(conv != '', conv, NA)
#   }
# }
# 
# for(row in 1:nrow(topcorr.kidney.M)){
#   print(row)
#   g <- topcorr.kidney.M[row,]$Gene
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Mouse
#     topcorr.kidney.M[row,]$Mouse<- ifelse(conv != '', conv, NA)
#   }
# }
# 
# for(row in 1:nrow(topcorr.kidney.M)){
#   print(row)
#   g <- topcorr.kidney.M[row,]$Gene
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Zebrafish
#     topcorr.kidney.M[row,]$Zebrafish<- ifelse(conv != '', conv, NA)
#   }
# }
# 
# EnhancedVolcano(topcorr.kidney.M,
#                 lab = topcorr.kidney.M$Human,
#                 x = 'cor_spear',
#                 y = 'padj_spear',
#                 selectLab = c('IRF4', 'RAG1', 'RAG2'),
#                 FCcutoff = 0.5,
#                 pCutoff = 10e-2) +
#   scale_x_continuous(breaks=seq(-1,1, 0.25))
# 
# ######
# topcorr.kidney.F <- as.data.frame(topcorr.kidney.F)
# topcorr.kidney.F$Gene <- rownames(topcorr.kidney.F)
# topcorr.kidney.F$Human <- NA
# topcorr.kidney.F$Mouse <- NA
# topcorr.kidney.F$Zebrafish <- NA
# 
# 
# for(row in 1:nrow(topcorr.kidney.F)){
#   g <- topcorr.kidney.F[row,]$Gene
#   print(g)
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Human
#     topcorr.kidney.F[row,]$Human<- ifelse(conv != '', conv, NA)
#   }
# }
# 
# for(row in 1:nrow(topcorr.kidney.F)){
#   print(row)
#   g <- topcorr.kidney.F[row,]$Gene
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Mouse
#     topcorr.kidney.F[row,]$Mouse<- ifelse(conv != '', conv, NA)
#   }
# }
# 
# for(row in 1:nrow(topcorr.kidney.F)){
#   print(row)
#   g <- topcorr.kidney.F[row,]$Gene
#   if(g %in% nfur.ortho$Killifish){
#     conv <- nfur.ortho[nfur.ortho$Killifish == g,]$Zebrafish
#     topcorr.kidney.F[row,]$Zebrafish<- ifelse(conv != '', conv, NA)
#   }
# }


# 
# ##### Volcano
# pdf(paste0(outdir,"Plots/Atlas_Kidney_Males_VolcanoPlot_TPMcutoff_240716.pdf"), width = 10, height = 7)
# EnhancedVolcano(topcorr.kidney.M,
#                 lab = topcorr.kidney.M$Human,
#                 x = 'cor_spear',
#                 y = 'padj_spear',
#                 selectLab = c('IRF4', 'RAG1', 'RAG2', 'PAX5', 'TCF7', 'CD247', 'CD79B'),
#                 max.overlaps = Inf,
#                 FCcutoff = 0.5,
#                 pCutoff = 10e-2,
#                 pointSize = 2.0,
#                 labSize = 4,
#                 labCol = 'black',
#                 labFace = 'bold',
#                 #boxedLabels = TRUE,
#                 colAlpha = 4/5,
#                 legendPosition = 'right',
#                 legendLabSize = 14,
#                 legendIconSize = 4.0,
#                 drawConnectors = TRUE,
#                 widthConnectors = 1.0,
#                 colConnectors = 'black') +
#   scale_x_continuous(breaks=seq(-1,1, 0.25))
# dev.off()
# 
# pdf(paste0(outdir,"Plots/Atlas_Kidney_Females_VolcanoPlot_TPMcutoff_240716.pdf"), width = 10, height = 7)
# EnhancedVolcano(topcorr.kidney.F,
#                 lab = topcorr.kidney.F$Human,
#                 x = 'cor_spear',
#                 y = 'padj_spear',
#                 selectLab = c('IRF4', 'RAG1', 'RAG2', 'PAX5', 'TCF7','CD247', 'CD79B'),
#                 max.overlaps = Inf,
#                 FCcutoff = 0.5,
#                 pCutoff = 10e-2,
#                 pointSize = 2.0,
#                 labSize = 4,
#                 labCol = 'black',
#                 labFace = 'bold',
#                 #boxedLabels = TRUE,
#                 colAlpha = 4/5,
#                 legendPosition = 'right',
#                 legendLabSize = 14,
#                 legendIconSize = 4.0,
#                 drawConnectors = TRUE,
#                 widthConnectors = 1.0,
#                 colConnectors = 'black') +
#   scale_x_continuous(breaks=seq(-1,1, 0.25))
# dev.off()
