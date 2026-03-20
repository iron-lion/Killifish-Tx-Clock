rm(list=ls())
# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory
setwd("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/")

#Sex: male, female 
colsex <-c('#068ec9','#ba1e2d')

#Bone, Brain, Eye, Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 

outdir = 'Output/'


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
#load the merged raw count matrix (finalized in 05_checksampleidentity.R)
countdata <- read.csv("Input/Counts_Atlas_allbatches_merged_v3.csv", row.names = 1)

# Load the edited experiment design matrix (finalized in 05_checksampleidentity.R)
sampleTable <- read.csv("Input/ExperimentDesign_allbatches_combined_v7.csv", row.names = 1)  

coldata <- DataFrame(sampleTable)
col_order <- rownames(coldata)
countdata <- countdata[,col_order]

table(colnames(countdata) == rownames(coldata)) # Make sure that the column names are identical

# Add a pseudocount of 1 to all the Count values - use this for normalized
# https://www.biostars.org/p/9510236/ - use of pseudocounts - only used before log transformation of ncounts to avoid log of zeros - do not do for DE testing
countdata_1 = countdata + 1

#load merged tpm matrix
tpm <- read.csv('Input/TPM_Atlas_allbatches_merged_v3.csv', row.names = 1)

#when we run the analyses, we want the flexibility to run Gonad-wide and sex-specific analyses
coldata$tissue <- gsub("Ovary", "Gonad", coldata$tissue) 
coldata$tissue <- gsub("Testis", "Gonad", coldata$tissue) 

#------------------------------------------------------------------
# Initial analysis - all samples and batches together
# ------------------------------------------------------------------
# Make dds object 
dds_all_tissue <- DESeqDataSetFromMatrix(countData = countdata,
                                         colData = coldata,
                                         design = ~ age_days + tissue)

#----------------Filter out not expressed genes ----------------
#rationale: Although DESeq2 does its own independent filtering, filtering out genes with few reads can save you time 
#source: https://support.bioconductor.org/p/65256/ 

### by raw count - ok to do for all tissue samples 
#how many samples have at least 677 sum count (# of samples)
table(rowSums(counts(dds_all_tissue)) > 0) #15 genes
#View(data.frame(rowSums(counts(dds_all_tissue))))

# only keep rows that have greater than 0 sum count - this gets rid of genes that have 0 counts
dds_all_tissue <- dds_all_tissue [ rowSums(counts(dds_all_tissue)) > 0, ] #now we have 25,107 genes

#save normcounts
dds_all_tissue <- estimateSizeFactors(dds_all_tissue)
ncounts <- counts(dds_all_tissue, normalized = T)
write.csv(ncounts, file = paste0(outdir,'CountsNormDESeq2_AllTissue_240708.csv'))

#---------------- Plot PCA ----------------
# Plot PCA using vst normalization
vsd_all <- vst(dds_all_tissue)
head(assay(vsd_all), 3)
p <- pca(assay(vsd_all), metadata = colData(dds_all_tissue))


### Fig 1B ###
pdf("Output/Plots/PCA/Atlas_PCA_alltissue_PC1PC2_240610.pdf", width = 10, height = 10)
biplot(p,
       x = 'PC1', y = 'PC2',
       lab = NA,
       colby = 'tissue',
       colkey = colTPS,
       shape = 'sex',
       legendPosition = 'top')
dev.off()

pdf("Output/Plots/PCA/Atlas_PCA_alltissue_seqbatch_PC1PC2_240610.pdf", width = 10, height = 10)
biplot(p,
       x = 'PC1', y = 'PC2',
       lab = NA,
       colby = 'cDNA_batch',
       legendPosition = 'top')
dev.off()

# ------------------------------------------------------------------
# Create the master list of DESeq2 objects
# ------------------------------------------------------------------
#This will create a list of DESEQ objects onto which we will add iteratively objects for each tissue individually
dds_TPS_list <- list(All=dds_all_tissue)

tissue.list <- as.character(unique(dds_all_tissue$tissue))

colage <- c('#395ea4','#395ea4','#395ea4','#3789c0','#3789c0','#3789c0','#dfeaf0','#dfeaf0','#f3d4ac','#f3d4ac', '#f7b77b','#f7b77b','#f7b77b','#a92b46','#a92b46')

#colage.rev = rev(colage)

#you may need to make subfolders for plots

#tissue = tissue.list[1] #debugging
#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective tissue
for (tissue in tissue.list) {
  print(tissue)
  #subset for current tissue and make DESeq object
  dds_TPS_list[[tissue]] <- dds_all_tissue[, dds_all_tissue$tissue == tissue]
  dds_TPS_list[[tissue]]$tissue <- droplevels(dds_TPS_list[[tissue]]$tissue)
  dds_TPS_list[[tissue]]$sex <- as.factor(dds_TPS_list[[tissue]]$sex)
  dds_TPS_list[[tissue]]$age_days <- as.factor(dds_TPS_list[[tissue]]$age_days)
  
  #save normcounts
  dds_TPS_list[[tissue]]<- estimateSizeFactors(dds_TPS_list[[tissue]])
  ncounts <- counts(dds_TPS_list[[tissue]], normalized = T)
  write.csv(ncounts, file = paste0(outdir,'CountsNormDESeq2_',tissue,'_240708.csv'))
  
  #this is where you run tissue by tissue PCA
  vsd.tissue <- vst(dds_TPS_list[[tissue]])
  head(assay(vsd.tissue), 3)
  p <- pca(assay(vsd.tissue), metadata = colData(dds_TPS_list[[tissue]]))

  pdf(file = paste0(outdir,"Plots/PCA/Atlas_PCA_",tissue,"_byage-sex_PC1PC2_240711.pdf", sep = ""), width = 10, height = 10)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'age_days',
               colkey = colage.rev,
               shape = 'sex',
               legendPosition = 'top'))
  dev.off()
  
  
  pdf(file = paste0("Plots/PCA/Atlas_PCA_",tissue,"_bysex_PC1PC2_240711.pdf", sep = ""), width = 10, height = 10)
  print(biplot(p,
               x = 'PC1', y = 'PC2',
               lab = NA,
               colby = 'sex',
               legendPosition = 'top'))
  dev.off()
  
}

save(dds_TPS_list, file = paste0('Output/robjects/dds_TPS_allsamples_Gonadcombo_240714.bin'))



#---------------- Single tissue PCAs ----------------
object.indir = "/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

outdir = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2/'
colage <- c('#395ea4','#395ea4','#395ea4','#3789c0','#3789c0','#3789c0','#dfeaf0','#dfeaf0','#f3d4ac','#f3d4ac', '#f7b77b','#f7b77b','#f7b77b','#a92b46','#a92b46')


dds_CA1_list <- dds_TPS_list
rm(dds_CA1_list)

tissue = 'Kidney'

#subset for current tissue and make DESeq object
dds_TPS_list[[tissue]] <- dds_all_tissue[, dds_all_tissue$tissue == tissue]
dds_TPS_list[[tissue]]$tissue <- droplevels(dds_TPS_list[[tissue]]$tissue)
dds_TPS_list[[tissue]]$sex <- as.factor(dds_TPS_list[[tissue]]$sex)
dds_TPS_list[[tissue]]$age_days <- as.factor(dds_TPS_list[[tissue]]$age_days)

#this is where you run tissue by tissue PCA
vsd.tissue <- vst(dds_TPS_list[[tissue]])
head(assay(vsd.tissue), 3)
p <- pca(assay(vsd.tissue), metadata = colData(dds_TPS_list[[tissue]]))


### Fig 4A ###
pdf(file = paste0(outdir,"Plots/PCA/Atlas_PCA_",tissue,"_byage-sex_PC1PC2_240920.pdf", sep = ""), width = 4, height = 4)
print(biplot(p,
             x = 'PC1', y = 'PC2',
             lab = NA,
             colby = 'age_days',
             colkey = colage,
             shape = 'sex',
             pointSize = 2,
             legendPosition = 'none'))
dev.off()



#---------------- DESeq2 Differential Expression Across Regions, with Binning ----------------
object.indir = "Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

tissue.list <- names(dds_TPS_list)[2:14]

dds_CA1_list <- dds_TPS_list
padj_cutoff = 0.05

results_list_CA1 <- list()

for (tissue in tissue.list) {
  print(tissue)
  #subset for current tissue and run Deseq2
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('47', '49', '52'), '1', NA)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('75', '77', '78'), '2', dds_CA1_list[[tissue]]$age_bin)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('102', '103'), '3', dds_CA1_list[[tissue]]$age_bin)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_CA1_list[[tissue]]$age_bin)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_CA1_list[[tissue]]$age_bin)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('147', '152', '155'), '5', dds_CA1_list[[tissue]]$age_bin)
  dds_CA1_list[[tissue]]$age_bin <- ifelse(dds_CA1_list[[tissue]]$age_days %in% c('161', '162'), '6', dds_CA1_list[[tissue]]$age_bin)
  
  dds_CA1_list[[tissue]]$age_bin <- as.factor(dds_CA1_list[[tissue]]$age_bin)
  design(dds_CA1_list[[tissue]]) <- ~age_bin + sex 
  
  #Run Deseq2 and store in the Deseq2 object list
  dds_CA1_list[[tissue]] <- DESeq(dds_CA1_list[[tissue]], fitType = 'local')

  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_CA1_list[[tissue]]
  results_list_tissue <- list()

  #Create a list of all pairwise comparisons
  comparison_list <- t(combn(unique(colData(dds_tissue_temp)$age_bin),2))

  setwd('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2/')
  #Iterate over all pairwise comparisons and extract the results table, store it in the results list. then flip the conditions (so that results for 3 vs 21 and 21 vs 3 months is stored)
  start_time <- Sys.time()
  for (row in 1:nrow(comparison_list)) {
    print(row)
    #Get the conditions/timepoints to test
    cond1 <- as.character(comparison_list[row,1])
    cond2 <- as.character(comparison_list[row,2])

    folder <- paste(cond1, cond2, sep = "_vs_")
    print(folder)

    aspect_to_test <- 'age_bin'

    results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
    results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
    resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
    resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
    dim(resSig)
    #Store the output tables in the results_list
    results_list_tissue[[folder]]$resall <- results_temp
    results_list_tissue[[folder]]$resOrdered <- resOrdered
    results_list_tissue[[folder]]$ressig <- resSig


    ##Flip conditions and re-run the results extraction

    cond1 <- as.character(comparison_list[row,2])
    cond2 <- as.character(comparison_list[row,1])
    folder <- paste(cond1, cond2, sep = "_vs_")
    print(folder)

    dds_tissue_temp <- dds_tissue_temp

    aspect_to_test <- 'age_bin'

    results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
    results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
    resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
    resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
    dim(resSig)
    #Store the output tables in the results_list
    results_list_tissue[[folder]]$resall <- results_temp
    results_list_tissue[[folder]]$resOrdered <- resOrdered
    results_list_tissue[[folder]]$ressig <- resSig

  }
  end_time <- Sys.time()
  end_time-start_time
  
  #Store tissue results in the major results list, then save the output and start the loop for the next tissue
  results_list_CA1[[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_CA1_list, results_list_CA1,  file='Robjects/dds_DESeq2_TPS_allsamples_agesbinned_240805.bin')

}
