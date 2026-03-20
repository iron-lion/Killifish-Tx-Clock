# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
#RStudio Version 2024.04.1+748 (2024.04.1+748)
#R Version 4.3.3

library(DESeq2)
library(PCAtools)

# Set wd to the current directory
setwd("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/")

# Load the DESeq2 object list
object.indir = "/Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

outdir = paste0(getwd(), '/Output/Tables/')
outdir.robj = paste0(getwd(), '/robjects/')

#------------------------------------------------------------------
# Male
# -----------------------------------------------------------------
#---------------- Determine the n per group, with binning, but sep by sex ----------------


#---------------- Make list of DESeq2 objs, male only ----------------
tissue.list <- names(dds_TPS_list)[2:14]

dds_TPS_male_list <- list()
dds_all_tissue <- dds_TPS_list[['All']]

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective tissue + sex
for (tissue in tissue.list) {
  print(tissue)
  
  #subset for current tissue and make DESeq object
  dds_TPS_male_list[[tissue]] <- dds_all_tissue[, dds_all_tissue$tissue == tissue & dds_all_tissue$sex == 'M']
  dds_TPS_male_list[[tissue]]$tissue <- droplevels(dds_TPS_male_list[[tissue]]$tissue)
  dds_TPS_male_list[[tissue]]$sex <- as.factor(dds_TPS_male_list[[tissue]]$sex)
  dds_TPS_male_list[[tissue]]$age_days <- as.factor(dds_TPS_male_list[[tissue]]$age_days)
  
  #save normcounts
  dds_TPS_male_list[[tissue]]<- estimateSizeFactors(dds_TPS_male_list[[tissue]])
  ncounts <- counts(dds_TPS_male_list[[tissue]], normalized = T)
  write.csv(ncounts, file = paste0(outdir,'CountsNormDESeq2_',tissue,'_MaleOnly_241106.csv'))
  
  #this is where you run tissue by tissue PCA
  vsd.tissue <- vst(dds_TPS_male_list[[tissue]])
  head(assay(vsd.tissue), 3)
  p <- pca(assay(vsd.tissue), metadata = colData(dds_TPS_male_list[[tissue]]))
  
}

save(dds_TPS_male_list, file = paste0(outdir.robj,'dds_TPS_alltissue_malesOnly_241106.bin'))

#---------------- DESeq2 Differential Expression Across Tissues, with Binning, male only ----------------
#We'll setup a results list in which we'll store the output of the following analysis
results_list_male <- list()

padj_cutoff = 0.05

for (tissue in tissue.list) {
  print(tissue)
  #subset for current tissue and run Deseq2
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('47', '49', '52'), '1', NA)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('75', '77', '78'), '2', dds_TPS_male_list[[tissue]]$age_bin)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('102', '103'), '3', dds_TPS_male_list[[tissue]]$age_bin)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_TPS_male_list[[tissue]]$age_bin)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_TPS_male_list[[tissue]]$age_bin)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('147', '152', '155'), '5', dds_TPS_male_list[[tissue]]$age_bin)
  dds_TPS_male_list[[tissue]]$age_bin <- ifelse(dds_TPS_male_list[[tissue]]$age_days %in% c('161', '162'), '6', dds_TPS_male_list[[tissue]]$age_bin)
  
  dds_TPS_male_list[[tissue]]$age_bin <- as.factor(dds_TPS_male_list[[tissue]]$age_bin)
  design(dds_TPS_male_list[[tissue]]) <- ~age_bin
  
  #Run Deseq2 and store in the Deseq2 object list
  dds_TPS_male_list[[tissue]] <- DESeq(dds_TPS_male_list[[tissue]], fitType = 'local')
  
  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_TPS_male_list[[tissue]]
  results_list_tissue <- list()
  
  #Create a list of all pairwise comparisons
  comparison_list <- t(combn(unique(colData(dds_tissue_temp)$age_bin),2))
  
  #Iterate over all pairwise comparisons and extract the results table, store it in the results list. then flip the conditions
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
  results_list_male[[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_TPS_male_list, results_list_male,  file='Output/robjects/dds_TPS_alltissue_malesOnly_DESeq2results_241106.bin')
  
}

#------------------------------------------------------------------
# Female
# -----------------------------------------------------------------
#---------------- Make list of DESeq2 objs, female only ----------------
dds_TPS_female_list <- list()

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective tissue + sex
for (tissue in tissue.list) {
  print(tissue)
  
  #subset for current tissue and make DESeq object
  dds_TPS_female_list[[tissue]] <- dds_all_tissue[, dds_all_tissue$tissue == tissue & dds_all_tissue$sex == 'F']
  dds_TPS_female_list[[tissue]]$tissue <- droplevels(dds_TPS_female_list[[tissue]]$tissue)
  dds_TPS_female_list[[tissue]]$sex <- as.factor(dds_TPS_female_list[[tissue]]$sex)
  dds_TPS_female_list[[tissue]]$age_days <- as.factor(dds_TPS_female_list[[tissue]]$age_days)
  
  #save normcounts
  dds_TPS_female_list[[tissue]]<- estimateSizeFactors(dds_TPS_female_list[[tissue]])
  ncounts <- counts(dds_TPS_female_list[[tissue]], normalized = T)
  write.csv(ncounts, file = paste0(outdir,'CountsNormDESeq2_',tissue,'_FemaleOnly_241106.csv'))
  
  #this is where you run tissue by tissue PCA
  vsd.tissue <- vst(dds_TPS_female_list[[tissue]])
  head(assay(vsd.tissue), 3)
  p <- pca(assay(vsd.tissue), metadata = colData(dds_TPS_female_list[[tissue]]))
  
}

save(dds_TPS_female_list, file = paste0(outdir.robj,'dds_TPS_alltissue_femalesOnly_241106.bin'))

#---------------- DESeq2 Differential Expression Across Tissues, with Binning, female only ----------------
#We'll setup a results list in which we'll store the output of the following analysis
results_list_female <- list()

padj_cutoff = 0.05

for (tissue in tissue.list) {
  print(tissue)
  #subset for current tissue and run Deseq2
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('47', '49', '52'), '1', NA)
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('75', '77', '78'), '2', dds_TPS_female_list[[tissue]]$age_bin)
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('102', '103'), '3', dds_TPS_female_list[[tissue]]$age_bin)
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_TPS_female_list[[tissue]]$age_bin)
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('133', '134'), '4', dds_TPS_female_list[[tissue]]$age_bin)
  dds_TPS_female_list[[tissue]]$age_bin <- ifelse(dds_TPS_female_list[[tissue]]$age_days %in% c('147', '152', '155'), '5', dds_TPS_female_list[[tissue]]$age_bin)

  dds_TPS_female_list[[tissue]]$age_bin <- as.factor(dds_TPS_female_list[[tissue]]$age_bin)
  design(dds_TPS_female_list[[tissue]]) <- ~age_bin
  
  #Run Deseq2 and store in the Deseq2 object list
  dds_TPS_female_list[[tissue]] <- DESeq(dds_TPS_female_list[[tissue]], fitType = 'local')
  
  #Now take the Deseq2 object and run the extraction of the pairwise comparisons
  dds_tissue_temp <- dds_TPS_female_list[[tissue]]
  results_list_tissue <- list()
  
  #Create a list of all pairwise comparisons
  comparison_list <- t(combn(unique(colData(dds_tissue_temp)$age_bin),2))
  
  #Iterate over all pairwise comparisons and extract the results table, store it in the results list. then flip the conditions
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
  results_list_female[[tissue]]  <- results_list_tissue
  #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
  save(dds_TPS_female_list, results_list_female,  file='Output/robjects/dds_TPS_alltissue_femalesOnly_DESeq2results_241106.bin')
  
}
