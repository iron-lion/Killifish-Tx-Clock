# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
#RStudio Version 2024.04.1+748 (2024.04.1+748)
#R Version 4.3.3

library(DESeq2)
library(PCAtools)
library(dplyr)
library(tidyr)
library(ggbreak)

# Set wd to the current directory
setwd("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/")

# Load the DESeq2 object list
object.indir = "Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

outdir = paste0(getwd(), '/Output/Tables/')
outdir.robj = paste0(getwd(), '/robjects/')

#------------------------------------------------------------------
# DESeq2
# -----------------------------------------------------------------

#list was saving as it was going. check where left off and continue from there

#---------------- Make list of DESeq2 objs, for each tissue - age bin combo bins 1-5 ----------------
tissue.list <- names(dds_TPS_list)[2:14]

dds_TPS_agebin_list <- list()
dds_all_tissue <- dds_TPS_list[['All']]

#Now we run a loop across all tissues/regions, where we subset the main Deseq2 object for samples of the respective tissue + sex
for (tissue in tissue.list) {
  print(tissue)
  
  #subset for current tissue and make DESeq object
  dds_temp <- dds_all_tissue[, dds_all_tissue$tissue == tissue]
  dds_temp$tissue <- droplevels(dds_temp$tissue)
  
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('47', '49', '52'), '1', NA)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('75', '77', '78'), '2', dds_temp$age_bin)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('102', '103'), '3', dds_temp$age_bin)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('133', '134'), '4', dds_temp$age_bin)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('133', '134'), '4', dds_temp$age_bin)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('147', '152', '155'), '5', dds_temp$age_bin)
  dds_temp$age_bin <- ifelse(dds_temp$age_days %in% c('161', '162'), '6', dds_temp$age_bin)
  dds_temp$age_bin <- as.factor(dds_temp$age_bin)
  dds_temp$sex <- as.factor(dds_temp$sex)
  
  for(i in c('1','2','3','4','5')){
    print(i)
    
    dds_TPS_agebin_list[[tissue]][[i]] <- dds_temp[, dds_temp$age_bin == i]
    dds_TPS_agebin_list[[tissue]][[i]]$age_bin <- droplevels(dds_TPS_agebin_list[[tissue]][[i]]$age_bin)

    #save normcounts
    dds_TPS_agebin_list[[tissue]][[i]]<- estimateSizeFactors(dds_TPS_agebin_list[[tissue]][[i]])
    
    #this is where you run tissue by tissue PCA
    vsd.tissue.age <- vst(dds_TPS_agebin_list[[tissue]][[i]])
  }
}

save(dds_TPS_agebin_list, file = paste0(outdir.robj,'dds_TPS_alltissue_byagebin_241121.bin'))

#---------------- DESeq2 Differential Expression Across Tissues, at agebin between sexes ----------------
#We'll setup a results list in which we'll store the output of the following analysis
load(file='Robjects/dds_TPS_alltissue_byagebin-sexdiffs_DESeq2results_241121.bin')

results_list <- list()
padj_cutoff = 0.05

for (tissue in tissue.list) {
  print(tissue)
  
  for(i in c('1','2','3','4','5')){
    print(i)
    
    #subset for current tissue and run Deseq2
    design(dds_TPS_agebin_list[[tissue]][[i]]) <- ~sex
    
    #Run Deseq2 and store in the Deseq2 object list
    dds_TPS_agebin_list[[tissue]][[i]] <- DESeq(dds_TPS_agebin_list[[tissue]][[i]], fitType = 'local')
    
    #Now take the Deseq2 object and run the extraction of the pairwise comparisons
    dds_tissue_temp <- dds_TPS_agebin_list[[tissue]][[i]]
    results_list_tissue.age <- list()
    
    #Create a list of all pairwise comparisons
    comparison_list <- t(combn(unique(colData(dds_tissue_temp)$sex),2))
    
    setwd("/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2_bysex/")
    #Iterate over all pairwise comparisons and extract the results table, store it in the results list. then flip the conditions
    start_time <- Sys.time()
    for (row in 1:nrow(comparison_list)) {
      print(row)
      #Get the conditions/timepoints to test
      cond1 <- as.character(comparison_list[row,1])
      cond2 <- as.character(comparison_list[row,2])
      
      folder <- paste(cond1, cond2, sep = "_vs_")
      print(folder)
      
      aspect_to_test <- 'sex'
      
      results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
      results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
      resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
      resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
      dim(resSig)
      #Store the output tables in the results_list
      results_list_tissue.age[[folder]]$resall <- results_temp
      results_list_tissue.age[[folder]]$resOrdered <- resOrdered
      results_list_tissue.age[[folder]]$ressig <- resSig
      
      
      ##Flip conditions and re-run the results extraction
      
      cond1 <- as.character(comparison_list[row,2])
      cond2 <- as.character(comparison_list[row,1])
      folder <- paste(cond1, cond2, sep = "_vs_")
      print(folder)
      
      dds_tissue_temp <- dds_tissue_temp
      
      aspect_to_test <- 'sex'
      
      results_temp <- results(dds_tissue_temp, contrast=c(aspect_to_test, cond1,cond2), cooksCutoff=T) #Set contrast and retrieve results
      results_temp$gene_symbol <- row.names(results_temp) #get the rownames as a separate column in the results report
      resOrdered <- as.data.frame(results_temp[order(results_temp$pvalue),]) #create a simple dataframe ordered by the padj
      resSig <- subset(resOrdered, padj < padj_cutoff) #create a simple dataframe containing only the significant results
      dim(resSig)
      #Store the output tables in the results_list
      results_list_tissue.age[[folder]]$resall <- results_temp
      results_list_tissue.age[[folder]]$resOrdered <- resOrdered
      results_list_tissue.age[[folder]]$ressig <- resSig
      
    }
    end_time <- Sys.time()
    end_time-start_time
    
    #Store tissue results in the major results list, then save the output and start the loop for the next tissue
    results_list[[tissue]][[i]]  <- results_list_tissue.age
    #Store both the results lists and the list of deseq2 objects, as these contain all the processed data
    save(dds_TPS_agebin_list, results_list,  file='Robjects/dds_TPS_alltissue_byagebin-sexdiffs_DESeq2results_241122.bin')
  }
}

# ------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------
setwd("/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2_bysex/")
outdir.deseq = 'Output/'

load('Robjects/dds_TPS_alltissue_byagebin-sexdiffs_DESeq2results_241122.bin')
tissue.list <- names(results_list) #no eye bc there wasnt enough sex balance
tissue.list <- tissue.list[-3]

#---------------- Explore the number of genes for each  ----------------
sex.terms <- c("M_vs_F", "F_vs_M")

sex.degs <- data.frame(tissue = NA)
sex.degs$term <- NA
sex.degs$pos.deg.number <- NA
sex.degs$neg.deg.number <- NA
sex.degs$age_bin <- NA

for(t in tissue.list){
  tissue.res <- results_list[[t]]
  for(i in c('1','2', '3', '4', '5')){
   tissue.age.res <- tissue.res[[i]] 
   for(int in sex.terms){
     res.sig = tissue.age.res[[int]]$ressig
     deg.num.down = length(which(res.sig$log2FoldChange < 0))
     deg.num.up = length(which(res.sig$log2FoldChange > 0))
     sex.degs[nrow(sex.degs)+1,] = c(t,int,deg.num.up,deg.num.down,i)
   }
  }
}

sex.degs = sex.degs[-1,]

#Bone, Brain, Fat, Gonad,Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550','#eee09b','#7962A3' ,'#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489' ,'#ab5673', '#f1a8a4','#ef9ac2') 

df.sex = sex.degs

df.sex$term <- factor(df.sex$term, levels = sex.terms)
df.sex$age_bin <- factor(df.sex$age_bin, levels = c('1', '2', '3', '4', '5'))
df.sex$pos.deg.number <- as.numeric(df.sex$pos.deg.number)
df.sex$neg.deg.number <- as.numeric(df.sex$neg.deg.number)

# Plotting
df.sex %>% filter(term == 'M_vs_F') %>%
ggplot(aes(x = age_bin, y = pos.deg.number, group = tissue, color = tissue)) +
  geom_line() +
  geom_point() +
  labs(x = "Age Bin", y = "Positive DEG Number", title = "Positive DEG Numbers by Age Bin for Each Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_color_manual(values = colTPS)

# Plotting
df.sex %>% filter(term == 'M_vs_F') %>%
ggplot(aes(x = age_bin, y = neg.deg.number, group = tissue, color = tissue)) +
  geom_line() +
  geom_point() +
  labs(x = "Age Bin", y = "Negative DEG Number", title = "Negative DEG Numbers by Age Bin for Each Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_color_manual(values = colTPS)


#express as % expressed 
outdir.cor = "/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Correlation/Output/"
num_Expressed <- read.csv(file = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/Jiachen_EngelhardtLab/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv')
#num_Expressed <- read.csv(paste0(outdir.cor,"Tables/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv")) #made in correlation_checkconfounds.R

df.sex = merge(df.sex, num_Expressed)
df.sex$tot.deg.number = df.sex$pos.deg.number + df.sex$neg.deg.number
df.sex$tot.deg.number.perc = ((df.sex$tot.deg.number) / (df.sex$total.expressed))*100

df.sex$pos.deg.number.perc = ((df.sex$pos.deg.number) / (df.sex$total.expressed))*100
df.sex$neg.deg.number.perc = ((df.sex$neg.deg.number) / (df.sex$total.expressed))*100

###TODO: Plot and expand the scale a bit, label the top 3 and bottom 3 tissues
#pdf(paste0(outdir.deseq, 'Plots/241218_totdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = tot.deg.number.perc, group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "DEG Number % of Transcriptome", title = "DEG by Age Bin for Each Tissue: M vs F") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,110)) +
  scale_y_break(c(40, 80)) 
#dev.off()

#pdf(paste0(outdir.deseq, 'Plots/241122_posdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
ggplot(aes(x = age_bin, y = pos.deg.number.perc, group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "Positive DEG Number % of Transcriptome", title = "Positive DEG Numbers by Age Bin for Each Tissue") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,55)) +
  scale_y_break(c(20, 40)) 
#dev.off()

# Plotting
#pdf(paste0(outdir.deseq, 'Plots/241122_negdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
ggplot(aes(x = age_bin, y = (neg.deg.number.perc), group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "Negative DEG Number % of Transcriptome", title = "Negative DEG Numbers by Age Bin for Each Tissue") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,55)) +
  scale_y_break(c(20, 40))  
#dev.off()

#save the output table
write.csv(df.sex, file = paste0(outdir.deseq, '241218_M_vs_F_DEGs_byage.csv'))

