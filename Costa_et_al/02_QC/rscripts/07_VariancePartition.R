#Variance Partition script

library(DESeq2)
library(variancePartition) #v1.33.11
library(ggplot2)
library(dplyr)
library(tidyr)

# ------------------------------------------------------------------
# Please run 000_merge_counts_TPM.R to get merged TPM matrices
# ------------------------------------------------------------------
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')

tpm <- read.csv('Output/TPM_Atlas_allbatches_merged_v3.csv', stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1)
metadata <- read.csv(file="Output/ExperimentDesign_allbatches_combined_v7.csv", stringsAsFactors = FALSE, header = TRUE, check.names = FALSE, row.names = 1 )

col_order <- rownames(metadata)
tpm <- tpm[,col_order]

table(rownames(metadata)==colnames(tpm)) #check order columns/rows in tpm and metadata file

# Load the DESeq2 object list
object.indir = "Path/to/DESeq2_Robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))


# ------------------------------------------------------------------
# Run Variance partition by tissue
# ------------------------------------------------------------------
#vignette: https://bioconductor.org/packages/release/bioc/vignettes/variancePartition/inst/doc/variancePartition.html
outdir = 'Output/'
subDir = 'Plots/'

#create subdir for plots
dir.create(file.path(outdir, subDir), showWarnings = FALSE)


tissue.list<- unique(metadata$tissue)
tissue.list <- tissue.list[-9] #run without gonad
tissue.list <- tissue.list[-10] #run without gonad
form1 <- c('Bone', 'Muscle', 'Fat') #all same RNA extractor
form2 <- c('Eye') #all same RNA extractor, all same RNA batch


for(t in tissue.list){
  metadata_tissue = metadata %>% filter(tissue == t)
  
  #add age_bins
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('47', '49', '52'), '1', NA)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('75', '77', '78'), '2', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('102', '103'), '3', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('147', '152', '155'), '5', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('161', '162'), '6', metadata_tissue$age_bin)

  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(metadata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue))
  print(paste0(t," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]

  ExpresedGenes <- rownames(expressed)

  # Variance Partition
  geneExpres <- expressed

  if(! t %in% c(form1, form2)){
    form = ~ (1|age_bin) + (1|sex) + (1|cohort) + (1|RNA_batch) + (1|RNA_extractor) + (1 | sex:age_bin)
  }
  else{
    if(t %in% form1){
      form = ~ (1|age_bin) + (1|sex) + (1|cohort) + (1|RNA_batch) + (1 | sex:age_bin)
    }
    else{
      form = ~ (1|age_bin) + (1|sex) + (1|cohort) + (1 | sex:age_bin)
    }
  }
  print(form)
  varPart <- fitExtractVarPartModel( geneExpres, form, metadata_tissue )
  vp <- sortCols(varPart )

  write.csv(vp, file = paste0(outdir,"Tables/241121_",t,"_variancePartition_results.csv"))

  # Compute Canonical Correlation Analysis (CCA)
  # between all pairs of variables
  # returns absolute correlation value
  C <- canCorPairs(form, metadata_tissue)

  # Plot correlation matrix
  # between all pairs of variables
  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_corMatrix.pdf"))
  print(plotCorrMatrix(C))
  dev.off()

  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_violinPlot.pdf"))
  print(plotVarPart(vp))
  dev.off()
}


#Model Testis and Ovary separately
remaining.tissue <- c('Testis', 'Ovary')

for(t in remaining.tissue){
  metadata_tissue = metadata %>% filter(tissue == t)
  
  #add age_bins
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('47', '49', '52'), '1', NA)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('75', '77', '78'), '2', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('102', '103'), '3', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('147', '152', '155'), '5', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('161', '162'), '6', metadata_tissue$age_bin)

  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(metadata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue))
  print(paste0(t," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]

  ExpresedGenes <- rownames(expressed)

  # Variance Partition
  geneExpres <- expressed

  form <- ~ (1 | age_bin) + (1|cohort) + (1|RNA_batch) + (1|RNA_extractor) #Testis, Ovary
  varPart <- fitExtractVarPartModel( geneExpres, form, metadata_tissue )
  vp <- sortCols(varPart )

  write.csv(vp, file = paste0(outdir,"Tables/241121_",t,"_variancePartition_results.csv"))

  # Compute Canonical Correlation Analysis (CCA)
  # between all pairs of variables
  # returns absolute correlation value
  C <- canCorPairs(form, metadata_tissue)

  # Plot correlation matrix
  # between all pairs of variables
  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_corMatrix.pdf"))
  print(plotCorrMatrix(C))
  dev.off()

  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_violinPlot.pdf"))
  print(plotVarPart(vp))
  dev.off()
}


#Model Testis and Ovary together 
remaining.tissue <- c('Gonad')

metadata$tissue <- gsub('Testis','Gonad', metadata$tissue)
metadata$tissue <- gsub('Ovary','Gonad', metadata$tissue)

for(t in remaining.tissue){
  metadata_tissue = metadata %>% filter(tissue == t)
  
  #add age_bins
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('47', '49', '52'), '1', NA)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('75', '77', '78'), '2', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('102', '103'), '3', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('133', '134'), '4', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('147', '152', '155'), '5', metadata_tissue$age_bin)
  metadata_tissue$age_bin <- ifelse(metadata_tissue$age_days %in% c('161', '162'), '6', metadata_tissue$age_bin)
  
  #subset tpm by tissue type and find genes for which TPM > 0.5 in 80% or more of samples
  tpm_tissue <- tpm[,which(colnames(tpm) %in% rownames(metadata_tissue))]
  pres = apply(tpm_tissue >=0.5,1,sum)
  Expressed_tpm = (pres >= 0.80*ncol(tpm_tissue))
  print(paste0(t," - Number of genes passing tpm filtering criteria: ",table(Expressed_tpm)[2]))
  expressed = tpm_tissue[Expressed_tpm,]
  
  ExpresedGenes <- rownames(expressed)
  
  # Variance Partition
  geneExpres <- expressed
  
  form <- ~ (1|age_bin) + (1|sex) + (1|cohort) + (1|RNA_batch) + (1|RNA_extractor) + (1 | sex:age_bin) #Gonad
  varPart <- fitExtractVarPartModel( geneExpres, form, metadata_tissue )
  vp <- sortCols(varPart )
  
  write.csv(vp, file = paste0(outdir,"Tables/241121_",t,"_variancePartition_results.csv"))
  
  # Compute Canonical Correlation Analysis (CCA)
  # between all pairs of variables
  # returns absolute correlation value
  C <- canCorPairs(form, metadata_tissue)
  
  # Plot correlation matrix
  # between all pairs of variables
  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_corMatrix.pdf"))
  print(plotCorrMatrix(C))
  dev.off()
  
  pdf(file = paste0(outdir,"Plots/241121_",t,"_variancePartition_violinPlot.pdf"))
  print(plotVarPart(vp))
  dev.off()
}


# ------------------------------------------------------------------
# Aggregate results
# ------------------------------------------------------------------

varPart.list = list()

for (t in c(tissue.list, remaining.tissue)){
  vp = read.csv(file = paste0(outdir,"Tables/241121_",t,"_variancePartition_results.csv"))
  vp$tissue = t
  if(! 'sex' %in% colnames(vp)){
    vp$sex = NA
  }
  if(! 'sex.age_bin' %in% colnames(vp)){
    vp$'sex.age_bin' = NA
  }
  if(! 'RNA_batch' %in% colnames(vp)){
    vp$'RNA_batch' = NA
  }
  if(! 'RNA_extractor' %in% colnames(vp)){
    vp$'RNA_extractor' = NA
  }
  varPart.list[[t]] = vp
}

##save aggregated list of variancePartition values
varPart.list = plyr::ldply(varPart.list, rbind)
varPart.list$.id = NULL
colnames(varPart.list)[1] = 'gene'
write.csv(varPart.list, file = paste0(outdir,"Tables/241121_Alltissue_variancePartition_results_Gonadcombo.csv"))

##summary df
data.summary <- data.frame(tissue = unique(varPart.list$tissue))
data.summary$age.median <- NA
data.summary$sex.median <- NA
data.summary$sex.age.median <- NA

for(i in 1:nrow(data.summary)){
  temp <- data.summary[i,]
  t <- temp$tissue
  df <- varPart.list %>% filter(tissue == t)
  data.summary[i,]$age.median = median(df$age_bin)*100
  data.summary[i,]$sex.median = median(df$sex)*100
  data.summary[i,]$sex.age.median = median(df$sex.age_bin)*100
}

# ------------------------------------------------------------------
# Plot results in median values
# ------------------------------------------------------------------
#Bone, Brain,Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Retina/RPE, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#C2B7D1','#b8b8c0','#ab5673', '#f1a8a4','#ef9ac2','#93cca8') 

##age.median barplot
#Retina, Ovary, Fat, Muscle, SpinalCord,  "Skin"       "Testis"     "Brain"      "Heart"      "Kidney"     "Spleen"     "Gut"        "Liver"      "Bone"    
#plotcols <- c('#b8b8c0', '#C2B7D1','#eee09b','#f4c489', '#f1a8a4','#ab5673','#93cca8', '#00a550','#f0932e','#fcd328','#ef9ac2', '#010101', '#6cc0ee','#d8413f' ) 

#Retina, Muscle, Fat, "Skin", Kidney, Gut, Spleen, Spinalcord, Heart, Gonad,Brain, Bone, Liver
plotcols <- c('#b8b8c0','#f4c489','#eee09b', '#ab5673', '#fcd328', '#010101', '#ef9ac2', '#f1a8a4', '#f0932e','#7962A3','#00a550', '#d8413f' , '#6cc0ee' ) 


data.summary <- data.summary[order(-data.summary$age.median), ]

tissue_levels <- data.summary$tissue
data.summary$tissue <- factor(data.summary$tissue, levels = tissue_levels)

pdf(file = paste0(outdir,"Plots/241121_Alltissue_variancePartition_Age-Median_barPlot.pdf"), width = 2, height = 2)
ggplot(data.summary, aes(x=tissue, y=age.median)) + 
  geom_bar(stat = 'identity', aes(fill=tissue)) + 
  scale_fill_manual(values=plotcols) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  scale_y_continuous(limits = c(0,20)) +
  ggtitle('Age - % Variance Explained') + 
  ylab('Median % Variance')
dev.off()

##sex.median barplot
#Gonad, Skin, Fat, Kidney, Heart, Bone, SpinalCord, Muscle, Liver, Gut, Spleen, Eye, Brain
plotcols2 <- c('#7962A3','#ab5673','#eee09b','#fcd328', '#f0932e','#d8413f', '#f1a8a4', '#f4c489', '#6cc0ee','#010101', '#ef9ac2', '#b8b8c0', '#00a550') 

data.summary$tissue <- droplevels(data.summary$tissue)
data.summary <- data.summary[order(-data.summary$sex.median), ]

tissue_levels <- data.summary$tissue
data.summary$tissue <- factor(data.summary$tissue, levels = tissue_levels)

pdf(file = paste0(outdir,"Plots/241121_Alltissue_variancePartition_Sex-Median_barPlot.pdf"), width = 2, height = 2)
ggplot(data.summary, aes(x=tissue, y=sex.median)) + #plot limited to stuff plottable (1 or greater)
  geom_bar(stat = 'identity', aes(fill=tissue)) + 
  scale_fill_manual(values=plotcols2) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  scale_y_continuous(limits = c(0,70)) +
  ggtitle('Sex - % Variance Explained') + 
  ylab('Median % Variance')
dev.off()


##sex.age. median barplot
#Liver, Bone, Skin, Gonad, Fat, Heart, Spinalcord, kidney, Muscle, Gut, Spleen, Eye, Brain
plotcols3 <- c('#6cc0ee','#d8413f','#ab5673','#7962A3','#eee09b','#f0932e', '#f1a8a4','#fcd328', '#f4c489', '#010101', '#ef9ac2', '#b8b8c0', '#00a550') 

data.summary$tissue <- droplevels(data.summary$tissue)
data.summary <- data.summary[order(-data.summary$sex.age.median), ]

tissue_levels <- data.summary$tissue
data.summary$tissue <- factor(data.summary$tissue, levels = tissue_levels)

pdf(file = paste0(outdir,"Plots/241121_Alltissue_variancePartition_SexAge-Median_barPlot.pdf"), width = 2, height = 2)
ggplot(data.summary, aes(x=tissue, y=sex.age.median)) + #plot limited to stuff plottable (1 or greater)
  geom_bar(stat = 'identity', aes(fill=tissue)) + 
  scale_fill_manual(values=plotcols3) +
  theme_classic() +
  theme(
    legend.position = 'none',
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank()
  ) + 
  scale_y_continuous(limits = c(0,20)) +
  ggtitle('Sex:Age - % Variance Explained') + 
  ylab('Median % Variance')
dev.off()

write.csv(data.summary, file = paste0(outdir,"Tables/241121_Alltissue_variancePartition_medians_gonadcombo.csv"))
