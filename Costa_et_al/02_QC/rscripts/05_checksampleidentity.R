rm(list=ls())
# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
library(ComplexHeatmap)
library('ggstatsplot')

# Set wd to the current directory
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')


#RStudio Version 2023.12.1+402 (2023.12.1+402)
#R Version 4.3.3

#Bone, Brain, Eye, Fat, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 

# ------------------------------------------------------------------
# Load data - edited in 000_merge_counts_TPM.R, checkreads, featurecounts
# ------------------------------------------------------------------
countdata <- read.csv("Output/Counts_Atlas_allbatches_merged_v3.csv", row.names = 1)
sampleTable <- read.csv('Output/ExperimentDesign_allbatches_combined_v5.csv', row.names = 1)

# ------------------------------------------------------------------
# Fixing data entry discrepancies
# ------------------------------------------------------------------
table(sampleTable$animalID) #it appears that some of the eye samples have slightly erroneous sampleIDs - let's standardize them
sampleTable$animalID <- gsub("C_03", "C03", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("D_05", "D05", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("I_03", "I03", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("J_04", "J04", as.character(sampleTable$animalID)) 
sampleTable$animalID <- gsub("-", "_", as.character(sampleTable$animalID)) #some of the names had dashes where they should have had underscores
View(table(sampleTable$animalID)) #looks good now

animal_sex <- paste0(sampleTable$animalID, "_", sampleTable$sex)
p <- table(animal_sex) #it appears that some of the eye samples have the sex erroneously labeled
p

p <- as.data.frame(p)
p <- p %>% filter(Freq == 1) %>% filter(!grepl('and',animal_sex)) %>% filter(!grepl('NA',animal_sex))
p$animalID <- sapply(strsplit(as.character(p$animal_sex), "_"), `[`, 1)
p$sex <- sapply(strsplit(as.character(p$animal_sex), "_"), `[`, 2)
p$sex_corr <- ifelse(p$sex == "M", "F", "M")
rownames(p) <- p$animalID

to_repl <- p$animalID
for(i in to_repl){
  idx = which(sampleTable$animalID == i & sampleTable$tissue == 'Eye')
  sampleTable[idx,]$sex <- p[i,]$sex_corr
}

animal_sex <- paste0(sampleTable$animalID, "_", sampleTable$sex)
View(table(animal_sex)) #looks good now



#save the edited Experiment Design file
write.csv(sampleTable, "Output/ExperimentDesign_allbatches_combined_v6.csv")

View(table(paste0(sampleTable$animalID, "_", sampleTable$age_days))) #looks like there is one sample from animal P_1B_1 erroneously annotated as being of age "7", instead of 77
idx = which(sampleTable$age_days == 7)
sampleTable$age_days[idx] = 77

#save the edited Experiment Design file
write.csv(sampleTable, "Output/ExperimentDesign_allbatches_combined_v7.csv")






# ------------------------------------------------------------------
# Check that you're ready to begin analysis
# ------------------------------------------------------------------
countdata <- read.csv("Output/Counts_Atlas_allbatches_merged_v3.csv", row.names = 1)
sampleTable <- read.csv('Output/ExperimentDesign_allbatches_combined_v6.csv', row.names = 1)


coldata <- DataFrame(sampleTable)
col_order <- rownames(coldata)
countdata <- countdata[,col_order]

table(colnames(countdata) == rownames(coldata)) # Make sure that the column names are identical

# ------------------------------------------------------------------
# Check normcounts of ncFem1 and other sex-biased genes
# ------------------------------------------------------------------
# Make dds object 
dds_all_tissue <- DESeqDataSetFromMatrix(countData = countdata,
                                         colData = coldata,
                                         design = ~ age_days + tissue)

# only keep rows that have at least 0 sum count
dds_all_tissue <- dds_all_tissue [ rowSums(counts(dds_all_tissue)) > 0, ]


#---------------- Plot PCA ----------------
# Plot PCA using vst normalization
vsd_all <- vst(dds_all_tissue)
head(assay(vsd_all), 3)

p <- pca(assay(vsd_all), metadata = colData(dds_all_tissue))

biplot(p,
       x = 'PC1', y = 'PC2',
       lab = NA,
       colby = 'tissue',
       colkey = colTPS,
       shape = 'sex',
       legendPosition = 'top')

#---------------- Normalize counts ----------------
# Run DEseq
dds_all_tissue <- DESeq(dds_all_tissue, fitType = 'local')

# Normalized counts
dds_all_tissue <- estimateSizeFactors(dds_all_tissue)
normcounts <- counts(dds_all_tissue, normalized=TRUE)
nrow(normcounts)
colnames(normcounts) <- rownames(sampleTable)
head(normcounts)


write.table(normcounts, file="Output/CountsNormDESeq2_AllTissue_240506.csv", sep= ",")
save(dds_all_tissue,file='Robjects/dds_TPS_allsamples_240506.bin')



####################
normcounts <- read.csv('Output/CountsNormDESeq2_AllTissue_240506.csv', row.names = 1)
load('Robjects/dds_TPS_allsamples_240506.bin')


#from Teefy et al in Cell Reports
#most frequently female biased genes: LOC107373896 (ncFem1), zp3, and vitellogenin-1-like - ncFem1 female bias holds in several tissues
#most frequently male biased genes: hpx, apoa1a, aldob

#vitellogenin-1-like not annotated as such in this assembly, so going to look at a few orthologs: LOC107392618, LOC107392704
#same for aldob: LOC107375611
#same for apoa1a: LOC107395221

fem.biased <- c('LOC107373896', 'LOC107392704','LOC107392618', 'LOC107381613', 'LOC107388674', 'LOC107388898')
masc.biased <- c('hpx', 'LOC107375611', 'LOC107395221')
sex.biased <- c(fem.biased, masc.biased)

fem.biased.cts <- subset(normcounts, rownames(normcounts) %in% fem.biased)
fem.biased.cts.t <- as.data.frame(t(fem.biased.cts))
fem.biased.cts.t$lib <- rownames(fem.biased.cts.t)

fem.biased.cts.t <- merge(fem.biased.cts.t, sampleTable, by = "lib")
plotCounts(dds_all_tissue, gene = 'LOC107373896', intgroup = 'tissue')

fem.biased.cts <- subset(normcounts, rownames(normcounts) %in% c('ppp2r5a'))
fem.biased.cts.t <- as.data.frame(t(fem.biased.cts))
fem.biased.cts.t$lib <- rownames(fem.biased.cts.t)
fem.biased.cts.t <- merge(fem.biased.cts.t, sampleTable, by = "lib")
plotCounts(dds_all_tissue, gene = 'ppp2r5a', intgroup = 'tissue')

ggbetweenstats(
  data = fem.biased.cts.t,
  x = age_days,
  y = 'ppp2r5a')

# Gathering the columns to have normalized counts to a single column
#fem.biased.cts <- as.data.frame(fem.biased.cts)
#gathered <- fem.biased.cts %>%
#  gather(colnames(fem.biased.cts), key = "sampleNames", value = "normalized_counts", rownames(fem.biased.cts))

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107388898') #zp3

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107373896') #ncFem1

ggplot(fem.biased.cts.t, aes(sex, LOC107373896)) + geom_jitter(width=0.2, alpha=0.5) + geom_boxplot(fill = NA, outlier.color = NA)

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107392704') #vitellogenin-1-like

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107392618') #vitellogenin-1-like

fem.biased.cts.t$tissue_sex <- paste(fem.biased.cts.t$tissue, fem.biased.cts.t$sex, sep = "_")

#separate by tissue as shown by Teefy et al (2023) Cell Reports, Fig 5F
#Liver only
liver.cts <- fem.biased.cts.t %>% filter(tissue == 'Liver')
ggbetweenstats(
  data = liver.cts,
  x = sex,
  y = 'LOC107388898') #zp3

ggbetweenstats(
  data = liver.cts,
  x = sex,
  y = 'LOC107373896') #ncfem1

ggplot(liver.cts, aes(sex, LOC107373896, fill = sex)) + 
  geom_jitter(width=0.2, alpha=0.5) + 
  geom_boxplot(fill = NA, outlier.color = NA, width = 0.3) +
  theme_classic()

ggbetweenstats(
  data = liver.cts,
  x = sex,
  y = 'LOC107392704') #vitellogenin-1-like

ggbetweenstats(
  data = liver.cts,
  x = sex,
  y = 'LOC107392618') #vitellogenin-1-like

#Spleen only
spleen.cts <- fem.biased.cts.t %>% filter(tissue == 'Spleen')
ggbetweenstats(
  data = spleen.cts,
  x = sex,
  y = 'LOC107388898') #zp3

ggbetweenstats(
  data = spleen.cts,
  x = sex,
  y = 'LOC107373896') #ncfem1

ggbetweenstats(
  data = spleen.cts,
  x = sex,
  y = 'LOC107392704') #vitellogenin-1-like

ggbetweenstats(
  data = spleen.cts,
  x = sex,
  y = 'LOC107392618') #vitellogenin-1-like

#Muscle only
muscle.cts <- fem.biased.cts.t %>% filter(tissue == 'Muscle')
ggbetweenstats(
  data = muscle.cts,
  x = sex,
  y = 'LOC107388898') #zp3

ggbetweenstats(
  data = muscle.cts,
  x = sex,
  y = 'LOC107373896') #ncfem1

ggbetweenstats(
  data = muscle.cts,
  x = sex,
  y = 'LOC107392704')  #vitellogenin-1-like

ggbetweenstats(
  data = muscle.cts,
  x = sex,
  y = 'LOC107392618')  #vitellogenin-1-like


# ------------------------------------------------------------------
# Check normcounts for ncFem1 (LOC107373896) in ALDR dataset
# ------------------------------------------------------------------
#normcounts interaction term
normcounts.aldr1 <- read.csv('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/DEseq2/Output/CountsNormDESeq2_Liver_ALDR_interaction_220331.csv', row.names = 1)
sampleTable.aldr <- read.csv('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishFeederPaper_AndrewMcKay/Revision/Code/RNAseq_Code_check/DEseq2/Input/DRexp20211219ExperimentDesign.csv', row.names = 1)

fem.biased.cts <- subset(normcounts.aldr1, rownames(normcounts.aldr1) %in% fem.biased)
fem.biased.cts.t <- as.data.frame(t(fem.biased.cts))
fem.biased.cts.t$lib <- rownames(fem.biased.cts.t)

fem.biased.cts.t <- merge(fem.biased.cts.t, sampleTable.aldr, by = "lib")

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107373896')

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107392704')

ggbetweenstats(
  data = fem.biased.cts.t,
  x = sex,
  y = 'LOC107392618')


write.csv(fem.biased.cts.t, 'Output/240429_Liver_ALDR_interaction_normcounts_femaleenrichgenes.csv')



# ------------------------------------------------------------------
# Check normcounts for tissue-specific genes as determined by Findallmarkers in the seurat pipeline
# ------------------------------------------------------------------
#Tissue color
tissue_labels <- c('Bone', 'Brain', 'Eye', 'Fat', 'Gut', 'Heart', 'Kidney', 'Liver', 'Muscle', 'Ovary', 'Skin', 'SpinalCord', 'Spleen', 'Testis')
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 

tissue_labels <- c('Bone', 'Brain', 'Fat', 'Gut', 'Heart', 'Kidney', 'Liver', 'Muscle', 'Ovary', 'Retina','Skin', 'SpinalCord', 'Spleen', 'Testis')
c('#d8413f', '#00a550', '#eee09b','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#c9bad4' ,'#b8b8c0','#ab5673', '#f1a8a4','#ef9ac2','#93cca8' ) 


#Sex: male, female 
sex_labels = c('M', "F")
colsex <-c('#068ec9','#ba1e2d')

#Ages (young to old): 
colage <- c('#a92b46', '#d15447', '#e98565', '#f7b77b', '#f3d4ac', '#dfeaf0', '#a2d3e9', '#6bb2dc', '#3789c0', '#395ea4')


#load the cluster markers
tissue.markers <- read.delim('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2/Output_tables/atlas_tissueMarkers_Deseq2.txt')

#1 - plot heatmap of top 10 markers from each tissue, across all samples
for (t in tissue_labels){
  if(exists('top.markers.alltissues')){
    top.markers.temp <- tissue.markers %>% 
      filter(cluster == paste(t)) %>%
      top_n(10, avg_log2FC)
    
    top.markers.alltissues <- rbind(top.markers.alltissues, top.markers.temp)
  }
  else {
    top.markers.alltissues <- tissue.markers %>% 
      filter(cluster == paste(t)) %>%
      top_n(10, avg_log2FC)
  }
}


ncts.alltissue <- subset(normcounts, rownames(normcounts) %in% top.markers.alltissues$gene)
ncts.alltissue.t <- as.data.frame(t(ncts.alltissue))
ncts.alltissue.t$lib <- rownames(ncts.alltissue.t) 
ncts.alltissue.t <- merge(ncts.alltissue.t, sampleTable, by = 'lib')

ncts.alltissue <- as.matrix(ncts.alltissue)


#create a dataframe to label your heatmap - inspo https://github.com/jokergoo/ComplexHeatmap/issues/1121 
heatmap_annotate <- sampleTable
heatmap_annotate$col <- plyr::mapvalues(x = sampleTable$tissue, from = tissue_labels, to = colTPS)
heatmap_annotate$sex_col <- plyr::mapvalues(x = sampleTable$sex, from = sex_labels, to = colsex)
col_list <- setNames(heatmap_annotate$col, heatmap_annotate$tissue)
sex_cols <- setNames(heatmap_annotate$sex_col, heatmap_annotate$sex)

column_ha = HeatmapAnnotation(df = heatmap_annotate[,c("tissue", "sex")],
                              col = list("tissue" = col_list, "sex" = sex_cols))
p <- Heatmap(ncts.alltissue, 
             bottom_annotation = column_ha,
             cluster_rows = T)
p <- draw(p)


#2 - check testis and ovary samples -  heatmap and normcounts
## By testis markers ##
top.testis.markers <- tissue.markers %>% 
  filter(cluster == 'Testis') %>%
  top_n(30, avg_log2FC)

ncts.testismarkers <- subset(normcounts, rownames(normcounts) %in% top.testis.markers$gene)
ncts.testis.t <- as.data.frame(t(ncts.testismarkers))
ncts.testis.t$lib <- rownames(ncts.testis.t)

ncts.testis.t <- merge(ncts.testis.t, sampleTable, by = "lib")

gonad.cts <- ncts.testis.t %>% filter(tissue %in% c('Ovary', 'Testis'))
ggbetweenstats(
  data = gonad.cts,
  x = sex,
  y = 'LOC107394067')


gonad.cts.subset <- gonad.cts[,c('LOC107394067', 'tissue', 'lib')]


ncts.testismarkers <- as.matrix(ncts.testismarkers)
Heatmap(ncts.testismarkers)

## By ovary markers ##
top.ovary.markers <- tissue.markers %>% 
  filter(cluster == 'Ovary') %>%
  top_n(30, avg_log2FC)

ncts.ovarymarkers <- subset(normcounts, rownames(normcounts) %in% top.ovary.markers$gene)
ncts.ovary.t <- as.data.frame(t(ncts.ovarymarkers))
ncts.ovary.t$lib <- rownames(ncts.ovary.t)

gonad.cts <- ncts.ovary.t %>% filter(tissue %in% c('Ovary', 'Testis'))
ggbetweenstats(
  data = gonad.cts,
  x = sex,
  y = 'LOC107394067') #top testis marker

ncts.ovarymarkers <- as.matrix(ncts.ovarymarkers)

#create a dataframe to label your heatmap - inspo https://github.com/jokergoo/ComplexHeatmap/issues/1121 
heatmap_annotate <- sampleTable
heatmap_annotate$col <- plyr::mapvalues(x = sampleTable$tissue, from = tissue_labels, to = colTPS)
col_list <- setNames(heatmap_annotate$col, heatmap_annotate$tissue)

column_ha = HeatmapAnnotation(tissue = heatmap_annotate[,c("tissue")],
                              col = list("tissue" = col_list))
p <- Heatmap(ncts.ovarymarkers, bottom_annotation = column_ha)
p <- draw(p)

#tips when saving heatmap -  100 inch is a good width





