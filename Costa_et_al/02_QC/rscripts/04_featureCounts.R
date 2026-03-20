setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')

library('tidyverse')
library('readxl')
library('base')
library(stringr)
library('ggstatsplot')

# Load experiment design and count matrices
sampleTable <- read.csv("Output/ExperimentDesign_allbatches_combined_v4.csv", row.names = 1)  
lib_name <- as.character(sampleTable$lib)

plate1.samples <- subset(sampleTable, cDNA_batch == 1)
plate2.samples <- subset(sampleTable, cDNA_batch == 2)

########## Load mapping summaries - original from multiqc ##########
# Input data description: multiqc run on featureCounts

# file 1: multiqc_featureCounts.txt
ft.ct1 <- read.delim('Input/multiqc_featureCounts_batch1.txt', header = T, sep = '\t')
ft.ct2 <- read.delim('Input/multiqc_featureCounts_batch2.txt', header = T, sep = '\t')
ft.ct3 <- read.delim('Input/multiqc_featureCounts_batch3.txt', header = T, sep = '\t')
ft.ct4 <- read.delim('Input/multiqc_featureCounts_batch4.txt', header = T, sep = '\t')


### Plate 1
ft.ct1$SampleID <- NA
ft.ct1$SampleID <- str_split_fixed(ft.ct1$Sample, "_",2)[,1]
ft.ct1.summary <- ft.ct1 %>% group_by(SampleID) %>%
  summarise(percent_assigned)
write.csv(ft.ct1.summary, "Output/featurecountsummary_batch1_240424.csv")

ft.ct3$SampleID <- NA
ft.ct3$SampleID <- str_split_fixed(ft.ct3$Sample, "_",2)[,1]
ft.ct3.summary <- ft.ct3 %>% group_by(SampleID) %>%
  summarise(percent_assigned)
write.csv(ft.ct3.summary, "Output/featurecountsummary_batch3_240424.csv")

#add readcounts together across seq batches
ft.ct.plate1 <- merge(ft.ct1.summary, ft.ct3.summary, by = 'SampleID')
colnames(ft.ct.plate1) <- c("sampleNames", "featureCount_batch1", "featureCount_batch3")
ft.ct.plate1 <- ft.ct.plate1 %>%
  mutate(Mean.Features.Assigned = rowMeans(across(c(featureCount_batch1, featureCount_batch3))))

#summarize by tissue type
plate1.all.data <- merge(ft.ct.plate1, plate1.samples, by='sampleNames')

plot1 <- ggbetweenstats(
  data = plate1.all.data,
  x = tissue,
  y = Mean.Features.Assigned)

plot1 <- plot1 + 
  # Add labels and title
  labs(
    x = "Tissue",
    y = "Mean % Features Assigned",
    title = "% Features Assigned across tissue types"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 20,
      face = "bold",
      color = "#2a475e"
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 15, 
      face = "bold",
      color="#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

plot1


### Plate 2
ft.ct2$SampleID <- NA
ft.ct2$SampleID <- str_split_fixed(ft.ct2$Sample, "_",2)[,1]
ft.ct2.summary <- ft.ct2 %>% group_by(SampleID) %>%
  summarise(percent_assigned)
write.csv(ft.ct2.summary, "Output/featurecountsummary_batch2_240424.csv")

ft.ct4$SampleID <- NA
ft.ct4$SampleID <- str_split_fixed(ft.ct4$Sample, "_",2)[,1]
ft.ct4.summary <- ft.ct4 %>% group_by(SampleID) %>%
  summarise(percent_assigned)
write.csv(ft.ct4.summary, "Output/featurecountsummary_batch4_240424.csv")

#add readcounts together across seq batches
ft.ct.plate2 <- merge(ft.ct2.summary, ft.ct4.summary, by = 'SampleID')
colnames(ft.ct.plate2) <- c("sampleNames", "featureCount_batch2", "featureCount_batch4")
ft.ct.plate2 <- ft.ct.plate2 %>%
  mutate(Mean.Features.Assigned = rowMeans(across(c(featureCount_batch2, featureCount_batch4))))


# edit the sampleNames in the plate2 to match, they have zero's in the single digit well names - better facilitates merging
plate2.samples$sampleNames <- gsub("0", "", as.character(plate2.samples$sampleNames)) 
to_repl <- which(substr(plate2.samples$lib, 2, 3) == '20')
for(i in to_repl){
  l = substring(plate2.samples[i,]$sampleNames,1,1)
  plate2.samples[i,]$sampleNames <- paste(l,"20", sep = "")
}
to_repl <- which(substr(plate2.samples$lib, 2, 3) == '10')
for(i in to_repl){
  l = substring(plate2.samples[i,]$sampleNames,1,1)
  plate2.samples[i,]$sampleNames <- paste(l,"10", sep = "")
}


#summarize by tissue type
plate2.all.data <- merge(ft.ct.plate2, plate2.samples, by='sampleNames')

plot1 <- ggbetweenstats(
  data = plate2.all.data,
  x = tissue,
  y = Mean.Features.Assigned)

plot1 <- plot1 + 
  # Add labels and title
  labs(
    x = "Tissue",
    y = "Mean % Features Assigned",
    title = "% Features Assigned across tissue types"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 8, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 20,
      face = "bold",
      color = "#2a475e"
    ),
    # Statistical annotations below the main title
    plot.subtitle = element_text(
      family = "Roboto", 
      size = 15, 
      face = "bold",
      color="#1b2838"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

plot1


### Plates 1 and 2 together
colnames(plate1.all.data)[2:4] <- c('PercAssigned_firstseq', 'PercAssigned_secondseq', 'Mean.Features.Assigned')
colnames(plate2.all.data)[2:4] <- c('PercAssigned_firstseq', 'PercAssigned_secondseq', 'Mean.Features.Assigned')
all.samples <- rbind(plate1.all.data, plate2.all.data)

#check here for additional params: https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html
plot1 <- ggbetweenstats(
  data = all.samples,
  x = tissue,
  y = Mean.Features.Assigned,
  results.subtitle = F,
  pairwise.display = "none",
  package = 'colorBlindness',
  palette = 'paletteMartin',
  centrality.plotting = F
  #violin.args = list(width = 0, linewidth = 0),
  #boxplot.args = list(width = 0)
)

plot1 <- plot1 + 
  # Add labels and title
  labs(
    x = "Tissue",
    y = "Mean % Features Assigned",
    title = "% Features Assigned across tissue types"
  ) + 
  # Customizations
  theme(
    # This is the new default font in the plot
    text = element_text(family = "Roboto", size = 5, color = "black"),
    plot.title = element_text(
      family = "Lobster Two", 
      size = 8,
      face = "bold",
      color = "#2a475e"
    ),
    plot.title.position = "plot", # slightly different from default
    axis.text = element_text(size = 10, color = "black"),
    axis.title = element_text(size = 12)
  )

plot1

#Figure out how to save these plots - right now I am just exporting them
ggsave(
  filename = here::here("img"),
  plot = plot1,
  width = 15,
  height = 10,
  device = pdf)



#TODO: compare batches - mapped reads across





#USE IF YOU WANT TO PLOT COUNTS OF UNIQUELY MAPPED READS
# file 2: mqc_star_alignment_plot_1.txt
# star.ct1 <- read.delim('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Sequencing_QC/batch_Z01_F001/mqc_star_alignment_plot_1_batch1.txt', header = T, sep = '\t')
# star.ct2 <- read.delim('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Sequencing_QC/batch_Z01_F002/mqc_star_alignment_plot_1_batch2.txt', header = T, sep = '\t')
# star.ct3 <- read.delim('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Sequencing_QC/batch_Z01_F003/mqc_star_alignment_plot_1_batch3.txt', header = T, sep = '\t')
# star.ct4 <- read.delim('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Sequencing_QC/batch_Z01_F004/mqc_star_alignment_plot_1_batch4.txt', header = T, sep = '\t')
# 

# ------------------------------------------------------------------
# Remove low perc mapped samples and save revised rawcount matrix and expt design
# ------------------------------------------------------------------
low.percmap.samples <- c('H19_1')
sampleTable <- sampleTable %>% filter(!lib %in% low.percmap.samples)

#save the edited Experiment Design file
write.csv(sampleTable, "Output/ExperimentDesign_allbatches_combined_v5.csv")

#Load the raw count matrix
countdata <- read.csv("Output/Counts_Atlas_allbatches_merged_v2.csv", row.names = 1)  
countdata <- countdata %>% select(order(colnames(countdata))) #alphabetize columns
countdata <- countdata[,!colnames(countdata) %in% low.percmap.samples]

#save the edited raw count matrix
write.csv(countdata, "Output/Counts_Atlas_allbatches_merged_v3.csv")

#Load the TPM count matrix
tpm_df <- read.csv("Output/TPM_Atlas_allbatches_merged_v2.csv", row.names = 1)
tpm_df <- tpm_df  %>% select(order(colnames(tpm_df))) #alphabetize columns
tpm_df <- tpm_df[, !colnames(tpm_df) %in% low.percmap.samples]

#save the edited TPM count matrix
write.csv(tpm_df, "Output/TPM_Atlas_allbatches_merged_v3.csv")

