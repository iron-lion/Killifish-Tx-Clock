setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')

library('tidyverse')
library('readxl')
library('base')
library(stringr)
library('ggstatsplot')

# Load experiment design and count matrices
sampleTable <- read.csv("Input/ExperimentDesign_allbatches_combined_v3.csv", row.names = 1)  #this will have been moved from the Output folder
lib_name <- as.character(sampleTable$lib)

plate1.samples <- subset(sampleTable, cDNA_batch == 1)
plate2.samples <- subset(sampleTable, cDNA_batch == 2)

# ------------------------------------------------------------------
# Load read count summaries - original from multiqc
# ------------------------------------------------------------------
# Input data description: multiqc run on fastqc from trimming?

# file 1: mqc_fastqc_sequence_counts_plot_1.txt
rd.ct1 <- read.delim("Input/mqc_fastqc_sequence_counts_plot_1_Batch1.txt", header = T, sep = '\t')
rd.ct2 <- read.delim('Input/mqc_fastqc_sequence_counts_plot_1_batch2.txt', header = T, sep = '\t')
rd.ct3 <- read.delim('Input/mqc_fastqc_sequence_counts_plot_1 _batch3.txt', header = T, sep = '\t')
rd.ct4 <- read.delim('Input/mqc_fastqc_sequence_counts_plot_1_batch4.txt', header = T, sep = '\t')


### Plate 1
rd.ct1$SampleID <- NA
rd.ct1$SampleID <- str_split_fixed(rd.ct1$Sample, "_",2)[,1]
rd.ct1$Total.Reads <- NA
rd.ct1 <- rd.ct1 %>%
  mutate(Total.Reads = rowSums(across(c(Unique.Reads, Duplicate.Reads))))
rd.ct1.summary <- rd.ct1 %>% group_by(SampleID) %>%
  summarise(sum(Total.Reads))
rd.ct1.summary <- rd.ct1.summary[-c(323),]
write.csv(rd.ct1.summary, "Output/readcountsummary_batch1_240425.csv")

rd.ct3$SampleID <- NA
rd.ct3$SampleID <- str_split_fixed(rd.ct3$Sample, "_",2)[,1]
rd.ct3$Total.Reads <- NA
rd.ct3 <- rd.ct3 %>%
  mutate(Total.Reads = rowSums(across(c(Unique.Reads, Duplicate.Reads))))
rd.ct3.summary <- rd.ct3 %>% group_by(SampleID) %>%
  summarise(sum(Total.Reads))
write.csv(rd.ct3.summary, "Output/readcountsummary_batch3_240422.csv")

#add readcounts together across seq batches
rd.ct.plate1 <- merge(rd.ct1.summary, rd.ct3.summary, by = 'SampleID')
colnames(rd.ct.plate1) <- c("sampleNames", "ReadCount_batch1", "ReadCount_batch3")
rd.ct.plate1 <- rd.ct.plate1 %>%
  mutate(Total.Sample.Reads = rowSums(across(c(ReadCount_batch1, ReadCount_batch3))))
write.csv(rd.ct.plate1, "Output/readcountsummary_plate1_allseqbatches_240425.csv")

#summarize by tissue type
plate1.all.data <- merge(rd.ct.plate1, plate1.samples, by='sampleNames')

plot1 <- ggbetweenstats(
  data = plate1.all.data,
  x = tissue,
  y = Total.Sample.Reads)

plot1 <- plot1 + 
  # Add labels and title
  labs(
    x = "Tissue",
    y = "# of Paired Reads per sample (R1 + R2)",
    title = "Distribution of reads counts (paired) across tissue types"
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
rd.ct2$SampleID <- NA
rd.ct2$SampleID <- str_split_fixed(rd.ct2$Sample, "_",2)[,1]
rd.ct2$Total.Reads <- NA
rd.ct2 <- rd.ct2 %>%
  mutate(Total.Reads = rowSums(across(c(Unique.Reads, Duplicate.Reads))))
rd.ct2.summary <- rd.ct2 %>% group_by(SampleID) %>%
  summarise(sum(Total.Reads))
write.csv(rd.ct2.summary, "Output/readcountsummary_batch2_240425.csv")

rd.ct4$SampleID <- NA
rd.ct4$SampleID <- str_split_fixed(rd.ct4$Sample, "_",2)[,1]
rd.ct4$Total.Reads <- NA
rd.ct4 <- rd.ct4 %>%
  mutate(Total.Reads = rowSums(across(c(Unique.Reads, Duplicate.Reads))))
rd.ct4.summary <- rd.ct4 %>% group_by(SampleID) %>%
  summarise(sum(Total.Reads))
write.csv(rd.ct4.summary, "Output/readcountsummary_batch4_240425.csv")



#add readcounts together across seq batches
rd.ct.plate2 <- merge(rd.ct2.summary, rd.ct4.summary, by = 'SampleID')
colnames(rd.ct.plate2) <- c("sampleNames", "ReadCount_batch2", "ReadCount_batch4")
rd.ct.plate2 <- rd.ct.plate2 %>%
  mutate(Total.Sample.Reads = rowSums(across(c(ReadCount_batch2, ReadCount_batch4))))
write.csv(rd.ct.plate2, "Output/readcountsummary_plate2_allseqbatches_240425.csv")

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
plate2.all.data <- merge(rd.ct.plate2, plate2.samples, by = 'sampleNames')

plot1 <- ggbetweenstats(
  data = plate2.all.data,
  x = tissue,
  y = Total.Sample.Reads)

plot1 <- plot1 + 
  # Add labels and title
  labs(
    x = "Tissue",
    y = "# of Paired Reads per sample (R1 + R2)",
    title = "Distribution of reads counts (paired) across tissue types"
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
colnames(plate1.all.data)[2:4] <- c('ReadCount_firstseq', 'ReadCount_secondseq', 'Total.Sample.Reads')
colnames(plate2.all.data)[2:4] <- c('ReadCount_firstseq', 'ReadCount_secondseq', 'Total.Sample.Reads')
all.samples <- rbind(plate1.all.data, plate2.all.data)

#check here for additional params: https://indrajeetpatil.github.io/ggstatsplot/reference/ggbetweenstats.html
plot1 <- ggbetweenstats(
  data = all.samples,
  x = tissue,
  y = Total.Sample.Reads,
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
    y = "# of Paired Reads per sample (R1 + R2)",
    title = "Distribution of reads counts (paired) across tissue types"
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


#TODO: compare batches - read counts across

# ------------------------------------------------------------------
# Remove low readcount samples and save revised rawcount matrix and expt design
# ------------------------------------------------------------------
low.readcount.samples <- c('L21_1', 'J6_1')
sampleTable <- sampleTable %>% filter(!lib %in% low.readcount.samples)

#save the edited Experiment Design file
write.csv(sampleTable, "Output/ExperimentDesign_allbatches_combined_v4.csv")

#Load the raw count matrix
countdata <- read.csv("Output/Counts_Atlas_allbatches_merged.csv", row.names = 1)  
countdata <- countdata %>% select(order(colnames(countdata))) #alphabetize columns
countdata <- countdata[,!colnames(countdata) %in% low.readcount.samples]

#save the edited raw count matrix
write.csv(countdata, "Output/Counts_Atlas_allbatches_merged_v2.csv")

#Load the TPM count matrix
tpm_df <- read.csv("Output/TPM_Atlas_allbatches_merged.csv", row.names = 1)
tpm_df <- tpm_df  %>% select(order(colnames(tpm_df))) #alphabetize columns
tpm_df <- tpm_df[, !colnames(tpm_df) %in% low.readcount.samples]

#save the edited TPM count matrix
write.csv(tpm_df, "Output/TPM_Atlas_allbatches_merged_v2.csv")
