# Set wd to the current directory
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')

# ------------------------------------------------------------------
# Merge Raw counts
# ------------------------------------------------------------------

# Data was sequenced in 2 orders, the first consisted of batch 1 and 2 which were 1 lane each of sequencing
# Batch 3 and 4 were for plates 01 and 02 respectively, but batch 3 was 1 lane of sequencing and batch 4 was 2 lanes (merged fastqs were used for mapping and featureCount)

### Merge count matrices
cts.batch3 <- read.csv("Input/Counts_Atlas_Plate01_Lane02_BatchF003_240407.csv")
cts.batch4 <- read.csv("Input/Counts_Atlas_Plate02_Lane02_BatchF004_240417.csv")
cts.batch.1.2 <- read.csv("Input/Counts_Atlas_Plate01-02_240319.csv")

colnames(cts.batch3) <- paste(colnames(cts.batch3), "1", sep = "_")
colnames(cts.batch4) <- paste(colnames(cts.batch4), "2", sep = "_")

#merge batch 3 and 4
cts.batch.3.4 <- merge(cts.batch3, cts.batch4, by = 'row.names', all = T)
rownames(cts.batch.3.4) <- cts.batch.3.4$Row.names
cts.batch.3.4 = subset(cts.batch.3.4, select = -c(Row.names))

cts.batch.1.2 <- cts.batch.1.2 %>% select(order(colnames(cts.batch.1.2)))
cts.batch.3.4 <- cts.batch.3.4 %>% select(order(colnames(cts.batch.3.4)))
colnames(cts.batch.1.2) == colnames(cts.batch.3.4)

row_order <- rownames(cts.batch.3.4)
cts.batch.1.2 <- cts.batch.1.2[row_order,]
table(rownames(cts.batch.1.2) == rownames(cts.batch.3.4))

#then merge the two big matrices
cts_df <- as.data.frame(as.matrix(cts.batch.1.2) + as.matrix(cts.batch.3.4))

# edit the sampleNames in the plate2 to match, they have zero's in the single digit well names - better facilitates merging
plate2.samples <- as.data.frame(colnames(cts_df))
plate2.samples$sampleNames = sapply(strsplit(as.character(plate2.samples$`colnames(cts_df)`), "_"), `[`, 1)
plate2.samples$plate = sapply(strsplit(as.character(plate2.samples$`colnames(cts_df)`), "_"), `[`, 2)
plate2.samples$lib <- plate2.samples$`colnames(cts_df)`

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
 
plate2.samples$lib_rename <- paste0(plate2.samples$sampleNames, "_", plate2.samples$plate)
colnames(cts_df) <- plate2.samples$lib_rename

#save the merged matrix
write.csv(cts_df, "Output/Counts_Atlas_allbatches_merged.csv")

# ------------------------------------------------------------------
# Merge TPM
# ------------------------------------------------------------------
# Data was sequenced in 2 orders, the first consisted of batch 1 and 2 which were 1 lane each of sequencing
# Batch 3 and 4 were for plates 01 and 02 respectively, but batch 3 was 1 lane of sequencing and batch 4 was 2 lanes (merged fastqs were used for mapping and featureCount)

### Merge TPM matrices
tpm.batch3 <- read.csv("Input/TPM_Atlas_Plate01_Lane02_BatchF003_240407.csv")
tpm.batch4 <- read.csv("Input/TPM_Atlas_Plate02_Lane03-4_BatchF004_240417.csv")
tpm.batch.1.2 <- read.csv("Input/TPM_Atlas_Plate01-02_240319.csv")

colnames(tpm.batch3) <- paste(colnames(tpm.batch3), "1", sep = "_")
colnames(tpm.batch4) <- paste(colnames(tpm.batch4), "2", sep = "_")

#merge batch 3 and 4
tpm.batch.3.4 <- merge(tpm.batch3, tpm.batch4, by = 'row.names', all = T)
rownames(tpm.batch.3.4) <- tpm.batch.3.4$Row.names
tpm.batch.3.4 = subset(tpm.batch.3.4, select = -c(Row.names))

tpm.batch.1.2 <- tpm.batch.1.2 %>% select(order(colnames(tpm.batch.1.2)))
tpm.batch.3.4 <- tpm.batch.3.4 %>% select(order(colnames(tpm.batch.3.4)))
colnames(tpm.batch.1.2) == colnames(tpm.batch.3.4)

row_order <- rownames(tpm.batch.3.4)
tpm.batch.1.2 <- tpm.batch.1.2[row_order,]
table(rownames(tpm.batch.1.2) == rownames(tpm.batch.3.4))

#then merge the two big matrices
tpm_df <- as.data.frame(as.matrix(tpm.batch.1.2) + as.matrix(tpm.batch.3.4))


# edit the sampleNames in the plate2 to match, they have zero's in the single digit well names - better facilitates merging
plate2.samples <- as.data.frame(colnames(tpm_df))
plate2.samples$sampleNames = sapply(strsplit(as.character(plate2.samples$`colnames(tpm_df)`), "_"), `[`, 1)
plate2.samples$plate = sapply(strsplit(as.character(plate2.samples$`colnames(tpm_df)`), "_"), `[`, 2)
plate2.samples$lib <- plate2.samples$`colnames(tpm_df)`

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

plate2.samples$lib_rename <- paste0(plate2.samples$sampleNames, "_", plate2.samples$plate)
colnames(tpm_df) <- plate2.samples$lib_rename


#save the merged matrix
write.csv(tpm_df, "Output/TPM_Atlas_allbatches_merged.csv")


# ------------------------------------------------------------------
# Experimental Design, Censoring animals
# ------------------------------------------------------------------

# Load experiment design and count matrices
sampleTable <- read.csv("Input/240430_alltissuescollected.csv")  

#remove the samples that were not included in sequenicng library
no.lib <- which(is.na(sampleTable$X))
sampleTable <- sampleTable[-c(no.lib),]

row.names(sampleTable) <- as.character(sampleTable$lib)
sampleTable <- subset(sampleTable, select = -c(X))

sampleTable$RNA_batch <- paste(sampleTable$tissue, sampleTable$RNA_batch, sep = '_') #Fix the RNA batch column
sampleTable <- sampleTable[order(sampleTable$lib),] #alphabetize rows
sampleTable <- sampleTable %>%
  rename(
    sampleNames = sampleName
  )

#Load the raw count matrix
countdata <- read.csv("Output/Counts_Atlas_allbatches_merged.csv", row.names = 1)  
countdata <- countdata %>% select(order(colnames(countdata))) #alphabetize columns

coldata <- DataFrame(sampleTable)
col_order <- rownames(coldata)

#remove the censored samples that were not in the sequencing lib
thing1 <- colnames(countdata)
thing2 <- rownames(coldata)
censor_list <- setdiff(union(thing1, thing2), intersect(thing1, thing2))

sampleTable <- sampleTable %>% filter(!lib %in% censor_list)
coldata <- DataFrame(sampleTable)
col_order <- rownames(coldata)

countdata <- countdata[,col_order]

table(colnames(countdata) == rownames(coldata)) # Make sure that the column names are identical

#save the edited Experiment Design file
write.csv(sampleTable, "Output/ExperimentDesign_allbatches_combined_v3.csv")
