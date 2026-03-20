# Title:  Get a large count table to be used as input for DEseq2 analysis
# Author: Emma Costa
# Date: code compiled on 20250108
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
# Set wd to the current directory
# setwd("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/01_GetCounts/")
setwd("/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/01_GetCounts/")

# Set home directory 
# homeDir = as.character("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/01_GetCounts")
homeDir = as.character("/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/01_GetCounts")



# Please create the following folder if not present: 'Output'

# ---------------  Input variables ---------------------------
genename = "gene_name"
sampleTable <- read.csv(paste0(homeDir, "/Input/ExperimentDesign_batch_Z01-F001_F002_combined_update.csv"), stringsAsFactors = F)

sampleTable$lib <- paste(sampleTable$sampleNames, sampleTable$seq_batch, sep = "_")
write.csv(sampleTable, paste0(homeDir, "/Input/ExperimentDesign_batch_Z01-F001_F002_combined_update.csv"))


rownames(sampleTable) <- sampleTable$lib


# ------------------------------------------------------------------
# Define functions for rpkm and tpm
# ------------------------------------------------------------------
## reference: https://gist.github.com/slowkow/c6ab0348747f86e2748b#file-counts_to_tpm-r-L44
## reference: https://www.biostars.org/p/171766/

rpkm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(counts) * 1e9
}

tpm <- function(counts, lengths) {
  rate <- counts / lengths
  rate / sum(rate) * 1e6
}


                    
# ------------------------------------------------------------------
# Loop to get the counts
# ------------------------------------------------------------------
for (i in 1:length(sampleTable$file)){
  
  outfilename <- file.path(homeDir, "/Output/", paste0(sampleTable$lib[i], "_counts.csv")) # Change this to the path of output file
  print (sampleTable$file[i])
  outfilename
  cts <- read.table(paste0(homeDir, "/featureCount/", sampleTable$file[i]), header = TRUE)
  colnames(cts)[7] <- paste("Counts")
  colnames(cts)[1] <- paste("GeneID")
  x = cts[c("GeneID", "Length", "Counts")]
  
  # get tpm, this is not normalized across samples
  x$TPM = tpm(x$Counts, x$Length)  
  write.table(x,file=outfilename,quote=FALSE,sep=",",row.names=FALSE)

  rm(outfilename, cts, x) # remove the temp variables for the next run

}

sessionInfo() 

