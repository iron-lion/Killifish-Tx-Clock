library(DESeq2)
library(PCAtools)
library(dplyr)
library(tidyr)
library(ggbreak)

# Set wd to the current directory
setwd()

# Load the DESeq2 object list
object.indir = "Path/to/DESeq2_Robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))

metadata <- as.data.frame(colData(dds_TPS_list$All))

metadata$age_bin <- ifelse(metadata$age_days %in% c('47', '49', '52'), '1', NA)
metadata$age_bin <- ifelse(metadata$age_days %in% c('75', '77', '78'), '2', metadata$age_bin)
metadata$age_bin <- ifelse(metadata$age_days %in% c('102', '103'), '3', metadata$age_bin)
metadata$age_bin <- ifelse(metadata$age_days %in% c('133', '134'), '4', metadata$age_bin)
metadata$age_bin <- ifelse(metadata$age_days %in% c('133', '134'), '4', metadata$age_bin)
metadata$age_bin <- ifelse(metadata$age_days %in% c('147', '152', '155'), '5', metadata$age_bin)
metadata$age_bin <- ifelse(metadata$age_days %in% c('161', '162'), '6', metadata$age_bin)
#metadata$age_bin <- as.factor(metadata$age_bin)

metadata$tissue_age_bin_sex <- paste(metadata$tissue, metadata$age_bin, metadata$sex, sep = '_')
metadata$tissue_sex <- paste(metadata$tissue,metadata$sex, sep = '_')


View(table(metadata$tissue_age_bin_sex))

View(table(metadata$tissue_sex))
