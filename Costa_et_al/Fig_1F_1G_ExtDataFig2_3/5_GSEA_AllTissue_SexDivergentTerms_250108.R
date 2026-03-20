# Title: Bubble plots for selected GO terms that show opposite directions of changes with age in males and females (sex-divergent terms)
# Author: Jingxun Chen
# Date: code compiled on 20241215
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025
# Associated figures: Extended Data Fig. 3


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
setwd("/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/GSEA/")
set.seed(1234)

library(dplyr)
library(ggplot2)
library(tidyr)

plot_outs = '/Users/jingxun/Dropbox/BrunetLab/Data/AtlasData/GSEA/Plots/'

# Input GSEA sex-divergent opposite terms
GSEA_all_sig <- read.csv('Output/Alltissues_sexsplit_corres_GOGSEA.csv', row.names = 1)

#------------------------------------------------------------------
# Select GO terms
# ------------------------------------------------------------------

IDs.byhallmarks <- c( 'GO:0019439','GO:0044272','GO:0001676','GO:0061620', 'GO:0006111',  # metabolism
                      'GO:0006986', 'GO:0045047','GO:0045862','GO:0006515', 'GO:0000502', #proteostasis
                      'GO:0016570', 'GO:0031124', 'GO:0034708', 'GO:0017053','GO:0070076','GO:0003712','GO:0016572', 'GO:0036260', #epigenetic alterations and mRNA expression
                      'GO:0022626','GO:0044391','GO:0030684','GO:0042254', 'GO:0006364', #ribosome
                      'GO:0005776','GO:0006914','GO:0005765','GO:0007033', #autophagy & lysosome
                      'GO:0042554', 'GO:0003954', 'GO:0022900', 'GO:0005746', 'GO:0019674', #mito dysfunction
                      'GO:0032570','GO:1902074','GO:1903305','GO:0030667', 'GO:0019838','GO:0030072','GO:0006865','GO:0006813','GO:0005261','GO:0046879', # transport & intercellular communication
                      'GO:0000819','GO:0045787','GO:0000922', #cell divisions
                      'GO:1903053', 'GO:0030198', 'GO:0030199', #ECM
                      'GO:0034340','GO:0032606','GO:0002444','GO:0070555', 'GO:0042119', 'GO:0036230'#inflammation
)


# Check if the selected GO terms are in the data
IDs.byhallmarks %in% GSEA_all_sig$ID 

# Make a dataframe with the selected GO terms
GSEA_all_sig_hallmark <- GSEA_all_sig %>% filter(ID %in% IDs.byhallmarks) 

# Set order of the dataframe 
order.df <- data.frame(IDs.byhallmarks) # set order based on hallmark
colnames(order.df) = c('GOBPID') # name the column 
rownames(order.df) = order.df$GOBPID # set row names to be the GOBPID

# Add a column called hallmark that denotes the different types of GO terms
order.df$hallmark = c('metabolism', 'metabolism','metabolism', 'metabolism','metabolism',
                      'proteostasis','proteostasis','proteostasis','proteostasis','proteostasis',
                      'epigenetics','epigenetics','epigenetics','epigenetics','epigenetics','epigenetics','epigenetics','epigenetics',
                      'ribosome','ribosome','ribosome','ribosome', 'ribosome',
                      'autophagy','autophagy','autophagy','autophagy',
                      'mitochondria', 'mitochondria','mitochondria','mitochondria','mitochondria',
                      'transport','transport','transport','transport','transport','transport','transport','transport','transport','transport',
                      'cellCycle','cellCycle','cellCycle',
                      'ECM','ECM','ECM',
                      'inflammation', 'inflammation', 'inflammation','inflammation', 'inflammation','inflammation')

# Add a new column in the GSEA data to keep track of hallmark types
GSEA_all_sig_hallmark$hallmark <- NA

# Loop to add the hallmark type to each GO term
for(i in 1:nrow(GSEA_all_sig_hallmark)){
  row.temp = GSEA_all_sig_hallmark[i,]   # temporary storage for each row
  ID = row.temp$ID   # find the ID of this row
  h = order.df[ID,]$hallmark   # find the hallmark type for this ID
  GSEA_all_sig_hallmark[i,]$hallmark = h   # record the hallmark type on the GSEA data
}


# Order the dataframe based on hallmark
GSEA_all_sig_hallmark <- GSEA_all_sig_hallmark[order(GSEA_all_sig_hallmark$hallmark), ]
GSEA_all_sig_hallmark$hallmark <- factor(GSEA_all_sig_hallmark$hallmark, levels = sort(unique(GSEA_all_sig_hallmark$hallmark )))
GSEA_all_sig_hallmark$Description <- factor(GSEA_all_sig_hallmark$Description, levels = unique(GSEA_all_sig_hallmark$Description))

#------------------------------------------------------------------
# Prepare to plot
# ------------------------------------------------------------------
GSEA_all_sig_hallmark$Condition = paste(GSEA_all_sig_hallmark$sex, GSEA_all_sig_hallmark$tissue, sep = "_")

# The code below is used to set the scale of dot size by including the lowest and highest p.adjust values for male and female
# *** Minimum p.adjust value = 4.859712e-09 from female
# *** Maximum p.adjust value = 1.0000000 from both male and female

# Minimum p.adjust value = 4.859712e-09. Create a row that include this info
min1 <- as.data.frame(c("BP", "GO", "maleScale", 0, 0, 0, 0, 4.859712e-09, 0, 0, 0, 'None', 'None', "Male", "MinScale", "Scale", "Male_Scale"))
min2 <- as.data.frame(c("BP", "GO", "maleScale", 0, 0, 0, 0, 4.859712e-09, 0, 0, 0, 'None', 'None', "Female", "MinScale", "Scale", "Female_Scale"))
min1.t <- t(min1) # transpose rows and columns
min2.t <- t(min2) # transpose rows and columns
rownames(min1.t) <- 1403 # set index name
rownames(min2.t) <- 1404 # set index name
colnames(min1.t) <- colnames(GSEA_all_sig_hallmark) # set the column names to be the same as the data
colnames(min2.t) <- colnames(GSEA_all_sig_hallmark) # set the column names to be the same as the data
min <- rbind(min1.t, min2.t)

# Maximum p.adjust value = 1.0000000. Create a row that include this info
max1 <- as.data.frame(c("BP", "GO", "maleScale", 0, 0, 0, 0, 1, 0, 0, 0, 'None', 'None', "Male", "MaxScale", "Scale", "Male_Scale"))
max2 <- as.data.frame(c("BP", "GO", "maleScale", 0, 0, 0, 0, 1, 0, 0, 0, 'None', 'None', "Female", "MaxScale", "Scale", "Female_Scale"))
max1.t <- t(max1)
max2.t <- t(max2)
rownames(max1.t) <- 1405
rownames(max2.t) <- 1406
colnames(max1.t) <- colnames(GSEA_all_sig_hallmark)
colnames(max2.t) <- colnames(GSEA_all_sig_hallmark)
max <- rbind(max1.t, max2.t)

# Create a data frame called scale for the min and max values 
scale <- as.data.frame(rbind(min, max))

# Set the NES and p.adjust terms to numeric for plotting
scale$NES <- as.numeric(scale$NES) 
scale$p.adjust <- as.numeric(scale$p.adjust)

# Combine the scale and GSEA data 
GSEA_all_sig_hallmark_scale <- rbind(GSEA_all_sig_hallmark, scale)

# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

# Results for females
pdf(paste0(plot_outs, '/Alltissue_GSEA_GO_F_oppositeGOterms_241229_1.pdf'), width = 15, height = 10)
ggplot(data=GSEA_all_sig_hallmark_scale %>% filter(sex == 'Female'), aes(x=Condition, y=Description)) +
  geom_point(aes(color=NES, size = -log10(p.adjust))) +
  scale_size(range = c(0, 5)) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  labs(title = 'GSEA by Aging Hallmark sex-opposite GO terms for females across tissue', x='Conditions', y = '',size='-log(FDR)', color='NES')
dev.off()

# Results for males
pdf(paste0(plot_outs, '/Alltissue_GSEA_GO_M_oppositeGOterms_241229_1.pdf'), width = 15, height = 10)
ggplot(data=GSEA_all_sig_hallmark_scale %>% filter(sex == 'Male'), aes(x=Condition, y=Description)) +
  geom_point(aes(color=NES, size = -log10(p.adjust))) +
  scale_size(range = c(0, 5)) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  labs(title = 'GSEA by Aging Hallmark GO terms for males across tissue', x='Conditions', y = '',size='-log(FDR)', color='NES')
dev.off()

# Note: In Extended Data 3, the male column and female column are placed side by side for each tissue. 
# Manual editing in Illustrator (by moving the female columns to proper location) is used to generate this side-by-side configuration. 


rm(list = ls())

# ------------------------------------------------------------------
sessionInfo() 
