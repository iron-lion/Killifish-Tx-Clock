rm(list=ls())
# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
library(ggplot2)
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/')

indir = '/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/02_QC/Input/'

# Set wd to the current directory
cts.batch1 <- read.csv(paste0(indir,'Counts_Atlas_Plate02_240319.csv'))
cts.batch3 <- read.csv(paste0(indir,'Counts_Atlas_Plate02_Lane02_BatchF004_240417.csv'))

col_order <- rownames(cts.batch1)
cts.batch3 <- cts.batch3[col_order,]

cts.batch3 <- cts.batch3[,colnames(cts.batch1)]
colnames(cts.batch1) == colnames(cts.batch3)


#make a dataframe with x = gene counts batch 1, y = gene counts batch 2, z = gene

cts.batch1$gene <- rownames(cts.batch1)
cts.batch3$gene <- rownames(cts.batch3)

 
df1 <- cbind(cts.batch1$gene, stack(cts.batch1[1:358]))
df1$gene <- df1$`cts.batch1$gene`

df3 <- cbind(cts.batch3$gene, stack(cts.batch3[1:358]))
df3$gene <- df3$`cts.batch3$gene`


df.combo <- data.frame(gene = df1$gene, x = df1$values, y = df3$values)

plot(x = df.combo$x,y = df.combo$y)

ggplot(df.combo, aes(x=x, y=y)) + 
  geom_point( color="#69b3a2") +
  xlim(0,4000000) + 
  ylim(0,4000000) + 
  geom_abline(aes(intercept = 0, slope =2))


ggplot(df.combo, aes(x=x, y=y)) + 
  geom_point( color="#69b3a2") +
  xlim(0,4000000) + 
  ylim(0,4000000)
