
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/')


#Sex: male, female 
colsex <-c('#068ec9','#ba1e2d')

#First, we load our Deseq2 object as it contains all the counts and metadata
object.indir = "/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/03_DESeq2/Output/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))


dds_CA1_list <- dds_TPS_list # there was an instance where the list of deseq objects was called dds_CA1_list. 

#Load if you are starting from here
load('/robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')



#### Get rho ####
# if you want to output the rho-value for a particular gene 
gene <- 'LOC107383908'
tissue <- 'Kidney'
cor_res <- as.data.frame(BulkSeq_Agecorrelation_results_list[[tissue]]$resall)
rownames(cor_res) <- cor_res$gene
rvalue <- cor_res[gene,]$cor_spear
print(rvalue)



### Fig 4G ###
##### plotting gene expression with age - irf4a
gene <- 'LOC107383908'
tissue <- 'Kidney'
#We extract the respective Deseq2 object
dds_tissue <- dds_CA1_list[[tissue]]

#We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
count_table <- counts(dds_tissue, normalized=T)
count_table <- as.data.frame(count_table)

temp <- t(count_table[gene,])
df <- cbind(temp, dds_tissue@colData)

#scatterplot
df$age_days <- as.numeric(as.character(df$age_days))

pdf(paste0(outdir,"Plots/LOC107383908-irf4_normcounts_Kidney_NEGcorrgene_240527.pdf"), width = 2, height = 2)
ggplot(df, aes(x=age_days, y=LOC107383908)) +
  geom_point(stroke = NA, shape=21, size = 2,
             aes(fill = factor(sex))) + 
  scale_fill_manual(values=rev(colsex))+
  ggtitle('LOC107383908 (norm counts) in Kidney')+
  #geom_smooth(method=lm , color = 'red',fill="#69b3a2",se=T) +
  scale_x_continuous(breaks = seq(40, 160, by = 20), limits = c(45,165)) +
  scale_y_continuous(breaks = seq(0,1000, by = 250), limits = c(0,1000)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.text=element_text(size=10))
dev.off()




### ExtData Fig 6F ###

##### plotting a gene trajectory with age - irf4b
gene <- 'irf4'
tissue <- 'Kidney'
#We extract the respective Deseq2 object
dds_tissue <- dds_CA1_list[[tissue]]

#We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
count_table <- counts(dds_tissue, normalized=T)
count_table <- as.data.frame(count_table)

temp <- t(count_table[gene,])
df <- cbind(temp, dds_tissue@colData)

#scatterplot
df$age_days <- as.numeric(as.character(df$age_days))

pdf(paste0(outdir,"Plots/irf4b_normcounts_Kidney_NEGcorrgene_241024.pdf"), width = 2, height = 2)
ggplot(df, aes(x=age_days, y=irf4)) +
  geom_point(stroke = NA, shape=21, size = 2,
             aes(fill = factor(sex))) + 
  scale_fill_manual(values=rev(colsex))+
  ggtitle('irf4b (norm counts) in Kidney')+
  #geom_smooth(method=lm , color = 'red',fill="#69b3a2",se=T) +
  scale_x_continuous(breaks = seq(40, 160, by = 20), limits = c(45,165)) +
  scale_y_continuous(breaks = seq(0,250, by = 50), limits = c(0,250)) +
  theme_classic() + 
  theme(legend.position = 'none', axis.text=element_text(size=10))
dev.off()