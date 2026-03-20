# ------------------------------------------------------------------
# Analysis
# ------------------------------------------------------------------
setwd("/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/DESeq2_bysex/")
outdir.deseq = 'Output/'

load('Robjects/dds_TPS_alltissue_byagebin-sexdiffs_DESeq2results_241122.bin')
tissue.list <- names(results_list) #no eye bc there wasnt enough sex balance
tissue.list <- tissue.list[-3]

#---------------- Explore the number of genes for each  ----------------
sex.terms <- c("M_vs_F", "F_vs_M")

sex.degs <- data.frame(tissue = NA)
sex.degs$term <- NA
sex.degs$pos.deg.number <- NA
sex.degs$neg.deg.number <- NA
sex.degs$age_bin <- NA

for(t in tissue.list){
  tissue.res <- results_list[[t]]
  for(i in c('1','2', '3', '4', '5')){
    tissue.age.res <- tissue.res[[i]] 
    for(int in sex.terms){
      res.sig = tissue.age.res[[int]]$ressig
      deg.num.down = length(which(res.sig$log2FoldChange < 0))
      deg.num.up = length(which(res.sig$log2FoldChange > 0))
      sex.degs[nrow(sex.degs)+1,] = c(t,int,deg.num.up,deg.num.down,i)
    }
  }
}

sex.degs = sex.degs[-1,]

#Bone, Brain, Fat, Gonad,Gut, Heart, Kidney, Liver, Muscle, Skin, SpinalCord, Spleen
colTPS <- c('#d8413f', '#00a550','#eee09b','#7962A3' ,'#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489' ,'#ab5673', '#f1a8a4','#ef9ac2') 

df.sex = sex.degs

df.sex$term <- factor(df.sex$term, levels = sex.terms)
df.sex$age_bin <- factor(df.sex$age_bin, levels = c('1', '2', '3', '4', '5'))
df.sex$pos.deg.number <- as.numeric(df.sex$pos.deg.number)
df.sex$neg.deg.number <- as.numeric(df.sex$neg.deg.number)

# Plotting
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = pos.deg.number, group = tissue, color = tissue)) +
  geom_line() +
  geom_point() +
  labs(x = "Age Bin", y = "Positive DEG Number", title = "Positive DEG Numbers by Age Bin for Each Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_color_manual(values = colTPS)

# Plotting
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = neg.deg.number, group = tissue, color = tissue)) +
  geom_line() +
  geom_point() +
  labs(x = "Age Bin", y = "Negative DEG Number", title = "Negative DEG Numbers by Age Bin for Each Tissue") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = 'none') +
  scale_color_manual(values = colTPS)


#express as % expressed 
outdir.cor = "/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Correlation/Output/"
num_Expressed <- read.csv(file = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/Jiachen_EngelhardtLab/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv')
#num_Expressed <- read.csv(paste0(outdir.cor,"Tables/240714_numexpressedgenes_alltissues_Gonadcombo_TPMcutoff.csv")) #made in correlation_checkconfounds.R

df.sex = merge(df.sex, num_Expressed)
df.sex$tot.deg.number = df.sex$pos.deg.number + df.sex$neg.deg.number
df.sex$tot.deg.number.perc = ((df.sex$tot.deg.number) / (df.sex$total.expressed))*100

df.sex$pos.deg.number.perc = ((df.sex$pos.deg.number) / (df.sex$total.expressed))*100
df.sex$neg.deg.number.perc = ((df.sex$neg.deg.number) / (df.sex$total.expressed))*100

df.sex.MvF <- df.sex %>% filter(term == 'M_vs_F')
df.sex.MvF$X = NULL


df.sex.gonad <- df.sex.MvF %>% filter(tissue == 'Gonad')
mean(df.sex.gonad$tot.deg.number.perc) #94.58%

df.sex.liver <- df.sex.MvF %>% filter(tissue == 'Liver')
mean(df.sex.liver$tot.deg.number.perc) #24.81%

df.sex.skin <- df.sex.MvF %>% filter(tissue == 'Skin')
mean(df.sex.skin$tot.deg.number.perc) #14.52%

df.sex.kidney <- df.sex.MvF %>% filter(tissue == 'Kidney')
mean(df.sex.kidney$tot.deg.number.perc) #14.31%

df.sex.fat <- df.sex.MvF %>% filter(tissue == 'Fat')
mean(df.sex.fat$tot.deg.number.perc) #5.74%

##### Fig 1E:
pdf(paste0(outdir.deseq, 'Plots/241227_totdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = tot.deg.number.perc, group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "DEG Number % of Transcriptome", title = "DEG by Age Bin for Each Tissue: M vs F") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,110)) +
  scale_y_break(c(40, 80)) 
dev.off()


##### Fig 1E: zoomed in
pdf(paste0(outdir.deseq, 'Plots/241227_totdegs-bytissue_percexpr_eachagebin_M_vs_F_zoomedIn.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = tot.deg.number.perc, group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "DEG Number % of Transcriptome", title = "DEG by Age Bin for Each Tissue: M vs F") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,40)) 
dev.off()


#pdf(paste0(outdir.deseq, 'Plots/241122_posdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = pos.deg.number.perc, group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "Positive DEG Number % of Transcriptome", title = "Positive DEG Numbers by Age Bin for Each Tissue") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,55)) +
  scale_y_break(c(20, 40)) 
#dev.off()

# Plotting
#pdf(paste0(outdir.deseq, 'Plots/241122_negdegs-bytissue_percexpr_eachagebin_M_vs_F.pdf'), width = 2, height = 2)
df.sex %>% filter(term == 'M_vs_F') %>%
  ggplot(aes(x = age_bin, y = (neg.deg.number.perc), group = tissue, color = tissue)) +
  geom_line(alpha = 0.5) +
  geom_point(alpha = 0.5, size = 0.5) +
  labs(x = "Age Bin", y = "Negative DEG Number % of Transcriptome", title = "Negative DEG Numbers by Age Bin for Each Tissue") +
  theme_classic() +
  theme(legend.position = 'none',axis.text.x = element_blank(), axis.title.y = element_blank(), title = element_blank(),
        axis.title.x = element_blank()) +
  scale_color_manual(values = colTPS) +
  ylim(c(0,55)) +
  scale_y_break(c(20, 40))  
#dev.off()

#save the output table
write.csv(df.sex, file = paste0(outdir.deseq, '241218_M_vs_F_DEGs_byage.csv'))

