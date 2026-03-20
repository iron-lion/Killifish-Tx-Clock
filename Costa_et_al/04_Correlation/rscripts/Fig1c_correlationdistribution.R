library(viridis)
library(dplyr)
library(ggplot2)

load('/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Robjects/dds_BulkSeq_Aging_CorrelationResults_allsamples_Gonadcombo_TPMcutoff_240714.bin')

outdir = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Correlation/Output/'

#Bone, Brain, Eye, Fat, Gonad, Gut, Heart, Kidney, Liver, Muscle, Ovary, Skin, SpinalCord, Spleen, Testis
colTPS <- c('#d8413f', '#00a550','#b8b8c0', '#eee09b','#7961a2','#010101','#f0932e', '#fcd328', '#6cc0ee','#f4c489','#ab5673', '#f1a8a4','#ef9ac2') 


tissue.list <- names(BulkSeq_Agecorrelation_results_list)

BulkSeq_Agecorrelation_results_list[[tissue.list[1]]]

#---------------- plot all genes, including sig and non sig---------------- 
all.corres <- BulkSeq_Agecorrelation_results_list[[tissue.list[1]]]$resall 
all.corres$tissue <- tissue.list[1]

for(t in tissue.list[2:13]){
  tissue.res <- BulkSeq_Agecorrelation_results_list[[t]]$resall 
  tissue.res$tissue <- t
  all.corres <- rbind(tissue.res,all.corres)
}

all.corres$abs_spear <- abs(all.corres$cor_spear)


#scatter distribution plot
all.corres$alpha <- ifelse(abs(all.corres$cor_spear) > 0.5, 1.0, 0.05 )

#use this one
pdf(paste0(outdir,'Plots/alltissue_sexcombined_cor-spear_scatterplot_bytissue_flippedaxis_smallest.pdf'), width = 2, height = 3)
all.corres %>%
  ggplot(aes(group = tissue, y = tissue, x = cor_spear)) +
  geom_point(position = 'jitter', size = 0.01, aes(color = tissue, alpha = alpha)) + 
  scale_x_continuous(limits = c(-1,1)) +
  geom_vline(xintercept=0.5, linetype='longdash', col = 'red') +
  geom_vline(xintercept=-0.5, linetype='longdash', col = 'red') +
  scale_color_manual(values = colTPS) + 
  geom_boxplot(outlier.shape = NA, aes(alpha = 0.1)) +
  theme_classic() +
  theme(legend.position = "none", 
        axis.title.x=element_blank(),
        axis.title.y=element_blank())
dev.off()


#----------------only plot significant genes----------------
# rm(all.corres)
# all.corres <- BulkSeq_Agecorrelation_results_list[[tissue.list[1]]]$resall %>% filter(padj_spear < 0.05)
# all.corres$tissue <- tissue.list[1]
# 
# for(t in tissue.list[2:13]){
#   tissue.res <- BulkSeq_Agecorrelation_results_list[[t]]$resall %>% filter(padj_spear < 0.05)
#   tissue.res$tissue <- t
#   all.corres <- rbind(tissue.res,all.corres)
# }
# 
# all.corres$abs_spear <- abs(all.corres$cor_spear)
# 
# 
# #scatter distribution plot
# all.corres$alpha <- ifelse(abs(all.corres$cor_spear) > 0.5, 1.0, 0.1 )
# 
# #pdf(paste0(outdir,'Plots/alltissue_sexcombined_cor-spear_scatterplot_bytissue_flippedaxis.pdf'), width = 3, height = 3)
# all.corres %>%
#   ggplot(aes(group = tissue, y = tissue, x = cor_spear)) +
#   geom_point(position = 'jitter', size = 0.2, aes(color = tissue, alpha = alpha)) + 
#   scale_x_continuous(limits = c(-1,1)) +
#   geom_vline(xintercept=0.5, linetype='longdash', col = 'red') +
#   geom_vline(xintercept=-0.5, linetype='longdash', col = 'red') +
#   scale_color_manual(values = colTPS) + 
#   geom_boxplot(outlier.shape = NA, aes(alpha = 0.1)) +
#   theme_classic() +
#   theme(legend.position = "none", 
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# #dev.off()
# 
# #negatively assoc.
# all.corres %>% filter(cor_spear < 0) %>%
#   ggplot(aes(group = tissue, y = tissue, x = cor_spear)) +
#   geom_point(position = 'jitter', size = 0.2, aes(color = tissue, alpha = alpha)) + 
#   scale_x_continuous(limits = c(-1,-0.25)) +
#   geom_vline(xintercept=-0.5, linetype='longdash', col = 'red') +
#   scale_color_manual(values = colTPS) + 
#   geom_boxplot(outlier.shape = NA, aes(alpha = 0.1)) +
#   theme_classic() +
#   theme(legend.position = "none", 
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# 
# #positively assoc.
# all.corres %>% filter(cor_spear > 0) %>%
#   ggplot(aes(group = tissue, y = tissue, x = cor_spear)) +
#   geom_point(position = 'jitter', size = 0.2, aes(color = tissue, alpha = alpha)) + 
#   scale_x_continuous(limits = c(0.25,1)) +
#   #geom_vline(xintercept=0.5, linetype='longdash', col = 'red') +
#   scale_color_manual(values = colTPS) + 
#   geom_boxplot(outlier.shape = NA, aes(alpha = 0.1)) +
#   theme_classic() +
#   theme(legend.position = "none", 
#         axis.title.x=element_blank(),
#         axis.title.y=element_blank())
# 
# #---------------- plot density plot ---------------- 
