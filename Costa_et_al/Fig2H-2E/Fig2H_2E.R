library(dplyr)

#First, we load our list of DeSeq2 objects as it contains all the counts and metadata
object.indir = "path/to/robjects/"
load(paste0(object.indir,'dds_TPS_allsamples_Gonadcombo_240714.bin'))
dds_CA1_list <- dds_TPS_list
rm(dds_TPS_list)

#instead of the list of DESeq2 object list you may also load a norm count data table and the experimental data table
#load norm counts here, and modify code as desired below if that's what you'd like to do

#colors, listed from youngest age bin to oldest age bin
colage <- c('#395ea4','#3789c0','#dfeaf0','#f3d4ac','#f7b77b', '#a92b46')

##### plotting a gene trajectory with age - for example IGF2BP3
gene <- 'LOC107383282'
tissue <- 'Gut'
#We extract the respective Deseq2 object 
dds_tissue <- dds_CA1_list[[tissue]]

#We use the DESE2-normalized data, which we extract by getting the counts from the Deseq2 object
count_table <- counts(dds_tissue, normalized=T)
count_table <- as.data.frame(count_table)

temp <- t(count_table[gene,])

dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('47', '49', '52'), '1', NA)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('75', '77', '78'), '2', dds_tissue@colData$age_bin)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('102', '103'), '3', dds_tissue@colData$age_bin)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('133', '134'), '4', dds_tissue@colData$age_bin)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('133', '134'), '4', dds_tissue@colData$age_bin)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('147', '152', '155'), '5', dds_tissue@colData$age_bin)
dds_tissue@colData$age_bin <- ifelse(dds_tissue@colData$age_days %in% c('161', '162'), '6', dds_tissue@colData$age_bin)

df <- cbind(temp, dds_tissue@colData)
df <- as.data.frame(df)

# Calculate mean and standard error by sex and age_bin
df_summary <- df %>%
  group_by(sex, age_bin) %>%
  summarize(
    LOC107383282_mean = mean(LOC107383282, na.rm = TRUE),
    LOC107383282_se = sd(LOC107383282, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

# Join the summary back to the original data frame
df <- left_join(df, df_summary, by = c("sex", "age_bin"))

###########

ggplot() + 
  # Use df_summary for the bars and error bars
  geom_col(data = df_summary, aes(x = factor(sex), y = LOC107383282_mean, 
                                  fill = factor(age_bin)), 
           width = 1, position = "dodge") +
  geom_errorbar(data = df_summary, 
                aes(x = factor(sex), 
                    y = LOC107383282_mean, 
                    ymin = LOC107383282_mean - LOC107383282_se, 
                    ymax = LOC107383282_mean + LOC107383282_se,
                    group = factor(age_bin)), 
                width = 0.2, 
                position = position_dodge(width = 1)) +
  # Use the original df for the dot plot
  geom_dotplot(data = df, 
               aes(x = factor(sex), 
                   y = LOC107383282, 
                   group = interaction(factor(sex), factor(age_bin))), 
               position = position_dodge(width = 1), 
               binaxis = "y", 
               dotsize = 0.25,
               stackdir = "center", 
               stackratio = 1) +
  scale_fill_manual(values = colage) +
  theme_classic() +
  theme(legend.position = 'none')