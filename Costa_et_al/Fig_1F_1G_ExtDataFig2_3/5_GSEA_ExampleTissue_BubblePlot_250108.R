# Title: Plot GSEA results (Figure 1F)
# Author: Jingxun Chen
# Date: code compiled on 20241215
# Related publication: Emma K. Costa and Jingxun Chen et al., Nature Aging, 2025


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/GSEA/")
set.seed(1234)

# Load packages
library('ggplot2')
library('dplyr')

path = "/Users/jingxunchen/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/GSEA/"


# ------------------------------------------------------------------
# Load data
# ------------------------------------------------------------------
# Define tissue types
# t = 'Liver'
# t = 'Brain'
t = 'Muscle'
# t = 'Gut'

# Load Male & Female GOGSEA csv list
Input_m <- read.csv(paste0(path, 'Output/241209_allCorrresults', t,'_maleOnly_GOGSEA.csv'), header = T)
Input_f <- read.csv(paste0(path, 'Output/241209_allCorrresults', t,'_femaleOnly_GOGSEA.csv'), header = T)


# ------------------------------------------------------------------
# Clean up, Filter, Sort the input data
# ------------------------------------------------------------------

# Extract only the NES and padj terms
Input_m <- Input_m[c('ID', 'Description', 'NES', 'p.adjust')]
Input_f <- Input_f[c('ID','NES', 'p.adjust')]

# Rename columns to keep track of Male vs. Female data
names(Input_m) <- c('ID', 'Description', 'Male_NES', 'Male_p.adjust')
names(Input_f) <- c('ID', 'Female_NES', 'Female_p.adjust')


# Make a large dataframe
data <- full_join(Input_m, Input_f, by = 'ID')


# ------------------------- Define GO terms of interest -----------------------------------------
# GO terms that are significantly altered by age in both males and females
# Pick the GO terms associated with each tissue to run the code with. For example, if you want to plot liver GO terms, select the first 'go' definition 

# Liver GO terms
go <- c('GO:0002253', 'GO:0006959', 'GO:0006364', 'GO:0042254', 'GO:0006260', 'GO:0051306', 'GO:0030198')
# GO:0002253 activation of immune response
# GO:0006959 humoral immune response
# GO:0006364 rRNA processing
# GO:0042254 ribosome biogenesis
# GO:0006260 DNA replication
# GO:0051306 mitotic sister chromatid separation
# GO:0030198 extracellular matrix organization

# Brain GO terms
go <- c('GO:0045087', 'GO:0006954','GO:0044391', 'GO:0006260', 'GO:0051306',  'GO:0007409', 'GO:0030900')
# GO:0045087 innate immune response
# GO:0006954 inflammatory response
# GO:0044391 ribosomal subunit
# GO:0006260 DNA replication
# GO:0051306 mitotic sister chromatid separation
# GO:0007409 axonogenesis
# GO:0030900 forebrain development

# Muscle GO terms
go <- c('GO:0004842', 'GO:0005776', 'GO:0016236', 'GO:0044391', 'GO:0001568', 'GO:0043062', 'GO:0001503', 'GO:0044786')
# GO:0004842 ubiquitin-protein transferase activity
# GO:0005776 autophagosome
# GO:0016236 macroautophagy
# GO:0044391 ribosomal subunit
# GO:0001568 blood vessel development
# GO:0043062 extracellular structure organization
# GO:0001503 ossification
# GO:0044786 cell cycle DNA replication

# Gut GO terms
go <- c('GO:0006111', 'GO:0043255', 'GO:0005759', 'GO:0098798', 'GO:0044843', 'GO:0098644', 'GO:0044391')
# GO:0006111 regulation of gluconeogenesis
# GO:0043255 regulation of carbohydrate biosynthetic process
# GO:0005759 mitochondrial matrix
# GO:0098798 mitochondrial protein-containing complex
# GO:0044843 cell cycle G1/S phase transition
# GO:0098644 complex of collagen trimers
# GO:0044391 ribosomal subunit


# ------------------------------------------------------------------
# Getting the data ready for plotting (common for all tissues)
# ------------------------------------------------------------------

# Keep only the terms of interest
data_go <- filter(data, data$ID %in% go)

# Make a long list of the data for plotting
data_go_f <- subset(data_go, select = -c(Male_NES, Male_p.adjust))
data_go_m <- subset(data_go, select = -c(Female_NES, Female_p.adjust))
names(data_go_f) <- c('ID', 'Description', 'NES', 'p.adjust')
names(data_go_m) <- c('ID', 'Description', 'NES', 'p.adjust')

# Keep track of the Condition
data_go_f['Condition']<- 'Female'
data_go_m['Condition']<- 'Male'

# Match the order so that ID = my GO list 
df <- data_go_f[match(go, data_go_f$ID), ]
dm <- data_go_m[match(go, data_go_m$ID), ]

# Make the final data frame for plotting
FinalData <- rbind(df, dm)

# Change the level of 'Description' based on the order I defined here.
order <- df$Description
FinalData$Description <- factor(FinalData$Description, levels=rev(order))

# ------------------------------------------------------------------
# Plot the bubble plot and save as PDF
# ------------------------------------------------------------------

pdf(paste0('Plots/', t,'_GSEA_GO_MvsF_selectedGOterms_241215.pdf'), width = 8, height = 10)
ggplot(data=FinalData, aes(x=Condition, y=Description)) +
  geom_point(aes(color=NES, size = -log10(p.adjust))) +
  scale_colour_gradient2(low = "blue", mid = "purple", high = "red", midpoint = 0, space = "Lab",na.value = "grey50", guide = "colourbar", aesthetics = "colour") +
  labs(title = paste0(t, ' GSEA by GO terms for males vs. females'), x='Conditions', y = '',size='-log(FDR)', color='NES')
dev.off()


sessionInfo() 
