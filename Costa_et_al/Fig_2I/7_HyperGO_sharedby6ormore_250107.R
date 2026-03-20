# Title: Atlas Hypergeometric GO analysis
# Author: Jingxun Chen
# Date: code compiled on 20240527
# Related publication: Emma K. Costa and Jingxun Chen, killifish atlas


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------

# Set wd to the current directory and seed (for reproducible results)
setwd("/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Analysis/HyperGO/")
set.seed(1234)

## This is to do GO enrichment analysis based on BH list from zebrafish/human
# *** IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####
library("GOstats")
library('dplyr')



# --------------------------------------------- 
# Load data
# ---------------------------------------------
# Pick a dataset to run the script Leave the rest labeled with '#' sign. 

# Load age positively correlated genes shared by at least 6 tissues
#genes <- read.csv("Input/Atlas_POScorrgene_intersections_sharedby6ormore_alltissue_TPMcutoff_240714.csv")
genes <- read.csv("Input/Atlas_POScorrgene_intersections_sharedby5ormore_alltissue_TPMcutoff_241113.csv")


# # Load age negatively correlated genes shared by at least 6 tissues
# genes <- read.csv("Input/Atlas_NEGcorrgene_intersections_sharedby6ormore_alltissue_TPMcutoff_240714.csv")
#genes <- read.csv("Input/Atlas_NEGcorrgene_intersections_sharedby5ormore_alltissue_TPMcutoff_241113.csv")

# --------------------- Process the data ------------------------
# clean up dataframe
colnames(genes) <- 'id'

# Add human name of the gene
hSymbols = read.table("Input/NCBI-Human-orthologs.txt", head = T, sep = "\t")
colnames(hSymbols) = c("Gene", "Human_name")

# --------------------------------------------- 
# Define background and genes of interest
# ---------------------------------------------
# List of universe genes (removed padj = NA). Background.

# all tissue universe
universe <- read.csv("Input/240527_expressedgeneset-corranalysisinput_alltissues_union.csv")

colnames(universe) <- 'id'


# ---------------------- Select Output files ----------------------
# Select the proper output file name given the specific data set, as well as upregulation or downregulation (condition of interest)  

# Pick one of these options:
# outfilename = "Output/GOBP_240526_Atlas_POScorrgene_AtLeast6Tissues.txt" #enrichment for the age positively correlated genes shared by at least 6 tissues
# outfilename = "Output/GOBP_240526_Atlas_NEGcorrgene_AtLeast6Tissues.txt" #enrichment for the age negatively correlated genes shared by at least 6 tissues
outfilename = "Output/GOBP_Atlas_POScorrgene_AtLeast5Tissues.txt" #enrichment for the age positively correlated genes shared by at least 5 tissues
#outfilename = "Output/GOBP_Atlas_NEGcorrgene_AtLeast5Tissues.txt" #enrichment for the age positively correlated genes shared by at least 5 tissues


# ---------------------- Select Output file with genes in GO terms ----------------------
# Select the proper output file name to record the genes associated with each GO term.

# Pick one of these options:
#gotermlist = "Output/GeneGOBP_240526_Atlas_POScorrgene_AtLeast6Tissues.txt" #enrichment for the age positively correlated genes shared by at least 6 tissues
#gotermlist = "Output/GeneGOBP_240526_Atlas_NEGcorrgene_AtLeast6Tissues.txt" #enrichment for the age negatively correlated genes shared by at least 6 tissues
gotermlist = "Output/GeneGOBP_Atlas_POScorrgene_AtLeast5Tissues.txt" #enrichment for the age positively correlated genes shared by at least 5 tissues
#gotermlist = "Output/GeneGOBP_Atlas_NEGcorrgene_AtLeast5Tissues.txt" #enrichment for the age positively correlated genes shared by at least 5 tissues


# --------------------------------------------- 
# Hypergeometric test for GO analysis
# ---------------------------------------------

# --------------------------------------------- INPUT FILES --------------------------------------------------------
# Go terms list. Use either zebrafish or human. Human works a bit better due to better annotations.
frame = read.table(file ="GO_terms/GO_killifish-human_best-hits.txt", header = T, colClasses=c(rep("factor",3)))
# ontology MF, BP, CC
ontolg = "BP"

# Minimum number of genes for a term to filter
mingenes = 5 # Bare minimum is 2. More will get more general terms, less more specific. 5-10 is a good number.
# Relative enrichment filter
relenrich = 0 # I generally use 0 to get all terms

# ------------------------------------------------------------------------------------------------------------------

# This is just to get the 3 column. I already have these so I don't need it, still.
goframeData = data.frame(frame$go_id, frame$Evidence, frame$gene_id)

# put your data into a GOFrame object
goFrame=GOFrame(goframeData,organism="Human")
head(goframeData)

# cast this object to a GOAllFrame object will tap into the GO.db package and populate this object with the implicated GO2All mappings for you
goAllFrame=GOAllFrame(goFrame)

# generate geneSetCollection objects
library("GSEABase")
gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())

# Process the universe list
universe = universe$id
universe = lapply(universe, as.character)
universe = unlist(universe)
head(universe)

# Process the gene list of interest
genes = genes$id
genes = lapply(genes, as.character)
genes = unlist(genes)
head(genes)

params <- GSEAGOHyperGParams(name="My Custom GSEA based annotation parameters", 
                             geneSetCollection=gsc, 
                             geneIds = genes, 
                             universeGeneIds = universe, 
                             ontology = ontolg,
                             pvalueCutoff = 1, # 1 will get all terms, and then we can filter later
                             conditional = F, # To consider GO DAG structure or not. Doesn't affect much, but we can always try
                             testDirection = "over") # Default hyper geometric test gives both over and under enriched terms. I am specifying the direction by "over". Fo depleted terms use under.

# call hyperGTest to do the test
Over <- hyperGTest(params)
head(summary(Over))

# calculate enrichment and add it to data frame.
# Relative enrichment factor (E-value) for a GO term = (count/size)/(size/universe)
enrichment = (summary(Over)[5]$Count / summary(Over)[6]$Size) / (summary(Over)[6]$Size / length(universe))

# create a new frame
SummaryOver = data.frame(summary(Over), enrichment)
head(SummaryOver)

# Filter the Over variable on parameters other than P-value
# Filter the summary of OVER with size of the term, at least 2 genes for a go term
FilteredSummaryOver = SummaryOver[(SummaryOver$Count >= mingenes & SummaryOver$enrichment >= relenrich),]
head(FilteredSummaryOver)

# adjust p value for multiple correction
padj = p.adjust(FilteredSummaryOver$Pvalue, "BH")

# Add padj to the data frame
FinalSummaryOver = data.frame(FilteredSummaryOver, padj)

# write to a file
write.table(FinalSummaryOver, outfilename, quote = F, row.names = F, sep = "\t") # Final result file


# --------------------- To get the genes for each GO terms ---------------------------------------------

# isolate indexes for the go terms in final results
ind.GO <- is.element(names(Over@goDag@nodeData@data), eval(parse(text=paste("FinalSummaryOver$", "GO",ontolg,"ID", sep=''))))
selected.GO <- Over@goDag@nodeData@data[which(ind.GO)]

# get a go terms and genes in a new variable for all the terms in the results of enrichment
goTerms <- lapply(selected.GO, 
                  function(x) x$geneIds)
names(goTerms) <- names(Over@goDag@nodeData@data)[ind.GO]

# This will create a new file that will have GO terms and gene names in each GO terms.
# Number of GO terms or lines should be equal to the enriched go terms as in the other file generated by this script.
# This needs to be further processed to generate the desired files.

for (i in 1:length(goTerms)){
  
  test = as.data.frame(do.call(rbind, goTerms[i]))
  write.table(test, file = gotermlist, quote = F, col.names = F, append = T) # append each line - so make sure the file is empty for each run, or renamed after each run
  rm(test)
}


rm(list = ls())

sessionInfo() 
