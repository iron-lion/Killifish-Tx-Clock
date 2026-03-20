# Title: Atlas Hypergeometric GO analysis
# Author: Jingxun Chen
# Date: code compiled on 20240527
# Related publication: Emma K. Costa and Jingxun Chen, killifish atlas


# ------------------------------------------------------------------
# Set up 
# ------------------------------------------------------------------
indir = '/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/Clustering_Hierarch/output/2024_12_31_outs/GO_enrichment/Input/'

# Set wd to the current directory and seed (for reproducible results)

set.seed(1234)
setwd('/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/05_HClustering_GeneTrajectories/')

## This is to do GO enrichment analysis based on BH list from zebrafish/human
# *** IMPORTANT: Run it one by one for each of BP, MF and CC in the function #####
library("GOstats")
library('dplyr')
library('AnnotationDbi')

indir.geneset = '/Users/emmacosta/Library/CloudStorage/Dropbox/KillifishAtlasPaper/FirstSubmission/Code_Check/EC/04_Correlation/Output/Tables/'


# --------------------------------------------- 
# Load data - Single Tissue HClust - Ex) Brain
# ---------------------------------------------
tissue.list <- c('Brain', 'Liver', 'Gut', 'Muscle')

#debugging
#tissue <- 'Brain'

for(tissue in tissue.list){
  # brain universe - 10,847 genes expressed in all tissues
  universe <- read.csv(file = paste0(indir.geneset,'240905_geneexpressedinalltissues.csv'), row.names = 1)
  colnames(universe) = 'id'
  
  # Load clusters brain hClust
  clusters <- read.csv(paste0(indir,'Atlas_hclust_loessfits_clusterdata_',tissue,'_confinterval_241231.csv'))
  clust.id <- unique(clusters$Cluster)
  
  for(clust in clust.id){
    select.cluster <- clusters %>% filter(Cluster == clust)
    genes <- as.data.frame(select.cluster[,c(1)])
    
    # --------------------- Process the data ------------------------
    # clean up dataframe
    colnames(genes) <- 'id'
    
    # Add human name of the gene
    hSymbols = read.table("Supporting_files/NCBI-Human-orthologs.txt", head = T, sep = "\t")
    colnames(hSymbols) = c("Gene", "Human_name")
    
    # --------------------------------------------- 
    # Define background and genes of interest
    # ---------------------------------------------
    # List of universe genes. Background.
    # See before the for loop for definition of the universe
    
    # ---------------------- Select Output files ----------------------
    # Select the proper output file name given the specific data set, as well as upregulation or downregulation (condition of interest)  
    
    # Pick one of these options:
    outfilename = paste0("Output/GOBP_250101_Atlas_",tissue,"hClust",clust,".txt")
    
    # ---------------------- Select Output file with genes in GO terms ----------------------
    # Select the proper output file name to record the genes associated with each GO term.
    gotermlist = paste0("Output/GeneGOBP_240101_Atlas_",tissue,"hClust",clust,".txt")
    
    # --------------------------------------------- 
    # Hypergeometric test for GO analysis
    # ---------------------------------------------
    
    # --------------------------------------------- INPUT FILES --------------------------------------------------------
    # Go terms list. Use either zebrafish or human. Human works a bit better due to better annotations.
    frame = read.table(file ="Supporting_files/GO_killifish-human_best-hits.txt", header = T, colClasses=c(rep("factor",3)))
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
    library("Category")
    gsc <- GeneSetCollection(goAllFrame, setType = GOCollection())
    
    # Process the universe list
    #universe = universe$id
    universe = lapply(universe, as.character)
    universe = unlist(universe)
    head(universe)
    
    # Process the gene list of interest
    genes = genes$id
    genes = lapply(genes, as.character)
    genes = unlist(genes)
    head(genes)
    
    library('GOstats')
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
    
    # adjust p value for multiple correctionBiocManager::install("Category")
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
    
  }
}


sessionInfo() 
