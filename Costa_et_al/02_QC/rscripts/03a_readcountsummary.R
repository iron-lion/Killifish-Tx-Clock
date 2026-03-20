indir = "/Users/emmacosta/Library/CloudStorage/Box-Box/Emma Costa Externally Shareable Files/TWC Lab/03_Collaborations/Atlas/01_2023 Full Atlas/Analysis/Analysis/01_DESeq2_allbatches/02_QC/Output/"

rd.ct.plate1 = read.csv(paste0(indir,"readcountsummary_plate1_allseqbatches_240425.csv"))
rd.ct.plate2 = read.csv(paste0(indir,"readcountsummary_plate2_allseqbatches_240425.csv"))

sum(rd.ct.plate1$Total.Sample.Reads > 30000000) / nrow(rd.ct.plate1) #92.55%
sum(rd.ct.plate2$Total.Sample.Reads > 30000000) / nrow(rd.ct.plate2) #96.65%

(sum(rd.ct.plate1$Total.Sample.Reads > 30000000) + sum(rd.ct.plate2$Total.Sample.Reads > 30000000))/(nrow(rd.ct.plate1) + nrow(rd.ct.plate2)) #94.71%
(sum(rd.ct.plate1$Total.Sample.Reads > 35000000) + sum(rd.ct.plate2$Total.Sample.Reads > 35000000))/(nrow(rd.ct.plate1) + nrow(rd.ct.plate2)) #83.68%
