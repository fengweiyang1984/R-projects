
setwd("C:/Users/z021w783/Desktop/GRA/Disertation/Real data example")

####Data Import, Cleaning and Preprocessing Steps 
####Import GTEx and TCGA count data
filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz")
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01))
#dim(filterGTEx01)
#[1] 19131     7
paraStudy = "GTEX" # Setting "GTEx" as the study of interest
paraPrimarySiteGTEx = "Breast" # Setting "Breast" as the primary site of interest
filterGTEx02 = subset(filterGTEx01,	study == paraStudy &  primarysite == paraPrimarySiteGTEx )
#dim(filterGTEx02)
#[1] 179   7

filterTCGA01 = fread("TCGA.BRCA.sampleMap_BRCA_clinicalMatrix")
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01))
#dim(filterTCGA01)
#[1] 1247     194

paraSampleType = "Primary Tumor" #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Breast" #Setting "Breast" as the primary site of interest.


filterTCGA02 = subset(filterTCGA01,	sampletype == paraSampleType &	primarysite == paraPrimarySiteTCGA) 
#dim(filterTCGA02)
#[1] 1101  194



##[Step 6-c] Merge GTEx and TCGA sample lists. Then pull expression profiles from .gz file by IDs on merged sample list.

filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample")
ExprSubsetBySamp =fread("TcgaTargetGtex_gene_expected_count.gz", select=filterExpr)

#dim(ExprSubsetBySamp)
#[1] 60498  1271



####[Main Text: Step 7] Subset expression data to include only protein coding genes

probemap =fread("zz_gencode.v23.annotation.csv", select =c(1,2))
exprALL = merge(probemap, ExprSubsetBySamp, by.x ="id", by.y = "sample")
genesPC = fread("zz_gene.protein.coding.csv")
#genesPC = read.csv("zz_gene.protein.coding.csv")
exprPC = subset(exprALL, gene %in% genesPC$Gene_Symbol)

# Remove duplicate gene symbols.
exprFinal = exprPC[!(duplicated(exprPC$gene)|
											duplicated(exprPC$gene, fromLast =TRUE)), ]

dim(exprFinal)
#[1] 18232  1272

