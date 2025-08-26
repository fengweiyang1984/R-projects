setwd("/Ovary Cancer/")

####Related to Data Import, Cleaning and Preprocessing Steps 1-10


##Load R packages:

library(UCSCXenaTools)
library(data.table)
library(R.utils)
library(dplyr)


data(XenaData)
write.csv(XenaData, "00_tblXenaHubInfo.csv")


####[Main Text: Step 5] Select then download target data sets from Xena Data Hubs.
##[Step 5-a] Target=RSEM expected counts provided by the UCSC toil Recompute Compendium

getOption('timeout')
options(timeout=10000)

GeneExpectedCnt_toil =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TcgaTargetGtex_gene_expected_count")
XenaQuery(GeneExpectedCnt_toil) %>%
	XenaDownload(destdir ="./")   ###, max_try = 1L)



##[Step 5-b] Target =TCGA Clinical data.

paraCohort = "TCGA Ovarian Cancer"
paraDatasets ="TCGA.OV.sampleMap/OV_clinicalMatrix"


Clin_TCGA = XenaGenerate(subset = XenaHostNames =="tcgaHub") %>%
	XenaFilter(filterCohorts = paraCohort) %>%
	XenaFilter(filterDatasets = paraDatasets)
XenaQuery(Clin_TCGA) %>%
	XenaDownload(destdir = "./")


##[Step 5-c] Target = TCGA Survival Data.

Surv_TCGA =XenaGenerate(subset = XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TCGA_survival_data")
XenaQuery(Surv_TCGA) %>%
 	XenaDownload(destdir ="./")



##[Step 5-d] Target = GTEx Phenotype data.

Pheno_GTEx = XenaGenerate(subset =XenaHostNames == "toilHub") %>%
	XenaFilter(filterCohorts = "TCGA TARGET GTEx") %>%
	XenaFilter(filterDatasets = "TcgaTargetGTEX_phenotype")
XenaQuery(Pheno_GTEx) %>%
	XenaDownload(destdir ="./")


#### [Main Text: Step 6] Subset expression data to include only desired samples.

## [Step 6-a] Retrive IDs for GTEx normal samples of desired tissue type(s).

filterGTEx01 = fread("TcgaTargetGTEX_phenotype.txt.gz")
names(filterGTEx01) = gsub("\\_", "", names(filterGTEx01))

#dim(filterGTEx01)
#[1] 19131     7

             sample detailedcategory primary disease or tissue primarysite
 1: TCGA-V4-A9EE-01   Uveal Melanoma            Uveal Melanoma         Eye
 2: TCGA-VD-AA8N-01   Uveal Melanoma            Uveal Melanoma         Eye
 3: TCGA-V4-A9EI-01   Uveal Melanoma            Uveal Melanoma         Eye
 4: TCGA-VD-AA8O-01   Uveal Melanoma            Uveal Melanoma         Eye
 5: TCGA-WC-A888-01   Uveal Melanoma            Uveal Melanoma         Eye
 6: TCGA-WC-A881-01   Uveal Melanoma            Uveal Melanoma         Eye
 7: TCGA-WC-A88A-01   Uveal Melanoma            Uveal Melanoma         Eye
 8: TCGA-YZ-A980-01   Uveal Melanoma            Uveal Melanoma         Eye
 9: TCGA-V4-A9EO-01   Uveal Melanoma            Uveal Melanoma         Eye
10: TCGA-WC-A87U-01   Uveal Melanoma            Uveal Melanoma         Eye
       sampletype gender study
 1: Primary Tumor   Male  TCGA
 2: Primary Tumor   Male  TCGA
 3: Primary Tumor   Male  TCGA
 4: Primary Tumor   Male  TCGA
 5: Primary Tumor   Male  TCGA
 6: Primary Tumor   Male  TCGA
 7: Primary Tumor   Male  TCGA
 8: Primary Tumor   Male  TCGA
 9: Primary Tumor   Male  TCGA
10: Primary Tumor   Male  TCGA


paraStudy = "GTEX" # Setting "GTEx" as the study of interest
paraPrimarySiteGTEx = "Ovary" # Setting "Ovary" as the primary site of interest

filterGTEx02 = subset(filterGTEx01,	study == paraStudy &  primarysite == paraPrimarySiteGTEx )

#dim(filterGTEx02)
#[1] 88  7



## [Step 6-b] Retrive IDs for TCGA primary tumor samples of desired histological type(s).

filterTCGA01 = fread(paraDatasets)
names(filterTCGA01) = gsub("\\_", "", names(filterTCGA01))

#dim(filterTCGA01)
#[1] 630 102
#which(colnames(filterTCGA01)=="sampletype")
#[1] 64
#which(colnames(filterTCGA01)=="primarysite")
#[1] 14
#which(colnames(filterTCGA01)=="histologicaltype")
#[1] 35


table(filterTCGA01$primarysite)
#Ovary 
#  630 


table(filterTCGA01$histologicaltype)
#                         Serous Cystadenocarcinoma 
#                      23                       607 


paraSampleType = "Primary Tumor" #Setting "Primary Tumor" as the sample type of interest.
paraPrimarySiteTCGA = "Ovary" #Setting "Ovary" as the primary site of interest.


filterTCGA02 = subset(filterTCGA01,	sampletype == paraSampleType &	primarysite == paraPrimarySiteTCGA) 
#dim(filterTCGA02)
#[1] 595 102



##[Step 6-c] Merge GTEx and TCGA sample lists. Then pull expression profiles from .gz file by IDs on merged sample list.

filterExpr = c(filterGTEx02$sample, filterTCGA02$sampleID, "sample")

ExprSubsetBySamp =fread("TcgaTargetGtex_gene_expected_count.gz", select=filterExpr)

#dim(ExprSubsetBySamp)
[1] 60498   508


################# Did not apply first. Use all genes to predict cases & controls
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


##################################################################################################
##### case-control classification analysis
##################################################################################################
##################################################################################################
#### single gene based DE analysis
#### 


sample_vec <- sapply(colnames(ExprSubsetBySamp), function(x) unlist(strsplit(x, "-"))[[1]][1]) 
table(sample_vec)
#sample_vec
#  GTEX sample   TCGA 
#    88      1    419 

# ExprSubsetBySamp[1:10, 503:508]
#    TCGA-61-2113-01 TCGA-OY-A56P-01 TCGA-OY-A56Q-01 TCGA-VG-A8LO-01
# 1:          1.0000          1.0000          0.0000          2.0000
# 2:          0.0000          0.0000          0.0000          0.0000
# 3:          2.8074          4.5236          0.0000          2.5850
# 4:         10.8002          9.5593         10.0221          9.9524
# 5:          0.0000          0.0000          0.0000          0.0000
# 6:          9.7549          8.6294          7.8765          8.9542
# 7:          8.6110          7.8329          7.1293         10.0169
# 8:          0.0000          0.0000          0.0000          0.0000
# 9:         13.5557         13.7352         13.6557         13.5191
#10:          5.9773          8.2527          5.7549          7.7415
#    TCGA-WR-A838-01             sample
# 1:          0.0000  ENSG00000242268.2
# 2:          0.0000  ENSG00000259041.1
# 3:          3.5850  ENSG00000270112.3
# 4:          9.5780 ENSG00000167578.16
# 5:          0.0000  ENSG00000278814.1
# 6:         10.8564  ENSG00000078237.5
# 7:          9.5736  ENSG00000269416.5
# 8:          0.0000  ENSG00000263642.1
# 9:         12.7722 ENSG00000146083.11
#10:          7.5622 ENSG00000158486.13


caseControlResponse <- rep(0, (dim(ExprSubsetBySamp)[2]-1))
caseControlResponse[which(sample_vec == "TCGA")] <- 1


ExprDD <- as.data.frame(ExprSubsetBySamp[,1:(dim(ExprSubsetBySamp)[2]-1)])
rownames(ExprDD) <- ExprSubsetBySamp$sample
dim(ExprDD)
[1] 60498   507


library(edgeR)



group <- caseControlResponse
y <- DGEList(counts=ExprDD,group=group)

keep <- filterByExpr(y)
y <- y[keep,,keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
design <- model.matrix(~group)
y <- estimateDisp(y,design)


fit <- glmQLFit(y,design)
qlf <- glmQLFTest(fit,coef=2)
topTags(qlf)


qlfT <- qlf$table
qlfT_sort <- qlfT[order(qlfT[,4]),]

qlfT_sort[1:10,]


Coefficient:  group 
                      logFC   logCPM        F        PValue           FDR
ENSG00000165215.6  2.576739 6.910046 3879.728 9.359050e-241 1.088270e-236
ENSG00000145103.12 3.200718 6.389086 3081.751 1.520418e-209 8.839712e-206
ENSG00000129354.11 2.460658 6.699637 3366.828 1.444872e-205 5.600323e-202
ENSG00000129455.15 2.945666 6.614035 2531.920 4.526949e-200 1.315984e-196
ENSG00000125850.10 3.463232 6.380057 2499.381 7.059642e-199 1.641790e-195
ENSG00000169035.11 2.429482 6.785169 2419.268 6.947620e-196 1.346449e-192
ENSG00000157765.11 2.124225 6.925950 2384.322 1.489970e-194 2.475052e-191
ENSG00000110195.11 2.655180 6.813235 2233.782 1.254092e-188 1.822823e-185
ENSG00000267795.5  2.686970 6.524976 2226.657 2.436423e-188 3.147859e-185
ENSG00000130768.14 2.457193 6.577795 3490.926 2.658978e-186 3.091859e-183




##########################################################################################
###########################################################################################
##### detect predictive gene networks using top 200 DE genes as hub genes
#####


first10000genes <- qlfT_sort[1:10000,]
mid <- match(rownames(first10000genes), rownames(ExprDD))
ExprDD10000 <- ExprDD[mid,]

dim(ExprDD10000)
#[1] 10000   507



source("library_CIS_imagingData_12222021.R")

geneName <- row.names(ExprDD10000)

X <- as.matrix(t(ExprDD10000))
Xsd <- apply(X, 2, sd)
length(which(Xsd==0))
#[1] 0

Xnew <- X #X[,-which(Xsd==0)]
geneName.new <- geneName[-which(Xsd==0)]

Xnew <- apply(Xnew, 2, function(x) (x-mean(x))/sd(x))
colnames(Xnew) <- geneName.new

### Need to add in covariates such as age, gender,etc.

Y <- as.vector(caseControlResponse)


system.time(try2  <- mLDA.pair(Xnew, Y, Xnew, Z=NULL, Z.new=NULL, pair=c(1,0), tau=50, alpha=0.9, nu=2000, d=2, nb=5))

pred <- as.numeric(try2$PredClass)

length(which(pred!=Y))








