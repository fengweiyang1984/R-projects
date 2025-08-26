#### The following codes are for creating a response variable of immune infiltration levels. 

### Reference paper: "Immune Inﬁltration Subtypes Characterization and Identiﬁcation of Prognosis-Related lncRNAs in Adenocarcinoma of the Esophagogastric Junction".

### Read in RNAseq expression data. I am using a GEO expression data. You should use the TCGA RNAseq expression data for this.
rm(list=ls())
setwd("C:/Users/z021w783/Desktop/GRA/Other paper/Breast Cancer/Datasets")
a <- read.table('GSE112996_merged_fpkm_table.txt.gz',
                header = T,
                row.names=1)
raw_data<- a[,-1]

#dim(raw_data)
#[1] 57773    45

#raw_data[1:5, 1:5]
#                X160003842T.3 X160003809T.2 X160003842T.2 X160003809T.3
#ENSG00000160293      2.008373      9.444742      3.712596      8.503404
#ENSG00000222293      0.000000      0.000000      0.000000      0.000000
#ENSG00000165410      1.831875      7.688678      3.402515      8.192393
#ENSG00000168803      2.313177      0.763094      3.128527      1.750566
#ENSG00000198890      0.000000      0.000000      4.431164      1.602184
#                X160003805T.3
#ENSG00000160293      9.323715
#ENSG00000222293      0.000000
#ENSG00000165410      8.090352
#ENSG00000168803      0.929858
#ENSG00000198890      0.000000


###Phenotype data. This is only for extracting subject information (IDs) to be matched with the subjects in the RNAseq data.
pheno <- read.csv(file = 'GSE112996_series_matrix.txt')
pheno <- data.frame(num1 = strsplit(as.character(pheno[42,]),split='\t')[[1]][-1],
                    num2 = gsub('patient: No.','P',strsplit(as.character(pheno[51,]),split='\t')[[1]][-1]))
                    
#pheno[1:5,]
#          num1 num2
#1 160003803T-1 P002
#2 160003803T-2 P002
#3 160003803T-3 P002
#4 160003803T-4 P002
#5 160003803T-5 P002


{
####Data filtering and normalization

data<- a[!apply(raw_data,1,sum)==0,]

data$median=apply(data[,-1],1,median)
data=data[order(data$GeneName,data$median,decreasing = T),]
data=data[!duplicated(data$GeneName),]
rownames(data)=data$GeneName

uni_matrix <- data[,grep('\\d+',colnames(data))]
uni_matrix <- log2(uni_matrix+1)
colnames(uni_matrix)<- gsub('X','',gsub('\\.','\\-',colnames(uni_matrix)))
uni_matrix<- uni_matrix[,order(colnames(uni_matrix))]
}
save(uni_matrix,pheno,file = 'uni_matrix.Rdata')




### You will need to install the following R packages. It is a bit hurdling when installing some of them. Please see the commends below.

{
library(genefilter)
library(GSVA)
library(Biobase)
library(stringr)
}

###
##if (!require("BiocManager", quietly = TRUE))
##    install.packages("BiocManager")

##BiocManager::install("genefilter")

##BiocManager::install("GSVA")

##install.packages("remotes")
##library(remotes)
##install_github("rcastelo/GSVA")

## load the Rdata
#load('uni_matrix.Rdata')




### Read in the marker genes for the immune cells
gene_set<- read.csv('mmc3.csv',header = T)
gene_set<-gene_set[, 1:2]
head(gene_set)
#  Metagene        Cell.type
#1   ADAM28 Activated B cell
#2    CD180 Activated B cell
#3    CD79B Activated B cell
#4      BLK Activated B cell
#5     CD19 Activated B cell
#6    MS4A1 Activated B cell



list<- split(as.matrix(gene_set)[,1], gene_set[,2])

#$`Activated B cell`
# [1] "ADAM28"   "CD180"    "CD79B"    "BLK"      "CD19"     "MS4A1"   
# [7] "TNFRSF17" "IGHM"     "GNG7"     "MICAL3"   "SPIB"     "HLA-DOB" 
#[13] "IGKC"     "PNOC"     "FCRL2"    "BACH2"    "CR2"      "TCL1A"   
#[19] "AKNA"     "ARHGAP25" "CCL21"    "CD27"     "CD38"     "CLEC17A" 
#[25] "CLEC9A"   "CLECL1"  
#
#$`Activated CD4 T cell`
# [1] "AIM2"   "BIRC3"  "BRIP1"  "CCL20"  "CCL4"   "CCL5"   "CCNB1"  "CCR7"  
# [9] "DUSP2"  "ESCO2"  "ETS1"   "EXO1"   "EXOC6"  "IARS"   "ITK"    "KIF11" 
#[17] "KNTC1"  "NUF2"   "PRC1"   "PSAT1"  "RGS1"   "RTKN2"  "SAMSN1" "SELL"  
#[25] "TRAT1" 


### Accessing Estimates GSVA enrichment scores.
gsva_matrix<- gsva(as.matrix(uni_matrix), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)

write.table(gsva_matrix, file="gsva_result_matrix.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

#gsva_matrix[1:5, 1:5]
#                               160003803T-1 160003803T-2 160003803T-3
#Activated B cell                  0.6737891    0.4632534    0.4661965
#Activated CD4 T cell              0.9779582    0.8309930    0.8877659
#Activated CD8 T cell              0.8953367    0.7369096    0.7792591
#Activated dendritic cell          0.7397454    0.6747484    0.6563200
#CD56bright natural killer cell    0.8970348    0.9174155    0.9041718
#                               160003803T-4 160003803T-5
#Activated B cell                  0.6571843    0.5014955
#Activated CD4 T cell              0.9650610    0.8940602
#Activated CD8 T cell              0.9124224    0.8674855
#Activated dendritic cell          0.7054462    0.6475732
#CD56bright natural killer cell    0.9289161    0.8922549



##########################################################
##########################################################
####

### Using ConsensusClusterPlus to cluster immune cell score into high-infiltration and low infiltraction levels

### BiocManager::install("ConsensusClusterPlus")


 library(ConsensusClusterPlus)
 
 d = sweep(gsva_matrix,1, apply(gsva_matrix,1,median,na.rm=T))

results = ConsensusClusterPlus(d,maxK=6,reps=50, pItem=0.8, pFeature=1, clusterAlg="hc",distance="pearson",seed=1262118388.71279,plot="png")

results[[2]][["consensusClass"]]

#160003803T-1 160003803T-2 160003803T-3 160003803T-4 160003803T-5 160003803T-6 
#           1            1            1            1            1            1 
#160003805T-2 160003805T-3 160003805T-4 160003805T-5 160003807T-2 160003807T-3 
#           2            2            2            2            1            1 
#160003807T-4 160003807T-5 160003807T-6 160003809T-1 160003809T-2 160003809T-3 
#           1            2            1            2            2            2 
#160003809T-4 160003809T-5 160003814T-1 160003814T-3 160003814T-4 160003814T-5 
#           2            2            2            2            2            2 
#160003818T-1 160003818T-2 160003818T-3 160003820T-1 160003820T-2 160003820T-3 
#           1            1            1            1            1            1 
#160003824T-1 160003824T-2 160003824T-3 160003828T-1 160003828T-2 160003828T-3 
#           2            2            2            2            2            2 
#160003830T-3 160003830T-4 160003830T-5 160003838T-1 160003838T-2 160003838T-3 
#           1            2            1            2            2            1 
#160003842T-1 160003842T-2 160003842T-3 
#           1            1            1 



