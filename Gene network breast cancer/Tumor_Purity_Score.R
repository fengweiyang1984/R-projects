### Reference paper: "Immune Inﬁltration Characteristics and a Gene Prognostic Signature Associated With the Immune Inﬁltration in Head and Neck Squamous Cell Carcinoma."


library(ggsci)
library(tidyr)
library(ggpubr)
library(utils) 
#rforge <- "http://r-forge.r-project.org"
#install.packages("estimate", repos=rforge, dependencies=TRUE)
library(estimate)
library(tidyverse)

#Readin the tumor patient expression file
expr <- read.table("GSE112996_merged_fpkm_table.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)

expr1 <- expr[!duplicated(expr[,1]),]
expr2 <- expr1[,-1]
rownames(expr2) <- expr1[,1]
write.table(expr2, file="GSE112996_fpkm_liyanming.txt", row.names=TRUE, col.names=TRUE, sep="\t", quote=FALSE)

#expr2[1:5, 1:5]
#         160003842T-3 160003809T-2 160003842T-2 160003809T-3 160003805T-3
#VAV2         2.008373     9.444742     3.712596     8.503404     9.323715
#RNU2-36P     0.000000     0.000000     0.000000     0.000000     0.000000
#CFL2         1.831875     7.688678     3.402515     8.192393     8.090352
#ADAL         2.313177     0.763094     3.128527     1.750566     0.929858
#PRMT6        0.000000     0.000000     4.431164     1.602184     0.000000



#calculate tumor scores
filterCommonGenes(input.f = "GSE112996_fpkm_liyanming.txt",
                  output.f = "GSE112996_fpkm_liyanming.gct",
     id = "GeneSymbol") 
     
estimateScore("GSE112996_fpkm_liyanming.gct",   
              "GSE112996_fpkm_liyanming_estimate_score.txt",   
              platform="affymetrix")   


### Output scores for each subject
result <- read.table("GSE112996_fpkm_liyanming_estimate_score.txt",sep = "\t",row.names = 1,check.names = F,stringsAsFactors = F,header = T)
result1 <- result[,-1]   
colnames(result1) <- result1[1,]
   
result2 <- as.data.frame(t(result1[-1,]))
colnames(result2) <- result[-1,1]


### We will probably use the categorized Immune score as another outcome variable.
result2[1:5,]
#                   StromalScore      ImmuneScore    ESTIMATEScore
#X160003842T.3 -462.638565061306 873.968521914912 411.329956853606
#X160003809T.2   223.47956817747 1033.08627889986 1256.56584707733
#X160003842T.2 -326.785429628925 879.273685640113 552.488256011188
#X160003809T.3   367.55140471736 1292.96469585506 1660.51610057242
#X160003805T.3  1321.60650099539 2622.19566582146 3943.80216681686
#                    TumorPurity
#X160003842T.3  0.78669109372014
#X160003809T.2 0.704245493642062
#X160003842T.2 0.773731464002248
#X160003809T.3 0.660935461213534
#X160003805T.3 0.377323741314622

