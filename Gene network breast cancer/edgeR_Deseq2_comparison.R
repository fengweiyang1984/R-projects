library(edgeR)
library( "DESeq2" )
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
qlfT_sort$fdr <- p.adjust(qlfT_sort[,4], method ="fdr")
edgeR.genes <- rownames(qlfT_sort)
length(edgeR.genes)
write.csv(as.data.frame(qlfT_sort), file="case_control_edgeR.csv")

metaData <- as.data.frame(colnames(ExprDD))
metaData$TC <- caseControlResponse
ExprDD1 <- cbind(rownames(ExprDD),ExprDD)
dds <- DESeqDataSetFromMatrix(countData=ExprDD1, colData=metaData, design=~TC, tidy = TRUE)
dds <- DESeq(dds)
resMFType <- results(dds)
#summary(resLRT) 
#resLRT_Sig = subset(resLRT, pvalue<0.05)
resLRTOrdered <- resMFType[order(resMFType$pvalue),]
write.csv(as.data.frame(resLRTOrdered), file="case_control_Deseq.csv")