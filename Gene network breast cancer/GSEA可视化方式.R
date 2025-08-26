library(clusterProfiler)
library(org.Hs.eg.db)
data(geneList, package="DOSE")

# enrichment
ego3 <- gseKEGG(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 100,
              maxGSSize    = 500,
              pvalueCutoff = 0.05,
              verbose      = FALSE)

temp <- as.data.frame(geneList)
