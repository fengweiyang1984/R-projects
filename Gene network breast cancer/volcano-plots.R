library(ggplot2)
#upload the table firstly
#change table names
de <- subtypes_edgeR_V2[,1:5]
colnames(de)
colnames(de) <- c("Geneid","log2FC", "logCPM","F","P_value")

#log2FC
de$P_value[which(de$P_value==0, arr.ind=TRUE)]=0.00000001
FC_f <- function(x){
  ifelse(x>0, 2^(x), -2^(-x))
}
round_f <- function(x){
  round(x,4)
}
de$FC <- lapply(de$log2FC, FC_f)
de$FC <- lapply(de$FC, round_f)

# add a column of NAs
de$diffexpressed <- "NO"
# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP" 
de$diffexpressed[de$FC > 1.2 & de$P_value < 0.05] <- "UP"
# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
de$diffexpressed[de$FC < -1.2 & de$P_value < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of genes differentially expressed (NA in case they are not)
#de$delabel <- NA
#de$delabel <- de$Geneid[abs(de$PFC.VTA.SAL)> 1.2& de$FDR < 0.05]

#generate volcano plots in ggplot2
de$expression = ifelse(de$P_value < 0.05 & abs(de$log2FC) >= 1, 
                      ifelse(de$log2FC> 1 ,'Up','Down'),
                      'Stable')
ggplot(data = de, 
       aes(x = de$log2FC, 
           y = -log10(de$P_value), 
           colour=expression,
           label = de$label)) +
  geom_point(alpha=0.4, size=3.5) +
  scale_color_manual(values=c("blue", "black","red"))+
  xlim(c(-1, 1)) +
  ylim(c(0,12))+
  geom_vline(xintercept=c(-1,1),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = 1.301,lty=4,col="black",lwd=0.8) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Differential expression")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank())

# Make a basic volcano plot
with(de, plot(log2FC, -log10(P_value), pch=20, main="Volcano plot", xlim=c(-0.6,0.6),ylim=c(0,50)))
# Add colored points: red if p<0.05&FC>1.2)
with(subset(de, P_value<=.05 & log2FC>0.2), points(log2FC, -log10(P_value), pch=20, col="red"))
with(subset(de, P_value<=.05 & log2FC< -0.2), points(log2FC, -log10(P_value), pch=20, col="blue"))
abline(v=0, , lty=2, lwd= 3)
