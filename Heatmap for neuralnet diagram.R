library(ggplot2)
library(plotly)
library(scico)
CDK_zscore_matrix <- read.csv("~/Rotation 2/Raw data/Neural network/CDK zscore matrix.csv")
rownames(CDK_zscore_matrix)= CDK_zscore_matrix$X
CDK_zscore_matrix= CDK_zscore_matrix[,-1]
CDK_zscore_matrix$AA5_zscore[is.na(CDK_zscore_matrix$AA5_zscore)]= 0
CDK_zscore_matrix$AA6_zscore[is.na(CDK_zscore_matrix$AA6_zscore)]= 0
CDK_zscore_matrix$AA7_zscore[is.na(CDK_zscore_matrix$AA7_zscore)]= 0
CDK_zscore_matrix$AA8_zscore[is.na(CDK_zscore_matrix$AA8_zscore)]= 0
x <- colnames(CDK_zscore_matrix)
y <- rownames(CDK_zscore_matrix)
data <- expand.grid(X=x, Y=y)
data$Z <- scale(CDK_zscore_matrix)[1:length(scale(CDK_zscore_matrix))]
v="scale"
v= ggplot(data, aes(X, Y, fill= Z)) + 
  geom_tile() +
  labs(x = "Category",y= "Amino acid",title ="",fill="Scale") +
  theme(plot.title = element_text(color="black", size=12,hjust = 0.5)) +
  theme(axis.text = element_text(color="black", size=12,hjust = 0.5)) +
  theme(axis.title = element_text(color="black", size=12,hjust = 0.5)) +
  theme( axis.line = element_line(colour = "black",size = 1, linetype = "solid")) +
  scale_fill_distiller(palette = "RdYlGn") +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))
v
heatmap(as.numeric(CDK_zscore_matrix$AA1_zscore))
