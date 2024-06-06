# Reading data
data <- read.delim("Sepsis_SMR_GTEx.smr", header = TRUE, stringsAsFactors = FALSE)
data <- read.delim("Sepsis_SMR_eQTLGen.smr", header = TRUE, stringsAsFactors = FALSE)

## Conversion if Gene is an ENSEMBL ID
library(clusterProfiler)
library(org.Hs.eg.db)
df=bitr(geneID = data$Gene,fromType = 'ENSEMBL',toType = 'SYMBOL',OrgDb = 'org.Hs.eg.db')
data=merge(data,df,by.x='Gene',by.y='ENSEMBL')
data=data[,-1]
colnames(data)[22]='Gene'
data=data[,c(1,2,22,3:21)]
# Output
write.csv(data, "SMR_results_eQTLGen.csv", row.names = FALSE)



