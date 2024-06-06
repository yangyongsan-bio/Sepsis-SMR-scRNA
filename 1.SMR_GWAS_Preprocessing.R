#The GWAS file input for SMR has 8 columns
##SNP：variant id
##A1：allele
##A2：another allele
##freq：the frequency of the allele，A1
##b：effect size，ES, beta
##se：standard error, SE
##p：p-value
##n：sample size
list=read.table("994.2_PheCode.v1.0.fastGWA.tsv",header = T)
listnew=list[,c(2,4,5,7,11,12,13,6)]
head(listnew)
colnames(listnew)<-c("SNP","A1","A2","freq","b","se","p","n")
write.table(listnew,"GWAS_Sepsis.txt",sep = "\t",row.names = FALSE,quote = FALSE)

  
