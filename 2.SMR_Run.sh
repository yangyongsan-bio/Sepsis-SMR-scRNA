#The native CMD runs
d:
cd D:\software\smr-1.3.1-win-x86_64\smr-1.3.1-win-x86_64

##Run SMR
##eQTLGen
smr-1.3.1-win.exe --bfile /g1000_eur/g1000_eur --gwas-summary /gwas/GWAS_Sepsis.txt --beqtl-summary /eQTLGen_cis_eQTL_hg19/eQTLGen --out Sepsis_SMR_eQTLGen --thread-num 10 
##GTEx
smr-1.3.1-win.exe --bfile /g1000_eur/g1000_eur --gwas-summary /gwas/GWAS_Sepsis.txt --beqtl-summary /GTEx_V8_Whole_Blood_eQTL_hg19/Whole_Blood/Whole_Blood --out Sepsis_SMR_GTEx --thread-num 10 
