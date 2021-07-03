#####################################MOLOC analysis##################################################
install.packages("devtools")
devtools::install_github("clagiamba/moloc")

install.packages("BiocManager")
BiocManager::install("snpStats")
#install.packages("truncnorm")
#install.packages("coloc")
install.packages("tidyverse")
#devtools::install_github("chr1swallace/coloc")


library(dplyr, lib.loc = "/Library/Frameworks/R.framework/Versions/3.6/Resources/library")


#library(snpStats)
#library(truncnorm)
#library(coloc)
library(moloc)
library(ggplot2)
#link to the test version
library(TwoSampleMR)
#ao<-available_outcomes()
#library(MRInstruments)
library("readxl")
library(tidyverse)
#library("optparse")
rm(list=ls(all=TRUE)) 

setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/pQTL-eQTL-combined-analysis/data/moloc-analysis/")

setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/Collaboration-projects/SARS-MR/COVID19-gwas/data/")

opt_all <- read.table("moloc-input/filename",header=F)  



i<-1
for (i in 1:nrow(opt_all)){

file <- paste0("moloc-input/",as.vector(opt_all$V1[i]))

df <- read.table(file,header=F)
#df <- read.table("moloc-input/temp.ENSG00000149564.rs6590090.ESAM.2981.9.3.ckqny.scz2snpres.gz.tab.assoc1",header=F)

for (j in 1:nrow(df)){
  if(is.na(df$V20[j])==TRUE){df$V20=df$V6}
  if(df$V20[j]==".") {df$V20=df$V6}
}
df <- df[complete.cases(df), ]
df <- df[df$V8!="Inf",]
df <- df[!duplicated(df$V1), ]

###tips from Mahsa
#No problem. So, below is what I did to run moloc for Atrial Fibrillation (AF) GWAS. The “final_data” dataframe is a merged file with all 3 datasets available after get overlapping SNPs (as you said the snp name is consistent across all data frames)
##make each dataset

#A for gene expression 
A <- df %>% select(1:10)
A$MAF <- NULL
for (j in 1:nrow(A)){
  if(A$V6[j]>0.5){A$MAF[j]=1-A$V6[j]} else {A$MAF[j]=A$V6[j]} 
}

colnames(A) <- c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","PVAL","N","MAF")

#B for protein expression
B <- df %>% select(1:3,11:17)
B$MAF <- NULL
for (j in 1:nrow(B)){
  if(B$V13[j]>0.5){B$MAF[j]=1-B$V13[j]} else {B$MAF[j]=B$V13[j]} 
}

colnames(B) <- c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","PVAL","N","MAF")


#C for disease outcomes
C <- df %>% select(1:3,18:24)
for (j in 1:nrow(C)){
  if(C$V20[j]>0.5){C$MAF[j]=1-C$V20[j]} else {C$MAF[j]=C$V20[j]} 
}

colnames(C) <- c("SNP","CHR","POS","A1","A2","EAF","BETA","SE","PVAL","N","MAF")

##check if trait is case control #already done for trait data (A)

test_genome <- list(eqtl=A, pqtl=B, trait=C)

results <- moloc_test(test_genome)
temp <- results$priors_lkl_ppa
n_snps <- results$nsnps
out <- cbind(as.vector(opt_all$V1[i]),n_snps,data.frame(matrix(unlist(temp$PPA), nrow=1, byrow=T)))

result_file <- paste0("moloc-results/",as.vector(opt_all$V2[i]))
if (exists("out")==TRUE){ write.table(out,file=result_file,sep="\t",col.names=F,row.names=F,quote=F)}
}
