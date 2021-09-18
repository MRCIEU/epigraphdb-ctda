#################################################eQTL-Gen-MR analysis-part1
# source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
#update the R package
#library(devtools)
#install_github("MRCIEU/TwoSampleMR")
#use the elder version of the package
#devtools::install_github("MRCIEU/TwoSampleMR@0.3.2")

library(ggplot2)
#link to the test version
library(TwoSampleMR)
#library(MRInstruments)
library("readxl")
#install.packages("MendelianRandomization")
#library(MendelianRandomization)
rm(list=ls(all=TRUE)) 
#toggle_api("test")

setwd("/newhome/epxjz/covid19/")

#ids<-as.character(unlist(read.table("outcome.id.hg19.txt",header=F)))

ids<-as.character(unlist(read.table("/newhome/epxjz/covid19/covid-19-2021-06/outcome.covid19.2021.txt",header=F)))


##read local files as exposure 
exposure_dat <- NULL
try(exposure_dat <- read.table("./instruments/eQTL-Gen-instruments-new.txt",header=T,stringsAsFactors=FALSE, colClasses = c("character")))

try(exposure_dat$exposure <- as.character(exposure_dat$exposure))
try(exposure_dat$SNP <- as.character(exposure_dat$SNP))
try(exposure_dat$effect_allele.exposure <- as.character(exposure_dat$effect_allele.exposure))
try(exposure_dat$other_allele.exposure <- as.character(exposure_dat$other_allele.exposure))
try(exposure_dat$pval_origin.exposure <- as.character(exposure_dat$pval_origin.exposure))
try(exposure_dat$id.exposure <- as.character(exposure_dat$id.exposure))
try(exposure_dat$pval.exposure <- as.numeric(exposure_dat$pval.exposure))
try(exposure_dat$beta.exposure <- as.numeric(exposure_dat$beta.exposure))
try(exposure_dat$se.exposure <- as.numeric(exposure_dat$se.exposure))


#i<-1

#for (i in 8:length(ids)){
#  print(i)
#  outcome_dat <- NULL
#  attempts <- 0
#  while(attempts<=10){
#
#    try(outcome_dat <- read_outcome_data(
#      snps = exposure_dat$SNP,
#      filename = ids[i],
#      sep = "\t",
#      snp_col = "snp",
#      beta_col = "beta",
#      se_col = "se",
#      effect_allele_col = "effect_allele",
#      other_allele_col = "other_allele",
#      eaf_col = "effect_allele_freq",
#      pval_col = "p"
#      ))
#    outcome_dat$outcome <- as.character(as.vector(strsplit(ids[i], "/")[[1]][6]))
#
#    if(is.null(outcome_dat)){
#      attempts<-attempts+1}
#    else{
#      break
#    }
#  }

#i<-4
for (i in 4:length(ids)){
  print(i)
  outcome_dat <- NULL
  attempts <- 0
  while(attempts<=10){

    #try(outcome_dat <- read_outcome_data(
    #  snps = exposure_dat$SNP,
    #  filename = ids[i],
    #  sep = "\t",
    #  snp_col = "rsid",
    #  beta_col = "all_inv_var_meta_beta",
    #  se_col = "all_inv_var_meta_sebeta",
    #  effect_allele_col = "ALT",
    #  other_allele_col = "REF",
    #  eaf_col = "all_meta_AF",
    #  pval_col = "all_inv_var_meta_p"
    #  ))

    try(outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      filename = ids[i],
      sep = "\t",
      snp_col = "snp",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "effect_allele_freq",
      pval_col = "p"
      ))

    outcome_dat$outcome <- as.character(as.vector(strsplit(ids[i], "/")[[1]][6]))

    if(is.null(outcome_dat)){
      attempts<-attempts+1}
    else{
      break
    }
  }

  dat <- NULL
  #try(dat <- harmonise_data(exposure_dat, outcome_dat))
  try(dat <- harmonise_data(exposure_dat, outcome_dat,action=1))

  mr_results <- NULL
  mr_hetero <- NULL
  mr_pleio <- NULL
  mr_single <- NULL
  try(mr_results <- mr(dat, method_list=c("mr_wald_ratio", "mr_ivw")))  
  mr_hetero <- mr_heterogeneity(dat)
  mr_pleio <- mr_pleiotropy_test(dat) 
  try(mr_single <- mr_singlesnp(dat))
  
  exposure <- NULL
  outcome <- NULL 
  exposure <- "eQTL-covid19-2021-06"
  outcome <- as.character(as.vector(strsplit(ids[i], "/")[[1]][6]))
  result_file0 <- paste0("./results/org-results/eQTL/",exposure,".",outcome,".harmonise.txt")
  result_file <- paste0("./results/org-results/eQTL/",exposure,".",outcome,".mr.txt")
  result_file2 <- paste0("./results/org-results/eQTL/",exposure,".",outcome,".mr_hetero.txt")
  result_file3 <- paste0("./results/org-results/eQTL/",exposure,".",outcome,".mr_pleio.txt")
  result_file4 <- paste0("./results/org-results/eQTL/",exposure,".",outcome,".mr_single.txt")
  if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}
}
