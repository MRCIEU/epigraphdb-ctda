#################################################MR analysis of ENTPD5
# source("https://bioconductor.org/biocLite.R")
#install.packages("devtools")
#update the R package
library(devtools)
install_github("MRCIEU/TwoSampleMR")
#use the elder version of the package
#devtools::install_github("MRCIEU/TwoSampleMR@0.4.26")
devtools::install_github("MRCIEU/MRInstruments")

library(ggplot2)
#link to the test version
library(TwoSampleMR)
#library(MRInstruments)
library("readxl")
#install.packages("MendelianRandomization")
#library(MendelianRandomization)
rm(list=ls(all=TRUE)) 
#toggle_api("test")
#library("RadialMR")

setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/Collaboration-projects/SARS-MR/COVID19-gwas/")

###step 1. ENTPD5 to selected outcomes
#ENTPD5 instruments
#pQTL
exposure_dat <- read_excel("ENTPD5-ins.xlsx",1)
#eQTL
exposure_dat <- read_excel("ENTPD5-ins.xlsx",2)

#COVID-19 severity
exposure_dat <- read_excel("ENTPD5-ins.xlsx",5)

exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_freq", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P",samplesize_col="N")
#exposure_dat <- clump_data(exposure_dat,clump_r2 = 0.001)

ao<-available_outcomes()

outcome_dat <- NULL

ids <- c("ukb-b-8909","ukb-b-13354","ukb-d-30710_irnt","ukb-d-30080_irnt")

#F8
ids <- c("prot-a-1009")
exposure_dat <- extract_instruments(outcomes=ids)

#D-dimer	
ids <- c("prot-a-1086")

outcome_dat<-extract_outcome_data(exposure_dat$SNP,ids,access_token=NULL) 

#HDL-C
ids <- c("ukb-d-30760_irnt")
outcome_dat<-extract_outcome_data(exposure_dat$SNP,ids,access_token=NULL) 

dat <- NULL
try(dat <- harmonise_data(exposure_dat, outcome_dat,action=2))

#mr_results <- mr(dat,method_list=c("mr_ivw_radial","mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
mr_results <- mr(dat,method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
mr_hetero <- mr_heterogeneity(dat)
mr_pleio <- mr_pleiotropy_test(dat) 
mr_single <- mr_singlesnp(dat)

exposure <- exposure_dat$exposure[2]
outcome <- outcome_dat$outcome[1]

exposure <- "ENTPD5-eQTL"
outcome <- "HDL" 

#exposure <- paste0(exposure_dat$exposure[1],".",outcome_dat$outcome[1])
##European
result_file0 <- paste0("./results/org-results/",exposure,"_",outcome,".harmonise.txt")
result_file <- paste0("./results/org-results/",exposure,"_",outcome,".mr.txt")
result_file2 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_hetero.txt")
result_file3 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_pleio.txt")
result_file4 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_single.txt")

if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}


###step 2. selected outcome to COVID19
###step 3. reverse MR of selected outcomes on ENTPD5
setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/Collaboration-projects/SARS-MR/COVID19-gwas/")
rm(list=ls(all=TRUE)) 

##selected outcomes
ids <- c("ukb-b-8909","ukb-b-13354","ukb-d-30710_irnt","ukb-d-30080_irnt")
exposure_dat <- extract_instruments(outcomes=ids)

##DD/PT/APTT
exposure_dat <- read_excel("data/meta-APTT-PT-DD-IV.xlsx",1)
exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_freq", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P",samplesize_col="N")

##HDL cholesterol - lipids
#univariate lipids MR instruments 
exposure_dat <- read_excel("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/SOST-George/kloto/klotho-MR-instruments-lipids.xlsx",1)

exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_frequency", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P",samplesize_col="N")

#multivariable MR instrument HDL-APOA1
exposure_dat <- read_excel("klotho-MR-instruments-lipids.xlsx",2)

#multivariable MR instrument LDL-TG-APOB
exposure_dat <- read_excel("klotho-MR-instruments-lipids.xlsx",3)

exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_frequency", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P",samplesize_col="N")



#CRP cis instruments

ids <- c("ukb-d-30710_irnt")	
exposure_dat <- extract_instruments(outcomes=ids)
if (exists("exposure_dat")==TRUE){write.table(exposure_dat,file="CPR.instruments.txt",sep="\t",col.names=T,row.names=F,quote=F)}

#F8 instruments
ids <- c("prot-a-1009")
exposure_dat <- extract_instruments(outcomes=ids)

#d-dimer
ids <- c("prot-a-1086")
exposure_dat <- extract_instruments(outcomes=ids)

#BMI
ids <- c("ieu-a-2")
exposure_dat <- extract_instruments(outcomes=ids)

#outcome instruments in ENTPD5
exposure_dat <- read_excel("ENTPD5-ins.xlsx",3)

exposure_dat <- read_excel("ENTPD5-ins.xlsx",4)
exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_freq", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P",samplesize_col="N")

outcome_dat <- NULL

attempts <- 0
while(attempts<=10){
  try(
    outcome_dat <- read_outcome_data(
      snps = exposure_dat$SNP,
      #COVID19 hosC-vs-non-hosC
      #filename = "./data/COVID19_HGI_ANA_B1_V2_20200701.b37.txt.tab",  
      #ENTPD5
      filename = "./data/ENTPD5.4437.56.3.tsv.gz.tab",  
      sep = "\t",
      snp_col = "snp",
      beta_col = "beta",
      se_col = "se",
      effect_allele_col = "effect_allele",
      other_allele_col = "other_allele",
      eaf_col = "effect_allele_freq",
      pval_col = "p",
      samplesize_col = "n")
  )
  if(is.null(outcome_dat)){
    attempts<-attempts+1}
  else{
    break
  }
}
#outcome_dat$outcome <- "COVID19"  
outcome_dat$outcome <- "ENTPD5"  

dat <- NULL
try(dat <- harmonise_data(exposure_dat, outcome_dat,action=2))

mr_results <- mr(dat,method_list=c("mr_wald_ratio","mr_ivw","mr_egger_regression","mr_weighted_median","mr_simple_mode","mr_weighted_mode"))
mr_hetero <- mr_heterogeneity(dat)
mr_pleio <- mr_pleiotropy_test(dat) 
mr_single <- mr_singlesnp(dat)

#exposure <- "selected_outcomes-ENTPD5"
#exposure <- "DD-APTT-PT"
exposure <- "lipids"
#exposure <- "CRP"
outcome <- "COVID19" 
#outcome <- "ENTPD5" 


#exposure <- paste0(exposure_dat$exposure[1],".",outcome_dat$outcome[1])
##European
result_file0 <- paste0("./results/org-results/",exposure,"_",outcome,".harmonise.txt")
result_file <- paste0("./results/org-results/",exposure,"_",outcome,".mr.txt")
result_file2 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_hetero.txt")
result_file3 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_pleio.txt")
result_file4 <- paste0("./results/org-results/",exposure,"_",outcome,".mr_single.txt")

if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}


###########MR PheWAS
cd /newhome/epxjz/covid19/
R

#library(ggplot2)
#link to the test version
library(TwoSampleMR)
library(MRInstruments)

library("readxl")
#install.packages("MendelianRandomization")
#library(MendelianRandomization)
rm(list=ls(all=TRUE))


#pQTL
exposure_dat <- read_excel("./instruments/ENTPD5-ins.xlsx",1)
#eQTL
exposure_dat <- read_excel("./instruments/ENTPD5-ins.xlsx",2)

#reformatting data
exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_freq", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P")

nrow(exposure_dat)

#exposure_dat <- clump_data(exposure_dat,clump_r2 = 0.01)
#if (exists("exposure_dat")==TRUE){write.table(exposure_dat,file="temp.ins",sep="\t",col.names=T,row.names=F,quote=F)}


##MRB phentoypes
ao<-available_outcomes(access_token=NULL)
ids<-as.character(unlist(read.table("/panfs/panasas01/sscm/epxjz/eQTL-Gen/outcome.id.txt",header=F)))

#ids<-ids[1:105]
#ids 111~123

ids <- ids[1:103]

outcome_dat <- NULL
attempts <- 0
while(attempts<=10){  
  try(outcome_dat<-extract_outcome_data(exposure_dat$SNP,ids,proxies = TRUE,access_token=NULL))    
  if(is.null(outcome_dat)){
    attempts<-attempts+1}
  else{
    break
  }
}

dat <- NULL
try(dat <- harmonise_data(exposure_dat, outcome_dat))

mr_results <- NULL
mr_hetero <- NULL
mr_pleio <- NULL
mr_single <- NULL
try(mr_results <- mr(dat))  
mr_hetero <- mr_heterogeneity(dat)
mr_pleio <- mr_pleiotropy_test(dat) 
try(mr_single <- mr_singlesnp(dat))

exposure <- NULL 
exposure <- "ENTPD5-pQTL-phewas"
result_file0 <- paste0("./results/org-results/",exposure,".harmonise.txt")
result_file <- paste0("./results/org-results/",exposure,".mr.txt")
result_file2 <- paste0("./results/org-results/",exposure,".mr_hetero.txt")
result_file3 <- paste0("./results/org-results/",exposure,".mr_pleio.txt")
result_file4 <- paste0("./results/org-results/",exposure,".mr_single.txt")
if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}

result_file0

##SAIGE phenotypes

library(ggplot2)
#link to the test version
library(TwoSampleMR)
library(MRInstruments)

library("readxl")
#install.packages("MendelianRandomization")
#library(MendelianRandomization)
rm(list=ls(all=TRUE))

#pQTL
exposure_dat <- read_excel("./instruments/ENTPD5-ins.xlsx",1)
#eQTL
exposure_dat <- read_excel("./instruments/ENTPD5-ins.xlsx",2)

#reformatting data
exposure_dat <-format_data(exposure_dat, type = "exposure", header = TRUE,
                           phenotype_col = "Phenotype", snp_col = "SNP", beta_col = "Beta",
                           se_col = "Se", eaf_col = "Effect_allele_freq", effect_allele_col = "Effect_allele",
                           other_allele_col = "Other_allele", pval_col = "P")

ids<-as.character(unlist(read.table("/panfs/panasas01/sscm/epxjz/eQTL-Gen/outcome.id.SAIGE.v2.txt",header=F)))

i<-1

for (i in 1:length(ids)){
  print(i)
  outcome_dat <- NULL
  attempts <- 0
  while(attempts<=10){
    
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
    #outcome_dat$outcome <- as.character(as.vector(strsplit(ids[i], "_")[[1]][2]))
    
    if(is.null(outcome_dat)){
      attempts<-attempts+1}
    else{
      break
    }
  }
  
  outcome_dat$id.outcome <- paste0("SAIGE-",as.character(as.vector(strsplit(ids[i], "_")[[1]][2])))
  
  
  dat <- NULL
  #try(dat <- harmonise_data(exposure_dat, outcome_dat))
  try(dat <- harmonise_data(exposure_dat, outcome_dat,action=2))
  
  mr_results <- NULL
  mr_hetero <- NULL
  mr_pleio <- NULL
  mr_single <- NULL
  try(mr_results <- mr(dat))  
  mr_hetero <- mr_heterogeneity(dat)
  mr_pleio <- mr_pleiotropy_test(dat) 
  try(mr_single <- mr_singlesnp(dat))
  
  exposure <- NULL 
  exposure <- as.character(as.vector(strsplit(ids[i], "_")[[1]][2]))
  
  result_file0 <- paste0("./results/org-results/",exposure,".harmonise.txt")
  result_file <- paste0("./results/org-results/",exposure,".mr.txt")
  result_file2 <- paste0("./results/org-results/",exposure,".mr_hetero.txt")
  result_file3 <- paste0("./results/org-results/",exposure,".mr_pleio.txt")
  result_file4 <- paste0("./results/org-results/",exposure,".mr_single.txt")
  
  
  #result_file0 <- paste0("./results/org-results/overall/",exposure,".harmonise.txt")
  #result_file <- paste0("./results/org-results/overall/",exposure,".mr.txt")
  #result_file2 <- paste0("./results/org-results/overall/",exposure,".mr_hetero.txt")
  #result_file3 <- paste0("./results/org-results/overall/",exposure,".mr_pleio.txt")
  #result_file4 <- paste0("./results/org-results/overall/",exposure,".mr_single.txt")
  
  if (exists("dat")==TRUE){ write.table(dat,file=result_file0,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_results")==TRUE){ write.table(mr_results,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_hetero")==TRUE){ write.table(mr_hetero,file=result_file2,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_pleio")==TRUE){write.table(mr_pleio,file=result_file3,sep="\t",col.names=T,row.names=F,quote=F)}
  if (exists("mr_single")==TRUE){write.table(mr_single,file=result_file4,sep="\t",col.names=T,row.names=F,quote=F)}
}


#####
devtools::install_github("MRCIEU/epigraphdb-r")
library("epigraphdb")

rm(list=ls(all=TRUE))

setwd("~/OneDrive - University of Bristol/Google_Drive/working_space/Post-doc-Bristol/projects/Collaboration-projects/SARS-MR/COVID19-gwas/")

ids<- c("AATF","NEU1","NPC2","NUP210","PIGS","PLOD2","PMPCA","RETREG3","SPART","SRP72","TLE3")

#ENSG00000141699	RETREG3
#ENSG00000133104	SPART

ids[i]<-"SPART"

i <- 8
for (i in 1:length(ids)){
  print(i)

  res <- query_epigraphdb(
  #route="/xqtl/single-snp-mr",
  route="/covid-19/ctda/single-snp-mr/gene",
  params=list(exposure_gene=ids[i], outcome_trait=NULL, variant=NULL, qtl_type="pQTL", pval_threshold=1),
  mode="table"
)

res  
  
exposure <- ids[i]  
result_file <- paste0("./results/org-results/",exposure,".lookup.results.txt")
if (exists("res")==TRUE){ write.table(res,file=result_file,sep="\t",col.names=T,row.names=F,quote=F)}
}
