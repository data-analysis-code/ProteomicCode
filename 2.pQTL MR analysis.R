pacman::p_load('data.table','devtools','MRInstruments','TwoSampleMR','data.table','devtools','MRInstruments','biomaRt','fdrtool',
               'MendelianRandomization','simex','MRPRESSO','stringr',"dplyr","remotes","coloc","readr","openxlsx","LDlinkR")
rm(list = ls())
dir="D:\\2.pQTL MR analysis\\IIBDGC-IBD"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  

out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-31', 
  proxies = FALSE, 
  maf_threshold = 0.01)

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


rm(list = ls())
dir="D:\\2.pQTL MR analysis\\FinnGen-IBD"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  
out <- read_tsv(str_c("FinnGen_IBD"))
out_dat <- format_data( 
  out, 
  snps = exp_dat$SNP,
  type='outcome',  
  snp_col = "SNP",                                                        
  chr_col = "CHR",
  pos_col = "BP",
  beta_col = "BETA",  
  se_col = "SE",  
  effect_allele_col ="A1", 
  other_allele_col = "A2", 
  eaf_col = "FRQ",
  pval_col = "P")
mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


IIBDGC_data <- load('D:\\2.pQTL MR analysis\\IIBDGC-IBD\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\2.pQTL MR analysis\\FinnGen-IBD\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\2.pQTL MR analysis\\meta-IBD.RData')



pacman::p_load('data.table','devtools','MRInstruments','TwoSampleMR','data.table','devtools','MRInstruments','biomaRt','fdrtool',
               'MendelianRandomization','simex','MRPRESSO','stringr',"dplyr","remotes","coloc","readr","openxlsx","LDlinkR")
rm(list = ls())
dir="D:\\2.pQTL MR analysis\\IIBDGC-CD"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  

out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-30', 
  proxies = FALSE, 
  maf_threshold = 0.01)

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


rm(list = ls())
dir="D:\\2.pQTL MR analysis\\FinnGen-CD"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  
out <- read_tsv(str_c("FinnGen_CD"))
out_dat <- format_data( 
  out, 
  snps = exp_dat$SNP,
  type='outcome',  
  snp_col = "SNP",                                                        
  chr_col = "CHR",
  pos_col = "BP",
  beta_col = "BETA",  
  se_col = "SE",  
  effect_allele_col ="A1", 
  other_allele_col = "A2", 
  eaf_col = "FRQ",
  pval_col = "P")
mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


IIBDGC_data <- load('D:\\2.pQTL MR analysis\\IIBDGC-CD\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\2.pQTL MR analysis\\FinnGen-CD\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\2.pQTL MR analysis\\meta-CD.RData')



pacman::p_load('data.table','devtools','MRInstruments','TwoSampleMR','data.table','devtools','MRInstruments','biomaRt','fdrtool',
               'MendelianRandomization','simex','MRPRESSO','stringr',"dplyr","remotes","coloc","readr","openxlsx","LDlinkR")
rm(list = ls())
dir="D:\\2.pQTL MR analysis\\IIBDGC-UC"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  

out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-32', 
  proxies = FALSE, 
  maf_threshold = 0.01)

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


rm(list = ls())
dir="D:\\2.pQTL MR analysis\\FinnGen-UC"   
setwd(dir)
exp_dat <- fread("exp_dat.txt"),header=T)  
out <- read_tsv(str_c("FinnGen_UC"))
out_dat <- format_data( 
  out, 
  snps = exp_dat$SNP,
  type='outcome',  
  snp_col = "SNP",                                                        
  chr_col = "CHR",
  pos_col = "BP",
  beta_col = "BETA",  
  se_col = "SE",  
  effect_allele_col ="A1", 
  other_allele_col = "A2", 
  eaf_col = "FRQ",
  pval_col = "P")
mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat)

steiger <- directionality_test(mydata)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata, mydata$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')
heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')


IIBDGC_data <- load('D:\\2.pQTL MR analysis\\IIBDGC-UC\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\2.pQTL MR analysis\\FinnGen-UC\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\2.pQTL MR analysis\\meta-UC.RData')
