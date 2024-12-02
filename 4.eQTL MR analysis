pacman::p_load('data.table','devtools','MRInstruments','TwoSampleMR','data.table','devtools','MRInstruments','biomaRt','fdrtool',
               'MendelianRandomization','simex','MRPRESSO','stringr',"dplyr","remotes","coloc","readr","openxlsx","LDlinkR")

dir="D:\\4.eQTL MR analysis\\IIBDGC_IBD"  
setwd(dir)
eqtl_gene_keep_rsid <- fread("eqtl_gene_keep_rsid.txt"),header=T)  
eqtl_gene_keep_rsid_sig <- subset(eqtl_gene_keep_rsid, pval < 5e-5)
exp_dat <- format_data( 
  eqtl_gene_keep_rsid_sig, 
  type='exposure',  
  snp_col = "rsid",
  chr_col = "CHR",
  pos_col = "POS",
  beta_col = "beta",  
  se_col = "se",  
  effect_allele_col ="effect_allele", 
  other_allele_col = "other_allele", 
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "unique_id",
  gene_col = "Gene",
  id_col = "unique_id")

out_name <- "IIBDGC_IBD"           
Ncase <- 12882                
Ncontrol <- 21770             
out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-31', 
  proxies = FALSE, 
  maf_threshold = 0.01)
out_dat$id.outcome <- out_name                                              
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                           "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

dir="D:\\4.eQTL MR analysis\\FinnGen_IBD"  
setwd(dir)
out_name <- "Finn_IBD"
Ncase <- 7625
Ncontrol <- 369652

out <- read_tsv(str_c("finngen_R9_K11_IBD_STRICT")) 
out <- dplyr::select(out, rsids,`#chrom`, pos, alt, ref, af_alt_controls,beta,sebeta,pval)
names(out) = c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","p")

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
out_dat$id.outcome <- out_name                                               
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

IIBDGC_data <- load('D:\\4.eQTL MR analysis\\IIBDGC-IBD\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\4.eQTL MR analysis\\FinnGen-IBD\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\4.eQTL MR analysis\\meta-IBD.RData')




dir="D:\\4.eQTL MR analysis\\IIBDGC_CD"  
setwd(dir)
eqtl_gene_keep_rsid <- fread("eqtl_gene_keep_rsid.txt"),header=T)  
eqtl_gene_keep_rsid_sig <- subset(eqtl_gene_keep_rsid, pval < 5e-5)
exp_dat <- format_data( 
  eqtl_gene_keep_rsid_sig, 
  type='exposure',  
  snp_col = "rsid",
  chr_col = "CHR",
  pos_col = "POS",
  beta_col = "beta",  
  se_col = "se",  
  effect_allele_col ="effect_allele", 
  other_allele_col = "other_allele", 
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "unique_id",
  gene_col = "Gene",
  id_col = "unique_id")

out_name <- "IIBDGC_CD"           
Ncase <- 5956                
Ncontrol <- 14927                    
out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-30', 
  proxies = FALSE, 
  maf_threshold = 0.01)
out_dat$id.outcome <- out_name                                              
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

dir="D:\\4.eQTL MR analysis\\FinnGen_CD"  
setwd(dir)
out_name <- "Finn_CD"
Ncase <- 1665
Ncontrol <- 375445

out <- read_tsv(str_c("finngen_R9_K11_CD_STRICT2")) 
out <- dplyr::select(out, rsids,`#chrom`, pos, alt, ref, af_alt_controls,beta,sebeta,pval)
names(out) = c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","p")

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
out_dat$id.outcome <- out_name                                               
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"  
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

IIBDGC_data <- load('D:\\4.eQTL MR analysis\\IIBDGC-CD\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\4.eQTL MR analysis\\FinnGen-CD\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\4.eQTL MR analysis\\meta-CD.RData')




dir="D:\\4.eQTL MR analysis\\IIBDGC_UC"  
setwd(dir)
eqtl_gene_keep_rsid <- fread("eqtl_gene_keep_rsid.txt"),header=T)  
eqtl_gene_keep_rsid_sig <- subset(eqtl_gene_keep_rsid, pval < 5e-5)
exp_dat <- format_data( 
  eqtl_gene_keep_rsid_sig, 
  type='exposure',  
  snp_col = "rsid",
  chr_col = "CHR",
  pos_col = "POS",
  beta_col = "beta",  
  se_col = "se",  
  effect_allele_col ="effect_allele", 
  other_allele_col = "other_allele", 
  pval_col = "pval",
  samplesize_col = "samplesize",
  phenotype_col = "unique_id",
  gene_col = "Gene",
  id_col = "unique_id")

out_name <- "IIBDGC_UC"           
Ncase <- 6968               
Ncontrol <- 20464              
out_dat <- extract_outcome_data( 
  snps=exp_dat$SNP, 
  outcomes='ieu-a-32', 
  proxies = FALSE, 
  maf_threshold = 0.01)
out_dat$id.outcome <- out_name                                              
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

dir="D:\\4.eQTL MR analysis\\FinnGen_UC"  
setwd(dir)
out_name <- "Finn_UC"
Ncase <- 5034
Ncontrol <- 371530

out <- read_tsv(str_c("finngen_R9_K11_UC_STRICT2")) 
out <- dplyr::select(out, rsids,`#chrom`, pos, alt, ref, af_alt_controls,beta,sebeta,pval)
names(out) = c("SNP","CHR","BP","A1","A2","FRQ","BETA","SE","p")

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
out_dat$id.outcome <- out_name                                               
out_dat$samplesize.outcome <- Ncase+Ncontrol                                

mydata <- harmonise_data(
  exposure_dat = exp_dat,
  outcome_dat = out_dat_subset)
head(mydata)    
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] >0.6] = "TRUE"   
mydata[,"mr_keep"][mydata[,"mr_keep"] == "FALSE" & mydata[,"eaf.outcome"] <0.4] = "TRUE"
mydata$mr_keep <- as.logical(mydata$mr_keep)

mydata_clump <- clump_data(mydata, clump_r2 = 0.001)
mydata_clump$id.exposure <- mydata_clump$exposure

steiger <- directionality_test(mydata_clump)
steiger_1 <- subset(steiger, correct_causal_direction == TRUE)             
mydata_steiger <- subset(mydata_clump, mydata_clump$exposure %in% steiger_1$exposure)

mr_mrbase <- mr(mydata_steiger,method_list=c("mr_wald_ratio",'mr_ivw_fe','mr_ivw_mre','mr_two_sample_ml','mr_egger_regression',
                                             "mr_weighted_median","mr_penalised_weighted_median")) 
save(mr_mrbase,file = 'mr_mrbase.RData')

heterogeneity <- mr_heterogeneity(mydata_steiger)
save(heterogeneity,file = 'heterogeneity.RData')
pleiotropy <- mr_pleiotropy_test(mydata_steiger)
save(pleiotropy,file = 'pleiotropy.RData')

IIBDGC_data <- load('D:\\4.eQTL MR analysis\\IIBDGC-UC\\clean_mr_mrbase.RData')
finn_data <-  load('D:\\4.eQTL MR analysis\\FinnGen-UC\\clean_mr_mrbase.RData')               
metadata <- rbind(IIBDGC_data, finn_data)
metamod <- rma(yi = metadata[,1], sei = metadata[,2], method = "FE")  
save(metamod,file = 'D:\\4.eQTL MR analysis\\meta-UC.RData')

