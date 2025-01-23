rm(list = ls())
pacman::p_load('data.table','devtools','MRInstruments','TwoSampleMR','data.table','devtools','MRInstruments','biomaRt','fdrtool','tidyr',
               'MendelianRandomization','simex','MRPRESSO','stringr',"dplyr","remotes","coloc","readr","openxlsx","LDlinkR")
dir="D:\\3.colocalization analysis"  
setwd(dir)

load("colocalization_data_IBD.RData")

for (i in 1:nrow(unique_exposures)) {
  dat1 <- dat0[which(dat0$exposure == unique_exposures[i,]), ]                  
  print(i)
  coloc_results_list <-
    coloc.abf(dataset1 = list( beta = dat1$beta.exposure,                       
                               varbeta = dat1$var.exposure,                    
                               N = dat1$samplesize.exposure,                   
                               type = "quant",
                               MAF = dat1$maf,
                               # sdY = 1,
                               snp = dat1$SNP),
              dataset2 = list( beta = dat1$beta.outcome,
                               varbeta = dat1$var.outcome,
                               type = "cc",                                     
                               N = dat1$samplesize.outcome,
                               s = cc_ratio,
                               snp = dat1$SNP)
    )
  coloc_result_summary <- cbind(t(data.frame(coloc_results_list$summary)), unique(dat1$exposure))
  coloc_result <- rbind(coloc_result, coloc_result_summary)
}
write.csv(coloc_result, str_c("coloc_result_IBD.csv"))

load("colocalization_data_CD.RData")

for (i in 1:nrow(unique_exposures)) {
  dat1 <- dat0[which(dat0$exposure == unique_exposures[i,]), ]                  
  print(i)
  coloc_results_list <-
    coloc.abf(dataset1 = list( beta = dat1$beta.exposure,                       
                               varbeta = dat1$var.exposure,                    
                               N = dat1$samplesize.exposure,                   
                               type = "quant",
                               MAF = dat1$maf,
                               # sdY = 1,
                               snp = dat1$SNP),
              dataset2 = list( beta = dat1$beta.outcome,
                               varbeta = dat1$var.outcome,
                               type = "cc",                                     
                               N = dat1$samplesize.outcome,
                               s = cc_ratio,
                               snp = dat1$SNP)
    )
  coloc_result_summary <- cbind(t(data.frame(coloc_results_list$summary)), unique(dat1$exposure))
  coloc_result <- rbind(coloc_result, coloc_result_summary)
}
write.csv(coloc_result, str_c("coloc_result_CD.csv"))

load("colocalization_data_UC.RData")

for (i in 1:nrow(unique_exposures)) {
  dat1 <- dat0[which(dat0$exposure == unique_exposures[i,]), ]                  
  print(i)
  coloc_results_list <-
    coloc.abf(dataset1 = list( beta = dat1$beta.exposure,                       
                               varbeta = dat1$var.exposure,                    
                               N = dat1$samplesize.exposure,                   
                               type = "quant",
                               MAF = dat1$maf,
                               # sdY = 1,
                               snp = dat1$SNP),
              dataset2 = list( beta = dat1$beta.outcome,
                               varbeta = dat1$var.outcome,
                               type = "cc",                                     
                               N = dat1$samplesize.outcome,
                               s = cc_ratio,
                               snp = dat1$SNP)
    )
  coloc_result_summary <- cbind(t(data.frame(coloc_results_list$summary)), unique(dat1$exposure))
  coloc_result <- rbind(coloc_result, coloc_result_summary)
}
write.csv(coloc_result, str_c("coloc_result_UC.csv"))
