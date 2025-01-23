pacman::p_load('dplyr','data.table','fst','plyr','tidyr','Hmisc','tableone','zoo','mgcv','splines',stringr,dplyr,progress,survival,openxlsx,plotRCS)
setwd('D:\\ProteomicCode\\1.Cox analysis')

load("1.data_protein_IBD.RData")

## Model 1  ####
cox_model_1 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+")))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_1 <- data.frame('Model1_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model1_HR' = HR,
                            'Model1_CI' = CI,
                            'Model1_P' = P)
  return(cox_model_1)
}

## Model 2  ####
cox_model_2 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_2 <- data.frame('Model2_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model2_HR' = HR,
                            'Model2_CI' = CI,
                            'Model2_P' = P)
  return(cox_model_2)
}

## Model 3  ####
cox_model_3 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI + smoking + alcohol + PA + diet + sleep + BMI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_3 <- data.frame('Model3_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model3_HR' = HR,
                            'Model3_CI' = CI,
                            'Model3_P' = P)
  return(cox_model_3)
}

results_1 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]   
  model_1 <- lapply(exp, cox_model_1)
  result_1 <- ldply(model_1,data.frame)
  results_1[[length(results_1) + 1]] <- result_1
}
write.xlsx(final_results_1,paste0("IBD_model1.xlsx"),row.names = FALSE)

results_2 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_2 <- lapply(exp, cox_model_2)
  result_2 <- ldply(model_2,data.frame)
  results_2[[length(results_2) + 1]] <- result_2
}
write.xlsx(final_results_2,paste0("IBD_model2.xlsx"),row.names = FALSE)

results_3 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_3 <- lapply(exp, cox_model_3)
  result_3 <- ldply(model_3,data.frame)
  results_3[[length(results_3) + 1]] <- result_3
}
write.xlsx(final_results_3,paste0("IBD_model3.xlsx"),row.names = FALSE)


load("2.data_protein_CD.RData")

## Model 1  ####
cox_model_1 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+")))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_1 <- data.frame('Model1_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model1_HR' = HR,
                            'Model1_CI' = CI,
                            'Model1_P' = P)
  return(cox_model_1)
}

## Model 2  ####
cox_model_2 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_2 <- data.frame('Model2_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model2_HR' = HR,
                            'Model2_CI' = CI,
                            'Model2_P' = P)
  return(cox_model_2)
}

## Model 3  ####
cox_model_3 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI + smoking + alcohol + PA + diet + sleep + BMI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_3 <- data.frame('Model3_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model3_HR' = HR,
                            'Model3_CI' = CI,
                            'Model3_P' = P)
  return(cox_model_3)
}

results_1 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]   
  model_1 <- lapply(exp, cox_model_1)
  result_1 <- ldply(model_1,data.frame)
  results_1[[length(results_1) + 1]] <- result_1
}
write.xlsx(final_results_1,paste0("CD_model1.xlsx"),row.names = FALSE)

results_2 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_2 <- lapply(exp, cox_model_2)
  result_2 <- ldply(model_2,data.frame)
  results_2[[length(results_2) + 1]] <- result_2
}
write.xlsx(final_results_2,paste0("CD_model2.xlsx"),row.names = FALSE)

results_3 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_3 <- lapply(exp, cox_model_3)
  result_3 <- ldply(model_3,data.frame)
  results_3[[length(results_3) + 1]] <- result_3
}
write.xlsx(final_results_3,paste0("CD_model3.xlsx"),row.names = FALSE)


load("3.data_protein_UC.RData")

## Model 1  ####
cox_model_1 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+")))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_1 <- data.frame('Model1_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model1_HR' = HR,
                            'Model1_CI' = CI,
                            'Model1_P' = P)
  return(cox_model_1)
}

## Model 2  ####
cox_model_2 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_2 <- data.frame('Model2_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model2_HR' = HR,
                            'Model2_CI' = CI,
                            'Model2_P' = P)
  return(cox_model_2)
}

## Model 3  ####
cox_model_3 <- function(x){
  FML <- as.formula(paste0("Surv(time, event==1) ~ ", paste(x, collapse = "+"), " + age + sex + 
   assessment center + education + employment  + income+ TDI + smoking + alcohol + PA + diet + sleep + BMI"))
  cox1 <- coxph(FML, data = mydata)
  cox2 <- summary(cox1)
  n <- length(grep(x, names(cox1$coefficients)))
  HR <- round(exp(coef(cox1)),3)[1:n]
  SE <- cox2$coefficients[,3][1:n]
  CI5 <- round(exp(coef(cox1)-1.96*SE),3)[1:n]
  CI95 <- round(exp(coef(cox1)+1.96*SE),3)[1:n]
  CI <- paste0(CI5,'-',CI95)
  P <- cox2$coefficients[,5][1:n]
  cox_model_3 <- data.frame('Model3_Characteristics' = names(cox1$coefficients)[1:n],
                            'Model3_HR' = HR,
                            'Model3_CI' = CI,
                            'Model3_P' = P)
  return(cox_model_3)
}

results_1 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]   
  model_1 <- lapply(exp, cox_model_1)
  result_1 <- ldply(model_1,data.frame)
  results_1[[length(results_1) + 1]] <- result_1
}
write.xlsx(final_results_1,paste0("UC_model1.xlsx"),row.names = FALSE)

results_2 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_2 <- lapply(exp, cox_model_2)
  result_2 <- ldply(model_2,data.frame)
  results_2[[length(results_2) + 1]] <- result_2
}
write.xlsx(final_results_2,paste0("UC_model2.xlsx"),row.names = FALSE)

results_3 <- list()  
for (exp in exps) {
  print(exp)
  table(complete.cases(data3[, ..exp]))
  mydata <- data3[complete.cases(data3[, ..exp]), ]
  model_3 <- lapply(exp, cox_model_3)
  result_3 <- ldply(model_3,data.frame)
  results_3[[length(results_3) + 1]] <- result_3
}
write.xlsx(final_results_3,paste0("UC_model3.xlsx"),row.names = FALSE)
