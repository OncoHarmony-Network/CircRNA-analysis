####OS######
rm(list=ls())
library(data.table)
library(tidyverse)
library(survival)
library(progressr)
result_file <- paste0("EGAD00001008549_ici.circRNA_OS.result.RData")
# Check if the file already exists
## load gene information
info <- read.delim('EGAD00001008549_info.tsv')
###Clin OS
rownames(info)<-info$Run
clin<-info[info$Treatment=="Atezolizumab",]

clin<-clin[order(rownames(clin)),]
clin$times<-as.numeric(clin$OS_Time)
clin$status<-as.numeric(clin$OS_Status)
clin<-clin[,c("times","status")]
###expreset
dat <- fread('EGAD00001008549_circRNA_Ensemble.tsv',data.table = F)
exprSet_processed = dat %>%
  column_to_rownames(var = "id")
expreset<-data.frame(t(exprSet_processed))
# colnames(expreset)<-rownames(exprSet_processed)
expreset<-expreset[order(rownames(expreset)),]
expreset<-expreset[rownames(expreset) %in% rownames(clin),]
clin<-clin[rownames(clin) %in% rownames(expreset),]
###expreset
expreset<-expreset[,!is.na(as.numeric(colSums(expreset)))]
expreset<-expreset[,(as.numeric(colSums(expreset))>0)]
###### Unvariate Cox anlysis: continuous variable###### 
# cox anlysis according the expression of your gene，lncRNA or signatures 
library(foreach)
library(doParallel)
library(survival)
continuous.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    surv <- as.formula(paste('Surv(times, status)~', i))
    cur_cox <- coxph(surv, data = data)
    x <- summary(cur_cox)
    HR <- x$coefficients[i, "exp(coef)"]
    HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
    HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
    CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
    p.value <- x$coef[i, "Pr(>|z|)"]
    data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
  }
  return(result_list)
}
######  Unvariate Cox anlysis: set continuous variable as Binary (median)###### 
### cox anlysis according the expression of your gene，lncRNA or signatures
library(foreach)
library(doParallel)
library(survival)
Binary.median.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    data[,i] <- ifelse(data[,i] >= median(data[,i]), 1, 0)
    surv <- as.formula(paste('Surv(times, status)~', i))
    cur_cox <- coxph(surv, data = data)
    x <- summary(cur_cox)
    HR <- x$coefficients[i, "exp(coef)"]
    HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
    HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
    CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
    p.value <- x$coef[i, "Pr(>|z|)"]
    data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
  }
  return(result_list)
}

###### Unvariate Cox anlysis: set continuous variable as Binary (best cutoff)#########
#cox anlysis according the expression of your gene，lncRNA or signatures 
library(foreach)
library(doParallel)
library(survival)
library(survminer)

Binary.best.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    tryCatch(
      {
        # Omit specifying minprop and maxprop to use default values
        res.cut <- surv_cutpoint(data, time = "times", event = "status", variables = i)
        data[, i] <- ifelse(data[, i] >= res.cut[["cutpoint"]][["cutpoint"]], 1, 0)
        surv <- as.formula(paste('Surv(times, status)~', i))
        cur_cox <- coxph(surv, data = data)
        x <- summary(cur_cox)
        HR <- x$coefficients[i, "exp(coef)"]
        HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
        HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
        CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
        p.value <- x$coef[i, "Pr(>|z|)"]
        data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
      },
      error = function(e) {
        cat("Error in processing variable", i, ":", conditionMessage(e), "\n")
        NULL  # Return NULL for failed iteration
      }
    )
  }
  return(result_list)
}
# Set the number of cores you want to use
num_cores <- 1  # You can adjust this based on your system
# Register parallel backend with doParallel
with_progress(cl <- makeCluster(detectCores() - 1))
registerDoParallel(cl)
# Ensure that the survival library is loaded in each worker
clusterEvalQ(cl, {
  library(survival)
  library(survminer)
})
# Export necessary objects to parallel workers
clusterExport(cl, c("clin", "expreset"))
# Apply the function to all columns
continuous.cox.result <- continuous.cox_parallel(1:ncol(expreset))
Binary.median.cox.result <- Binary.median.cox_parallel(1:ncol(expreset))
Binary.best.cox.result <- Binary.best.cox_parallel(1:ncol(expreset))
# Stop the parallel backend
stopCluster(cl)
# Save the results only if the file doesn't exist
save(Binary.best.cox.result,continuous.cox.result,
     Binary.median.cox.result,
     file =result_file) 

result_file <- paste0("EGAD00001008549_nonici.circRNA_OS.result.RData")
# Check if the file already exists
## load gene information
info <- read.delim('EGAD00001008549_info.tsv')
###Clin OS
rownames(info)<-info$Run
clin<-info[info$Treatment=="Docetaxel",]

clin<-clin[order(rownames(clin)),]
clin$times<-as.numeric(clin$OS_Time)
clin$status<-as.numeric(clin$OS_Status)
clin<-clin[,c("times","status")]
###expreset
dat <- fread('EGAD00001008549_circRNA_Ensemble.tsv',data.table = F)
exprSet_processed = dat %>%
  column_to_rownames(var = "id")
expreset<-data.frame(t(exprSet_processed))
# colnames(expreset)<-rownames(exprSet_processed)
expreset<-expreset[order(rownames(expreset)),]
expreset<-expreset[rownames(expreset) %in% rownames(clin),]
clin<-clin[rownames(clin) %in% rownames(expreset),]
###expreset
expreset<-expreset[,!is.na(as.numeric(colSums(expreset)))]
expreset<-expreset[,(as.numeric(colSums(expreset))>0)]
###### Unvariate Cox anlysis: continuous variable###### 
# cox anlysis according the expression of your gene，lncRNA or signatures 
library(foreach)
library(doParallel)
library(survival)
continuous.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    surv <- as.formula(paste('Surv(times, status)~', i))
    cur_cox <- coxph(surv, data = data)
    x <- summary(cur_cox)
    HR <- x$coefficients[i, "exp(coef)"]
    HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
    HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
    CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
    p.value <- x$coef[i, "Pr(>|z|)"]
    data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
  }
  return(result_list)
}
######  Unvariate Cox anlysis: set continuous variable as Binary (median)###### 
### cox anlysis according the expression of your gene，lncRNA or signatures
library(foreach)
library(doParallel)
library(survival)
Binary.median.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    data[,i] <- ifelse(data[,i] >= median(data[,i]), 1, 0)
    surv <- as.formula(paste('Surv(times, status)~', i))
    cur_cox <- coxph(surv, data = data)
    x <- summary(cur_cox)
    HR <- x$coefficients[i, "exp(coef)"]
    HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
    HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
    CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
    p.value <- x$coef[i, "Pr(>|z|)"]
    data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
  }
  return(result_list)
}

###### Unvariate Cox anlysis: set continuous variable as Binary (best cutoff)#########
#cox anlysis according the expression of your gene，lncRNA or signatures 
library(foreach)
library(doParallel)
library(survival)
library(survminer)

Binary.best.cox_parallel <- function(num) {
  data = cbind(clin, expreset[, num])
  UniNames <- colnames(data)[-c(1:2)]
  result_list <- foreach(i = UniNames, .combine = rbind, .export = "num") %dopar% {
    tryCatch(
      {
        # Omit specifying minprop and maxprop to use default values
        res.cut <- surv_cutpoint(data, time = "times", event = "status", variables = i)
        data[, i] <- ifelse(data[, i] >= res.cut[["cutpoint"]][["cutpoint"]], 1, 0)
        surv <- as.formula(paste('Surv(times, status)~', i))
        cur_cox <- coxph(surv, data = data)
        x <- summary(cur_cox)
        HR <- x$coefficients[i, "exp(coef)"]
        HR.confint.lower <- signif(x$conf.int[i, "lower .95"], 3)
        HR.confint.upper <- signif(x$conf.int[i, "upper .95"], 3)
        CI <- paste0("(", HR.confint.lower, "-", HR.confint.upper, ")")
        p.value <- x$coef[i, "Pr(>|z|)"]
        data.frame(gene = i, HR = HR, CI = CI, p.value = p.value)
      },
      error = function(e) {
        cat("Error in processing variable", i, ":", conditionMessage(e), "\n")
        NULL  # Return NULL for failed iteration
      }
    )
  }
  return(result_list)
}
# Set the number of cores you want to use
num_cores <- 1  # You can adjust this based on your system
# Register parallel backend with doParallel
with_progress(cl <- makeCluster(detectCores() - 1))
registerDoParallel(cl)
# Ensure that the survival library is loaded in each worker
clusterEvalQ(cl, {
  library(survival)
  library(survminer)
})
# Export necessary objects to parallel workers
clusterExport(cl, c("clin", "expreset"))
# Apply the function to all columns
continuous.cox.result <- continuous.cox_parallel(1:ncol(expreset))
Binary.median.cox.result <- Binary.median.cox_parallel(1:ncol(expreset))
Binary.best.cox.result <- Binary.best.cox_parallel(1:ncol(expreset))
# Stop the parallel backend
stopCluster(cl)
# Save the results only if the file doesn't exist
save(Binary.best.cox.result,continuous.cox.result,
     Binary.median.cox.result,
     file =result_file)