rm(list=ls())
library("glmnet")
library('purrr')
library('dplyr')
library("survival")
library('pROC')
library("survminer")
library('tidyr')
library('survivalROC')
library('caret')
df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_info.tsv')
df_rna <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_circRNA_Ensemble.tsv')
df_info2 <- read.delim('~/circRNA-lung cancer/data/EGAD00001008548_info.tsv')
df_rna2 <- read.delim('~/circRNA-lung cancer/data/EGAD00001008548_circRNA_Ensemble.tsv')
id <- data.frame(lapply(df_rna$id, function(x) {  
  gsub("[^a-zA-Z0-9]", ".", x)  
}))  
df_rna$id <- t(data.frame(id))
id <- data.frame(lapply(df_rna2$id, function(x) {  
  gsub("[^a-zA-Z0-9]", ".", x)  
}))  
df_rna2$id <- t(data.frame(id))
colnames(df_rna2)[1] <- 'id'

p_value <- NULL
model <- list()
for (r in 1:100000) {
  
  load('train.RData')
  load('test.RData')
  clin_dd<-train[-c(3,19)]
  validation<-test[-c(3,19)]
  names <- rownames(clin_dd)                # 设置种子，以便后续结果重复
  set.seed(46)
  samp <- createDataPartition(Surv(clin_dd$times, clin_dd$status), p = 0.7, list = FALSE)    # 随机选择60%的数据作为模拟训练集，实际操作用训练集应该是一个来源的数据
  train <- clin_dd[samp,]
  test <- clin_dd[-samp,]
  ##scale
  test1<-test
  dd1<-scale(test[,-c(1:2)])
  sd(dd1[,1])
  test<-data.frame(test[,c(1:2)],dd1)
  train1<-train
  dd1<-scale(train[,-c(1:2)])
  sd(dd1[,1])
  train<-data.frame(train[,c(1:2)],dd1)
  validation1<-validation
  dd1<-scale(validation[,-c(1:2)])
  sd(dd1[,1])
  validation<-data.frame(validation[,c(1:2)],dd1)
  
  a<-train[-c(1:2)]
  b<-train[c(1:2)]
  
  cvfit = cv.glmnet(data.matrix(a), Surv(b$times,b$status), 
                    family = "cox",nfolds = 10,type.measure = "deviance"
  ) 
  
  gene_name <- gsub("\\....*$", "", colnames(a), perl = TRUE)
  fit <- glmnet(data.matrix(a), Surv(b$times,b$status), alpha=1,
                family = "cox") 
  c <- coef(fit,s = cvfit$lambda.min)
  rna <- rownames(c)[which(c!=0)]
  #用包自带的函数画图
  #########
  if(length(rna)<2 | length(rna)>15){
    next
  }
  
  expr <- train[,colnames(train) %in% rna]
  expr$status <- train$status
  expr$times <- train$times
  new_formula <- as.formula(paste("Surv(times, status) ~ ", paste(rna, collapse = " + "))) 
  fit_cox <- coxph(new_formula, data = expr)
  summary(fit_cox)
  
  train$Pre <- predict(fit_cox,type='risk') # 建立预测变量
  expr <- test[,colnames(test) %in% rna]
  expr$status <- test$status
  test$Pre <- predict(fit_cox,newdata = expr,type='risk') # 建立预测变量
  expr <- validation[,colnames(validation) %in% rna]
  expr$status <- validation$status
  validation$Pre <- predict(fit_cox,newdata = expr,type='risk') # 建立预测变量
  #####AUC_plot#####
  set_time <- c(10,12,14,16,18,20,22,24)
  auc_plot <- data.frame(
    time = set_time,
    train = numeric(length(set_time)),
    test = numeric(length(set_time)),
    validation = numeric(length(set_time)),
    stringsAsFactors = FALSE
  )
  #####train#####
  ttrain<-train
  res.cut<-median(train$Pre)
  risk<-as.vector(ifelse(train$Pre>=res.cut,"high","low"))
  ttrain$risk<-risk
  ttrain2 <- ttrain
  ttrain2$times <- ttrain$times/365*12
  survivalROC_helper <- function(t,data) {
    survivalROC(Stime=data$times, status=data$status, marker = data$Pre, 
                predict.time =t, method="NNE",span = 0.25*nrow(data)^(-0.20))
  }
  survivalROC_data <- data_frame(t = set_time) %>%
    mutate(survivalROC = map(t, data=ttrain2, survivalROC_helper),
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  for (i in 1:length(set_time)){
    num <- 10+(i-1)*2
    auc_plot$train[i] <- survivalROC_data$auc[survivalROC_data$t == num][1]
  }
  
  surdiff_os <- survdiff(Surv(times, status) ~ risk, data = ttrain)
  p_os <- 1 - pchisq(surdiff_os$chisq, length(surdiff_os$n) - 1)
  
  fit <- survfit(Surv(times, status) ~ risk,  
                 data = ttrain)
  fit
  
  #####test#####
  ttest<-test
  res.cut<-median(test$Pre)
  ttest.risk<-as.vector(ifelse(test$Pre>=median(test$Pre),"high","low"))
  ttest$risk<-ttest.risk
  #b 测试集不同时间点的ROC曲线
  ttest2 <- ttest
  ttest2$times <- ttest$times/365*12
  survivalROC_data <- data_frame(t = set_time) %>%
    mutate(survivalROC = map(t, data=ttest2,survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  for (i in 1:length(set_time)){
    num <- 10+(i-1)*2
    auc_plot$test[i] <- survivalROC_data$auc[survivalROC_data$t == num][1]
  }
  
  
  data<-validation
  res.cut<-median(data$Pre)
  data.risk<-as.vector(ifelse(data$Pre>=res.cut,"high","low"))
  data$risk<-data.risk
  #b 测试集不同时间点的ROC曲线
  data2 <- data
  data2$times <- data$times/365*12
  
  survivalROC_data <- data_frame(t = set_time) %>%
    mutate(survivalROC = map(t, data=data2,survivalROC_helper),
           ## Extract scalar AUC
           auc = map_dbl(survivalROC, magrittr::extract2, "AUC"),
           ## Put cut off dependent values in a data_frame
           df_survivalROC = map(survivalROC, function(obj) {
             as_data_frame(obj[c("cut.values","TP","FP")])
           })) %>%
    dplyr::select(-survivalROC) %>%
    unnest() %>%
    arrange(t, FP, TP)
  for (i in 1:length(set_time)){
    num <- 10+(i-1)*2
    auc_plot$validation[i] <- survivalROC_data$auc[survivalROC_data$t == num][1]
  }
  
  expr <- df_rna[df_rna$id %in% rna,]
  results <- list()  
  
  # 遍历df中唯一的id值  
  unique_ids <- unique(expr$id)  
  for (id in unique_ids) {  
    # 提取当前id的所有行  
    group_df <- expr[expr$id == id, ]  
    
    # 计算第9列到第707列的和  
    sums <- rowSums(group_df[, 8:706], na.rm = TRUE)  
    
    # 找出和最大的行的索引  
    max_index <- which.max(sums)  
    
    # 将和最大的行添加到结果列表中  
    results[[length(results) + 1]] <- group_df[max_index, ]  
  }  
  
  # 将结果列表合并为一个data.frame  
  expr <- data.frame(do.call(rbind, results))
  rownames(expr) <- expr$id
  expr <- expr[8:706]
  new_df <- data.frame(t(expr))
  
  new_df<-data.frame(scale(new_df))
  
  rownames(new_df) <- df_info$Patient_ID
  new_df$OS_time <- df_info$OS_Time
  new_df$OS_Status <- df_info$OS_Status
  new_df$treatment <- df_info$Treatment
  new_df$Pre <- predict(fit_cox,new_df,type='risk') # 建立预测变量
  
  
  cut<-median(new_df$Pre)
  new_df$num <- ifelse(new_df$Pre > cut, 'High', 'Low')
  
  fit <- survfit(Surv(OS_time,OS_Status) ~ num + treatment,  
                 data = new_df)
  fit
  
  if(length(unique(new_df$num))<2){
    next
  }
  
  #计算p值
  surdiff_hi <- survdiff(Surv(OS_time, OS_Status) ~ treatment, data = subset(new_df, num == "High"))
  surdiff_lo <- survdiff(Surv(OS_time, OS_Status) ~ treatment, data = subset(new_df, num == "Low"))
  p_hi <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
  p_lo <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)
  
  expr <- df_rna2[df_rna2$id %in% rna,]
  results <- list()  
  
  # 遍历df中唯一的id值  
  unique_ids <- unique(expr$id)  
  for (id in unique_ids) {  
    # 提取当前id的所有行  
    group_df <- expr[expr$id == id, ]  
    
    # 计算第9列到第707列的和  
    sums <- rowSums(group_df[, 8:199], na.rm = TRUE)  
    
    # 找出和最大的行的索引  
    max_index <- which.max(sums)  
    
    # 将和最大的行添加到结果列表中  
    results[[length(results) + 1]] <- group_df[max_index, ]  
  }  
  
  # 将结果列表合并为一个data.frame  
  expr <- data.frame(do.call(rbind, results))
  rownames(expr) <- expr$id
  expr <- expr[8:199]
  new_df <- data.frame(t(expr))
  
  new_df<-data.frame(scale(new_df))
  
  rownames(new_df) <- df_info2$Patient_ID
  new_df$OS_time <- df_info2$OS_Time
  new_df$OS_Status <- df_info2$OS_Status
  new_df$treatment <- df_info2$Treatment
  new_df$Pre <- predict(fit_cox,new_df,type='risk') # 建立预测变量
  
  
  cut<-median(new_df$Pre)
  new_df$num <- ifelse(new_df$Pre > cut, 'High', 'Low')
  
  fit <- survfit(Surv(OS_time,OS_Status) ~ num + treatment,  
                 data = new_df)
  fit
  
  #计算p值
  surdiff_hi <- survdiff(Surv(OS_time, OS_Status) ~ treatment, data = subset(new_df, num == "High"))
  surdiff_lo <- survdiff(Surv(OS_time, OS_Status) ~ treatment, data = subset(new_df, num == "Low"))
  p_hi_pop <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
  p_lo_pop <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)
  
  rna <- paste(rna, collapse = ",")
  p_value1 <- data.frame(time=r,p_hi,p_lo,p_hi_pop,p_lo_pop,mean(auc_plot$train),mean(auc_plot$test),mean(auc_plot$validation),rna)
  p_value <- rbind(p_value,p_value1)
  model[[r]] <- fit_cox
}
save(p_value,model,file='result.RData')

