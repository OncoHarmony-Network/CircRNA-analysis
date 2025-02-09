rm(list=ls())

input_oak <- read.csv('~/circRNA-lung cancer/figure 5/anonymized_OAK-TPMs2.csv')
gene_list <- c('RPL31',
               'RPS16',
               'EIF3G',
               'PLOD2',
               'FRAT2',
               'MYD88' ,
               'TNFSF4',
               'VNN1',
               'SLC12A8',
               'SNX10',
               'CHI3L1',
               'BCL2A1',
               'CFLAR')
input_oak <- input_oak[input_oak$id %in% gene_list, ]
input_oak <- data.frame(t(input_oak))
colnames(input_oak) <- input_oak[rownames(input_oak)=='id',];input_oak <- input_oak[-1,]
input_poplar <- read.csv('~/circRNA-lung cancer/figure 5/anonymized_POPLAR-TPMs2.csv')
input_poplar <- input_poplar[input_poplar$id %in% gene_list, ]
input_poplar <- data.frame(t(input_poplar))
colnames(input_poplar) <- input_poplar[rownames(input_poplar)=='id',];input_poplar <- input_poplar[-1,]

load('~/circRNA-lung cancer/figure 3/final_data/train.RData')
load('~/circRNA-lung cancer/figure 3/final_data/test.RData')
rownames(train) <- gsub("-", ".", rownames(train)) 
rownames(test) <- gsub("-", ".", rownames(test)) 
train$times <- train$times/365*12
test$times <- test$times/365*12

clin_dd<-train
validation<-test
names <- rownames(clin_dd)
set.seed(46)                  # 设置种子，以便后续结果重复
samp <- createDataPartition(Surv(clin_dd$times, clin_dd$status), p = 0.7, list = FALSE)  
train <- clin_dd[samp,]
test <- clin_dd[-samp,]

# 尝试将列转换为数值型  
convertToNumeric <- function(df) {  
  for (col in names(df)) {  
    df[[col]] <- as.numeric(as.character(df[[col]]))  
  }  
  return(df)  
}

train1 <- input_oak[rownames(input_oak) %in% rownames(train), ]
train <- cbind(train[c(1,2)],scale(convertToNumeric(train1)))
test1 <- input_oak[rownames(input_oak) %in% rownames(test), ]
test <- cbind(test[c(1,2)],scale(convertToNumeric(test1)))
validation1 <- input_poplar[rownames(input_poplar) %in% rownames(validation), ]
validation <- cbind(validation[c(1,2)],scale(convertToNumeric(validation1)))

# risk score
#calculate_risk_score <- function(data) {
#  risk_score <- (data$RPL31*0.010127) +
#                 (data$RPS16*0.007807)+ 
#                 (data$EIF3G*0.021535 )+
#                 (data$PLOD2*-0.02265)+ 
#                 (data$FRAT2*0.220046)+ 
#                 (data$MYD88*-0.02697)+ 
#                ( data$TNFSF4*-0.10415 )+
#                 (data$VNN1*-0.01664)+ 
#                 (data$SLC12A8*-0.03837)+ 
#                 (data$SNX10*-0.06819)+ 
#                 (data$CHI3L1*-0.03997)+ 
#                 (data$BCL2A1*-0.24106)+ 
#                 (data$CFLAR*-0.04284 )
#  return(risk_score)
#}

new_formula <- as.formula(paste("Surv(times, status) ~ ", paste(gene_list, collapse = " + "))) 
fit_cox <- coxph(new_formula, data = train)
summary(fit_cox)

train$risk_score <- predict(fit_cox,type='risk')
test$risk_score <- predict(fit_cox,newdata = test,type='risk')
validation$risk_score <- predict(fit_cox,newdata = validation,type='risk')

library(timeROC)
calculate_auc <- function(time, status, risk_score, time_points) {
  
  roc_results <- timeROC(
    T = time,                     # 生存时间
    delta = status,               # 状态（1: 事件发生, 0: 截尾）
    marker = risk_score,          # 风险评分
    cause = 1,                    # 关注的事件（默认是1）
    weighting = "marginal",       # 权重方式
    times = time_points,          # 特定时间点
    iid = TRUE                    # 是否计算置信区间
  )
  
  return(roc_results)
}

# 定义时间点
time_points <- c(6, 9, 12, 15, 18)

train_auc <- calculate_auc(
  time = train$times,
  status = train$status,
  risk_score = train$risk_score,
  time_points = time_points
)

test_auc <- calculate_auc(
  time = test$times,
  status = test$status,
  risk_score = test$risk_score,
  time_points = time_points
)

validation_auc <- calculate_auc(
  time = validation$times,
  status = validation$status,
  risk_score = validation$risk_score,
  time_points = time_points
)

all <- data.frame(train_auc$AUC, test_auc$AUC, validation_auc$AUC)
save(all,file='all.RData')
