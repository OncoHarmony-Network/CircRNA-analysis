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
library('timeROC')
load('train.RData')
load('test.RData')

clin_dd<-train
validation<-test
names <- rownames(clin_dd)
set.seed(46)                  # 设置种子，以便后续结果重复
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

load('rna.RData')
expr <- train[,colnames(train) %in% rna]
expr <- expr*train$times
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

#####train#####
ttrain<-train
res.cut<-median(ttrain$Pre)
risk<-as.vector(ifelse(train$Pre>=res.cut,"high","low"))
ttrain$risk<-risk
ttrain2 <- ttrain
ttrain2$times <- ttrain$times/365*12

plot_km_os <- function(ttrain2){
  cox_model <- coxph(Surv(times, status) ~ risk, data = ttrain2)
  summ <- summary(cox_model)
  hr_values <- exp(coef(cox_model))  
  p_values <- summ$coefficients[, "Pr(>|z|)"]  
  lower_ci <- signif(summ$conf.int[3]) 
  upper_ci <- signif(summ$conf.int[4])  
  
  fit <- survfit(Surv(times, status) ~ risk,  
                 data = ttrain2)
  fit
  #计算p值
  
  p <- ggsurvplot(fit, risk.table = T,
                  pval = F,pval.method = F,conf.int = F, risk.table.title="",
                  surv.median.line = "hv", risk.table.col="strata", ylab="OS(%)", break.x.by = 4,
                  xlab="Time(months)", legend="none", palette = c("blue", "red"))
  p$plot<-p$plot + annotate("text", x = Inf, y = 1, label = paste("log-rank p=", sprintf("%.3f",p_values)),  
                            parse = FALSE, hjust = 1, vjust = 1.5,)+
    annotate("text", x = Inf, y = 0.9, label = paste("HR ", sprintf("%.3f",hr_values), '[95%CI:',sprintf("%.2f",lower_ci),'-',sprintf("%.2f",upper_ci),']'),  
             parse = FALSE, hjust = 1, vjust = 1.5,)
  p
}

plot_km_pfs <- function(ttrain){
  cox_model <- coxph(Surv(PFS_Time, PFS_Status) ~ risk, data = ttrain)
  summ <- summary(cox_model)
  hr_values <- exp(coef(cox_model))  
  p_values <- summ$coefficients[, "Pr(>|z|)"]  
  lower_ci <- signif(summ$conf.int[3]) 
  upper_ci <- signif(summ$conf.int[4]) 
  
  fit <- survfit(Surv(PFS_Time, PFS_Status) ~ risk,  
                 data = ttrain)
  fit
  #计算p值
  
  p <- ggsurvplot(fit, risk.table = T,
                  pval = F,pval.method = F,conf.int = F, risk.table.title="",
                  surv.median.line = "hv", risk.table.col="strata", ylab="PFS(%)", break.x.by = 4,
                  xlab="Time(months)", legend="none", palette = c("blue", "red"))
  p$plot<-p$plot + annotate("text", x = Inf, y = 1, label = paste("p_value=", sprintf("%.3f",p_values)),  
                            parse = FALSE, hjust = 1, vjust = 1.5,)+
    annotate("text", x = Inf, y = 0.9, label = paste("HR ", sprintf("%.3f",hr_values), '[95%CI:',sprintf("%.2f",lower_ci),'-',sprintf("%.2f",upper_ci),']'),  
             parse = FALSE, hjust = 1, vjust = 1.5,)
  p
}

p <- plot_km_os(ttrain2)
pdf("KMcurve_OS.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_info.tsv')
c<-df_info[df_info$Patient_ID %in% rownames(ttrain),]
ttrain$PFS_Time <- c$PFS_Time/365*12
ttrain$PFS_Status <- c$PFS_Status

p <- plot_km_pfs(ttrain)
pdf("KMcurve_PFS.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

#####test#####
ttest<-test
res.cut<-median(test$Pre)
ttest.risk<-as.vector(ifelse(test$Pre>=res.cut,"high","low"))
ttest$risk<-ttest.risk
#b 测试集不同时间点的ROC曲线
ttest2 <- ttest
ttest2$times <- ttest$times/365*12

p <- plot_km_os(ttest2)
pdf("KMcurve_OS_test.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

c<-df_info[df_info$Patient_ID %in% rownames(ttest),]
ttest$PFS_Time <- c$PFS_Time/365*12
ttest$PFS_Status <- c$PFS_Status

p <- plot_km_pfs(ttest)
pdf("KMcurve_PFS_test.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

data<-validation
res.cut<-median(data$Pre)
data.risk<-as.vector(ifelse(data$Pre>=res.cut,"high","low"))
data$risk<-data.risk
#b 测试集不同时间点的ROC曲线
data2 <- data
data2$times <- data$times/365*12

a<-survdiff(Surv(times, status) ~ risk,  
            data = data2)
p <- plot_km_os(data2)
pdf("KMcurve_OS_validation.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008548_info.tsv')
c<-df_info[df_info$Patient_ID %in% rownames(data2),]
data2$PFS_Time <- c$PFS_Time/365*12
data2$PFS_Status <- c$PFS_Status

p <- plot_km_pfs(data2)
pdf("KMcurve_PFS_validation.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

save(fit_cox,file='model.RData')

ROC.a <- timeROC(T=ttrain2$times, 
                 delta=ttrain2$status, marker=ttrain2$Pre,
                 cause=1,
                 weighting="marginal",
                 times=c(0,4,8,12,16,20,24),
                 iid=TRUE)

ROC.b <- timeROC(T=ttest2$times, 
                 delta=ttest2$status, marker=ttest2$Pre,
                 cause=1,
                 weighting="marginal",
                 times=c(0,4,8,12,16,20,24),
                 iid=TRUE)

ROC.c <- timeROC(T=data2$times, delta=data2$status,marker=data2$Pre,
                 cause=1,weighting="marginal",
                 times=c(0,4,8,12,16,20,24),
                 iid=TRUE)

#pdf("timeROC.pdf", 6, 5)
#plotAUCcurve(ROC.a, conf.int=FALSE, col="red")
#plotAUCcurve(ROC.b, conf.int=FALSE, col="darkblue", add=TRUE)
#plotAUCcurve(ROC.c, conf.int=FALSE, col="darkgreen", add=TRUE)

# 图例设置
#legend("topright", c("Train","Internal Validation","External Validation"),
#       col=c("red","darkblue","darkgreen"),
#       bty='n', lty=1, lwd=2, cex=0.8)
#dev.off()

train <- data.frame(
  time <- ROC.a$times,
  auc <- ROC.a$AUC
)
colnames(train) <- c('time','auc')
train$group <- 'train'
train <- train[-1,]

test <- data.frame(
  time <- ROC.b$times,
  auc <- ROC.b$AUC
)
colnames(test) <- c('time','auc')
test$group <- 'test'
test <- test[-1,]

validation <- data.frame(
  time <- ROC.c$times,
  auc <- ROC.c$AUC
)
colnames(validation) <- c('time','auc')
validation$group <- 'validation'
validation <- validation[-1,]

all <- rbind(train,test,validation)

library(scales)
mycols <- c('train' = "pink",'test' = "lightblue",'validation' = "lightgreen")
mycols1 <- c('train' = "red",'test' = "darkblue",'validation' = "darkgreen")

p <- ggplot(all, aes(x=all$time)) + 
  # 画柱形图
  geom_bar(aes(y = all$auc , fill = all$group), position=position_dodge(),stat = "identity", size = 10, width = 3) +
  scale_fill_manual(values = mycols) + #自定义bar的颜色
  
  # 画折线图
  #geom_point(aes(y = all$auc), shape = 21) +
  geom_line(aes(y = all$auc, group = all$group, color=all$group),size = 1.5 ) +
  scale_color_manual(values = mycols1) + 
  xlab('Time') +
  scale_x_continuous(breaks = c(4,8,12,16,20,24))+
  
  labs(x="") + 
  theme_bw() + #去除背景色
  theme(panel.grid = element_blank())  #去除网格线
pdf("timeROC.pdf", 6, 5)
print(p)
dev.off()
