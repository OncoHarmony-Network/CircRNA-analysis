rm(list=ls())
library(survival)
library(survminer)
library("glmnet")

df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_info.tsv')
df_rna <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_circRNA_Ensemble.tsv')
load('model.RData')
id <- data.frame(lapply(df_rna$id, function(x) {  
  gsub("[^a-zA-Z0-9]", ".", x)  
}))  
df_rna$id <- t(data.frame(id))

load('~/circRNA-lung cancer/figure 3/final_data/11.4test/rna.RData')
genelist <- rna

expr <- df_rna[df_rna$id %in% genelist,]
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
new_df$Pre <- predict(fit_cox, newdata = new_df, proximity = T)$predicted  #RFSRC
new_df$PFS_time <- df_info$PFS_Time
new_df$PFS_Status <- df_info$PFS_Status
new_df$treatment <- df_info$Treatment
new_df$response <- df_info$Response


cut<-median(new_df$Pre)
#cut <- 1.118855
new_df$num <- ifelse(new_df$Pre > cut, 'High', 'Low')

new_df$OS_time <- new_df$OS_time/365*12
new_df$PFS_time <- new_df$PFS_time/365*12

fit <- survfit(Surv(OS_time,OS_Status) ~ num + treatment,  
               data = new_df)
fit

#计算p值
surdiff_hi <- survdiff(Surv(OS_time,OS_Status) ~ treatment, data = subset(new_df, num == "High"))
surdiff_lo <- survdiff(Surv(OS_time,OS_Status) ~ treatment, data = subset(new_df, num == "Low"))
p_hi <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
p_lo <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)

p1 <- ggsurvplot(fit, 
                 surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
                 pval = T,                          # 显示统计检验p值
                 pval.method = T,
                 conf.int = F,                      # 显示置信区间
                 risk.table = T,# 显示风险表格，展示每个时间点的风险数量
                 test.for.trend = F,
                 risk.table.title="",
                 risk.table.col = "strata",
                 xlab = "Time(months)",               # x轴标签
                 ylab = "OS%",     # y轴标签
                 legend.title = "Group",              # 图例标题
                 legend.labs = c("hi-Ate", "hi-Dtx",'lo-Ate','lo-Dtx'), # 图例标签
                 break.x.by = 4,                  # x轴刻度间隔
                 palette = c("blue1","lightblue","red1","pink"),   # 自定义颜色
                 fontsize = 3,
                 pval.size = 3,
                 tables.y.text = F,
                 tables.height = 0.3)

p1$plot<-p1$plot + annotate("text", x = Inf, y = 1, label = paste("hiAte-hiDtx p=", sprintf("%.3f",p_hi)),  
                            parse = FALSE, hjust = 1.1, vjust = 0.5)+
  annotate("text", x = Inf, y = 0.85, label = paste("loAte-loDtx p=", sprintf("%.3f",p_lo)),  
           parse = FALSE, hjust = 1.1, vjust = 0.5)
p1
pdf("KMcurve_OS_OAK.pdf",width = 10, height = 6)
print(p1,newpage = FALSE)
dev.off()

fit <- survfit(Surv(PFS_time,PFS_Status) ~ num + treatment,  
               data = new_df)
fit

#计算p值
surdiff_hi <- survdiff(Surv(PFS_time,PFS_Status) ~ treatment, data = subset(new_df, num == "High"))
surdiff_lo <- survdiff(Surv(PFS_time,PFS_Status) ~ treatment, data = subset(new_df, num == "Low"))
p_hi <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
p_lo <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)

p1 <- ggsurvplot(fit, 
                 surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
                 pval = T,                          # 显示统计检验p值
                 pval.method = T,
                 conf.int = F,                      # 显示置信区间
                 risk.table = T,# 显示风险表格，展示每个时间点的风险数量
                 test.for.trend = F,
                 risk.table.title="",
                 risk.table.col = "strata",         
                 xlab = "Time(months)",               # x轴标签
                 ylab = "PFS%",     # y轴标签
                 legend.title = "Group",              # 图例标题
                 legend.labs = c("hi-Ate", "hi-Dtx",'lo-Ate','lo-Dtx'), # 图例标签
                 break.x.by = 4,                  # x轴刻度间隔
                 palette = c("blue1","lightblue","red1","pink"),   # 自定义颜色
                 fontsize = 3,
                 pval.size = 3,
                 tables.y.text = FALSE,
                 tables.height = 0.3)

p1$plot<-p1$plot + annotate("text", x = Inf, y = 1, label = paste("hiAte-hiDtx p=", sprintf("%.3f",p_hi)),  
                            parse = FALSE, hjust = 1.1, vjust = 0.5)+
  annotate("text", x = Inf, y = 0.85, label = paste("loAte-loDtx p=", sprintf("%.3f",p_lo)),  
           parse = FALSE, hjust = 1.1, vjust = 0.5)
p1
pdf("KMcurve_PFS_OAK.pdf",width = 10, height = 6)
print(p1,newpage = FALSE)
dev.off()

df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008548_info.tsv')
df_rna <- read.delim('~/circRNA-lung cancer/data/EGAD00001008548_circRNA_Ensemble.tsv')
id <- data.frame(lapply(df_rna$id, function(x) {  
  gsub("[^a-zA-Z0-9]", ".", x)  
}))  
df_rna$id <- t(data.frame(id))

expr <- df_rna[df_rna$id %in% genelist,]
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

rownames(new_df) <- df_info$Patient_ID
new_df$OS_time <- df_info$OS_Time
new_df$OS_Status <- df_info$OS_Status
new_df$Pre <- predict(fit_cox, newdata = new_df, proximity = T)$predicted  #RFSRC
new_df$PFS_time <- df_info$PFS_Time
new_df$PFS_Status <- df_info$PFS_Status
new_df$treatment <- df_info$Treatment
new_df$response <- df_info$Response

cut<-median(new_df$Pre)
#cut <- 1.118855
new_df$num <- ifelse(new_df$Pre > cut, 'High', 'Low')

new_df$OS_time <- new_df$OS_time/365*12
new_df$PFS_time <- new_df$PFS_time/365*12

fit <- survfit(Surv(OS_time,OS_Status) ~ num + treatment,  
               data = new_df)
fit

#计算p值
surdiff_hi <- survdiff(Surv(OS_time,OS_Status) ~ treatment, data = subset(new_df, num == "High"))
surdiff_lo <- survdiff(Surv(OS_time,OS_Status) ~ treatment, data = subset(new_df, num == "Low"))
p_hi <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
p_lo <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)

p1 <- ggsurvplot(fit, 
                 surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
                 pval = T,                          # 显示统计检验p值
                 pval.method = T,
                 conf.int = F,                      # 显示置信区间
                 risk.table = T,# 显示风险表格，展示每个时间点的风险数量
                 test.for.trend = F,
                 risk.table.col = "strata",         
                 xlab = "Time(months)",               # x轴标签
                 ylab = "OS%",     # y轴标签
                 legend.title = "Group",              # 图例标题
                 legend.labs = c("hi-Ate", "hi-Dtx",'lo-Ate','lo-Dtx'), 
                 risk.table.title="",
                 break.x.by = 4,                  # x轴刻度间隔
                 palette = c("blue1","lightblue","red1","pink"),   # 自定义颜色
                 fontsize = 3,
                 pval.size = 3,
                 tables.y.text = FALSE,
                 tables.height = 0.3)

p1$plot<-p1$plot + annotate("text", x = Inf, y = 1, label = paste("hiAte-hiDtx p=", sprintf("%.3f",p_hi)),  
                            parse = FALSE, hjust = 1.1, vjust = 0.5)+
  annotate("text", x = Inf, y = 0.85, label = paste("loAte-loDtx p=", sprintf("%.3f",p_lo)),  
           parse = FALSE, hjust = 1.1, vjust = 0.5)
p1
pdf("KMcurve_OS_POPLAR.pdf",width = 10, height = 6)
print(p1,newpage = FALSE)
dev.off()

fit <- survfit(Surv(PFS_time,PFS_Status) ~ num + treatment,  
               data = new_df)
fit

#计算p值
surdiff_hi <- survdiff(Surv(PFS_time,PFS_Status) ~ treatment, data = subset(new_df, num == "High"))
surdiff_lo <- survdiff(Surv(PFS_time,PFS_Status) ~ treatment, data = subset(new_df, num == "Low"))
p_hi <- 1 - pchisq(surdiff_hi$chisq, length(surdiff_hi$n) - 1)
p_lo <- 1 - pchisq(surdiff_lo$chisq, length(surdiff_lo$n) - 1)

p1 <- ggsurvplot(fit, 
                 surv.median.line = "hv",           # 同时显示中位生存时间的垂直和水平参考线
                 pval = T,                          # 显示统计检验p值
                 pval.method = T,
                 conf.int = F,                      # 显示置信区间
                 risk.table = T,# 显示风险表格，展示每个时间点的风险数量
                 test.for.trend = F,
                 risk.table.col = "strata",         
                 xlab = "Time(months)",               # x轴标签
                 ylab = "PFS%",     # y轴标签
                 legend.title = "Group",              # 图例标题
                 legend.labs = c("hi-Ate", "hi-Dtx",'lo-Ate','lo-Dtx'), # 图例标签
                 risk.table.title="",
                 break.x.by = 4,                  # x轴刻度间隔
                 palette = c("blue1","lightblue","red1","pink"),   # 自定义颜色
                 fontsize = 3,
                 pval.size = 3,
                 tables.y.text = FALSE,
                 tables.height = 0.3)

p1$plot<-p1$plot + annotate("text", x = Inf, y = 1, label = paste("hiAte-hiDtx p=", sprintf("%.3f",p_hi)),  
                            parse = FALSE, hjust = 1.1, vjust = 0.5)+
  annotate("text", x = Inf, y = 0.85, label = paste("loAte-loDtx p=", sprintf("%.3f",p_lo)),  
           parse = FALSE, hjust = 1.1, vjust = 0.5)
p1
pdf("KMcurve_PFS_POPLAR.pdf",width = 10, height = 6)
print(p1,newpage = FALSE)
dev.off()
