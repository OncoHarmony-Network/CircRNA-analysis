rm(list=ls())
# 加载必要的库
library(tools)  # 用于 file_path_sans_ext 函数

# 设置主目录路径，即包含多个子文件夹的目录
main_directory <- "~/circRNA-lung cancer/figure 5"

# 获取主目录下的所有子文件夹
subfolders <- list.dirs(main_directory, recursive = FALSE, full.names = FALSE)

# 初始化空的dataframe，分别用于存储train, test, validation的数据
train <- data.frame(t6 = numeric(0), t12 = numeric(0), t18 = numeric(0))
test <- data.frame(t6 = numeric(0), t12 = numeric(0), t18 = numeric(0))
validation <- data.frame(t6 = numeric(0), t12 = numeric(0), t18 = numeric(0))

# 遍历每个子文件夹
for (folder in subfolders) {
  folder_path <- file.path(main_directory, folder)
  
  # 构建 all.RData 文件的完整路径
  rdata_path <- file.path(folder_path, "all.RData")
  
  # 检查文件是否存在
  if (file.exists(rdata_path)) {
    # 加载 all.RData 文件
    load(rdata_path)
    
    # 从 all 中提取对应的 auc 值
    new_row <- data.frame(
      t6 = all$train_auc.AUC[1],     
      t9 = all$train_auc.AUC[2],
      t12 = all$train_auc.AUC[3],
      t15 = all$train_auc.AUC[4],
      t18 = all$train_auc.AUC[3]     
    )
    
    # 将该行数据添加到 train dataframe
    train <- rbind(train, new_row)
    
    new_row <- data.frame(
      t6 = all$test_auc.AUC[1],     
      t9 = all$test_auc.AUC[2],
      t12 = all$test_auc.AUC[3],
      t15 = all$test_auc.AUC[4],
      t18 = all$test_auc.AUC[3]     
    )
    
    # 将该行数据添加到 test dataframe
    test <- rbind(test, new_row)
    
    new_row <- data.frame(
      t6 = all$validation_auc.AUC[1],     
      t9 = all$validation_auc.AUC[2],
      t12 = all$validation_auc.AUC[3],
      t15 = all$validation_auc.AUC[4],
      t18 = all$validation_auc.AUC[3]     
    )
    
    # 将该行数据添加到 validation dataframe
    validation <- rbind(validation, new_row)
    
    # 设置行名为文件夹名
    rownames(train)[nrow(train)] <- folder
    rownames(test)[nrow(test)] <- folder
    rownames(validation)[nrow(validation)] <- folder
  } else {
    # 如果没有all.RData文件，跳过该文件夹
    cat("Skipping folder:", folder, "- all.RData not found.\n")
  }
}

train$Model <- rownames(train)
test$Model <- rownames(test)
validation$Model <- rownames(validation)

library(dplyr)
train$Model <- factor(train$Model, levels = c("Cox","Cox_Binary","Time-cox","Time-cox_Binary","RF","RF_Binary","SVM","SVM_Binary",
                                              "XGBoost","XGBoost_Binary","checkpoint","CSIS", "CYT","Davoli","expanded immune", "IFN-γ",
                                              "IPRES","Roh immune","T effector"))
test$Model <- factor(test$Model, levels = c("Cox","Cox_Binary","Time-cox","Time-cox_Binary","RF","RF_Binary","SVM","SVM_Binary",
                                              "XGBoost","XGBoost_Binary","checkpoint","CSIS", "CYT","Davoli","expanded immune", "IFN-γ",
                                              "IPRES","Roh immune","T effector"))
validation$Model <- factor(validation$Model, levels = c("Cox","Cox_Binary","Time-cox","Time-cox_Binary","RF","RF_Binary","SVM","SVM_Binary",
                                              "XGBoost","XGBoost_Binary","checkpoint","CSIS", "CYT","Davoli","expanded immune", "IFN-γ",
                                              "IPRES","Roh immune","T effector"))
train <- train %>% arrange(Model)
test <- test %>% arrange(Model)
validation <- validation %>% arrange(Model)

# 加载必要的包
library(ggplot2)
library(tidyr)
library(dplyr)

# 转换为长格式数据
auc_long <- train %>%
  gather(key = "Time", value = "AUC", t6, t9, t12, t15, t18)

p1 = ggplot(auc_long, aes(x = Model, y = AUC)) +
  geom_boxplot() +   # 绘制箱线图
  geom_jitter(size = 1, width = 0.2, alpha = 0.6) +  # 添加数据点，避免重叠
  labs(x = "", y = "AUC", title = "AUC for Different Models") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签
        axis.text = element_text(size = 8),  # 调整坐标轴文字大小
        axis.title = element_text(size = 12),  # 调整坐标轴标题大小
        plot.title = element_text(size = 14, face = "bold")) +  # 调整标题大小
  coord_cartesian(ylim = c(0, 1))  # 限制y轴范围为0到1

# 转换为长格式数据
auc_long <- test %>%
  gather(key = "Time", value = "AUC", t6, t9, t12, t15, t18)

p2 = ggplot(auc_long, aes(x = Model, y = AUC)) +
  geom_boxplot() +   # 绘制箱线图
  geom_jitter(size = 1, width = 0.2, alpha = 0.6) +  # 添加数据点，避免重叠
  labs(x = "", y = "AUC", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签
        axis.text = element_text(size = 8),  # 调整坐标轴文字大小
        axis.title = element_text(size = 12)) +  # 调整坐标轴标题大小
  coord_cartesian(ylim = c(0, 1))  # 限制y轴范围为0到1

# 转换为长格式数据
auc_long <- validation %>%
  gather(key = "Time", value = "AUC", t6, t9, t12, t15, t18)

p3 = ggplot(auc_long, aes(x = Model, y = AUC)) +
  geom_boxplot() +   # 绘制箱线图
  geom_jitter(size = 1, width = 0.2, alpha = 0.6) +  # 添加数据点，避免重叠
  labs(x = "Model", y = "AUC", title = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),  # 旋转x轴标签
        axis.text = element_text(size = 8),  # 调整坐标轴文字大小
        axis.title = element_text(size = 12)) +  # 调整坐标轴标题大小
  coord_cartesian(ylim = c(0, 1))  # 限制y轴范围为0到1

# 保存为PDF文件，调整图形尺寸
pdf("Horizontal_comparison_train.pdf", width = 8, height = 3)
print(p1)
dev.off()

pdf("Horizontal_comparison_test.pdf", width = 8, height = 3)
print(p2)
dev.off()

pdf("Horizontal_comparison_validation.pdf", width = 8, height = 3)
print(p3)
dev.off()
