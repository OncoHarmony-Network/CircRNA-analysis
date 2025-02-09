rm(list=ls())
library(ggplot2)
library(plyr)
library(dplyr)
library(tools)  # 用于file_path_sans_ext函数

# 设置主目录路径，即包含多个子文件夹的目录
main_directory <- "~/circRNA-lung cancer/figure 3/final_data/11.4test/plot"

# 获取主目录下的所有子文件夹
subfolders <- list.dirs(main_directory, recursive = FALSE, full.names = FALSE)

# 遍历每个子文件夹
for (folder in subfolders) {
  folder_path <- file.path(main_directory, folder)
  
  # 构建all.RData文件的完整路径
  rdata_path <- file.path(folder_path, "time.RData")
  load(rdata_path)
  all$method <- folder
  split_dfs <- split(all, all$group)
  if (folder == subfolders[1]){
    train <- split_dfs$train
    test <- split_dfs$test
    validation <- split_dfs$validation
  } else {
    train <- rbind(train,split_dfs$train)
    test <- rbind(test,split_dfs$test)
    validation <- rbind(validation,split_dfs$validation)
  }
}
# 构建示例数据
labs <- unique(train$method)

# 示例 AUC 数据
result <- train %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))
train_auc <- result$mean_auc
result <- test %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))
test_auc <- result$mean_auc
result <- validation %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))
validation_auc <- result$mean_auc

# 创建数据框
df <- data.frame(
  labs = rep(labs, 3),
  values = c(-train_auc, -test_auc, validation_auc),  # 左侧顺序：Test AUC -> Train Increment，右侧为正值
  group = rep(c("Train AUC", "Test AUC", "Validation AUC"), each = length(train_auc)),
  side = c(rep("Left", 2 * length(test_auc)), rep("Right", length(validation_auc)))
)

# 设置因子排序
df$group <- factor(df$group, levels = c("Train AUC", "Test AUC", "Validation AUC"))  # 修改图例顺序

# 绘制金字塔图
p = ggplot(df, aes(x = labs, y = values, fill = group)) +
  geom_bar(data = subset(df, side == "Left" & group == "Train AUC"), aes(y = values), stat = "identity") +  # Train AUC
  geom_bar(data = subset(df, side == "Left" & group == "Test AUC"), aes(y = values), stat = "identity") +  # Test AUC
  geom_bar(data = subset(df, side == "Right"), aes(y = values), stat = "identity") +  # Validation AUC
  geom_hline(yintercept = 0, colour = "white", lwd = 2) +
  theme_bw() +
  coord_flip() +
  scale_fill_manual(values = c("red", "pink", "lightblue")) +  # 修改颜色
  scale_y_continuous(breaks = seq(-1, 1, 0.2), labels = c(1, 0.8, 0.6, 0.4, 0.2, 0, 0.2, 0.4, 0.6, 0.8, 1),limits = c(-1, 1)) +
  scale_x_discrete(limits = labs) +  # 设置X轴显示的标签顺序
  labs(y = "AUC", x = "Method", fill = "Group") +  # 图例标题
 # ggtitle("Test + Train (Left)   Validation (Right)") +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),        # 去除网格线
    axis.text = element_text(size = 8),   # 调整坐标轴文字大小
    axis.title = element_text(size = 10)   # 调整坐标轴标题文字大小
  ) 
pdf("pyramid_plot.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()
