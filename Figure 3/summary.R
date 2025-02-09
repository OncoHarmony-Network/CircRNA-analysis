rm(list=ls())
# 加载必要的库
library(ggplot2)

# 加载必要的库
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


p = ggplot(train, aes(x = time, y = auc, group = method, color = method)) +
  geom_line(linewidth = 1.5) +
  # 主题设置
  theme_classic() +
  # 轴设置
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 22)) +
  coord_cartesian(xlim = c(3.5, 24.5), ylim = c(0,1)) +
  ylab("AUC") +
  xlab("Time(m)") +
  scale_x_continuous(breaks = seq(0, 28, by = 4)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  # 颜色
  scale_color_manual(values = c('red', 'pink', 'blue', 'lightblue', 'darkgreen', 'lightgreen', 'purple', 'thistle','#A52A2A', '#D3B091')) +
  # 图例设置
  theme(
    legend.position = "top",  # 将图例放置在上方
    legend.direction = "horizontal"  # 横向排列图例
  )

pdf("method_train.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

p = ggplot(test, aes(x = time, y = auc, group = method, color = method)) +
  geom_line(linewidth = 1.5) +
  # 主题设置
  theme_classic() +
  # 轴设置
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 22)) +
  coord_cartesian(xlim = c(3.5, 24.5), ylim = c(0,1)) +
  ylab("AUC") +
  xlab("Time(m)") +
  scale_x_continuous(breaks = seq(0, 28, by = 4)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  # 颜色
  scale_color_manual(values = c('red', 'pink', 'blue', 'lightblue', 'darkgreen', 'lightgreen', 'purple', 'thistle','#A52A2A', '#D3B091')) +
  # 图例设置
  theme(
    legend.position = "top",  # 将图例放置在上方
    legend.direction = "horizontal"  # 横向排列图例
  )

pdf("method_test.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()

p = ggplot(validation, aes(x = time, y = auc, group = method, color = method)) +
  geom_line(linewidth = 1.5) +
  # 主题设置
  theme_classic() +
  # 轴设置
  theme(axis.text = element_text(size = 20)) +
  theme(axis.title = element_text(size = 22)) +
  coord_cartesian(xlim = c(3.5, 24.5), ylim = c(0,1)) +
  ylab("AUC") +
  xlab("Time(m)") +
  scale_x_continuous(breaks = seq(0, 28, by = 4)) + 
  scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
  # 颜色
  scale_color_manual(values = c('red', 'pink', 'blue', 'lightblue', 'darkgreen', 'lightgreen', 'purple', 'thistle','#A52A2A', '#D3B091')) +
  # 图例设置
  theme(
    legend.position = "top",  # 将图例放置在上方
    legend.direction = "horizontal"  # 横向排列图例
  )

pdf("method_validation.pdf",width = 10, height = 6)
print(p,newpage = FALSE)
dev.off()


result_train <- train %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))

result_test <- test %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))

result_validation <- validation %>%
  group_by(method) %>%
  summarise(mean_auc = mean(auc, na.rm = TRUE))
