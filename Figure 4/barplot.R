rm(list=ls())
library(ggplot2)

tsv_files <- list.files(pattern = "\\.tsv$")  
# 初始化一个空列表来存储数据框  
data_frames <- list()  
# 使用for循环读取每个文件  
for(file in tsv_files) {  
  # 读取每个文件到数据框  
  data_frames[[file]] <- read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE) 
}
high <- data_frames[[1]]
low <- data_frames[[2]]
all <- rbind(high,low)
all$group <- ifelse(all$NES >0, "high", "low")


x <- data.frame(
  all$NAME,
  all$NES,
  all$NOM.p.val,
  all$group
)
colnames(x) <- c('ID','NES','p_value','group')


#按照score排序
sortdf <- x[order(x$NES),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)
head(sortdf)


ggplot(sortdf, aes(ID, NES, fill = group)) + geom_bar(stat = 'identity') + 
  coord_flip() + 
  scale_fill_manual(values = c('#00DAE2', '#FF9389'), guide = FALSE) + 
  
  #写label
  geom_text(data = subset(x, NES < 0),
            aes(x=ID, y= -2, label= ID, color = group),#bar跟坐标轴间留出间隙
            size = 2, #字的大小
            hjust = 0 ) +  #字的对齐方式
  geom_text(data = subset(x, NES > 0),
            aes(x=ID, y= -2, label=ID, color = group),
            size = 2, hjust = 0) +
  geom_text(data = subset(x, NES < 0),
            aes(x=ID, y= 0, label= p_value, color = group),#bar跟坐标轴间留出间隙
            size = 3, #字的大小
            hjust = 'outward' ) +  #字的对齐方式
  geom_text(data = subset(x, NES > 0),
            aes(x=ID, y= 0, label=p_value, color = group),
            size = 3, hjust = 'outward') +  
  scale_colour_manual(values = c("black","black"), guide = FALSE) +
  xlab("") +ylab("NES")+
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank()) + #去除网格线
  theme(  
    panel.grid = element_blank(), # 去除网格线  
    panel.border = element_blank(), # 去除面板边框（但这里已经是默认行为）  
    axis.line.y = element_blank(), # 去除y轴  
    axis.ticks.y = element_blank(), # 去除y轴刻度  
    axis.text.y = element_blank(), # 去除y轴文本（如果您确实不想显示y轴标签）  
    axis.line.x = element_line(color = "black", size = 0.5) # 保留x轴并设置其样式  
  )

ggsave("gsva.pdf", width = 6, height = 8)
