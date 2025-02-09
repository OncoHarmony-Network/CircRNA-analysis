rm(list=ls())
# 加载R包，没有安装请先安装  install.packages("包名") 
library(ggplot2)
library(reshape2)
library(readxl)
library(RColorBrewer) 
library(ClassDiscovery)
library(cowplot)
library(grid)
library(dplyr)

# 读取簇状、堆积、填充柱形图数据文件
df <- read.delim('data_bar.tsv')
df_heat <- read.delim('data_heat.tsv')
annCol <- df_heat[,c(2,3)]

df_summarised <- df %>%  
  group_by(X) %>%  
  summarise(sum_value = sum(value), .groups = 'drop') 
df_summarised <- df_summarised %>% arrange(X);df_heat <- df_heat %>% arrange(X)
df_summarised$treatment <- df_heat$treatment
df_summarised$histology <- df_heat$histology
df_summarised <- df_summarised %>%  
  group_by(treatment) %>%  # 根据B列的值进行分组  
  arrange(sum_value)
df_summarised <- df_summarised %>% arrange(treatment)
df_heat <- df_summarised[-2]

df <- df %>% arrange(X)
order_index <- match(df_summarised$X, df$X)  
df1_sorted1 <- df[order_index, ]
df1_sorted1$sort <- 1:699
order_index <- order_index+1
df1_sorted2 <- df[order_index, ]
df1_sorted2$sort <- 1:699
order_index <- order_index+1
df1_sorted3 <- df[order_index, ]
df1_sorted3$sort <- 1:699
order_index <- order_index+1
df1_sorted4 <- df[order_index, ]
df1_sorted4$sort <- 1:699
df <- rbind(df1_sorted1,df1_sorted2,df1_sorted3,df1_sorted4)
df <- df%>%arrange(sort)

# 绘图

My_Theme1 <- theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_rect(size = 1)) + 
  theme(legend.position="right") +
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 12,angle = 90),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
  )+
  theme(plot.margin = margin(0.1,0.1,0,0.1, "cm"))

My_Theme2 = theme_minimal()+
  theme(legend.position="right") +
  theme(legend.margin = margin(0.5, 0.1, 0.5, 0.1),legend.text = element_text(size = 7)) +
  theme(legend.key.width = unit(1, "lines"), legend.key.height = unit(1, "lines"))+ # 调整图例大小
  theme(legend.title=element_blank())+
  theme(axis.title.x=element_blank(),
        axis.title.y=element_text(size = 6,angle = 0),
        axis.text = element_blank(),
        axis.ticks = element_blank())+
  theme (panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(plot.margin = margin(0.1,0.1,0,0.1, "cm"))

p = ggplot(df, aes(x=factor(X,levels =unique(X)),  # 转化为因子，目的是显示顺序与文件顺序相同，否则按照字母顺序排序
                   y=value, 
                   fill=factor(variable,levels = unique(variable)), 
))+
  labs(
    x="",   # 调整x轴名称
    y="",   # 调整y轴名称
    fill="" # 调整图例名称
  )
blues_palette <- brewer.pal(n = 4, name = "Blues")[rev(seq_along(brewer.pal(n = 4, name = "Blues")))]  #调色板中提取配色
p1 <- p +  geom_bar(
  position="stack",
  stat="identity",
  width = 1,
) + 
  scale_fill_manual(values = blues_palette)+  #配色
  My_Theme1

df_heat$X <- 1:nrow(df_heat)
p2 <- ggplot(df_heat,aes(X,1))+
  geom_tile(aes(fill = treatment))+
  My_Theme2+
  labs(y = "Treatment")+
  scale_fill_manual(values = c("#74C065","#EFF882"),labels = c('Atezolizumab','Docetaxel')) +
  scale_x_continuous(expand = c(0,0))#不留空

p3 <- ggplot(df_heat,aes(X,1))+
  geom_tile(aes(fill = histology))+
  My_Theme2+
  labs(y = "Histology",)+
  scale_fill_manual(values = c("#E00F0A","#3D6197")) +
  scale_x_continuous(expand = c(0,0)) #不留空

legend_a <- get_legend(p2+theme(legend.position = "top"))
legend_b <- get_legend(p3+theme(legend.position = "top"))
p <- plot_grid(p3,p2,p1,
               align = "v",axis = "lr",
               ncol = 1, rel_heights = c(1,1,8),
               legend_a,legend_b)
p
ggsave("Barplot.pdf",width=10,height=5)
