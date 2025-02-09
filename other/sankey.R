rm(list=ls())
library(ggalluvial)
df_info <- read.delim('~/circRNA-lung cancer/data/EGAD00001008549_info.tsv')
df <- df_info[c(6,7,5)]
#定义足够多的颜色，后面从这里选颜色
mycol <- rep(c("#E0367A","#D20A13","#FFD121","#5D90BA","#11AA4D","#58CDD9","#7A142C","#029149","#431A3D","#91612D","#6E568C","#D8D155","#64495D","#7CC767"),2)

#格式转换
UCB_lodes <- to_lodes_form(df[,1:ncol(df)],
                           axes = 1:ncol(df),
                           id = "Cohort")



ggplot(UCB_lodes,
       aes(x = x, stratum = stratum, alluvium = Cohort,
           fill = stratum, label = stratum)) +
  scale_x_discrete(expand = c(0, 0)) + 
  
  #用aes.flow参数控制线从哪边来，颜色就跟哪边一致。
  #默认是forward，此处用backward。
  geom_flow(width = 1/8,aes.flow = "backward") +
  
  coord_flip() + #旋转90度
  geom_stratum(alpha = .9,width = 1/10) +
  geom_text(stat = "stratum", size = 3,color="black") +
  
  #如果分组少，可以用scale_fill_brewer。修改type和palette两个参数选择配色方案，更多配色方案参考下图，分别对应type参数的“Seq”、“qual”、“Div”
  #scale_fill_brewer(type = "Div", palette = "BrBG") +
  
  #如果分组太多，就要用前面自己写的配色方案
  scale_fill_manual(values = mycol) +
  
  xlab("") + ylab("") +
  theme_bw() + 
  theme(panel.grid =element_blank()) + 
  theme(panel.border = element_blank()) + 
  theme(axis.line.x = element_blank(),axis.ticks = element_blank(),axis.text.x = element_blank()) + #显示分组名字
  ggtitle("") +
  guides(fill = FALSE) 

ggsave("sankey.pdf")