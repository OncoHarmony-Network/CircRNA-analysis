rm(list=ls())
library(ggplot2)
library(ggrepel)
library(ggthemes)
#载入数据
#load('~/circRNA-lung cancer/figure 4/GSEA/volcano_plot/B_cell_Ma J-15/median/all.RData')
#B_cell <- all
#load('~/circRNA-lung cancer/figure 4/GSEA/volcano_plot/T_cell_Chu, Y/median/all.RData')
#T_cell <- all
#load('~/circRNA-lung cancer/figure 4/GSEA/volcano_plot/PC_Ma J-10/median/all.RData')
#PC <- all
#all <- rbind(B_cell,T_cell,PC)
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

x <- data.frame(
  all$NAME,
  all$NES,
  all$NOM.p.val
)
colnames(x) <- c('ID','NES','p_value')
x$label <- x$ID
x$p_value <- replace(x$p_value, x$p_value < 0.001, 0.001)

# 突出展示感兴趣的基因
selectgenes <- x$ID[x$p_value<0.05]

#plot_mode <- "classic" #经典版
plot_mode <- "advanced" #酷炫版

logFCcut <- 0 #log2-foldchange/NES
pvalCut <- 0.05 #P.value
adjPcut <- 0.05 #adj.P.value

#for advanced mode
#logFCcut2 <- NA
#logFCcut3 <- 5
#pvalCut2 <- 0.0001
#pvalCut3 <- 0.00001

#置x，y軸的最大最小位置
xmin <- -2
xmax <- 2
ymin <- 0
ymax <- 3

# 基因名的颜色，需大于等于pathway的数量，这里自定义了足够多的颜色
mycol <- c("darkgreen","chocolate4","blueviolet","#223D6C","#D20A13","#088247","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

if (plot_mode == "classic"){
  # 簡單的setting for color
  x$color_transparent <- ifelse((x$p_value < pvalCut & x$NES > logFCcut), "red", ifelse((x$p_value < pvalCut & x$NES < -logFCcut), "blue","grey30"))
  # 簡單的setting for size
  size <- ifelse((x$p_value < pvalCut & abs(x$NES) > logFCcut), 4, 2)
  
} else if (plot_mode == "advanced") {
  # 複雜的的setting for color
  n1 <- length(x[, 1])
  cols <- rep("grey30", n1)
  names(cols)<- rownames(x)
  
  #不同阈值的点的颜色
  #cols[x$p_value > pvalCut & x$NES >logFCcut]<- "#FB9A99"
  cols[x$p_value < pvalCut & x$NES > logFCcut]<- "#00DAE2"
  #cols[x$p_value > pvalCut & x$NES < logFCcut]<- "#B2DF8A"
  cols[x$p_value < pvalCut & x$NES < logFCcut]<- "#FF9389"
  color_transparent <- adjustcolor(cols, alpha.f = 0.5)
  x$color_transparent <- color_transparent
  
  # 複雜的的setting for size
  n1 <- length(x[, 1])
  size <- rep(1, n1)
  
  #不同阈值的点的大小
  size[x$p_value > pvalCut & x$NES >logFCcut]<- 2
  size[x$p_value < pvalCut & x$NES > logFCcut]<- 4
  size[x$p_value > pvalCut & x$NES < logFCcut]<- 2
  size[x$p_value < pvalCut & x$NES < logFCcut]<- 4
  
} else {
  stop("Unsupport mode")
}

# Construct the plot object
p1 <- ggplot(data=x, aes(NES, -log10(p_value), label = label, color = pathway)) +
  geom_point(alpha = 0.6, size = size, colour = x$color_transparent) +
  
  labs(x=bquote(~"NES"), y=bquote(~-Log[10]~italic("P-value")), title="") + 
  ylim(c(ymin,ymax)) + 
  xlim(c(xmin, xmax)) + 
  
  #画阈值分界线
  geom_vline(xintercept = c(-logFCcut, logFCcut), color="grey40", 
             linetype="longdash", lwd = 0.5) + #虚线的形状和粗细
  geom_hline(yintercept = -log10(pvalCut), color="grey40", 
             linetype="longdash", lwd = 0.5) +
  
  theme_bw(base_size = 12#, base_family = "Times" #修改字体
  ) +
  theme(panel.grid=element_blank())

#if (plot_mode == "advanced") {
#  p1 <- p1 + 
#    geom_vline(xintercept = c(-logFCcut2, logFCcut2), color="grey40", 
#               linetype="longdash", lwd = 0.5) +
#    geom_hline(yintercept = -log10(pvalCut2), color="grey40", 
#               linetype="longdash", lwd = 0.5)
#}
p1

# 显示 logFC > n 的基因的基因名
n = 0.05
genelist <- c(
'PID_CD8_TCR_PATHWAY',
'REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR',
'BIOCARTA_IL12_PATHWAY',
'WP_IL1_AND_MEGAKARYOCYTES_IN_OBESITY',
'PID_IL12_2PATHWAY',
'PID_IL27_PATHWAY',
'KEGG_MEDICUS_REFERENCE_TYPE_I_IFN_SIGNALING_PATHWAY',
'KEGG_NATURAL_KILLER_CELL_MEDIATED_CYTOTOXICITY',
'REACTOME_ACTIVATED_NTRK2_SIGNALS_THROUGH_CDK5',
'BIOCARTA_CDC25_PATHWAY'
)
p1 + geom_text_repel(aes(label = ifelse(ID %in% genelist, ID,"")),
                     colour="darkred", size = 1.8, box.padding = unit(0.4, "lines"), 
                     point.padding = unit(0.5, "lines"),
                     max.overlaps = 500)

# 突出显示候选基因
#p2 <- p1 + 
  # 在感兴趣的基因外面画个黑色圈
#  geom_point(data = selectgenes, alpha = 1, size = 4.6, shape = 1, 
#             stroke = 1, #圈粗细
#             color = "black") +
  
  # 显示感兴趣的基因的基因名
#  scale_color_manual(values = mycol) + 
#  geom_text_repel(data = selectgenes, 
#                  show.legend = FALSE, #不显示图例
#                  size = 5, box.padding = unit(0.35, "lines"), 
#                  point.padding = unit(0.3, "lines")) +
#  guides(color=guide_legend(title = NULL)) 

#p2

# 保存到PDF文件
ggsave("Volcano_advanced.pdf",width=8,height=6)
