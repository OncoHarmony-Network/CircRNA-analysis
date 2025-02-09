rm(list=ls())
library('pheatmap')
library('ggplot2')
library(org.Hs.eg.db)
library(clusterProfiler)

load('LM22_gene.RData')

# 将列名设置为列表的名称  
geneSet <- result_list
input <- read.csv('~/circRNA-lung cancer/figure 4/ssGSEA/anonymized_OAK-TPMs2.csv')

convertToNumeric <- function(df) {  
  # 遍历DataFrame的每一列  
  for (col in names(df)) {  
    # 尝试将列转换为数值型  
    # 注意：如果列中包含非数字字符，这些值将被转换为NA  
    df[[col]] <- as.numeric(as.character(df[[col]]))  
  }  
  # 返回转换后的DataFrame  
  return(df)  
}
input[-1] <- convertToNumeric(input[-1])

colnames(input)[1] <- 'SYMBOL'
final_df <- data.frame(
  symbol <- input$SYMBOL,
  gene <- input$EA.5e194ff13c
)
colnames(final_df) <- c('SYMBOL','Cells')
gene <- final_df$SYMBOL
#开始ID转换，会有丢失
gene=bitr(gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db") 
gene <- dplyr::distinct(gene,SYMBOL,.keep_all=TRUE)
data_all <- input %>% 
  inner_join(gene,by="SYMBOL")
rownames(data_all) <- data_all$ENTREZID;data_all <- data_all[-701]

library(GSVA)

gsvaPar <- ssgseaParam(as.matrix.data.frame(data_all[-1]),geneSet)
ssgsea <- gsva(gsvaPar)
## Estimating ssGSEA scores for 22 gene sets.
##   |=======================================================================================================================| 100%
load('~/circRNA-lung cancer/figure 4/ssGSEA/group.RDATA')
input <- data.frame(t(input))

a <- ssgsea %>% t() %>% as.data.frame()
rownames(oak) <- gsub("-", ".", rownames(oak)) 
dataframe_sorted <- oak[match(rownames(a), rownames(oak)), ]
a$group <- dataframe_sorted$num

write.table(a, "ssGSEA.txt", sep = "\t", row.names = T, col.names = T, quote = F)
######更改数字
a <- a[order(a$group), ]
heat_a <- data.frame(t(a))

convertToNumeric <- function(df) {  
  # 遍历DataFrame的每一列  
  for (col in names(df)) {  
    # 尝试将列转换为数值型  
    # 注意：如果列中包含非数字符，这些值将被转换为NA  
    df[[col]] <- as.numeric(as.character(df[[col]]))  
  }  
  # 返回转换后的DataFrame  
  return(df)  
}
heat_a <- convertToNumeric(heat_a)
group <- as.factor(a$group)
group <- data.frame(group)
rownames(group) <- colnames(heat_a)
annotation_colors = c(Low = "blue", High = "red1")
heat_a <- heat_a[-nrow(heat_a),]

pdf("heatmap.pdf", width=8, height=6)
pheatmap(heat_a,
         show_colnames = F,
         cluster_rows = F,
         cluster_cols = F,
         fontsize=5,  
         annotation = group)#画热图
dev.off()

