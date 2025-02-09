library(ggplot2)

df_gene <- read.delim('data_gene.tsv')
df_num <- read.delim('data_num.tsv')
df_gene$X <- factor(df_gene$X,levels=c('1','2','3','4','5','6','7','8','9','10','>10'))

ggplot(df_num, aes(x = X, y = value)) +  
  geom_bar(stat = "identity", fill = "steelblue") +  
  labs(title = "", x = "Samples", y = "Number of circRNAs") +  
  theme_minimal()
ggsave("Barplot_num.pdf",width=10,height=5)

ggplot(df_gene, aes(x = X, y = value)) +  
  geom_bar(stat = "identity", fill = "steelblue") +   
  labs(title = "", x = "Number of circRNAs", y = "Number of genes") +  
  theme_minimal()
ggsave("Barplot_gene.pdf",width=10,height=5)