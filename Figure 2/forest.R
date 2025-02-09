rm(list=ls())
library(tidyverse)
library(ggplot2)
library(ggstatsplot)
library(survival)
library(stringr)
library(viridis)
library(scales)
library(stringr)

load('EGAD00001008549_ici.circRNA_OS.result.Rdata')
ici_bi_best <- Binary.best.cox.result
ici_bi_median <- Binary.median.cox.result
ici_continuous <- continuous.cox.result
load('EGAD00001008549_nonici.circRNA_OS.result.Rdata')
nonici_bi_best <- Binary.best.cox.result
nonici_bi_median <- Binary.median.cox.result
nonici_continuous <- continuous.cox.result


df_info1 <- read.delim('EGAD00001008549_info.tsv')
df_rna1 <- read.delim('EGAD00001008549_circRNA_Ensemble.tsv')
names(df_rna1) <- gsub(".*?_EA_([^_]+).*", "EA-\\1", names(df_rna1))
df_info <- df_info1[df_info1$Treatment %in% 'Atezolizumab' & df_info1$Histology %in% 'SQUAMOUS', ]
pattern <- paste(df_info$Patient_ID, collapse = "|") 
cols_to_keep <- grep(pattern, colnames(df_rna1))
df_rna <- df_rna1[, cols_to_keep]
df_rna$id <- df_rna1$id #Ensemble

Sinfo <- data.frame(
  OS.time = df_info$OS_Time,
  OS = df_info$OS_Status
)
rownames(Sinfo) <- df_info$Patient_ID

gene = gsub("[^0-9]", "", plotCoxoutput$gene)
df_rna$sym = gsub("[^0-9]", "", df_rna$id)
expr <- df_rna[df_rna$sym %in% gene, ]
rownames(expr) <- expr$id
expr <- expr[1:nrow(df_info)]
colnames(expr) <- df_info$Patient_ID
expr <- expr[apply(expr, 1, function(x) {sum(x > 0) > 0.05*ncol(expr)}),]
##train
train <- cbind(Sinfo,data.frame(t(expr)))

##select gene
names(train)[c(1, 2)] <- c("times", "status") 

#test
df_info1 <- read.delim('EGAD00001008548_info.tsv')
df_rna1 <- read.delim('EGAD00001008548_circRNA_Ensemble.tsv')
names(df_rna1) <- gsub(".*?_EA_([^_]+).*", "EA-\\1", names(df_rna1))
df_info <- df_info1[df_info1$Treatment %in% 'Atezolizumab' & df_info1$Histology %in% 'SQUAMOUS', ]
pattern <- paste(df_info$Patient_ID, collapse = "|") 
cols_to_keep <- grep(pattern, colnames(df_rna1))
df_rna <- df_rna1[, cols_to_keep]
df_rna$id <- df_rna1$id #Ensemble

Sinfo <- data.frame(
  OS.time = df_info$OS_Time,
  OS = df_info$OS_Status
)
rownames(Sinfo) <- df_info$Patient_ID

gene = gsub("[^0-9]", "", plotCoxoutput$gene)
df_rna$sym = gsub("[^0-9]", "", df_rna$id)
expr <- df_rna[df_rna$sym %in% gene, ]
rownames(expr) <- expr$id
expr <- expr[1:nrow(df_info)]
colnames(expr) <- df_info$Patient_ID
expr <- expr[apply(expr, 1, function(x) {sum(x > 0) > 0.05*ncol(expr)}),]

##test
test <- cbind(Sinfo,data.frame(t(expr)))
names(test)[c(1, 2)] <- c("times", "status") 

matching_ids <- colnames(train) %in% colnames(test)
train <- train[matching_ids] 
matching_ids <- colnames(test) %in% colnames(train)
test <- test[matching_ids]
matching_ids <- plotCoxoutput$gene %in% colnames(train[-c(1,2)])  
plotCoxoutput <- plotCoxoutput[matching_ids, ] 
save(train,file='train.RData')
save(test,file='test.RData')

#提取CI值的上下限
numbers_list <- str_extract_all(plotCoxoutput$CI.x, "\\(([0-9.]+)-([0-9.]+)\\)")  
numbers_vector <- unlist(lapply(numbers_list, function(x) {  
  # 移除括号并分割字符串  
  gsub("[()]", "", x) %>%  
    str_split("-") %>%  
    unlist()  
})) 
upper <- list()  
lower <- list()  
# 遍历列表的索引  
for (i in seq_along(numbers_vector)) {  
  if (i %% 2 == 1) {  
    upper[[length(upper) + 1]] <- numbers_vector[[i]]  
  } else {  
    lower[[length(lower) + 1]] <- numbers_vector[[i]]  
  }  
}  
upper <- sapply(upper, as.numeric)
lower <- sapply(lower, as.numeric)

#####plot####
#让误差线跟着变色
ggplot(data=plotCoxoutput,aes(x=HR.x,y=gene,color=p.value.x))+
  geom_errorbarh(aes(xmax=upper,xmin=lower,color = p.value.x),height=0,size=1)+
  geom_point(aes(x=HR.x,y=gene),size=2,shape=16)+   #画成菱形
  geom_vline(xintercept = 1,linetype='dashed',size=0.5)+
  scale_x_continuous(breaks = c(0.5,1,1.5))+
  coord_trans(x='log2')+ 
  ylab("RNA")+  #标签
  xlab("Hazard ratios")+
  labs(color="P value",title ="Cox anlysis" )+
  scale_color_viridis()+ 
  theme_ggstatsplot()+  #好看的主题，同原文一致
  theme_bw(base_size = 12)+ 
  theme(axis.text.y = element_text(face="bold",  color="black", size=9),
        axis.title.x = element_text(face="bold", color="black", size=11),
        axis.title.y = element_text(face="bold",color="black", size=11),
        legend.text= element_text(face="bold", color="black", size=9),
        legend.title = element_text(face="bold", color="black", size=11),
        panel.border = element_rect(colour = 'black',size=0.5)) #去除网格线
ggsave('forest plot.pdf',width = 8,height = 9)
save(plotCoxoutput,file='significance.RData')
