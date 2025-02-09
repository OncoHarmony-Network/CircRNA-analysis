rm(list=ls())

load('corr_input.RData')
df <- input
head(df)

colnames(df) <- gsub("\\.\\.\\..*", "", colnames(df))  
par(bg = "#fdfdfd")

# 左下角
panel.raters <- function (x, y, corr = NULL, ...) {
  if (!is.null(corr)) 
    return()
  plot.xy(xy.coords(x, y), type = "p", 
          pch = 20, #点形状
          cex = .5, #点大小
          ...)
  abline(lm(y ~ x), lwd = 2) #画拟合线
  box(col = "black", lwd = 2) #黑色粗边框
}

# 对角线
textPanel <- function (x = 0.5, y = 0.5, txt, cex, font) {
  text(x, y, txt, cex = cex, font = font)
  box(col = "black", lwd = 2)
}

# 右上角
panel.fill.cor <- function (x, y, corr = NULL, ...) 
{
  # 计算相关系数
  corr <- round(cor(x, y, use = "pairwise", method = "pearson"),2) # 可以换成"kendall"或 "spearman"
  
  # 自定义背景颜色
  ncol <- 14
  col.regions <- colorRampPalette(c('darkslateblue', 'navy', 'white', 'firebrick3', 'red'))
  pal <- col.regions(ncol)
  col.ind <- as.numeric(cut(corr, breaks = seq(from = -1, to = 1, length.out = ncol + 1), include.lowest = TRUE))
  
  # 画背景
  par(new=TRUE)
  plot(0, type='n', xlim=c(-1,1), ylim=c(-1,1), axes=FALSE, asp=1)
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[4], col = pal[col.ind], 
       border = NA)
  
  # 写相关系数
  text(0, 0, labels = corr, cex = 2.5, col = ifelse(corr > 0, "black", "white"))
  box(col = "black") #黑色窄边框
}

# 画图并保存到pdf文件
pdf("corrgram0.pdf",13,13)
pairs(df[1:12], 
      gap = .5, #小图之间的空隙
      text.panel = textPanel, #对角线
      lower.panel = panel.raters, #左下角
      upper.panel = panel.fill.cor) #右上角

dev.off()
