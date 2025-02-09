library(VennDiagram)
library(colorfulVennPlot)
library(ggplot2)
library(dplyr)
library(magrittr)
library(readr)
library(purrr)
library(RColorBrewer)
library(grDevices)
library(Cairo)
library(stringr)
library(tibble)
library(tidyr)

#导入数据
name_file <- c('circexplorer2', 'circRNA_finder','CIRI','find_circ')
filenames <- paste0("input_", name_file, ".txt")
data_ls <- filenames %>% map(., ~read_table(.x, col_names = FALSE)) %>% map(., ~.x$X1)
names(data_ls) <- paste0("G", 1:length(data_ls)) # 重命名
str(data_ls)

# 子区域数量
number_area <- 2^length(data_ls) - 1 # 

# 子区域编号
## 自定义一个函数，将整数转换成二进制字符串
intToBin <- function(x){
  if (x == 1)
    1
  else if (x == 0)
    NULL
  else {
    mod <- x %% 2
    c(intToBin((x-mod) %/% 2), mod)
  }
}

x_area <- seq(number_area) %>% map(., ~intToBin(.x)) %>%  # 转换成二进制字符串
  map_chr(., ~paste0(.x, collapse = "")) %>% 
  map_chr(., ~str_pad(.x, width = length(data_ls), side = "left", pad = "0"))

# 取data_ls中的向量取并集
G_union <- data_ls$G1 %>% union(data_ls$G2) %>% 
  union(data_ls$G3) %>% union(data_ls$G4)

# 自定义一个函数，计算子区域中元素的数量
area_calculate <- function(data_ls, character_area){
  character_num <- 1:4 %>% map_chr(., ~substr(character_area, .x, .x)) %>% 
    as.integer() %>% as.logical()
  
  element_alone <- G_union
  for (i in 1:4) {
    element_alone <- 
      if (character_num[i]) {
        intersect(element_alone, data_ls[[i]])
      } else {
        setdiff(element_alone, data_ls[[i]])
      }
  }
  return(element_alone)
}

# 调用自定义的函数，求各个子区域的元素
element_ls <- map(x_area, ~area_calculate(data_ls = data_ls, character_area = .x))

# 计算各个子区域中元素的数量：
quantity_area <- map_int(element_ls, length)

# 计算各个子区域中元素数量的百分比
percent_area <- (quantity_area / sum(quantity_area)) %>% round(3) # 四舍五入保留3位小数
percent_area <- (percent_area * 100) %>% paste0("%")

# 色板长度
length_pallete <- max(quantity_area) - min(quantity_area) + 1
color_area <- colorRampPalette(brewer.pal(n = 7, name = "YlGn"))(length_pallete) #其中YlGn可以换成其他颜色组合，例如Blues BuGn BuPu GnBu Greens Greys Oranges OrRd PuBu PuBuGn PuRd Purples RdPu Reds YlGn YlGnBu YlOrBr YlOrRd
color_tb <- tibble(quantity = scolor_tb <- tibble(quantity = scolor_tb <- tibble(quantity = seq(min(quantity_area), max(quantity_area), by = 1),
                   color = color_area))) 

nest1 <- tibble(quantity = quantity_area, percent = percent_area, area = x_area) %>% 
  group_by(quantity) %>% nest() %>% 
  left_join(color_tb$quantity$quantity, by = "quantity") %>% 
  arrange(quantity) %>% 
  unnest()
  
  # Venn图中显示数量
regions <- nest1$quantity
names(regions) <- nest1$area

CairoPDF(file = "venn_num.pdf", width = 8, height = 6)
plot.new()
plotVenn4d(regions, Colors = c('#C59DAF', '#AABEAD', '#BCB2C3', '#AABBA3','#E0B3B1','#C5A0B9','#C4D1A5','#C5DABC', 
                               '#BDA2A7','#FEE9C6','#C8ADAD','#C9AAA4','#BEA0A3','#E0BCC8','#BCBBD9'), Title = "", 
           labels = name_file,rot=45,shrink = 0.75) #从左至右与输入文件顺序对应
dev.off() # 关闭绘图以保存

# Venn图中显示百分比图
regions <- nest1$percent
names(regions) <- nest1$area

CairoPDF(file = "venn_percent.pdf", width = 8, height = 6)
plot.new()
plotVenn4d(regions, Colors = c('#C59DAF', '#AABEAD', '#BCB2C3', '#AABBA3','#E0B3B1','#C5A0B9','#C4D1A5','#C5DABC', 
                               '#BDA2A7','#FEE9C6','#C8ADAD','#C9AAA4','#BEA0A3','#E0BCC8','#BCBBD9'), Title = "", 
           labels = name_file,shrink = 0.75) #从左至右与输入文件顺序对应
dev.off() # 关闭绘图以保存
  
