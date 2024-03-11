
# gene correlation of SOX9 with others
# JL


library(dplyr)
library(ggplot2)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(Rcpp)
library(biomaRt)
library(stringr)
library(reshape)
library(data.table)
library(ggbreak)
library(patchwork)
library(devtools)
library(easyGgplot2)
library(grid)
library(gridExtra)


setwd("C:/")
options(stringsAsFactors = FALSE)

tcga_sox9_1 = read.table('Colon and Rectal Cancer_RNAseq.txt')
tcga_sox9_1 = read.table('HiSeqV2.txt')


tcga_sox9 = read.csv('TCGA_colonRectal_RNAseq.csv')


tcga_sox9 = transpose(tcga_sox9_1)

# rownames(tcga_sox9) = tcga_sox9[,1]
# tcga_sox9 = tcga_sox9[,-1]
colnames(tcga_sox9) = tcga_sox9[1,]
tcga_sox9 = tcga_sox9[-1,]
colnames(tcga_sox9)[1] = "sampleID"

write.csv(tcga_sox9, "TCGA_colonRectal_RNAseq.csv")


###########################################################

rm(list = ls(all = TRUE))


tcga_sox9 = read.csv('TCGA_colonRectal_RNAseq.csv')


other_genes = dplyr::select(tcga_sox9, "sampleID","KRAS", "SOX9", "APC","PTEN",  'CTNNB1', 'BRAF'
                            ,'AXIN1', 'AXIN2','TCF7L2', "TP53",'FBXW7', 'NRAS', "SOX2" )
#sox9 = read.csv("spliteall.csv")
sox9 = read.csv("all_tcga_surv_phen.csv")
sox9 = sox9[,-1]
sox9 = dplyr::select(sox9, "sampleID")
other_genes = merge(other_genes, sox9, by = "sampleID")
 

 



cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$APC), method="pearson")


result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$APC), method="pearson")

r2 = round(result1$estimate ^ 2,3)
r = round(result1$estimate,3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes( as.numeric(other_genes$APC),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm,  linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 APC: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("APC") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  


plot(as.numeric(other_genes$APC), as.numeric(other_genes$SOX9), xlab = "SOX9", ylab ="APC", main = title )
model <- lm(as.numeric(other_genes$SOX9) ~ as.numeric(other_genes$APC))
abline(model)


# TP53

cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$TP53), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$TP53), method="pearson")

r2 = round(result1$estimate ^ 2,3)
r = round(result1$estimate,3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$TP53),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 TP53: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("TP53") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=32),
  axis.text.y = element_text(color="black", size=32),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  



plot(as.numeric(other_genes$SOX9), as.numeric(other_genes$TP53), xlab = "SOX9", ylab ="TP53", main = title )
model <- lm(as.numeric(other_genes$SOX9) ~ as.numeric(other_genes$TP53))

abline(model)


# PTEN
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$PTEN), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$PTEN), method="pearson")

r2 = round(result1$estimate ^ 2,2)
r = round(result1$estimate,2)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes( as.numeric(other_genes$PTEN),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 PTEN: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("PTEN") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  



# KRAS
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$KRAS), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$KRAS), method="pearson")

r2 = round(result1$estimate^ 2,2)
r = round(result1$estimate, 2)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$KRAS),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 KRAS: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("KRAS") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  



#CTNNB1  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$CTNNB1 ), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$CTNNB1 ), method="pearson")

r2 = round(result1$estimate^ 2,2)
r = round(result1$estimate, 2)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$CTNNB1 ),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 CTNNB1: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("CTNNB1") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  


# BRAF  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$BRAF), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$BRAF), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$BRAF),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 BRAF: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("BRAF") + ylab("SOX9")
p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  

# AXIN1  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$AXIN1), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$AXIN1), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$AXIN1),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 AXIN1: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("AXIN1") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  


# AXIN2  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$AXIN2), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$AXIN2), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$AXIN2),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 AXIN2: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("AXIN2") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  


# TCF7L2  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$TCF7L2), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$TCF7L2), method="pearson")

r2 = round(result1$estimate^ 2,2)
r = round(result1$estimate, 2)
p1 = printf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$TCF7L2),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 TCF7L2: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("TCF7L2") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  

# FBXW7  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$FBXW7), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$FBXW7), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$FBXW7),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 FBXW7: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("FBXW7") + ylab("SOX9")

p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=22, face="bold.italic", hjust = 0.5),
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  


# NRAS  
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$NRAS ), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$NRAS ), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$NRAS ),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 NRAS: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("NRAS") + ylab("SOX9")
p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0.5),
    
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  + coord_cartesian(xlim =c(9, 12.5) )
# + coord_cartesian(xlim =c(9, 12.5), ylim = c(8, 14))

# SOX2 
cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$SOX2 ), method="pearson")
result1 = cor.test(as.numeric(other_genes$SOX9), as.numeric(other_genes$SOX2 ), method="pearson")

r2 = round(result1$estimate^ 2,3)
r = round(result1$estimate, 3)
p1 = sprintf("%0.3g",result1$p.value)

p = ggplot(other_genes, aes(as.numeric(other_genes$SOX2 ),as.numeric(other_genes$SOX9))) + 
  geom_point(shape = 18, color = "blue") +
  geom_smooth(method = lm, linetype = "dashed",
              color = "darkred", se=FALSE) + ggtitle(paste0('SOX9 SOX2: p = ', p1, "; Cor Coefficient = ", r)) +
  xlab("SOX2") + ylab("SOX9")
p + theme( 
  axis.line = element_line(colour = "black"),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  panel.background = element_blank(),
  axis.line.x = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.line.y = element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks=element_line(colour = 'black', size=1.2, linetype='solid'),
  axis.ticks.length=unit(0.5, "lines"), 
  axis.text.x = element_text(color="black", size=36),
  axis.text.y = element_text(color="black", size=36),
  plot.title = element_text(color="black", size=20, face="bold.italic", hjust = 0.5),
  
  axis.title.x = element_text(color="black", size=36, face="bold"),
  axis.title.y = element_text(color="black", size=36, face="bold", angle = 90)
)  # + coord_cartesian(xlim =c(9, 12.5) )
# + coord_cartesian(xlim =c(9, 12.5), ylim = c(8, 14))
