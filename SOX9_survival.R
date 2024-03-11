
# 03/21/2021
# JL
# SOX9 survival


library(ggplot2)
library(dplyr)
library(survival)
library(surviminer)
library(readxl)
library(reshape2)
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(Rcpp)
library(stringr)
library(reshape)
library(data.table)
library(ggsurvplot)
 

if(!require(devtools)) install.packages("devtools")
devtools::install_github("kassambara/survminer")

install.packages("survminer")
library("survminer")


options(stringsAsFactors = FALSE)
setwd("C:/")


# read the data
tcga_sox9 = read.table("Colon and Rectal Cancer_RNAseq.txt")
tcga_sox9_1 = t(tcga_sox9)
rownames(tcga_sox9_1) = colnames(tcga_sox9)
colnames(tcga_sox9_1) = rownames(tcga_sox9)
tcga_sox9_1 = data.table(tcga_sox9_1)
setnames(tcga_sox9_1, as.character(tcga_sox9_1[1,]))
tcga_sox9_1 = tcga_sox9_1[-1,]
colnames(tcga_sox9_1)[1] = "sampleID"

write.csv(tcga_sox9_1, 'tcga_colonal_rnaseq.csv')


tcga_surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
tcga_phen = read.table("Colon and rectal_clinicalmatrix.txt", header = T, sep = "\t")
colnames(tcga_surv)[1] <- "sampleID"


# select the primary tumor for RNAseq data
primary_tumor = filter(tcga_phen, tcga_phen$sample_type == "Primary Tumor")
primary_tumor = dplyr::select(primary_tumor, c("sampleID","sample_type"))
 

# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sampleID, SOX9))


tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- merge(tcga_survsox9, tcga_surv, by = "sampleID")

# normal tissue sox9
tumor_sox9 = select(tcga_survsox9, "sampleID", "SOX9")
sox9 = dplyr::select(tcga_sox9_2, c(sampleID, SOX9))
 
normal_sox9 = setdiff(sox9, tumor_sox9)
normal_sox9 = mutate(normal_sox9, sampleType = 1)


write.csv(normal_sox9, "normal_sox9.csv")
write.csv(tumor_sox9, "tumor_sox9.csv")
write.csv(tcga_survsox9, "tcga_sox9_surv.csv")


kruskal.test()

result = t.test(as.numeric(tumor_sox9$SOX9), as.numeric(normal_sox9$SOX9))
p.adjust(result$p.value, method = 'bonferroni', n = 460)

# p.adjust(p, method = p.adjust.methods, n = length(p))
# p.adjust.methods
# c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
 
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

# subset of first quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[96:380,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 1 )
tcga_survsox9b = tcga_survsox9[1:95,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 0)

# subset of bottom quartile of patients vs the rest of patients
tcga_survsox9a = tcga_survsox9[286:380,]  
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[1:285,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)
 
# subset of second quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[1:95,]  
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[189:285,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of third quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[1:95,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[95:190,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# for the 20% percentage survival
# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

# subset of 20 and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[305:380,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 1)
tcga_survsox9b = tcga_survsox9[1:76,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 0)

# subset of 40 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9  [1:76,]
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[228:305,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of 60 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9 [1:76,]
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[152:227,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of 80 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9 [1:76,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[76:152,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# sort by SOX9 RNAseq data

tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]


# for the 10% percentage survival
# sort by SOX9 RNAseq data

tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

tcga_survsox9 = read.csv('tcga_sox9_phen_surv2.csv')
#  
tcga_survsox9a = tcga_survsox9[39:380,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 1)
tcga_survsox9b = tcga_survsox9[1:38,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 0)

#  
tcga_survsox9a = tcga_survsox9 [1:38,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[305:342,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9 [1:38,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[266:304,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# 
tcga_survsox9a = tcga_survsox9 [1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[227:265,] 
tcga_survsox9b = mutate(tcga_survsox9b , pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9 [1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[190:227,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[152:190,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[114:152,] 
tcga_survsox9b = mutate(tcga_survsox9b , pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[76:114,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[38:76,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)


# combine data of RNAseq with survival
tcga_survsox9 = rbind(tcga_survsox9a, tcga_survsox9b)
 
for (i in 1:length(tcga_survsox9$DFI.time)) 
{
  if(tcga_survsox9$DFI.time[!is.na(tcga_survsox9$DFI.time)][[i]] > 1500 )
  {
    tcga_survsox9$DFI[i] = 0
  }
}

os <- Surv(time = as.numeric(tcga_survsox9$OS.time/30), event = as.numeric(tcga_survsox9$OS))
fit_km <- survfit(os ~ pecentage, data = tcga_survsox9)

# summary(fit_km)
plot(fit_km)

os_p <- ggsurvplot(
  fit_km,                     # survfit object with calculated statistics.
  data = tcga_survsox9,       # data used to fit survival curves.
  risk.table = TRUE,          # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 5,
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,150),          # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",     # customize X axis label.
  break.time.by = 20,        # break X axis in time intervals by 500.
  ggtheme = theme_light(),   # customize plot and risk table with a theme.
  # risk.table.y.text.col = T,# colour risk table text annotations.
  # risk.table.height = 0.25, # the height of the risk table
  # risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "",
  legend.labs = c("Low expression of SOX9", "High expression of sox9") # change legend labels.
   
)

# jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300 )
os_p
# dev.off()



# ------------ dfi survival ----------------

options(stringsAsFactors = FALSE)
setwd("C:/Users/lijinju/downloads/projects/Ying_project/")

# read the data
tcga_sox9 = read.table("Colon and Rectal Cancer_RNAseq.txt")
tcga_sox9_1 = t(tcga_sox9)
rownames(tcga_sox9_1) = colnames(tcga_sox9)
colnames(tcga_sox9_1) = rownames(tcga_sox9)
tcga_sox9_1 = data.table(tcga_sox9_1)
setnames(tcga_sox9_1, as.character(tcga_sox9_1[1,]))
tcga_sox9_1 = tcga_sox9_1[-1,]


tcga_surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
tcga_phen = read.table("Colon and rectal_clinicalmatrix.txt", header = T, sep = "\t")
colnames(tcga_surv)[1] <- "sampleID"
tcga_surv = subset( tcga_surv, select = c(sampleID, DFI, DFI.time))

# tcga_surv = mutate(tcga_surv, id = rownames(tcga_surv))
# select the primary tumor for RNAseq data

primary_tumor = filter(tcga_phen, tcga_phen$sample_type == "Primary Tumor")
primary_tumor = dplyr::select(primary_tumor, c("sampleID","sample_type"))
normal = filter(tcga_phen, tcga_phen$sample_type == "Solid Tissue Normal")
normal = dplyr::select(normal, c("sampleID","sample_type"))

# sort by SOX9 RNAseq data

tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9_tumor <- left_join(primary_tumor, tcga_survsox9, by = 'sampleID')
tcga_survsox9_tumor <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

tcga_survsox9_normal <- left_join(normal, tcga_survsox9, by = 'sampleID')
tcga_survsox9_normal <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]


# subset of first quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[286:380,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 1)
tcga_survsox9b = tcga_survsox9[1:95,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 0)

# subset of first quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[286:380,]  
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[1:285,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of first quartile and the 2nd quartile of patients
tcga_survsox9a = tcga_survsox9[1:95,]  
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[189:285,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of first quartile and the third quartile of patients
tcga_survsox9a = tcga_survsox9[1:95,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[95:190,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)


# combine data of RNAseq with survival
tcga_survsox9 = rbind(tcga_survsox9a, tcga_survsox9b)
tcga_survsox9 <- merge(tcga_survsox9, tcga_surv, by = "sampleID")

for (i in 1:length(tcga_survsox9$OS.time)) 
{
  if (tcga_survsox9$OS.time[i] > 800)
  {
    tcga_survsox9$OS[i] = 0
  }
}

os <- Surv(time = tcga_survsox9$DFI.time/30, event = tcga_survsox9$DFI)
fit_km <- survfit(os ~ pecentage, data = tcga_survsox9)

# summary(fit_km)
plot(fit_km)

os_p <- ggsurvplot(
  fit_km,                     # survfit object with calculated statistics.
  data = tcga_survsox9,       # data used to fit survival curves.
  risk.table = TRUE,          # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 5,
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,150),          # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",     # customize X axis label.
  break.time.by = 20,        # break X axis in time intervals by 500.
  ggtheme = theme_light(),   # customize plot and risk table with a theme.
  # risk.table.y.text.col = T,# colour risk table text annotations.
  # risk.table.height = 0.25, # the height of the risk table
  # risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,      # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "",
  legend.labs = c("Low expression of SOX9", "High expression of sox9") # change legend labels.
  
)

#jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300 )
os_p
#dev.off()


## for the 20% percentage survival
# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

# subset of 20 and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[305:380,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[1:76,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of 40 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9  [1:76,]
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[228:305,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of 60 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9 [1:76,]
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[152:227,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# subset of 80 quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9 [1:76,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[76:152,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]


## for the 10% percentage survival

# sort by SOX9 RNAseq data
tcga_survsox9 = read.csv('tcga_sur_phen_sox9all.csv')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]


#  
tcga_survsox9a = tcga_survsox9[1:38,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[39:376,]
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  9
tcga_survsox9a = tcga_survsox9 [339:376,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 1)
tcga_survsox9b = tcga_survsox9[1:38,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 0)

#  8
tcga_survsox9a = tcga_survsox9 [1:38,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[302:338,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  7
tcga_survsox9a = tcga_survsox9 [1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[265:301,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  6
tcga_survsox9a = tcga_survsox9 [1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[227:264,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  5
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[189:226,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  4
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[151:188,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  3
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[113:150,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  2
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[76:112,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)

#  1
tcga_survsox9a = tcga_survsox9[1:38,]   
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[38:75,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)


# combine data of RNAseq with survival
tcga_survsox9 = rbind(tcga_survsox9a, tcga_survsox9b)


#for(i in 1:length(tcga_survsox9$OS.time)) 
#{
#  if(tcga_survsox9$OS.time[i] > 800)
#  {
#    tcga_survsox9$OS[i] = 0
#  }
#}
os <- Surv(time = as.numeric(tcga_survsox9$DSS.time)/30, event = tcga_survsox9$DSS)
os <- Surv(time = as.numeric(tcga_survsox9$DFI.time)/30, event = tcga_survsox9$DFI)

os <- Surv(time = as.numeric(tcga_survsox9$OS.time)/30, event = tcga_survsox9$OS)
fit_km <- survfit(os ~ pecentage, data = tcga_survsox9)

# summary(fit_km)
plot(fit_km)

# png(filename="0 vs.OS_1.png", width = 6.5, height = 6.5, units = 'in', res=150) 

os_p <- ggsurvplot(
  fit_km,                     # survfit object with calculated statistics.
  size = 1.25,
  linetype = 1,
  data = tcga_survsox9,       # data used to fit survival curves.
  risk.table = FALSE,         # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 9,
  #pval.text = element_text(color="black", size=24,  angle = 60),
  #palette = c("#E7B800", "#2E9FDF"),
  palette = c("red", "blue"),
  xlim = c(0,150),           # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",     # customize X axis label.
  break.time.by = 30,        # break X axis in time intervals by 500.
  #ggtheme = theme_light(),  # customize plot and risk table with a theme.
  risk.table.pos = "out",
  risk.table.col = "black",
  risk.table.y.text = FALSE,
  #risk.table.y.text.col=TRUE,
  tables.theme = theme_cleantable(),
  font.tickslab = c(16),
  font.labs = "black",
  ggtheme = theme_classic2(base_size = 16, base_family = "Arial"),
  font.family = "Arial",
  #ggtheme = theme_classic(),
  # risk.table.y.text.col = T,# colour risk table text annotations.
  # risk.table.height = 0.25, # the height of the risk table
  # risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,    # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  #surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "",
  legend.labs = c("Low SOX9", "High SOX9"), # change legend labels.
  legend.coord = c(3,3),
  axis.ticks = element_line(size = 2),
  axis.line = element_line(size = 2),
  axis.labs = "black",
  axis.text.labs =  "black",
  cumevents.text.col = "black"
  
) 

 

#jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300 )
os_p
dev.off()




# --------------- Lymphovascular invasion --------------

options(stringsAsFactors = FALSE)
setwd("C:/Users/lijinju/downloads/projects/Ying_project/")


# read the data
tcga_sox9 = read.table("Colon and Rectal Cancer_RNAseq.txt")
tcga_sox9_1 = t(tcga_sox9)
rownames(tcga_sox9_1) = colnames(tcga_sox9)
colnames(tcga_sox9_1) = rownames(tcga_sox9)
tcga_sox9_1 = data.table(tcga_sox9_1)
setnames(tcga_sox9_1, as.character(tcga_sox9_1[1,]))
tcga_sox9_1 = tcga_sox9_1[-1,]


tcga_surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
tcga_phen = read.table("Colon and rectal_clinicalmatrix.txt", header = T, sep = "\t")
colnames(tcga_surv)[1] <- "sampleID"
tcga_surv = subset( tcga_surv, select = c(sampleID, PFI, PFI.time))
# tcga_surv = mutate(tcga_surv, id = rownames(tcga_surv))

tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))
colnames(tcga_survsox9)[1] = "sampleID"


# select the primary tumor for RNAseq data
primary_tumor = filter(tcga_phen, tcga_phen$sample_type == "Primary Tumor")
primary_tumor = dplyr::select(primary_tumor, c("sampleID","sample_type"))


# Lymphovascular invasion indicator
colnames(tcga_phen)[1] = "sampleID"
colnames(tcga_phen)[63] = "LymVInvasion"
tcga_phen = subset (tcga_phen, select = c(sampleID, LymVInvasion))


# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 = merge(tcga_survsox9, tcga_phen, by = "sampleID")

# lymphavasuclar invasion for "yes"
tcga_survsox9_y = filter(tcga_survsox9, tcga_survsox9$LymVInvasion == "YES")
tcga_survsox9_y = mutate(tcga_survsox9_y, pecentage = 0)

# lymphavasuclar invasion for "no"
tcga_survsox9_n = filter(tcga_survsox9, tcga_survsox9$LymVInvasion == "NO")
tcga_survsox9_n = mutate(tcga_survsox9_n, pecentage = 1)


tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]

# subset of 20 and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[11:102,] 
tcga_survsox9a = mutate(tcga_survsox9a, pecentage = 0)
tcga_survsox9b = tcga_survsox9[1:10,] 
tcga_survsox9b = mutate(tcga_survsox9b, pecentage = 1)


# combine data of RNAseq with survival
tcga_survsox9 = rbind(tcga_survsox9_y, tcga_survsox9_n)
tcga_survsox9 <- merge(tcga_survsox9, tcga_surv, by = "sampleID")

 

os <- Surv(time = tcga_survsox9$PFI.time/30, event = tcga_survsox9$PFI)
fit_km <- survfit(os ~ pecentage, data = tcga_survsox9)

# summary(fit_km)
plot(fit_km)

os_p <- ggsurvplot(
  fit_km,                     # survfit object with calculated statistics.
  size = 1.5,
  data = tcga_survsox9,       # data used to fit survival curves.
  risk.table = TRUE,          # show risk table.
  pval = TRUE,                # show p-value of log-rank test.
  conf.int = FALSE,           # show confidence intervals for 
  pval.coord = c(0.15,0.15),
  pval.size = 5,
  # palette = c("#E7B800", "#2E9FDF"),
  xlim = c(0,150),            # present narrower X axis, but not affect
  # survival estimates.
  xlab = "Time (Month)",     # customize X axis label.
  break.time.by = 20,        # break X axis in time intervals by 500.
  #ggtheme = theme_light(),   # customize plot and risk table with a theme.
  risk.table.pos="out",
  risk.table.col="black",
  risk.table.y.text.col=FALSE,
  tables.theme = theme_cleantable(),
  font.tickslab = c(12),
  ggtheme = theme_classic2(base_size=12, base_family = "Arial"),
  font.family = "Arial",
  # ggtheme = theme_classic(),
  # risk.table.y.text.col = T,# colour risk table text annotations.
  # risk.table.height = 0.25, # the height of the risk table
  # risk.table.y.text = FALSE,# show bars instead of names in text annotations
  # in legend of risk table.
  # ncensor.plot = TRUE,    # plot the number of censored subjects at time t
  # ncensor.plot.height = 0.25,
  conf.int.style = "step",  # customize style of confidence intervals
  surv.median.line = "hv",  # add the median survival pointer.
  legend.title = "",
  legend.labs = c("Low expression of SOX9", "High expression of sox9") # change legend labels.
  
)

# jpeg('km_aneu.jpeg',width = 10, height = 7.5, units = 'in',res=300 )
os_p
#dev.off()

### Cox PH models

# read the data
tcga_sox9 = read.table("Colon and Rectal Cancer_RNAseq.txt")
tcga_sox9_1 = t(tcga_sox9)
rownames(tcga_sox9_1) = colnames(tcga_sox9)
colnames(tcga_sox9_1) = rownames(tcga_sox9)
tcga_sox9_1 = data.table(tcga_sox9_1)
setnames(tcga_sox9_1, as.character(tcga_sox9_1[1,]))
tcga_sox9_1 = tcga_sox9_1[-1,]

# sort by SOX9 RNAseq data
tcga_sox9 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_sox9 = dplyr::select(tcga_sox9, c(sample, SOX9))
colnames(tcga_sox9)[1] = "sampleID"

tcga_phen = read.table("Colon and rectal_clinicalmatrix.txt", header = T, sep = "\t")
tcga_surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
tcga_phen = filter(tcga_phen, tcga_phen$sample_type == "Primary Tumor")

# select the primary tumor for RNAseq data

colnames(tcga_phen)[1] = "sampleID"
colnames(tcga_phen)[22] = "Age"
colnames(tcga_phen)[42] = "Gender"
colnames(tcga_phen)[63] = "Lymphatic_invasion"

tcga_phen = filter(tcga_phen, Lymphatic_invasion == "YES" | Lymphatic_invasion == "NO")
tcga_phen = select(tcga_phen, "sampleID", "Age","Gender", "Lymphatic_invasion")
colnames(tcga_surv)[1] = "sampleID"

tcga_cox = merge(tcga_surv, tcga_phen, by = "sampleID")
tcga_cox = merge(tcga_sox9, tcga_cox, by = "sampleID")

for(i in 1:length(tcga_cox$sampleID ))
{
  round(as.numeric(tcga_cox$SOX9[i]), digits = 1)
}

write.csv(tcga_cox, "tcga_cox_sox9_1.csv")
tcga_cox = read.csv("tcga_cox_sox9_1.csv")

write.csv(tcga_cox, "tcga_cox_sox9.csv")
tcga_cox = read.csv("tcga_cox_sox9.csv")


tcga_cox = read.csv('TCGA_phen_tumor_surv.csv')
tcga_cox = filter(tcga_cox, Lymphatic_invasion == "YES" | Lymphatic_invasion == "NO")
tcga_cox = filter(tcga_cox, Venous_invasion == "YES" | Venous_invasion == "NO")
tcga_cox = filter(tcga_cox, Perineural_invasion == "YES" | Perineural_invasion == "NO")
 

#tcga_cox  <- lapply( tcga_cox, function(x)  (as.character(x)))
os <- Surv(time = tcga_cox$OS.time/30,event = tcga_cox$OS)

cox <- coxph(os ~ `Age`  + `Gender` + `Lymphatic_invasion` + `Venous_invasion`,   data = tcga_cox)


cox <- coxph(os ~ `Age` + `Gender` + `Lymphatic_invasion` + `Perineural_invasion` + `Venous_invasion`, na.action = na.exclude,  data = tcga_cox)
cox <- coxph(os ~ `Age` + `Gender`  + `Venous_invasion`,  data = tcga_cox)
cox <- coxph(os ~ `Age` + `Gender`  + `Lymphatic_invasion`,  data = tcga_cox)
cox <- coxph(os ~ `Age` + `Gender`  + `Perineural_invasion`,  data = tcga_cox)
cox <- coxph(os ~ `Age` + `Gender`  + `Pathologic_stage`,  data = tcga_cox)

install.packages("devtools")
devtools::install_github("NightingaleHealth/ggforestplot")
library(ggforestplot)

summary(cox)
cox.zph(cox)
# plot(cox.zph(cox))
# ggforestplot::forestplot(cox,data=tcga_cox,fontsize = 1.8, noDigits = 4)
                     
cox_plot <-  ggforest(cox,data=tcga_cox,fontsize = 2,
                     noDigits = 4)
forest_model(cox)
# forestplot(cox)
# ('cox_aneu_hyper_rrmm_bm2.#',width = 10, height = 7.5)
cox_plot
# dev.off()

# translocation
# cox <- coxph(os ~ `High Aneuploidy` + `t(4;14)` + `t(11;14)`,
#             data = surv_aneu)
# summary(cox)
# cox.zph(cox)
# cox_plot <- ggforest(cox,data=surv_aneu,fontsize = 1.2,noDigits = 4)
# ('cox_aneu_trans_rrmm_bm2.#',width = 10, height = 7.5)
# cox_plot
# dev.off()
# write.csv(cox_df,"maxpro_aneusurv.csv")
# chi-sq test for lymphovascular invasion indicator 
# read the data
tcga_sox9 = read.csv("tcga_colonal_rnaseq.csv")

tcga_surv = read_excel("TCGA Colon and Rectal Cancer_survival.xlsx")
tcga_phen = read.table("Colon and rectal_clinicalmatrix.txt", header = T, sep = "\t")
colnames(tcga_surv)[1] <- "sampleID"
 


# select the primary tumor for RNAseq data
primary_tumor = filter(tcga_phen, tcga_phen$sample_type == "Primary Tumor")
primary_tumor = dplyr::select(primary_tumor, c("sampleID","sample_type"))


# sort by SOX9 RNAseq data
tcga_sox9_2 <- tcga_sox9_1[order(as.numeric(tcga_sox9_1$SOX9)),]
tcga_survsox9 = dplyr::select(tcga_sox9_2, c(sample, SOX9))

colnames(tcga_survsox9)[1] = "sampleID"

tcga_survsox9 <- merge(tcga_survsox9, primary_tumor, by = 'sampleID')
tcga_survsox9 <- tcga_survsox9[order(as.numeric(tcga_survsox9$SOX9)),]


# subset of first quartile and the bottom quartile of patients
tcga_survsox9a = tcga_survsox9[96:380,] 
tcga_survsox9b = tcga_survsox9[1:95,] 
 
# subset of 10%  and the bottom patients
tcga_survsox9a = tcga_survsox9[39:380,] 
tcga_survsox9b = tcga_survsox9[1:38,] 
 


colnames(tcga_phen)[1] = "sampleID"
colnames(tcga_phen)[22] = "Age"
colnames(tcga_phen)[42] = "Gender"
colnames(tcga_phen)[63] = "Lymphatic_invasion"
colnames(tcga_phen)[101] = "Tissue_Retrospective_Indicator"
colnames(tcga_phen)[56] = "Kras_Mutation"
colnames(tcga_phen)[93] = "sample_type"
colnames(tcga_phen)[94] = "sample_typeID"
colnames(tcga_phen)[83] = "perineural_invasion_present"

tcga_phen = select(tcga_phen, "sampleID","Age","Gender","Lymphatic_invasion","Tissue_Retrospective_Indicator","Kras_Mutation","sample_type","sample_typeID","perineural_invasion_present")
tcga_phen = filter(tcga_phen, Lymphatic_invasion == "YES" | Lymphatic_invasion == "NO")

tcga_sox9_phen = merge(tcga_survsox9, tcga_phen, by = "sampleID")
write.csv(tcga_sox9_phen, "tcga_sox9_phen_surv2.csv")

tcga_sursox_a_phen = merge(tcga_survsox9a, tcga_phen)
tcga_sursox_b_phen = merge(tcga_survsox9b, tcga_phen)

tcga_sursox_a_phen %>% count(Lymphatic_invasion)
tcga_sursox_b_phen %>% count(Lymphatic_invasion)

# sox9 divivded into 10 groups for cox 

tcga_sox9_cox = read.csv('tcga_sox9_phen_surv_RNAseq.csv')
colnames(tcga_sox9_cox)[2] = "sox9"
tcga_sox9_cox <- tcga_sox9_cox[order(as.numeric(tcga_sox9_cox$sox9)),]

# subset of 20 and the bottom quartile of patients

tcga_sox9_cox_10 = tcga_sox9_cox[1:37,] 
tcga_sox9_cox_10 = mutate(tcga_sox9_cox_10, SOX9 = '0')
tcga_sox9_cox_9 = tcga_sox9_cox[38:75,] 
tcga_sox9_cox_9 = mutate(tcga_sox9_cox_9, SOX9 = '1')
tcga_sox9_cox_8 = tcga_sox9_cox[76:112,] 
tcga_sox9_cox_8 = mutate(tcga_sox9_cox_8, SOX9 = '2')
tcga_sox9_cox_7 = tcga_sox9_cox[113:150,] 
tcga_sox9_cox_7 = mutate(tcga_sox9_cox_7, SOX9 = '3')
tcga_sox9_cox_6 = tcga_sox9_cox[151:188,] 
tcga_sox9_cox_6 = mutate(tcga_sox9_cox_6, SOX9 = '4')
tcga_sox9_cox_5 = tcga_sox9_cox[189:226,] 
tcga_sox9_cox_5 = mutate(tcga_sox9_cox_5, SOX9 = '5')
tcga_sox9_cox_4 = tcga_sox9_cox[227:264,] 
tcga_sox9_cox_4 = mutate(tcga_sox9_cox_4, SOX9 = '6')
tcga_sox9_cox_3 = tcga_sox9_cox[265:301,] 
tcga_sox9_cox_3 = mutate(tcga_sox9_cox_3, SOX9 = '7')
tcga_sox9_cox_2 = tcga_sox9_cox[302:338,] 
tcga_sox9_cox_2 = mutate(tcga_sox9_cox_2, SOX9 = '8')
tcga_sox9_cox_1 = tcga_sox9_cox[339:376,] 
tcga_sox9_cox_1 = mutate(tcga_sox9_cox_1, SOX9 = '9')


tcga_sox9_cox_10 = tcga_sox9_cox[1:37,] 
tcga_sox9_cox_10 = mutate(tcga_sox9_cox_10, SOX9 = "Low")
tcga_sox9_cox_1 = tcga_sox9_cox[38:376,] 
tcga_sox9_cox_1 = mutate(tcga_sox9_cox_1, SOX9 = "High")



tcga_survcox = rbind(tcga_sox9_cox_10,tcga_sox9_cox_9,tcga_sox9_cox_8,tcga_sox9_cox_7,tcga_sox9_cox_6,
                     tcga_sox9_cox_5,tcga_sox9_cox_4,tcga_sox9_cox_3,tcga_sox9_cox_2,tcga_sox9_cox_1)


tcga_survcox = rbind(tcga_sox9_cox_10,tcga_sox9_cox_1)


os <- Surv(time = as.numeric(tcga_survcox$OS.time)/30,event = tcga_survcox$OS)
os <- Surv(time = as.numeric(tcga_survcox$DSS.time)/30,event = tcga_survcox$DSS)
os <- Surv(time = as.numeric(tcga_survcox$DFI.time)/30,event = tcga_survcox$DFI)

cox <- coxph(os ~ `SOX9`, data = tcga_survcox)

summary(cox)

cox.zph(cox)

cox_plot <- ggforest(cox,data=tcga_survcox,fontsize = 1.9,
                     noDigits = 4)

cox_plot <- ggforest (cox, data=tcga_survcox, main = c("Hazard ratio"),
                      cpositions=c(0.02, 0.22, 0.4),
                      fontsize = 2.2, refLabel = "", noDigits=2) 

forest_model(cox)  
  
#('cox_aneu_hyper_rrmm_bm2.#',width = 10, height = 7.5)
cox_plot
#dev.off()

               

mat = rbind(c(17, 85), c(18,210))
result = fisher.test(mat)
result = chisq.test(mat) 
pval = signif(result$p.value,3)
mat
result

# X-squared = 4.8316, df = 1, p-value = 0.02794

tcga_sursox_a_phen %>% count(Tissue_Retrospective_Indicator)
tcga_sursox_b_phen %>% count(Tissue_Retrospective_Indicator)

mat = rbind(c(19, 177), c(19,144))
result = fisher.test(mat)
result = chisq.test(mat) 
pval = signif(result$p.value,3)
mat
result

# X-squared = 0.18448, df = 1, p-value = 0.6676

tcga_sursox_a_phen %>% count(Kras_Mutation)
tcga_sursox_b_phen %>% count(Kras_Mutation)

mat = rbind(c(20, 1), c(5,0), c(1,0))
result = fisher.test(mat)
result = chisq.test(mat) 
pval = signif(result$p.value,3)
mat
result

# X-squared = 0.2967, df = 2, p-value = 0.8621


tcga_sursox_a_phen %>% count(perineural_invasion_present)
tcga_sursox_b_phen %>% count(perineural_invasion_present)

mat = rbind(c(153, 13), c(52,6))
result = fisher.test(mat)
result = chisq.test(mat) 
pval = signif(result$p.value,3)
mat
result

# X-squared = 0.10095, df = 1, p-value = 0.7507


tcga_sursox_a_phen %>% count(Age)
tcga_sursox_b_phen %>% count(Age)
mean(tcga_sursox_a_phen$Age, na.rm = TRUE)
sd(tcga_sursox_a_phen$Age, na.rm = TRUE)

mean(tcga_sursox_b_phen$Age, na.rm = TRUE)
sd(tcga_sursox_b_phen$Age, na.rm = TRUE)

t.test(tcga_sursox_a_phen$Age, tcga_sursox_b_phen$Age, conf.level = 0.95)

# mean of x mean of y 
# 65.42899  56.60526  t = 3.9485, df = 45.339, p-value = 0.0002716


tcga_sursox_a_phen %>% count(Gender)
tcga_sursox_b_phen %>% count(Gender)

mat = rbind(c(149, 20), c(189,18))
result = fisher.test(mat)
result = chisq.test(mat) 
pval = signif(result$p.value,3)
mat
result

# X-squared = 0.69297, df = 1, p-value = 0.4052


 
