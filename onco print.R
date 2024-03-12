# 08252022
# JL


rm(list = ls(all = TRUE))

library(dplyr)
library(ComplexHeatmap)
library(tidyverse)
library(reshape)
library(reshape2)
library(data.table)


setwd("C:/")



# set up gene of dataframe
# wnt pathway: APC, CTNNB1, RNF43, RSPO2, ZNRF3, 

df = read.csv("df2.csv")
df = df[,-1]

sox10 = read.csv('tcga_sur_phen_10.csv')
sox90 = read.csv('tcga_sur_phen_90.csv')
sox9 = read.csv('tcga_sox9_RNAseq_tumor_normal.csv')

onco = read.csv("onco_unique.csv")

soxonco = left_join(df, onco, by = "name")
soxonco = distinct(soxonco)
soxonco = data.table(soxonco)

onco6 = reshape2::dcast(soxonco,sample.x~gene.x, value.var = 'effect' ) 

colnames(onco6)[1] = "sampleID"

sox10onco = merge(sox10, onco6, by = "sampleID" )
sox90onco = merge(sox90, onco6, by = "sampleID" )

# a = sox10onco[,c('APC','CTNNB1','AXIN1','TCF7L2','AXIN2','FBXW7',"RNF43",'RSPO2','ZNRF3','SOX9')]

#################################################################

library(tidyverse)

sox10onco = sox10onco[,-(2:21)]
rownames(sox10onco) = sox10onco[, 1]
sox10onco = sox10onco[,-1]
sox10onco = sox10onco %>% dplyr::rename(SOX9 = SOX9.y)
  
sox10onco1 = t(sox10onco)

sox10onco1 = subset (sox10onco1, rownames(sox10onco1) %in% c('TP53','APC','CTNNB1','AXIN1','AXIN2',"TP53",'BRAF', 
               'NOTCH3', 'NCOR1','NCOR2', 'RPS6KA2', 'CTBP1', 'CTBP1', 'AR', 'ERBB3', 'FN1', 'MAPK3K4', 'MAPK14',
               'MAPK6', 'MAPK7', 'MMP16', 'NFKB1', 'SMAD3','TGFBR2', 'VEGFA',
               'TCF7L2','FBXW7',"RNF43",'RSPO2','ZNRF3','SOX9','KRAS','PTEN', 'MTOR', 'NOTCH3', 'BRAF', 'EP300',
               'PIK3R1', 'RPS6KA2', 'AKT1', 'AKT2', 'AR', 'E2F5', 'E2F8','NFKB1','RB1', 'PIK3CA','SMAD2', 'SMAD3', 'SMAD4',
               'TGFBR2','TSC2', 'BRCA2', 'NF1', 'NRAS', 'STAT1', 'NCOR1', 'NCOR2', 'GSK3B', 'TGFBR2', 'NFKB1' 
               ))

sox90onco = sox90onco[,-(2:21)]
rownames(sox90onco) = sox90onco[, 1]
sox90onco = sox90onco[,-1]

sox90onco = sox90onco %>% dplyr::rename( SOX9 = SOX9.y)
  
sox90onco1 = t(sox90onco)

sox90onco1 = subset (sox90onco1, rownames(sox90onco1) %in% c('TP53','APC','CTNNB1','AXIN1','AXIN2',"TP53",'BRAF', 
                                                             'NOTCH3', 'NCOR1','NCOR2', 'RPS6KA2', 'CTBP1', 'CTBP1', 'AR', 'ERBB3', 'FN1', 'MAPK3K4', 'MAPK14',
                                                             'MAPK6', 'MAPK7', 'MMP16', 'NFKB1', 'SMAD3','TGFBR2', 'VEGFA',
                                                             'TCF7L2','FBXW7',"RNF43",'RSPO2','ZNRF3','SOX9','KRAS','PTEN', 'MTOR', 'NOTCH3', 'BRAF', 'EP300',
                                                             'PIK3R1', 'RPS6KA2', 'AKT1', 'AKT2', 'AR', 'E2F5', 'E2F8','NFKB1','RB1', 'PIK3CA','SMAD2', 'SMAD3', 'SMAD4',
                                                             'TGFBR2','TSC2', 'BRCA2', 'NF1', 'NRAS', 'STAT1', 'NCOR1', 'NCOR2', 'GSK3B', 'TGFBR2', 'NFKB1' 
))


# (onco,(onco$effect == "Missense_Mutation" | onco$effect == "Frame_Shift_Ins"| onco$effect == "Frame_Shift_Del"|
#       onco$effect == "In_Frame_Del" | onco$effect == "Splice_Site" | onco$effect == "Translation_Start_Site" | onco$effect == "large deletion" |
#       onco$effect == "In_Frame_Ins" | onco$effect == "Nonstop_Mutation" ))

for (i in 1:dim(sox10onco1)[1])
{
  for (j in 1:dim(sox10onco1)[2])
    
  if(!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Missense_Mutation')
  {
    sox10onco1[i,j] = 'MUT'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Frame_Shift_Ins')
  {
    sox10onco1[i,j] = 'FRA'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Frame_Shift_Del')
  {
    sox10onco1[i,j] = 'FRA'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'In_Frame_Del')
  {
    sox10onco1[i,j] = 'INDEL'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Splice_Site')
  {
    sox10onco1[i,j] = 'SPL'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Translation_Start_Site')
  {
    sox10onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'large deletion')
  {
    sox10onco1[i,j] = 'DEL'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'In_Frame_Ins')
  {
    sox10onco1[i,j] = 'INDEL'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Nonstop_Mutation')
  {
    sox10onco1[i,j] = 'NonSense'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Nonsense_Mutation')
  {
    sox10onco1[i,j] = 'NonSense'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Silent')
  {
    sox10onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == '3Flank')
  {
    sox10onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == '5UTR')
  {
    sox10onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == '3UTR')
  {
    sox10onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'Intron')
  {
    sox10onco1[i,j] = 'Intron'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == 'RNA')
  {
    sox10onco1[i,j] = 'MUT'
  }
  else if (!is.na(sox10onco1[i,j]) & sox10onco1[i,j] == '5Flank')
  {
    sox10onco1[i,j] = 'OTHER'
  }
}

for (i in 1:dim(sox90onco1)[1])
{
  for (j in 1:dim(sox90onco1)[2])
    
  if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Missense_Mutation')
  {
     sox90onco1[i,j] = 'MUT'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Frame_Shift_Ins')
  {
    sox90onco1[i,j] = 'FRA'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Frame_Shift_Del')
  {
    sox90onco1[i,j] = 'FRA'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'In_Frame_Del')
  {
    sox90onco1[i,j] = 'INDEL'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Splice_Site')
  {
    sox90onco1[i,j] = 'SPL'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Translation_Start_Site')
  {
    sox90onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'large deletion')
  {
    sox90onco1[i,j] = 'DEL'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'In_Frame_Ins')
  {
    sox90onco1[i,j] = 'INDEL'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Nonstop_Mutation')
  {
    sox90onco1[i,j] = 'NonSense'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Nonsense_Mutation')
  {
    sox90onco1[i,j] = 'NonSense'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Silent')
  {
    sox90onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == '3Flank')
  {
    sox90onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == '5UTR')
  {
    sox90onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == '3UTR')
  {
    sox90onco1[i,j] = 'OTHER'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'Intron')
  {
    sox90onco1[i,j] = 'Intron'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == 'RNA')
  {
    sox90onco1[i,j] = 'MUT'
  }
  else if (!is.na(sox90onco1[i,j]) & sox90onco1[i,j] == '5Flank')
  {
    sox90onco1[i,j] = 'OTHER'
  }
}

# In order to clearly visualize the data, genes for which less than 10% patients have mutations and patients which have 
# mutations in less than 5% of the genes are removed.
# l1 = apply(sox10onco1, 1, function(x) sum(!grepl("^\\s*$", x))/length(x) > 0.1)
# l2 = apply(sox10onco1, 2, function(x) sum(!grepl("^\\s*$", x))/length(x) > 0.5)
# sox10onco1 = sox10onco1[l1, l2]

anno_df = read.csv("tcga_sur_phen_10.csv")
Gender = anno_df[, "Gender"]
yearstobirth = as.numeric(anno_df[, "Age"])
pathologicstage = anno_df[, "Pathologic_stage"]
Lymphatic_INV = anno_df[, "Lymphatic_invasion"]

ha = HeatmapAnnotation(Gender = Gender, stage = pathologicstage, Lymphatic_INV =  Lymphatic_INV, 
                       Age = anno_points(yearstobirth, ylim = c(0, max(yearstobirth, na.rm = TRUE)), axis = TRUE),
                       col = list(Gender = c("MALE" = "gray", "FEMALE" = "red"),  
                                  stage = c("Stage I" = 'darkolivegreen1', "Stage IA" = "darkolivegreen2", 
                                            "Stage II" = "darkorange", "Stage IIB" = "darkorange1", "Stage IIC" = "darkorange2",'Stage IIA' = 'darkorange3',
                                            'Stage III' = 'darkorchid', "Stage IIIA" = "darkorchid1", "Stage IIIB" = "darkorchid2","Stage IIIC" = "darkorchid3",
                                            "Stage IV" = "blue","Stage IVA" = "blue3", "Stage IVB" = "blue4", 
                                            '[Discrepancy]' = "cornsilk2"),
                                  Lymphatic_INV = c("YES" = "red","NO" = "gray")),
                       annotation_height = unit(c(5, 5, 5,15), "mm"),
                       annotation_legend_param = list(Gender = list(title = "Gender"),stage = list(title = "Stage")))



anno_df1 = read.csv("tcga_sur_phen_90.csv")
Gender1 = anno_df1[, "Gender"]
yearstobirth1 = as.numeric(anno_df1[, "Age"])
pathologicstage1 = anno_df1[, "Pathologic_stage"]
Lymphatic_INV1 = anno_df1[, "Lymphatic_invasion"]
ha1 = HeatmapAnnotation(Gender = Gender1, stage = pathologicstage1,Lymphatic_INV =  Lymphatic_INV1, 
                        Age = anno_points(yearstobirth1, ylim = c(0, max(yearstobirth1, na.rm = TRUE)), axis = TRUE),
                        col = list(Gender = c("MALE" = "gray", "FEMALE" = "red"), 
                                    stage = c("Stage I" = 'darkolivegreen1', "Stage IA" = "darkolivegreen2", 
                                             "Stage II" = "darkorange", "Stage IIB" = "darkorange1", "Stage IIC" = "darkorange2",'Stage IIA' = 'darkorange3',
                                             'Stage III' = 'darkorchid', "Stage IIIA" = "darkorchid1", "Stage IIIB" = "darkorchid2","Stage IIIC" = "darkorchid3",
                                             "Stage IV" = "blue","Stage IVA" = "blue3", "Stage IVB" = "blue4", 
                                             '[Discrepancy]' = "cornsilk2"),
                                   Lymphatic_INV = c("YES" = "red","NO" = "gray")),
                        annotation_height = unit(c(5, 5, 5,15), "mm"),
                        annotation_legend_param = list(Gender = list(title = "Gender"),stage = list(title = "Stage")))




col = c("DEL" = "blue", "INDEL" = "red", 'FRA' = 'orange',"MUT" = "dark green", "SPL" = "dark red" ,"Intron" = "light green",  "NonSense" = "yellow", "OTHER" = "purple")
alter_fun = list(background = function(x, y, w, h)
  {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["DEL"], col = NA))
  },
  # big red
  INDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["INDEL"], col = NA))
  },
  FRA = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["FRA"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  },
  SPL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["SPL"], col = NA))
  },
  NonSense = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["NonSense"], col = NA))
  },
  Intron = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["Intron"], col = NA))
  },
  OTHER = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["OTHER"], col = NA))
  }
  
)

alter_fun = list(
  background = alter_graphic("rect", fill = "#CCCCCC"),   
  DEL = alter_graphic("rect", fill = col["DEL"]),
  INDEL = alter_graphic("rect", fill = col["INDEL"]),
  FRA = alter_graphic("rect", fill = col["FRA"]),
  MUT = alter_graphic("rect", fill = col["MUT"]),
  SPL = alter_graphic("rect", height = 0.33, fill = col["SPL"]),
  NonSense = alter_graphic("rect", height = 0.33, fill = col["NonSense"]),
  Intron = alter_graphic("rect", height = 0.33, fill = col["Intron"]),
  OTHER = alter_graphic("rect", height = 0.33, fill = col["OTHER"]))

# row_order = 1:nrow(sox10onco1)

column_title = "Low SOX9"
heatmap_legend_param = list(title = "Alternations", at = c("DEL", "INDEL", "FRA", "MUT", "SPL", 'NonSense', "Intron", "OTHER"), 
                            labels = c("Large Deletion", "INDEL", "Frame Shift","Missense_Mutation",  "Splice Site","Non_Sense", "Intron", "Other" ))
ht= oncoPrint(sox10onco1, alter_fun = alter_fun, col = col, alter_fun_is_vectorized = FALSE,
              column_title = column_title, heatmap_legend_param = heatmap_legend_param, bottom_annotation = ha)
       
column_title = "High SOX9"
heatmap_legend_param = list(title = "Alternations", at = c("DEL", "INDEL", "FRA", "MUT", "SPL", 'NonSense', "Intron", "OTHER"), 
                            labels = c("Large Deletion", "INDEL", "Frame Shift","Missense_Mutation",  "Splice Site","Non_Sense", "Intron", "Other" ))
ht1 = oncoPrint(sox90onco1,alter_fun = alter_fun, col = col, alter_fun_is_vectorized = FALSE,  
                column_title = column_title, heatmap_legend_param = heatmap_legend_param, bottom_annotation = ha1)

ht + ht1



