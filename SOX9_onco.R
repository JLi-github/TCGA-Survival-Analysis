# 08252022
# JL

 
library(dplyr)
library(ComplexHeatmap)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ComplexHeatmap")


# For mat, we can define the function as:
mat = read.table(textConnection(
"s1,s2,s3
g1,snv;indel,snv,indel
g2,,snv;indel,snv
g3,snv,,indel;snv"), row.names = 1, header = TRUE, sep = ",", stringsAsFactors = FALSE)
mat = as.matrix(mat)
mat

# Colors for different alterations are defined in col. It should be a named vector for which names correspond
# to alteration types. It is used to generate the barplots.

get_type_fun = function(x) strsplit(x, ";")[[1]]
get_type_fun(mat[1, 1])

get_type_fun(mat[1, 2])

col = c(snv = "red", indel = "blue")
oncoPrint(mat,
          alter_fun = list(
          snv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["snv"], col = NA)),
          indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["indel"], col = NA))
          ), col = col)


# You can see the order in barplots also correspond to the order defined in alter_fun. The grahpics in legend are based 
# on the functions defined in alter_fun.
oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = list(
          snv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = "red", col = NA)),
          indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = "blue", col = NA))
          ), col = c(snv = "red", indel = "blue"))


# input data is a list of matrix for which each matrix contains binary value representing whether the alteration is 
# absent or present.
mat_list = list(snv = matrix(c(1, 0, 1, 1, 1, 0, 0, 1, 1), nrow = 3),indel = matrix(c(1, 0, 0, 0, 1, 0, 1, 0, 0), nrow = 3))
                
rownames(mat_list$snv) = rownames(mat_list$indel) = c("g1", "g2", "g3")
colnames(mat_list$snv) = colnames(mat_list$indel) = c("s1", "s2", "s3")
mat_list

oncoPrint(mat_list,
          alter_fun = list(
          snv = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = col["snv"], col = NA)),
          indel = function(x, y, w, h) grid.rect(x, y, w*0.9, h*0.4, gp = gpar(fill = col["indel"], col = NA))
          ), col = col)



mat = read.table(paste0(system.file("extdata", package = "ComplexHeatmap"), 
                 "/tcga_lung_adenocarcinoma_provisional_ras_raf_mek_jnk_signalling.txt"), 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t")
mat[is.na(mat)] = ""
rownames(mat) = mat[, 1]
mat = mat[, -1]
mat=  mat[, -ncol(mat)]
mat = t(as.matrix(mat))
mat[1:3, 1:3]



col = c("HOMDEL" = "blue", "AMP" = "red", "MUT" = "#008000")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  HOMDEL = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["HOMDEL"], col = NA))
  },
  # big red
  AMP = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["AMP"], col = NA))
  },
  # small green
  MUT = function(x, y, w, h) {
    grid.rect(x, y, w-unit(2, "pt"), h*0.33, 
              gp = gpar(fill = col["MUT"], col = NA))
  }
)


col = c("MUT" = "#008000", "AMP" = "red", "HOMDEL" = "blue")

# Make the oncoPrint and adjust heatmap components such as the title and the legend

oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("AMP", "HOMDEL", "MUT"), 
                                      labels = c("Amplification", "Deep deletion", "Mutation")))


# Apply to cBioPortal dataset

setwd("C:/Users/lijinju/downloads/")
mat = read.table("COADREAD_mc3_gene_level.txt", 
                 header = TRUE,stringsAsFactors=FALSE, sep = "\t", fill = TRUE)
                
 
new_mat = mat[!duplicated(mat$xena_sample), ]
new_mat = new_mat[-7224,]
rownames(new_mat)= new_mat[,1]
new_mat = new_mat[,-1]


alter_fun = list(
  background = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  # big blue
  Mutation = function(x, y, w, h) {
  grid.rect(x, y, w-unit(2, "pt"), h-unit(2, "pt"), 
              gp = gpar(fill = col["1"], col = NA))
  })


oncoPrint(mat, get_type = function(x) strsplit(x, ";")[[1]],
          alter_fun = alter_fun, col = col, 
          column_title = "OncoPrint for TCGA Lung Adenocarcinoma, genes in Ras Raf MEK JNK signalling",
          heatmap_legend_param = list(title = "Alternations", at = c("Mutation"), 
                                      labels = c("Mutation")))

