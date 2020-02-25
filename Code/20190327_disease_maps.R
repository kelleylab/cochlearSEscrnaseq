library(Seurat)
library(dplyr)
library(matrixStats)
library(reshape)
library(qlcMatrix)
library(ggplot2)
library(pheatmap)

genes <- read.delim('/data/kollal2/hearing_loss_genes.txt')
gwas <- read.delim("/data/kollal2/gwas_hearing_loss_genes.txt")

p1 <- read.delim("/data/kollal2/scenic/p1_epi_final.Rds")
p1 <- SetIdent(p1, cells.use = WhichCells(p1, ident = c("earlyIHC", "earlyOHC", "OHC", "IHC")), ident.use = "HC")
p1 <- SetIdent(p1, cells.use = WhichCells(p1, ident = c("Deiters_1st/2nd_row", "Deiters_3rd_row")), ident.use = "Deiters")
p1 <- SetIdent(p1, cells.use = WhichCells(p1, ident = c("Lateral_GER_1", "Lateral_GER_2", "Lateral_GER_3")), ident.use = "Lateral_GER")

p7 <- read.delim("/data/kollal2/scenic/p7_ls_epi_final.Rds")
p7 <- SetIdent(p7, cells.use = WhichCells(p7, ident = c("IHC", "OHC")), ident.use = "HC")

e14 <- read.delim("/data/kollal2/scenic/e14_epi_final.Rds")

e16 <- read.delim("/data/kollal2/scenic/e16_epi_final.Rds")
e16 <- SetIdent(e16, cells.use = WhichCells(e16, ident = c("earlyOHC", "OHC", "IHC")), ident.use = "HC")
e16 <- SetIdent(e16, cells.use = WhichCells(e16, ident = c("GER")), ident.use = "Medial_GER")

merged <- read.delim("/data/kollal2/scenic/merged_20190305.Rds")
merged <- SetAllIdent(merged, "timepoint")

#making the large genelist
#dominant non syndromic hearing loss
dns <- genes$Gene[genes$Disease == "Autosomal Dominant Non-syndromic Hearing Loss"]
dns <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% dns]

#recessive nonsyndromic hearing loss
rns <- genes$Gene[genes$Disease == "Autosomal Recessive Non-syndromic Hearing Loss"]
rns <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% rns]


#x-linked nonsyndromic hearing loss
xlink <- genes$Gene[genes$Disease == "X-Linked Nonsyndromic Hearing Loss Genes"]
xlink <-rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% xlink]


##Syndromic hearing loss
shl <- genes$Gene[!(genes$Disease %in% c("X-Linked Nonsyndromic Hearing Loss Genes", "Autosomal Recessive Non-syndromic Hearing Loss",
                                         "Autosomal Dominant Non-syndromic Hearing Loss"))]
shl <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% shl]

#noise induced hearing loss
ni <- gwas$Gene[gwas$Disease == "Noise_induced"]
ni <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% ni]

#age related hearing loss
ar <- gwas$Gene[gwas$Disease == "Age_related"]
ar <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% ar]

#cisplatin_response
cr <- gwas$Gene[gwas$Disease == "Cisplatin_reponse"]
cr <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% cr]

#ototoxicity
ot <- gwas$Gene[gwas$Disease == "Ototoxicity"]
ot <- rownames(merged@raw.data)[toupper(rownames(merged@raw.data)) %in% ot]

#making a merged list
new <- c(dns, rns, xlink, shl, ar, ni, cr)
new_labels <- c(rep("DFNA", length(dns)), rep("DFNB", length(rns)), rep("DFNX", length(xlink)), rep("Syndromic", length(shl)), 
                rep("ARHL", length(ar)), rep("NIHL", length(ni)), rep("CRHL", length(cr)))
genes <- as.data.frame(cbind(new, new_labels))
genes_merged <- aggregate(new_labels ~ new, data = genes, paste, collapse = "/")
rownames(genes_merged) <- genes_merged$new
all_genes <- as.character(genes_merged$new)

#functions: 
colMax <- function(data) sapply(data, max, na.rm = TRUE)
disease_maps_2 <- function(object, genes, defined_threshold) {
  genes <- genes[genes %in% rownames(object@raw.data)]
  temp_data <- rowSums(as.data.frame(as.matrix(object@raw.data[genes,])))
  genes <- names(temp_data[temp_data>10])
  object_data <- as.data.frame(as.matrix(object@data[genes,]))
  object_mean <- rowMeans(object_data)
  object_sd <- rowSds(as.matrix(object_data))
  
  new_data <- c()
  for(i in 1:nrow(object_data)){
    object_data[i, ] <- (object_data[i, ] - object_mean[i]) / object_sd[i]
    temp <- c()
    temp$ID <- as.matrix(object@ident)
    temp$Value <- object_data[i,]
    temp <- aggregate(unlist(temp$Value) ~ temp$ID, temp, mean)
    temp$`unlist(temp$Value)` <- (temp$`unlist(temp$Value)` - mean(temp$`unlist(temp$Value)`))/ sd(temp$`unlist(temp$Value)`)
    new_data <- cbind(new_data, temp$`unlist(temp$Value)`)
    temp <- c()
  }
  
  colnames(new_data) <- rownames(object_data) 
  rownames(new_data) <- levels(factor(object@ident))
  max <- colMax(as.data.frame(new_data))
  new_data <- new_data[, max>=defined_threshold]
  return(new_data)
}


#running the functions 
p1_all_genes <- disease_maps_2(p1, all_genes,2)
saveRDS(p1_all_genes, "/data/kollal2/Routput/20190327_p1_disease.Rds")

p7_all_genes <- disease_maps_2(p7, all_genes,2)
saveRDS(p7_all_genes, "/data/kollal2/Routput/20190327_p7_disease.Rds")

e14_all_genes <- disease_maps_2(e14, all_genes,2)
saveRDS(e14_all_genes, "/data/kollal2/Routput/20190327_e14_disease.Rds")

e16_all_genes <- disease_maps_2(e16, all_genes,2)
saveRDS(e16_all_genes, "/data/kollal2/Routput/20190327_e16_disease.Rds")

merged_all_genes <- disease_maps_2(merged, all_genes,1)
saveRDS(merged_all_genes, "/data/kollal2/Routput/20190327_merged_disease.Rds")


#for non-biowulf purposes
p1_all_genes <- disease_maps_2(p1, all_genes,2)
p1_ls_all_genes <- disease_maps_2(p1, all_genes,2)
saveRDS(p1_all_genes, "~/Desktop/RDS/20190327_p1_disease.Rds")
saveRDS(p1_ls_all_genes, "~/Desktop/RDS/20190327_p1_ls_disease.Rds")

p7_all_genes <- disease_maps_2(p7, all_genes,2)
saveRDS(p7_all_genes, "~/Desktop/RDS/20190327_p7_disease.Rds")

e14_all_genes <- disease_maps_2(e14, all_genes,2)
saveRDS(e14_all_genes, "~/Desktop/RDS/20190327_e14_disease.Rds")

e16_all_genes <- disease_maps_2(e16, all_genes,2)
saveRDS(e16_all_genes, "~/Desktop/RDS/20190327_e16_disease.Rds")
e16_ls_all_genes <- disease_maps_2(e16, all_genes,2)
saveRDS(e16_ls_all_genes, "~/Desktop/RDS/20190327_ls_e16_disease.Rds")

merged_all_genes <- disease_maps_2(merged, all_genes,1)
saveRDS(merged_all_genes, "~/Desktop/RDS/20190327_merged_disease.Rds")

saveRDS(genes_merged, "~/Desktop/RDS/20190327_genes_merged.Rds")


#making the heatmaps
p1_disease <- readRDS("~/Desktop/RDS/20190327_p1_disease.Rds")
p1_ls_disease <- readRDS( "~/Desktop/RDS/20190327_p1_ls_disease.Rds")
e14_disease <- readRDS("~/Desktop/RDS/20190327_e14_disease.Rds")
e16_disease <- readRDS("~/Desktop/RDS/20190327_e16_disease.Rds")
e16_ls_disease <- readRDS("~/Desktop/RDS/20190327_ls_e16_disease.Rds")
p7_disease <- readRDS("~/Desktop/RDS/20190327_p7_disease.Rds")
merged_disease <- readRDS("~/Desktop/RDS/20190327_merged_disease.Rds")
genes_merged <- readRDS("~/Desktop/RDS/20190327_genes_merged.Rds")

library(pheatmap)
creat_ann_by_disease <- function(disease_data, genes_merged){
  disease_data <- disease_data[,duplicated(colnames(disease_data))==FALSE]
  ann_col <- genes_merged[genes_merged$new %in% colnames(disease_data), ]
  ann_col <- ann_col %>% group_by(new_labels) %>% arrange(desc(new_labels), .by_group = TRUE)
  genes <- ann_col$new
  ann_col <- as.data.frame(ann_col$new_labels)
  colnames(ann_col) <- c("Disease")
  rownames(ann_col) <- genes
  return(ann_col)
}



make_disease_map <- function(disease_data,genes_merged, path){
  ann <- creat_ann_by_disease(disease_data, genes_merged)
  disease_data <- disease_data[, rownames(ann)]
  pdf(path, width = 20, height = 5)
  pheatmap(disease_data, cluster_rows = TRUE, cluster_cols = FALSE,annotation_col = ann, treeheight_row = 0, treeheight_col = 0)
  dev.off()
}

make_disease_map(p1_disease, genes_merged, "~/Desktop/disease_p1.pdf")
make_disease_map(p1_ls_disease, genes_merged, "~/Desktop/disease_p1_ls.pdf")
make_disease_map(p7_disease, genes_merged, "~/Desktop/disease_p7.pdf")
make_disease_map(e14_disease, genes_merged, "~/Desktop/disease_e14.pdf")
make_disease_map(e16_disease, genes_merged, "~/Desktop/disease_e16.pdf")
make_disease_map(e16_ls_disease, genes_merged, "~/Desktop/disease_e16_ls.pdf")
make_disease_map(merged_disease, genes_merged, "~/Desktop/disease_merged.pdf")

make_disease_map_by_group <- function(disease_data,genes_merged, path){
  ann <- creat_ann_by_disease(disease_data, genes_merged)
  ann$Gene <- rownames(ann)
  ann <- ann[colnames(disease_data),]
  genes <- ann$Gene
  ann <- as.data.frame(ann$Disease)
  colnames(ann) <- "Disease"
  rownames(ann) <- genes
  pdf(path, width = 20, height = 5)
  pheatmap(disease_data, cluster_rows = FALSE, cluster_cols = TRUE, treeheight_row = 0, treeheight_col = 0, annotation_col = ann)
  dev.off()
}
ann_merged <- creat_ann_by_disease(merged_disease, genes_merged)
ann_merged <- as.data.frame(ann_merged[colnames(merged_disease),])
ann_merged <- make_disease_map_by_group(merged_disease, genes_merged, "~/Desktop/disease_merged.pdf")

e14_disease <- e14_disease[, rownames(ann_col_e14)]
table(factor(colnames(e14_disease) == rownames(ann_col_e14)))
pdf("~/Desktop/lala.pdf", width = 30, height = 5)
pheatmap(e14_disease, cluster_rows = TRUE, cluster_cols = FALSE,annotation_col = ann_col_e14)
dev.off()

pdf("~/Desktop/lala.pdf", width = 30, height = 5)
pheatmap(p1_disease, treeheight_row = 0, treeheight_col = 0, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = ann_col)
dev.off()

p1_all_genes <- p1_all_genes[,duplicated(colnames(p1_all_genes))==FALSE]
ann_col <- genes_merged[genes_merged$new %in% colnames(p1_all_genes), ]
ann_col <- genes_merged[colnames(p1_all_genes), ]
ann_col_2 <- as.data.frame(ann_col$new_labels)
colnames(ann_col_2) <- "Disease"
rownames(ann_col_2) <- rownames(ann_col)

temp_ann_col <- ann_col
temp_ann_col <- temp_ann_col %>% group_by(new_labels) %>% arrange(desc(new_labels), .by_group = TRUE)
rownames(temp_ann_col) <- temp_ann_col$new

temp_p1_all_genes <- p1_all_genes
temp_p1_all_genes <- temp_p1_all_genes[, rownames(temp_ann_col)]
table(factor(colnames(temp_p1_all_genes) == rownames(temp_ann_col)))

ann_col_3 <- as.data.frame(temp_ann_col$new_labels)
colnames(ann_col_3) <- "Disease"
rownames(ann_col_3) <- rownames(temp_ann_col)

png("~/Desktop/lala.png")
pheatmap(temp_p1_all_genes, treeheight_row = 0, treeheight_col = 0, cluster_rows = TRUE, cluster_cols = FALSE, annotation_col = ann_col_3)
dev.off()
