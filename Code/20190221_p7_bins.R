#script for SCENIC.. 

library(Seurat)
library(SCENIC)
library(GENIE3)
library(RcisTarget)
library(SingleCellExperiment)
library(scales)
library(dplyr)
library(AUCell)


#need to add merged to the biowulf account and create a new scenic file.. 

#opening the merged and merged metadata files
P7 <- readRDS("/data/kollal2/scenic/p7_ls_epi_final.Rds")
P7@meta.data$clus.ident <- P7@ident

timepoint <- as.data.frame(as.matrix(P7@meta.data$timepoint))
clus.ident <- as.data.frame(as.matrix(P7@meta.data$clus.ident))
stim <- as.data.frame(as.matrix(P7@meta.data$stim))
metadata <- cbind(timepoint, clus.ident, stim)
colnames(metadata) <- c("timepoint", "ident", "batch")
rownames(metadata) <- colnames(P7@data) 


#get the raw counts data and metadata from Seurat

raw.data <- as.matrix(P7@raw.data)
data <- raw.data[rownames(raw.data) %in% rownames(as.matrix(P7@data)),
                 colnames(raw.data) %in% colnames(as.matrix(P7@data))]
data <- data[rownames(P7@data), colnames(P7@data)]
mydata <- as.data.frame(data)


#make the rownames of the metadata  == colnames of the data
table(factor(rownames(metadata) == colnames(mydata)))

##POOLING/BINNING THE AVERAGE COUNTS DATA FOR EVERY 20 CELLS OF THE SAME CELL-TYPE TOGETHER.. need a nested for loop where the inner loop creates a new dataframe for each cell-type in the data.. 
pooled_data <- c()
pooled_metadata <- c()
celltypes <- levels(factor(metadata$ident))

for(cell in celltypes){
  new_data <- mydata[, metadata$ident == cell]
  if(length(new_data) < 20){
    #SAMPLE FOR EVERY 5
    new_data <- new_data[,sample(ncol(new_data), floor(length(new_data)/5)*5)] #to check if the colnames get scrambled everytime, make vectors with colnames for each random sampling and see if the colnames equal each other
    #find the average counts of every 10 cells.. 
    new_data <- t(apply(new_data, 1, tapply, gl(ncol(new_data)/5, 5), mean)) #https://stackoverflow.com/questions/40365902/how-to-average-every-12-columns-for-each-row-in-r?rq=1
    
  }
  else if(length(new_data) < 100){
    #SAMPLE FOR EVERY 10
    new_data <- new_data[,sample(ncol(new_data), floor(length(new_data)/10)*10)] #to check if the colnames get scrambled everytime, make vectors with colnames for each random sampling and see if the colnames equal each other
    #find the average counts of every 10 cells..
    new_data <- t(apply(new_data, 1, tapply, gl(ncol(new_data)/10, 10), mean))
  }
  else{
    #SAMPLE FOR EVERY 20
    new_data <- new_data[,sample(ncol(new_data), floor(length(new_data)/20)*20)]#get a multiple of 20 for each dataset...
    new_data <- t(apply(new_data, 1, tapply, gl(ncol(new_data)/20, 20), mean)) #https://stackoverflow.com/questions/40365902/how-to-average-every-12-columns-for-each-row-in-r?rq=1
  }
  pooled_metadata <- c(pooled_metadata, rep(cell, ncol(new_data))) 
  pooled_data <-cbind(pooled_data, new_data)
  new_data <- c()
}

#colnames(pooled_data) <- pooled_metadata
colnames(pooled_data) <- paste0(pooled_metadata, "-", 1:ncol(pooled_data))
pooled_metadata <- cbind(pooled_metadata, rep("P7", length(pooled_metadata)))
colnames(pooled_metadata) <- c("ident","timepoint")
rownames(pooled_metadata) <- colnames(pooled_data)




###RUNNING SCENIC:
##create directory 
dir.create("/data/kollal2/scenic/20190221_p7_bins")
setwd("/data/kollal2/scenic/20190221_p7_bins")
dir.create("int")

#set up the scenic object.. preprocessing
exprMat <- as.matrix(pooled_data)
cellInfo <- as.data.frame(pooled_metadata)

cellInfo$CellType <- cellInfo$ident
cellInfo$nGene <- colSums(exprMat>0)
cellInfo <- data.frame(cellInfo)
cellInfo <- cellInfo[,-2:-4]

saveRDS(cellInfo, file = "int/cellInfo.Rds")

p1_epi_levels <- levels(factor(cellInfo))
p1_epi_colors <-  hue_pal()(length(p1_epi_levels))

colVars <- list(CellType=setNames(p1_epi_colors, p1_epi_levels))
saveRDS(colVars, file="int/colVars.Rds")

plot.new(); legend(0,1, fill=colVars$CellType, legend=names(colVars$CellType))

##setting up the scenic object continued
org="mgi" # or hgnc, or dmel
dbDir="/data/kollal2/scenic/databases" # RcisTarget databases location
myDatasetTitle="Pooled P7 20190221 scenic" # choose a name for your analysis
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, datasetTitle=myDatasetTitle, nCores=20) 

scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"

saveRDS(scenicOptions, file="int/scenicOptions.Rds") 

nCellsPerGene <- apply(exprMat, 1, function(x) sum(x>0))
nCountsPerGene <- apply(exprMat, 1, sum)

summary(nCellsPerGene)
summary(nCountsPerGene)

max(exprMat)
sum(exprMat>0) / sum(exprMat==0)

minReads <- 3*.01*ncol(exprMat)
genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minReads)]
length(genesLeft_minReads)

minSamples <- ncol(exprMat)*.01
nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
length(genesLeft_minCells)


motifRankings <- importRankings(getDatabases(scenicOptions)[[1]]) 
genesInDatabase <- colnames(getRanking(motifRankings))

genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
length(genesLeft_minCells_inDatabases)

#interestingGenes <- c("Neurod6", "Pou4f3", "Ccer2", "Atoh1","Pou4f3", "Sall1", "Sall2", "Sall3" ,"Sall4", "Sox2", "Lrrn1", "Lrrn3",  "Fam19a3", "Grp")
#interestingGenes[which(!interestingGenes %in% genesLeft_minCells_inDatabases)]

genesKept <- genesLeft_minCells_inDatabases
saveRDS(genesKept, file=getIntName(scenicOptions, "genesKept"))

exprMat_filtered <- exprMat[genesKept, ]

rm(exprMat)

#correlation 
corrMat <- cor(t(exprMat_filtered), method="spearman")
saveRDS(corrMat, file=getIntName(scenicOptions, "corrMat"))

##run GENIE3

exprMat_filtered <- log2(exprMat_filtered+1) 
runGenie3(exprMat_filtered, scenicOptions)

#scenicOptions@settings$modules


##
#scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 20
scenicOptions@settings$seed <- 123

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, exprMat_filtered)
# 



##FOLLOWUP ANALYSIS.. SHANNON ENTROPY FOR EACH CELL TYPE.. 
regulon_auc <- readRDS("int/3.4_regulonAUC.Rds")
scenicOptions <- readRDS("int/scenicOptions.Rds")
cellInfo <- readRDS("int/cellInfo.Rds")
levels <- levels(factor(cellInfo))

auc_scores <- getAUC(regulon_auc)

#making normalized auc matrix 
normalized_regulons <- auc_scores/rowSums(auc_scores)

#making an "identity" matrix and then normalizing it
identity <- c()
for(i in 1:length(levels(cellInfo))){
  identity <- rbind(identity, as.character(cellInfo))
}

for(i in 1:nrow(identity)){
  for(n in 1:ncol(identity)){
    if(identity[i,n] == levels[i]){
      identity[i,n] = 1
    }
    else{
      identity[i,n] = 0
    }
  }
}

identity <- as.data.frame(identity)
identity <- data.frame(sapply(identity, function(x) as.numeric(as.character(x))))
normalized_identity <- identity/rowSums(identity)
rownames(normalized_identity) <- levels
colnames(normalized_identity) <- colnames(auc_scores)

#calculating the JSD
h_func <- function(p){
  p <- p*log(p)
  p[is.na(p)] <- 0
  p <- -sum(p)
  return(p)
}

returned_matrix <- c()
for(i in 1:length(levels)){
  rss_matrix <- c()
  for(n in 1:nrow(normalized_regulons)){
    h_rc <- h_func((normalized_identity[i,] + normalized_regulons[n,])/2)
    h_r <- h_func(normalized_regulons[n,])
    h_c <- h_func(normalized_identity[i,])
    rss <- 1-sqrt(h_rc - ((h_c + h_r)/2))
    rss_matrix <- c(rss_matrix, rss)
  }
  returned_matrix <- cbind(returned_matrix, rss_matrix)
}

colnames(returned_matrix) <- levels
rownames(returned_matrix) <- rownames(normalized_regulons)

saveRDS(returned_matrix,"int/cell_type_specificity_score.Rds")

for(i in 1:nrow(normalized_identity)){print(paste0("min for ", levels[i], " is ", min(returned_matrix[,i]), "and the max is ", max(returned_matrix[,i])))}






