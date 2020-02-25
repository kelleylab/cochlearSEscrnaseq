library(URD)
#library(rgl)

merged <- readRDS("/data/kollal2/scenic/merged_20190305.Rds")
merged@meta.data$stage.nice <- c(rep(1, 5253), rep(2, 7961), rep(3, 14043))
merged@meta.data$time_batch <- paste0(merged@meta.data$timepoint, "-", merged@meta.data$stim)
merged@meta.data$clus.ident <- merged@meta.data$ident
merged <- SetAllIdent(merged, "time_batch")
merged <- SubsetData(merged, ident.remove = c("P1-se_1", "P1-se_2", "P1-se_3"))
merged <- SetAllIdent(merged, "clus.ident")

urd_creation <- function(object, identities){
  object <- SubsetData(object, ident.use = identities)
  raw.data <- as.matrix(object@raw.data)
  data <- raw.data[,colnames(raw.data) %in% colnames(as.matrix(object@data))]
  data <- data[, colnames(object@data)]
  timepoint <- as.data.frame(as.matrix(object@meta.data$timepoint))
  clus.ident <- as.data.frame(as.matrix(object@meta.data$ident))
  stim <- as.data.frame(as.matrix(object@meta.data$stim))
  stage.nice <- as.data.frame(as.matrix(object@meta.data$stage.nice))
  metadata <- cbind(timepoint, clus.ident, stim, stage.nice)
  rownames(metadata) <- colnames(object@data)
  colnames(metadata) <- c("timepoint", "clus.ident", "stim", "stage.nice")
  urd_object <- createURD(count.data = data, meta = metadata, min.cells = 3, min.counts = 3)
  urd_object@meta$timepoint <- metadata$timepoint
  urd_object@meta$cell.ident <- metadata$clus.ident
  urd_object@meta$stim <- metadata$stim
  urd_object@meta$stage.nice <- metadata$stage.nice
  urd_object@group.ids$stage <- as.character(urd_object@meta[rownames(urd_object@group.ids), "stage.nice"])
  return(urd_object)
}


urd_analysis <- function(object){
  stages <- sort(unique(object@group.ids$stage))
  var.by.stage <- lapply(1:length(stages), function(n) {
    findVariableGenes(object, cells.fit=cellsInCluster(object, "stage", stages[n]), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=paste0("Stage:", stages[n]), do.plot=F)
  })
  var.genes <- sort(unique(unlist(var.by.stage)))
  object@var.genes <- var.genes
  object <- calcPCA(object, mp.factor = 2)
  set.seed(19)
  object <- calcTsne(object = object)
  object <- calcDM(object, knn = sqrt(length(object@group.ids$init)), sigma=NULL)
  root.cells <- cellsInCluster(object, "stage", "1")
  floods <- floodPseudotime(object, root.cells = root.cells, n=50, minimum.cells.flooded = 2, verbose=F)
  object <- floodPseudotimeProcess(object, floods, floods.name="pseudotime")
  object_tips <- urdSubset(object, cells.keep=cellsInCluster(object, "stage", "3"))
  tips <- object_tips@meta$clus.ident
  object@group.ids[rownames(object_tips@group.ids), "tip.clusters"] <- tips
  object.ptlogistic <- pseudotimeDetermineLogistic(object, "pseudotime", optimal.cells.forward=20, max.cells.back=40, do.plot = T)
  object.biased.tm <- as.matrix(pseudotimeWeightTransitionMatrix(object, "pseudotime", logistic.params=object.ptlogistic))
  object.walks <- simulateRandomWalksFromTips(object, tip.group.id="tip.clusters", root.cells=root.cells, transition.matrix = object.biased.tm, n.per.tip = 25000, root.visits = 1, max.steps = 5000, verbose = F)
  object <- processRandomWalksFromTips(object, object.walks, verbose = F)
  return(object)
}

hc_pro_urd <- urd_creation(merged, c("IHC", "earlyIHC", "earlyOHC", "OHC", "Lateral Prosensory", "Medial_Prosensory"))
saveRDS(hc_pro_urd,"~/Desktop/RDS/hc_pro_urd_creation.Rds")

hc_pro_urd <- urd_analysis(hc_pro_urd)
saveRDS(hc_pro_urd, "/data/kollal2/2019_hc_prosensory_urd.Rds")

everything_urd <- urd_creation(merged, c("Medial_Prosensory", "IHC", "IPhC", "Lateral Prosensory", 
                                     "OHC", "Deiters_1st/2nd_row", "Deiters_3rd_row", "OPC", "IPC", "earlyOHC", "earlyIHC"))
saveRDS(everything_urd,"~/Desktop/RDS/everything_urd_creation.Rds")
everything_urd <- urd_analysis(everything_urd)
saveRDS(everything_urd, "/data/kollal2/2019_se_4_everything_urd.Rds")

lateral_pc_urd <- urd_creation(merged, c("Lateral Prosensory", "OHC",
                                         "Deiters_1st/2nd_row", "Deiters_3rd_row", "OPC", "IPC", "earlyOHC"))
saveRDS(lateral_pc_urd,"~/Desktop/RDS/lateral_pc_urd_creation.Rds")
lateral_pc_urd  <- urd_analysis(lateral_pc_urd )
saveRDS(lateral_urd, "/data/kollal2/2019_se_4_lateral_pc_urd.Rds")

medial_pc_urd <- urd_creation(merged, c("IHC", "earlyIHC", "Medial_Prosensory", "IPhC", "IPC"))


