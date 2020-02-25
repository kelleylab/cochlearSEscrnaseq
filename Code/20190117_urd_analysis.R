library(URD)
library(rgl)

merged <- readRDS("/data/kollal2/merged_20190117.Rds")
merged@meta.data$stage.nice <- c(rep(1, 5253), rep(2, 7961), rep(3, 14043))


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
  urd_object@meta$cell.ident <- metadata$cell.ident
  urd_object@meta$stim <- metadata$stim
  urd_object@meta$stage.nice <- metadata$stage.nice
  urd_object@group.ids$stage <- as.character(urd_object@meta[rownames(urd_object@group.ids), "stage.nice"])
  return(urd_object)
}


urd_analysis <- function(object){
  stages <- sort(unique(object@group.ids$stage))
  var.by.stage <- lapply(1:length(stages), function(n) {
    findVariableGenes(object, cells.fit=cellsInCluster(object, "stage", stages[n]), set.object.var.genes=F, diffCV.cutoff=0.3, mean.min=.005, mean.max=100, main.use=paste0("Stage:", stages[n]), do.plot=T)
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

hc_urd <- urd_creation(merged, c("IHC", "earlyIHC", "earlyOHC", "OHC"))
hc_urd <- urd_analysis(hc_urd)

saveRDS(hc_urd, "/data/kollal2/2019_hc_urd.Rds")

ihc_urd <- urd_creation(merged, c("IHC", "earlyIHC"))
ihc_urd <- urd_analysis(ihc_urd)
saveRDS(ihc_urd, "/data/kollal2/2019_ihc_urd.Rds")

just_ihc_urd <- urd_creation(merged, c("IHC"))
just_ihc_urd <- urd_analysis(just_ihc_urd)
saveRDS(just_ihc_urd, "/data/kollal2/2019_just_ihc_urd.Rds")

ohc_urd <- urd_creation(merged, c("OHC", "earlyOHC"))
ohc_urd <- urd_analysis(ohc_urd)
saveRDS(ohc_urd, "/data/kollal2/2019_ohc_urd.Rds")

just_ohc_urd <- urd_creation(merged, c("OHC"))
just_ohc_urd <- urd_analysis(just_ohc_urd)
saveRDS(just_ohc_urd, "/data/kollal2/2019_just_ohc_urd.Rds")

medial_urd <- urd_creation(merged, c("Medial_Prosensory", "IHC", "IPhC"))
medial_urd <- urd_analysis(medial_urd)
saveRDS(medial_urd, "/data/kollal2/2019_medial_urd.Rds")

lateral_urd <- urd_creation(merged, c("Lateral Prosensory", "OHC", "Deiters_1st/2nd_row", "Deiters_3rd_row"))
lateral_urd <- urd_analysis(lateral_urd)
saveRDS(lateral_urd, "/data/kollal2/2019_lateral_urd.Rds")









