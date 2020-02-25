library(URD)

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

hc_pro_urd <- readRDS("/data/kollal2/hc_pro_urd_creation.Rds")
hc_pro_urd <- urd_analysis(hc_pro_urd)
saveRDS(hc_pro_urd, "/data/kollal2/2019_hc_prosensory_urd.Rds")

everything_urd <- readRDS("/data/kollal2/everything_urd_creation.Rds")
everything_urd <- urd_analysis(everything_urd)
saveRDS(everything_urd, "/data/kollal2/2019_se_4_everything_urd.Rds")

lateral_pc_urd<- readRDS("/data/kollal2/lateral_pc_urd_creation.Rds")
lateral_pc_urd  <- urd_analysis(lateral_pc_urd )
saveRDS(lateral_pc_urd, "/data/kollal2/2019_se_4_lateral_pc_urd.Rds")


