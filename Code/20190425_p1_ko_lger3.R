library(DEsingle)
#library(BiocParallel)
library(Seurat)

##P1 LATERAL VS. MEDIAL GER
ko <- readRDS("/data/kollal2/scenic/20190418_P1_KO.Rds")
#ko <- SetIdent(ko, cells.use = WhichCells(ko, ident = c()), ident.use = "Prosensory")
idents <- levels(factor(ko@ident))
ko <- SetIdent(ko, cells.use = WhichCells(ko, ident = idents[-5]), ident.use = "AOC")


#for DESingle, you need a counts matrix and you need another vector with the metadata.. aka just the gorup names.. 
#and there can only be two levels.

object <- SubsetData(ko, ident.use = c( "Lateral_GER_3", "AOC"))
raw.data <- as.matrix(object@raw.data)
data <- raw.data[,colnames(raw.data) %in% colnames(as.matrix(object@data))]
data <- data[, colnames(object@data)]
idents <- as.data.frame(factor(object@ident))

colnames(idents) <- "ident"

group <- factor(idents$ident)
counts <- data

# 
# param <- MulticoreParam(workers = 18, progressbar = TRUE)
# register(param)

results <- DEsingle(counts = counts, group = group)
saveRDS(results, "/data/kollal2/desingle/20190425_ko_lger3_aoc_results.Rds")

results.classified <- DEtype(results = results, threshold = 0.05)
saveRDS(results.classified, "/data/kollal2/desingle/20190425_ko_lger3_aoc_results_classified.Rds")

results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]
saveRDS(results.sig, "/data/kollal2/desingle/20190425_ko_lger3_aoc_results_sig.Rds")



# results.DEs <- results.sig[results.sig$Type == "DEs", ]
# results.DEa <- results.sig[results.sig$Type == "DEa", ]
# results.DEg <- results.sig[results.sig$Type == "DEg", ]


