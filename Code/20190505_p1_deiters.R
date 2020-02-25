library(DEsingle)
#library(BiocParallel)
library(Seurat)

p1 <- readRDS("/data/kollal2/scenic/p1_epi_final.Rds")

#for DESingle, you need a counts matrix and you need another vector with the metadata.. aka just the gorup names.. 
#and there can only be two levels.


object <- SubsetData(p1, ident.use = c("Deiters_1st/2nd_row", "Deiters_3rd_row"))
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
saveRDS(results, "/data/kollal2/desingle/20190505_p1_deiters_results.Rds")

results.classified <- DEtype(results = results, threshold = 0.05)
saveRDS(results.classified, "/data/kollal2/desingle/20190505_p1_deiters_results_classified.Rds")

results.sig <- results.classified[results.classified$pvalue.adj.FDR < 0.05, ]
saveRDS(results.sig, "/data/kollal2/desingle/20190505_p1_deiters_results_sig.Rds")



# results.DEs <- results.sig[results.sig$Type == "DEs", ]
# results.DEa <- results.sig[results.sig$Type == "DEa", ]
# results.DEg <- results.sig[results.sig$Type == "DEg", ]


