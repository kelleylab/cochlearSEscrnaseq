library(RcisTarget)
library(RcisTarget.mm9.motifDatabases.20k)

p1_list <- readRDS("/data/kollal2/scenic/20190104_p1_se_4_bins_regulon_genelists.Rds")
e14_list <- readRDS("/data/kollal2/scenic/20190104_e14_bins_regulon_genelists.Rds")
e16_list <- readRDS("/data/kollal2/scenic/20190104_e16_bins_regulon_genelists.Rds")

data(motifAnnotations_mgi)
motifRankings <- importRankings("/data/kollal2/scenic/databases/mm9-500bp-upstream-7species.mc9nr.feather")


p1_motifEnrichmentTable_wGenes <- cisTarget(p1_list, motifRankings,
                                            motifAnnot=motifAnnotations_mgi)
p1_motifEnrichmentTable_wGenes_wLogo <- addLogo(p1_motifEnrichmentTable_wGenes)
saveRDS(p1_motifEnrichmentTable_wGenes, "/data/kollal2/scenic/20190204_p1_motifEnrichmentTable_wGenes.Rds")
saveRDS(p1_motifEnrichmentTable_wGenes_wLogo, "/data/kollal2/scenic/20190204_p1_motifEnrichmentTable_wGenes_wLogo.Rds")


e14_motifEnrichmentTable_wGenes <- cisTarget(e14_list, motifRankings,
                                            motifAnnot=motifAnnotations_mgi)
e14_motifEnrichmentTable_wGenes_wLogo <- addLogo(e14_motifEnrichmentTable_wGenes)
saveRDS(e14_motifEnrichmentTable_wGenes, "/data/kollal2/scenic/20190204_e14_motifEnrichmentTable_wGenes.Rds")
saveRDS(e14_motifEnrichmentTable_wGenes_wLogo, "/data/kollal2/scenic/20190204_e14_motifEnrichmentTable_wGenes_wLogo.Rds")


e16_motifEnrichmentTable_wGenes <- cisTarget(e16_list, motifRankings,
                                            motifAnnot=motifAnnotations_mgi)
e16_motifEnrichmentTable_wGenes_wLogo <- addLogo(e16_motifEnrichmentTable_wGenes)
saveRDS(e16_motifEnrichmentTable_wGenes, "/data/kollal2/scenic/20190204_e16_motifEnrichmentTable_wGenes.Rds")
saveRDS(e16_motifEnrichmentTable_wGenes_wLogo, "/data/kollal2/scenic/20190204_e16_motifEnrichmentTable_wGenes_wLogo.Rds")


