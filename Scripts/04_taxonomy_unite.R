library(dada2)

source("scripts/00_config.R")

seqtab.nochim <- readRDS(file.path(out_dirs$rds, "seqtab_nochim.rds"))
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = params$multithread, tryRC = TRUE)

saveRDS(taxa, file.path(out_dirs$rds, "taxa_assignments.rds"))
