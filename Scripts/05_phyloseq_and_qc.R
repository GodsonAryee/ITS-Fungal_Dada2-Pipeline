library(phyloseq)
library(Biostrings)
library(data.table)
library(dplyr)

source("scripts/00_config.R")

seqtab.nochim <- readRDS(file.path(out_dirs$rds, "seqtab_nochim.rds"))
taxa <- readRDS(file.path(out_dirs$rds, "taxa_assignments.rds"))

Metadata <- read.delim(metadata_file, header = TRUE, row.names = 1, sep = "\t")

# Ensure rownames match (your original code does cleaning)
# Optional: clean sample IDs from seqtab
sample_ids_clean <- sub("_.*", "", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- sample_ids_clean

# Check overlap
missing <- setdiff(rownames(seqtab.nochim), rownames(Metadata))
if(length(missing) > 0) {
  warning("Some seqtab samples missing from Metadata:\n", paste(missing, collapse=", "))
}

ps1 <- phyloseq(
  otu_table(seqtab.nochim, taxa_are_rows = FALSE),
  sample_data(Metadata),
  tax_table(taxa)
)

# Replace ASV sequences with ASV IDs
dna <- DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))

# Filter non-target taxa (same logic as your code)
PS_clean <- ps1 %>%
  subset_taxa(Kingdom != "Eukaryota" & Kingdom != "NA" &
                Order != "Chloroplast" & Family != "Mitochondria")

Reads_per_sample <- data.table(
  as(sample_data(PS_clean), "data.frame"),
  TotalReads = sample_sums(PS_clean), keep.rownames = TRUE
)
setnames(Reads_per_sample, "rn", "SampleID")

write.csv(Reads_per_sample, file.path(out_dirs$tables, "Reads_per_sample_Water.csv"), row.names = FALSE)
saveRDS(PS_clean, file.path(out_dirs$rds, "ps_water.rds"))
