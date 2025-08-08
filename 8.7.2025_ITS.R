# -----------------------------
# Load Required Libraries
# -----------------------------
library(dada2); packageVersion("dada2")
library(ShortRead); packageVersion("ShortRead")
library(Biostrings); packageVersion("Biostrings")

# -----------------------------
# Set Paths
# -----------------------------
path <- "/Users/godsonaryee/Documents/Collaborators_Data/Juddy/ITS PRESSED WATER DATA/"  # CHANGE THIS IF NEEDED
setwd(path)

# -----------------------------
# List and Sort Input Files
# -----------------------------
all_files <- list.files(path, pattern = ".fastq.gz$", full.names = TRUE)
fnFs <- sort(all_files[grepl("_R1.fastq.gz$", all_files)])
fnRs <- sort(all_files[grepl("_R2.fastq.gz$", all_files)])

# Sanity check
if(length(fnFs) != length(fnRs)) stop("Forward and reverse reads do not match!")

# -----------------------------
# Extract Sample Names
# -----------------------------
get.sample.name <- function(fname) {
  sub(" \\(paired\\)_ITS1_R[12]\\.fastq\\.gz", "", basename(fname))
}
sample.names <- unname(sapply(fnFs, get.sample.name))

# Confirm correct pairing
stopifnot(all(get.sample.name(fnFs) == get.sample.name(fnRs)))

# -----------------------------
# Pre-Filter Reads with Ns
# -----------------------------
filtN_dir <- file.path(path, "filtN")
dir.create(filtN_dir, showWarnings = FALSE)

fnFs.filtN <- file.path(filtN_dir, basename(fnFs))
fnRs.filtN <- file.path(filtN_dir, basename(fnRs))

out <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN,
                     maxN = 0, maxEE = c(2,2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out)

# -----------------------------
# Primer Definitions
# -----------------------------
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"

# All Orientations
allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna,
               Complement = complement(dna),
               Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Primer Hit Checker
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

# Check Hits (optional)
rbind(
  FWD.ForwardReads  = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
  FWD.ReverseReads  = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
  REV.ForwardReads  = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
  REV.ReverseReads  = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
)

# -----------------------------
# Run Cutadapt
# -----------------------------
#cutadapt <- "/usr/local/bin/cutadapt"  # CHANGE THIS TO YOUR CUTADAPT PATH
cutadapt <-"/Users/godsonaryee/miniconda3/bin/cutadapt"
stopifnot(file.exists(cutadapt))

system2(cutadapt, args = "--version")  # Check it's working

path.cut <- file.path(path, "cutadapt")
dir.create(path.cut, showWarnings = FALSE)

fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

#for(i in seq_along(fnFs)) {
#system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2,
#"-o", fnFs.cut[i], "-p", fnRs.cut[i],
#fnFs.filtN[i], fnRs.filtN[i]))
#}
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(
    R1.flags, R2.flags, "-n", "2",
    "-o", shQuote(fnFs.cut[i]), 
    "-p", shQuote(fnRs.cut[i]),
    shQuote(fnFs.filtN[i]), 
    shQuote(fnRs.filtN[i])
  ))
}

# Optional check after cutadapt
rbind(
  FWD.ForwardReads  = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
  FWD.ReverseReads  = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
  REV.ForwardReads  = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
  REV.ReverseReads  = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])
)

# -----------------------------
# Filter Again After Primer Removal
# -----------------------------
filtered_dir <- file.path(path.cut, "filtered")
dir.create(filtered_dir, showWarnings = FALSE)

filtFs <- file.path(filtered_dir, basename(fnFs.cut))
filtRs <- file.path(filtered_dir, basename(fnRs.cut))

out <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs,
                     maxN = 0, maxEE = c(2,2), truncQ = 2,
                     minLen = 50, rm.phix = TRUE,
                     compress = TRUE, multithread = TRUE)
head(out)

# -----------------------------
# Learn Error Rates
# -----------------------------
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)
plotErrors(errF, nominalQ = TRUE)

# -----------------------------
# Dereplication & Denoising
# -----------------------------
dadaFs <- dada(filtFs, err = errF, multithread = TRUE)
dadaRs <- dada(filtRs, err = errR, multithread = TRUE)

# Merge Paired Reads
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)

# -----------------------------
# Sequence Table & Chimera Removal
# -----------------------------
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = TRUE, verbose = TRUE)
table(nchar(getSequences(seqtab.nochim)))

# -----------------------------
# Track Read Counts Through Pipeline
# -----------------------------
getN <- function(x) sum(getUniques(x))
track <- cbind(out,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

# -----------------------------
# Taxonomic Assignment
# -----------------------------
unite.ref <- "/Users/godsonaryee/Documents/Collaborators_Data/Juddy/ITS PRESSED WATER DATA/sh_general_release_dynamic_19.02.2025.fasta"  # CHANGE IF NEEDED
taxa <- assignTaxonomy(seqtab.nochim, unite.ref, multithread = TRUE, tryRC = TRUE)
# library(Biostrings)
# 
# # Load the full UNITE reference used for assignTaxonomy
# unite_ref <- readDNAStringSet("/Users/godsonaryee/Documents/Collaborators_Data/Juddy/ITS PRESSED PULP DATA/sh_general_release_dynamic_19.02.2025.fasta")
# 
# # Extract only genus and species from header
# get_species_name <- function(x) {
#   tx <- strsplit(x, ";")[[1]]
#   sp <- grep("s__", tx, value = TRUE)
#   if (length(sp) > 0) {
#     sp_name <- sub("s__", "", sp)
#     sp_name <- gsub("_", " ", sp_name)
#     return(sp_name)
#   } else {
#     return(NA)
#   }
# }
# 
# species_names <- sapply(names(unite_ref), get_species_name)
# valid <- !is.na(species_names)
# species_ref <- unite_ref[valid]
# names(species_ref) <- species_names[valid]
# 
# # Save new species-level FASTA
# writeXStringSet(species_ref, "unite_species_assign.fa")
# 
# # Assign species-level taxonomy
# taxa <- assignSpecies(taxa, "unite_species_assign.fa", tryRC = TRUE, verbose = TRUE)
# 
# #taxa <- assignSpecies(taxa, unite.ref, tryRC = TRUE, verbose = TRUE)
# # -----------------------------
# # Save Outputs
# -----------------------------
saveRDS(seqtab.nochim, file = "seqtab_nochim.rds")
saveRDS(taxa, file = "taxa_assignments.rds")
write.csv(track, file = "read_tracking.csv")
###Making Phyloseq Object###
Metadata <- read.delim("METADATA.txt", header=T, row.names = 1)
head(Metadata); dim(Metadata); class(Metadata)
dim(track)

sample_names(seqtab.nochim)
rownames(Metadata)
sample_names(taxa)
dim(seqtab.nochim)
head(seqtab.nochim)
nrow(seqtab.nochim)
nrow(Metadata)
length(rownames(Metadata))

rownames(seqtab.nochim) <- rownames(Metadata)
sample_ids_seqtab <- rownames(seqtab.nochim)
# Extract sample IDs from seqtab.nochim rownames
sample_ids_clean <- sub("_.*", "", rownames(seqtab.nochim))
rownames(seqtab.nochim) <- sample_ids_clean
all(rownames(seqtab.nochim) %in% rownames(Metadata))
setdiff(rownames(seqtab.nochim), rownames(Metadata))

##create phyloseq object - ASV table
ps1 <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
                sample_data(Metadata),
                tax_table(taxa))
ps1

##replace sequences with ASV numbers
dna <- Biostrings::DNAStringSet(taxa_names(ps1))
names(dna) <- taxa_names(ps1)
ps1 <- merge_phyloseq(ps1, dna)
taxa_names(ps1) <- paste0("ASV", seq(ntaxa(ps1)))
ps1

##Remove chloroplasts, mitochondria, eukaryota, and ASVs without a domain.
PS_clean <- ps1 %>%
  subset_taxa(Kingdom != "Eukaryota" & Kingdom != "NA" & Order   != "Chloroplast" & Family  != "Mitochondria" &            
                Order   != "Chloroplast")                
PS_clean
##Summarizing the number of reads per sample
Reads_per_sample = data.table(as(sample_data(PS_clean), "data.frame"),
                              TotalReads = sample_sums(PS_clean), keep.rownames = TRUE)
setnames(Reads_per_sample, "rn", "SampleID")
write.csv(Reads_per_sample, file="Reads_per_sample_Water.csv")
saveRDS(PS_clean, "ps_water.rds")
#########Genus level Abundances##
# Load necessary packages
library(phyloseq)
library(ggplot2)
library(dplyr)

# Assume `PS_clean` is your phyloseq object
# Agglomerate taxa to genus level
PS_genus <- tax_glom(PS_clean, taxrank = "Genus")

# Transform counts to relative abundances
PS_genus_rel <- transform_sample_counts(PS_genus, function(x) x / sum(x))

# Convert phyloseq object to a data frame for plotting
df_genus <- psmelt(PS_genus_rel)

# Create the stacked bar plot
ggplot(df_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Relative Abundance", title = "Genus-Level Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 1))  # Adjust legend layout if needed

# Load necessary packages
library(phyloseq)
library(ggplot2)
library(dplyr)

# Step 1: Agglomerate to Genus level
PS_genus <- tax_glom(PS_clean, taxrank = "Genus")

# Step 2: Transform counts to relative abundances
PS_genus_rel <- transform_sample_counts(PS_genus, function(x) x / sum(x))

# Step 3: Convert phyloseq object to data frame for plotting
df_genus <- psmelt(PS_genus_rel)

# Step 4: Filter to keep only the top 15 genera by mean relative abundance
top15_genera <- df_genus %>%
  group_by(Genus) %>%
  summarise(mean_abundance = mean(Abundance, na.rm = TRUE)) %>%
  top_n(15, mean_abundance) %>%
  pull(Genus)

df_top15 <- df_genus %>%
  filter(Genus %in% top15_genera)

# Step 5: Plot with ggplot2, faceting by SampleType and grouping by Treatment
ggplot(df_top15, aes(x = Location, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Location", y = "Relative Abundance", title = "Top 15 Genera-Level Relative Abundance by Sample Type") +
  facet_wrap(~ Processing.Period, scales = "free_x") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 1))  # Adjust legend layout if needed
write.csv(df_top15, file = "top15_genera.csv")
####RELATIVE ABUNDANCE FOR EACH SAMPLE TYPE###
##Analyzing only a subset of samples (e.g. Nasal samples)
PS_Nasal_swab <- subset_samples(ps_reads_1000, Sample_type == "Nasal swab")
#Remove samples with fewer than 1000 sequences
ps_reads_1000 = prune_samples(sample_sums(PS_clean) > 1000, PS_clean)

#Create genus level summaries				 
ps_phylum <- tax_glom(PS_clean, taxrank = 'Genus',NArm = FALSE,bad_empty=c(NA, "", " ", "\t"))
ps_phylum_rel = transform_sample_counts(ps_phylum, function(x) x/sum(x)*100)
dat <- psmelt(ps_phylum_rel)
dat$Genus<- as.character(dat$Genus)
phylum_abundance <- aggregate(Abundance~Sample+Genus, dat, FUN=sum)
phylum_abundance <- cast(phylum_abundance, Sample ~ Genus)


# Make sure the columns to use for the merge are identical
identical( sort(Metadata$Description), sort(phylum_abundance$Sample))
# [1] TRUE
genus_abundance_with_metadata <-  merge(Metadata, genus_abundance, by.x = "Description", by.y = "Sample") 
write.csv(phylum_abundance, file="Genus_summary_with_metadata.csv")
