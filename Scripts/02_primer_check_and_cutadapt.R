library(dada2)
library(ShortRead)
library(Biostrings)

source("scripts/00_config.R")
paths <- readRDS(file.path(out_dirs$rds, "paths_step01.rds"))

fnFs.filtN <- paths$fnFs.filtN
fnRs.filtN <- paths$fnRs.filtN

allOrients <- function(primer) {
  dna <- DNAString(primer)
  orients <- c(Forward = dna,
               Complement = complement(dna),
               Reverse = reverse(dna),
               RevComp = reverseComplement(dna))
  sapply(orients, toString)
}
primerHits <- function(primer, fn) {
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  sum(nhits > 0)
}

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

# Optional primer check on first sample
primer_check_before <- rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]])
)
write.csv(primer_check_before, file.path(out_dirs$tables, "primer_hits_before_cutadapt.csv"))

stopifnot(file.exists(cutadapt))
system2(cutadapt, args="--version")

dir.create(out_dirs$cut_dir, showWarnings = FALSE, recursive = TRUE)
fnFs.cut <- file.path(out_dirs$cut_dir, basename(paths$fnFs))
fnRs.cut <- file.path(out_dirs$cut_dir, basename(paths$fnRs))

FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)

for(i in seq_along(fnFs.filtN)) {
  system2(cutadapt, args = c(
    R1.flags, R2.flags, "-n", "2",
    "-o", shQuote(fnFs.cut[i]),
    "-p", shQuote(fnRs.cut[i]),
    shQuote(fnFs.filtN[i]),
    shQuote(fnRs.filtN[i])
  ))
}

# Optional primer check after
primer_check_after <- rbind(
  FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
  FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
  REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
  REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])
)
write.csv(primer_check_after, file.path(out_dirs$tables, "primer_hits_after_cutadapt.csv"))

saveRDS(list(fnFs.cut=fnFs.cut, fnRs.cut=fnRs.cut),
        file.path(out_dirs$rds, "paths_step02.rds"))
