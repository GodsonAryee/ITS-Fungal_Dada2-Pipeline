library(dada2)
library(ShortRead)
library(Biostrings)

source("scripts/00_config.R")
setwd(path)

all_files <- list.files(path, pattern = ".fastq.gz$", full.names = TRUE)
fnFs <- sort(all_files[grepl("_R1.fastq.gz$", all_files)])
fnRs <- sort(all_files[grepl("_R2.fastq.gz$", all_files)])
if(length(fnFs) != length(fnRs)) stop("Forward and reverse reads do not match!")

get.sample.name <- function(fname) {
  sub(" \\(paired\\)_ITS1_R[12]\\.fastq\\.gz", "", basename(fname))
}
sample.names <- unname(sapply(fnFs, get.sample.name))
stopifnot(all(get.sample.name(fnFs) == get.sample.name(fnRs)))

fnFs.filtN <- file.path(out_dirs$filtN_dir, basename(fnFs))
fnRs.filtN <- file.path(out_dirs$filtN_dir, basename(fnRs))

out_prefilter <- filterAndTrim(
  fnFs, fnFs.filtN, fnRs, fnRs.filtN,
  maxN = params$maxN, maxEE = params$maxEE, truncQ = params$truncQ,
  minLen = params$minLen, rm.phix = params$rm.phix,
  compress = params$compress, multithread = params$multithread
)

write.csv(out_prefilter, file.path(out_dirs$tables, "prefilter_filterN_summary.csv"))
saveRDS(list(fnFs=fnFs, fnRs=fnRs, fnFs.filtN=fnFs.filtN, fnRs.filtN=fnRs.filtN,
             sample.names=sample.names),
        file.path(out_dirs$rds, "paths_step01.rds"))
