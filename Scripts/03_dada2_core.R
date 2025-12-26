library(dada2)

source("scripts/00_config.R")
paths1 <- readRDS(file.path(out_dirs$rds, "paths_step01.rds"))
paths2 <- readRDS(file.path(out_dirs$rds, "paths_step02.rds"))

fnFs.cut <- paths2$fnFs.cut
fnRs.cut <- paths2$fnRs.cut

filtFs <- file.path(out_dirs$filt_dir, basename(fnFs.cut))
filtRs <- file.path(out_dirs$filt_dir, basename(fnRs.cut))

out_filter2 <- filterAndTrim(
  fnFs.cut, filtFs, fnRs.cut, filtRs,
  maxN = params$maxN, maxEE = params$maxEE, truncQ = params$truncQ,
  minLen = params$minLen, rm.phix = params$rm.phix,
  compress = params$compress, multithread = params$multithread
)
write.csv(out_filter2, file.path(out_dirs$tables, "post_cutadapt_filter_summary.csv"))

errF <- learnErrors(filtFs, multithread = params$multithread)
errR <- learnErrors(filtRs, multithread = params$multithread)
pdf(file.path(out_dirs$figures, "error_profile_forward.pdf"))
plotErrors(errF, nominalQ = TRUE)
dev.off()

dadaFs <- dada(filtFs, err = errF, multithread = params$multithread)
dadaRs <- dada(filtRs, err = errR, multithread = params$multithread)

mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
seqtab <- makeSequenceTable(mergers)
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus",
                                    multithread = params$multithread, verbose = TRUE)

getN <- function(x) sum(getUniques(x))
track <- cbind(out_filter2,
               sapply(dadaFs, getN),
               sapply(dadaRs, getN),
               sapply(mergers, getN),
               rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- paths1$sample.names

write.csv(track, file.path(out_dirs$tables, "read_tracking.csv"))

saveRDS(seqtab.nochim, file.path(out_dirs$rds, "seqtab_nochim.rds"))
saveRDS(track, file.path(out_dirs$rds, "track_table.rds"))
