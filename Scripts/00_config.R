# -----------------------------
# CONFIG: paths + parameters
# -----------------------------

# INPUTS (edit these)
path <- "/Users/godsonaryee/Documents/Collaborators_Data/Juddy/ITS PRESSED WATER DATA/"
metadata_file <- file.path(path, "METADATA.txt")

cutadapt <- "/Users/godsonaryee/miniconda3/bin/cutadapt"
unite.ref <- file.path(path, "sh_general_release_dynamic_19.02.2025.fasta")

# Primer sequences
FWD <- "CTTGGTCATTTAGAGGAAGTAA"
REV <- "GCTGCGTTCTTCATCGATGC"

# Filtering parameters
params <- list(
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  minLen = 50,
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

# OUTPUT directories
out_dirs <- list(
  filtN_dir = file.path(path, "filtN"),
  cut_dir   = file.path(path, "cutadapt"),
  filt_dir  = file.path(path, "cutadapt", "filtered"),
  results   = file.path(path, "results"),
  tables    = file.path(path, "results", "tables"),
  figures   = file.path(path, "results", "figures"),
  rds       = file.path(path, "results", "rds"),
  logs      = file.path(path, "results", "logs")
)

invisible(lapply(out_dirs, dir.create, showWarnings = FALSE, recursive = TRUE))
