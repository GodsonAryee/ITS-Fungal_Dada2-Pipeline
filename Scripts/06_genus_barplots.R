library(phyloseq)
library(ggplot2)
library(dplyr)

source("scripts/00_config.R")

PS_clean <- readRDS(file.path(out_dirs$rds, "ps_water.rds"))

PS_genus <- tax_glom(PS_clean, taxrank = "Genus")
PS_genus_rel <- transform_sample_counts(PS_genus, function(x) x / sum(x))
df_genus <- psmelt(PS_genus_rel)

p <- ggplot(df_genus, aes(x = Sample, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Sample", y = "Relative Abundance", title = "Genus-Level Relative Abundance") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
  guides(fill = guide_legend(ncol = 1))

ggsave(filename = file.path(out_dirs$figures, "genus_relative_abundance.png"),
       plot = p, width = 14, height = 6, dpi = 300)
sink(file.path(out_dirs$logs, "sessionInfo.txt"))
cat("Run date:", as.character(Sys.time()), "\n\n")
sessionInfo()
sink()
