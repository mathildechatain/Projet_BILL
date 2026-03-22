library(ggplot2)

infile <- "../../Files/plot_files/snp_sc_counts.csv"
outfile <- "../../Plots/count_snp_by_echantillon_gen.png"

d <- read.csv(infile, stringsAsFactors = FALSE)
d$phase <- factor(d$phase, levels = c("control_P15_P25", "stressed_P25_P27", "control_P30_P50"))
d$generation <- factor(d$generation, levels = c("P15", "P25", "P27", "P30", "P50"))
d$lineage_group <- factor(d$lineage_group, levels = c("lineages_1_5", "lineages_6_10"))

p <- ggplot(d, aes(x = generation, y = n_SNP, group = replicate, color = factor(replicate))) +
  geom_line(linewidth = 0.85, alpha = 0.85) +
  geom_point(size = 2.0) +
  facet_grid(phase ~ lineage_group, scales = "free_x") +
  labs(
    title = "SNP totaux: controle P15->P25 vs stress P25->P27 vs controle P30->P50",
    x = "Passage",
    y = "Nombre de SNP",
    color = "Lignee"
  ) +
  theme_minimal(base_size = 12)


ggsave(outfile, p, width = 10, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")
