library(ggplot2)

# infile = fichier CSV récapitulatif des SV avec annotation ORF
infile <- "../../Files/plot_files/summary_table_sv_with_ORF.csv"
# out_counts = fichier CSV de sortie pour les statistiques de nombre moyen de SV par génération
out_counts <- "../../Plots/sv_per_generation_summary.csv"
# out_plot = chemin du plot de sortie pour le plot
out_plot <- "../../Plots/sv_per_generation.png"

#générations d'intérêt
gen_order <- c("P15", "P25", "P27", "P30", "P50", "P65")

# Lecture du fichier CSV
d <- read.csv(infile, stringsAsFactors = FALSE)
d$replicate <- suppressWarnings(as.numeric(d$replicate)) 
d <- d[!is.na(d$replicate) & d$replicate >= 1 & d$replicate <= 10, ]
d <- d[d$generation %in% gen_order, ] # Filtrer pour ne garder que les générations d'intérêt

# Compter le nombre de SV par échantillon (génération + replicate)
sample_counts <- aggregate(
  position ~ generation + replicate,
  data = d,
  FUN = length
)
colnames(sample_counts)[colnames(sample_counts) == "position"] <- "n_sv"

summary_mean <- aggregate(n_sv ~ generation, data = sample_counts, FUN = mean)
summary_sd <- aggregate(n_sv ~ generation, data = sample_counts, FUN = sd)
summary_n <- aggregate(n_sv ~ generation, data = sample_counts, FUN = length)

summary_df <- Reduce(
  function(x, y) merge(x, y, by = "generation"),
  list(summary_mean, summary_sd, summary_n)
)
colnames(summary_df) <- c("generation", "mean_n_sv", "sd_n_sv", "n_samples")
summary_df$sem_n_sv <- summary_df$sd_n_sv / sqrt(summary_df$n_samples)
summary_df$generation <- factor(summary_df$generation, levels = gen_order)
summary_df <- summary_df[order(summary_df$generation), ]

write.csv(summary_df, out_counts, row.names = FALSE)

p <- ggplot(summary_df, aes(x = generation, y = mean_n_sv, group = 1)) +
  geom_line(linewidth = 1.0, color = "#2B6CB0") +
  geom_point(size = 2.5, color = "#2B6CB0") +
  geom_errorbar(
    aes(ymin = mean_n_sv - sem_n_sv, ymax = mean_n_sv + sem_n_sv),
    width = 0.15,
    linewidth = 0.8,
    color = "#2B6CB0"
  ) +
  labs(
    title = "Evolution du nombre moyen de mutations SV par generation",
    x = "Generation",
    y = "Nombre moyen de SV par echantillon"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank())

ggsave(out_plot, p, width = 8, height = 5, dpi = 300)
cat("Saved:", out_counts, "\n")
cat("Saved:", out_plot, "\n")
