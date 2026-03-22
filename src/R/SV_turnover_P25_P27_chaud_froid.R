library(ggplot2)

args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))

infile <- file.path(script_dir, "../../Files/plot_files/summary_table_sv_with_ORF.csv")
outfile <- file.path(script_dir, "../../Plots/sv_turnover_P25_P27_chaud_froid.png")
summary_outfile <- file.path(script_dir, "../../Files/plot_files/sv_turnover_P25_P27_chaud_froid.csv")

d <- read.csv(infile, stringsAsFactors = FALSE)
d$replicate <- as.character(d$replicate)
d <- d[d$generation %in% c("P25", "P27"), ]

common_replicates <- sort(intersect(
  unique(d$replicate[d$generation == "P25"]),
  unique(d$replicate[d$generation == "P27"])
))

d <- d[d$replicate %in% common_replicates, ]
d$variant_id <- paste(d$chrom, d$position, d$svtype, d$ref, d$alt, d$svlen, sep = "|")

result_list <- list()
idx <- 1

for (rep in common_replicates) {
  rep_data <- d[d$replicate == rep, ]
  p25_ids <- unique(rep_data$variant_id[rep_data$generation == "P25"])
  p27_ids <- unique(rep_data$variant_id[rep_data$generation == "P27"])

  gained <- length(setdiff(p27_ids, p25_ids))
  lost <- length(setdiff(p25_ids, p27_ids))
  condition <- ifelse(as.numeric(rep) <= 5, "Froid", "Chaud")

  result_list[[idx]] <- data.frame(
    replicate = rep,
    condition = condition,
    category = "SV acquis",
    n_sv = gained,
    stringsAsFactors = FALSE
  )
  idx <- idx + 1

  result_list[[idx]] <- data.frame(
    replicate = rep,
    condition = condition,
    category = "SV perdus",
    n_sv = lost,
    stringsAsFactors = FALSE
  )
  idx <- idx + 1
}

plot_data <- do.call(rbind, result_list)
plot_data$condition <- factor(plot_data$condition, levels = c("Froid", "Chaud"))
plot_data$category <- factor(plot_data$category, levels = c("SV perdus", "SV acquis"))

write.csv(plot_data, summary_outfile, row.names = FALSE)

summary_mean <- aggregate(n_sv ~ condition + category, data = plot_data, FUN = mean)
summary_sd <- aggregate(n_sv ~ condition + category, data = plot_data, FUN = sd)
summary_n <- aggregate(n_sv ~ condition + category, data = plot_data, FUN = length)

summary_df <- Reduce(
  function(x, y) merge(x, y, by = c("condition", "category")),
  list(summary_mean, summary_sd, summary_n)
)
colnames(summary_df) <- c("condition", "category", "mean_n_sv", "sd_n_sv", "n_samples")
summary_df$sem_n_sv <- summary_df$sd_n_sv / sqrt(summary_df$n_samples)

p <- ggplot(summary_df, aes(x = condition, y = mean_n_sv, fill = category)) +
  geom_col(
    position = position_dodge(width = 0.75),
    width = 0.65,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(ymin = mean_n_sv - sem_n_sv, ymax = mean_n_sv + sem_n_sv),
    position = position_dodge(width = 0.75),
    width = 0.15,
    linewidth = 0.8
  ) +
  geom_jitter(
    data = plot_data,
    aes(x = condition, y = n_sv, color = category),
    size = 1.8,
    alpha = 0.7,
    position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.75),
    show.legend = FALSE
  ) +
  scale_fill_manual(values = c("SV perdus" = "#d95f02", "SV acquis" = "#1b9e77")) +
  scale_color_manual(values = c("SV perdus" = "#d95f02", "SV acquis" = "#1b9e77")) +
  labs(
    title = "Gains et pertes de SV entre P25 et P27",
    subtitle = "Comparaison des groupes froid et chaud",
    x = "Condition thermique",
    y = "Nombre moyen de SV par échantillon",
    fill = "Catégorie"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor = element_blank(),
    legend.position = "top"
  )

ggsave(outfile, p, width = 7.5, height = 5.5, dpi = 300)
cat("Saved:", summary_outfile, "\n")
cat("Saved:", outfile, "\n")
