library(ggplot2)

infile <- "../../Files/plot_files/sv_sc_turnover.csv"
outfile <- "../../Plots/sv_selection_par_passage.png"

d <- read.csv(infile, stringsAsFactors = FALSE)
d$phase <- factor(d$phase, levels = c("control_P15_P25", "stressed_P25_P27", "control_P30_P50"))
d$category <- factor(d$category, levels = c("shared", "gained_end", "lost_from_start"))

p <- ggplot(d, aes(x = factor(replicate), y = n_variants, fill = category)) +
  geom_col(width = 0.75) +
  facet_grid(~phase) +
  labs(
    title = "Selection SV: passages controles et phase stress",
    x = "Echantillon (replicate)",
    y = "Nombre de variants",
    fill = "Categorie"
  ) +
  theme_minimal(base_size = 12)

ggsave(outfile, p, width = 10, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")
