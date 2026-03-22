library(ggplot2)

infile <- "../../Files/plot_files/snp_turnover.csv"
outfile <- "../../Plots/snp_selection_par_passage.png"

d <- read.csv(infile, stringsAsFactors = FALSE)
d$phase <- factor(d$phase, levels = c("control_P15_P25", "stressed_P25_P27", "control_P30_P50"))
d$category <- factor(d$category, levels = c("shared", "gained", "lost"))

p <- ggplot(d, aes(x = factor(replicate), y = n_variants, fill = category)) +
  geom_col(width = 0.75) +
  facet_grid(~phase) +
  labs(
    title = "Dinamique des SNP entre passage",
    x = "Echantillon (replicate)",
    y = "Nombre de variants",
    fill = "Categorie"
  ) +
  theme_minimal(base_size = 12)

ggsave(outfile, p, width = 10, height = 6, dpi = 300)
cat("Saved:", outfile, "\n")


### Partie stress thermique
library(dplyr)

outfile_stress <- "../../Plots/snp_selection_par_passage.png"

d2 <- read.csv(infile, stringsAsFactors = FALSE)

# Préparation des données
d_sub <- d2 %>%
  filter(phase == "stressed_P25_P27") %>%
  mutate(
    condition = ifelse(replicate <= 5, "stress froid", "stress chaud"),
    category = factor(category, levels = c("shared", "gained", "lost"))
  )

# Étape 1 : total par replicate
d_rep <- d_sub %>%
  group_by(replicate, condition, category) %>%
  summarise(n_variants = sum(n_variants), .groups = "drop")

# Étape 2 : moyenne + erreur standard
d_plot <- d_rep %>%
  group_by(condition, category) %>%
  summarise(
    mean_variants = mean(n_variants),
    sd = sd(n_variants),
    se = sd / sqrt(n()),
    .groups = "drop"
  )

# Plot
p2 <- ggplot(d_plot, aes(x = condition, y = mean_variants, fill = category)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_errorbar(
    aes(ymin = mean_variants - se, ymax = mean_variants + se),
    position = position_dodge(width = 0.8),
    width = 0.2
  ) +
  labs(
    title = "Comparaison du nombre de SNP en fonction de leur categorie entre la condition chaud et froid entre les passages P25 et P27",
    x = "Condition",
    y = "Nombre moyen de variants",
    fill = "Categorie"
  ) +
  theme_minimal(base_size = 12)

ggsave(outfile_stress, p2, width = 8, height = 5, dpi = 300)
cat("Saved:", outfile, "\n")
