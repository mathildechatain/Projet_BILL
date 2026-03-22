# Charger les librairies nécessaires
library(tidyverse)

infile <- "../../Files/plot_files/snp_sc_turnover.csv"
outfile <- "../../Plots/snp_selection_par_passage.png"

# === 1. Charger le fichier ===
data <- read.csv(infile, stringsAsFactors = FALSE)


# === 3. Séparer les ORF multiples ===
data_long <- data %>%
  separate_rows(!!sym("orf"), sep = ";")

# === 4. Compter le nombre de variants par ORF ===
orf_counts <- data_long %>%
  group_by(!!sym("orf")) %>%
  summarise(nb_variants = n()) %>%
  arrange(desc(nb_variants))

top20_orf <- data_long %>%
  group_by(!!sym("orf")) %>%
  summarise(nb_variants = n(), .groups = "drop") %>%
  slice_max(order_by = nb_variants, n = 20)

top20_orf_cat <- data_long %>%
  filter(!!sym("orf") %in% top20_orf[["orf"]]) %>%
  group_by(!!sym("orf"), category) %>%
  summarise(nb_variants = n(), .groups = "drop")

p <- ggplot(top20_orf_cat,
       aes(x = reorder(!!sym("orf"), nb_variants),
           y = nb_variants,
           fill = category)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  labs(
    title = "nombre de snp par orf en fonction de leur categorie",
    x = "ORF",
    y = "Nombre de variants",
    fill = "Catégorie"
  ) +
  theme_minimal()

ggsave(outfile, p, width = 8, height = 5, dpi = 300)
cat("Saved:", outfile, "\n")
