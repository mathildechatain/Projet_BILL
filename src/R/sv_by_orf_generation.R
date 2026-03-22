library(ggplot2) 
# Récupération du chemin du script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))  # Dossier du script

# Fichier d'entrée : tableau des SV par ORF et génération 
infile <- file.path(script_dir, "../../Files/plot_files/sv_by_orf_generation_wide.csv")
# Fichier de sortie : heatmap
out_plot <- file.path(script_dir, "../../Plots/sv_by_orf_generation.png")
# Ordre des générations
gen_order <- c("P15", "P25", "P27", "P30", "P50", "P65")

# Nombre d'ORFs à afficher (top 30 ici)
top_n_orfs <- 30

# Lecture du CSV 
d <- read.csv(infile, stringsAsFactors = FALSE, check.names = FALSE)
# Sélection des ORFs les plus importants (top N)
d <- head(d, top_n_orfs)


plot_long <- reshape(
  d[, c("orf", gen_order)],  # colonnes utilisées
  varying = gen_order,       # colonnes à transformer
  v.names = "n_unique_sv",   # nom de la variable créée
  timevar = "generation",    # nom de la variable génération
  times = gen_order,         # valeurs de génération
  direction = "long"         # transformation vers format long
)

# Conversion en facteur ordonné pour garder l'ordre des générations
plot_long$generation <- factor(plot_long$generation, levels = gen_order, ordered = TRUE)
# Ordre des ORFs inversé pour affichage (du plus important en haut)
plot_long$orf <- factor(plot_long$orf, levels = rev(d$orf))

# Création de la heatmap
p <- ggplot(plot_long, aes(x = generation, y = orf, fill = n_unique_sv)) +
  geom_tile(color = "white", linewidth = 0.3) +  # cases de la heatmap
  geom_text(aes(label = n_unique_sv), size = 2.5) +  # valeurs affichées dans les cases
  scale_fill_gradient(low = "#fff7bc", high = "#d95f0e") +  # couleur foncé = bcp de variants et clair - peu de variants
  labs(
    title = paste("Top", nrow(d), "ORF par nombre de SV uniques et par generation"),
    x = "Generation",
    y = "ORF",
    fill = "Nombre de SV uniques"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),  
    axis.text.y = element_text(size = 8)  
  )
height <- max(8, 0.3 * nrow(d) + 2)

# Sauvegarde du graphique
ggsave(out_plot, p, width = 9, height = height, dpi = 300)

# Message de confirmation
cat("Sauvegardé:", out_plot, "\n")