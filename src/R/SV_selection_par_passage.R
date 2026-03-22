library(ggplot2)  
# Récupération du chemin du script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))  # Dossier contenant le script

# Fichier d'entrée
infile <- file.path(script_dir, "../../Files/plot_files/sv_sc_turnover.csv")
# Fichier de sortie : plot enregistré dans le dossier Plots
outfile <- file.path(script_dir, "../../Plots/sv_selection_par_passage.png")
# Lecture du fichier CSV
d <- read.csv(infile, stringsAsFactors = FALSE)

# Définition de l'ordre des passages 
d$passage <- factor(
  d$passage,
  levels = c("control_P15_P25", "stressed_P25_P27", "control_P30_P50")
)

# Définition de l'ordre des catégories de variants, ceux gagnés, perdus ou conservés dans la generation d'arrivée
d$category <- factor(
  d$category,
  levels = c("shared", "gained_end", "lost_from_start")
)


# Construction du plot
p <- ggplot(d, aes(x = factor(replicate), y = n_variants, fill = category)) +
  geom_col(width = 0.75) +  
  facet_grid(~passage) +  # Un graphique par passage 
    labs(
    title = "Dynamique des variants structuraux entre passages",
    x = "Echantillon (replicate)",
    y = "Nombre de variants",
    fill = "Categorie"
  ) +
  theme_minimal(base_size = 12)  

# Sauvegarde du graphique en fichier PNG
ggsave(outfile, p, width = 10, height = 6, dpi = 300)
# Message de confirmation
cat("Sauvegardé:", outfile, "\n")