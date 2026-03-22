library(ggplot2) # librairie pour créer des graphiques (plots)

# Récupération des arguments de la ligne de commande
args <- commandArgs(trailingOnly = FALSE)

# Définition du motif pour trouver le chemin du script
file_arg <- "--file="
# Extraction du chemin du script 
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
# Dossier contenant le script
script_dir <- dirname(normalizePath(script_path))

# Fichier d'entrée (csv contenant les données des SVs)
infile <- file.path(script_dir, "../../Files/plot_files/sv_sc_counts.csv") 
# Fichier de sortie (image du plot) dans le dossier Plots
outfile <- file.path(script_dir, "../../Plots/count_sv_by_echantillon_gen.png") 
# Lecture du fichier d'entrée
d <- read.csv(infile, stringsAsFactors = FALSE)

# Transformation de la variable "passage" en facteur 
# permet de contrôler l'ordre d'affichage dans le plot
d$passage <- factor(d$passage, levels = c("control_P15_P25", "stressed_P25_P27", "control_P30_P50"))

# Transformation de la variable "generation" en facteur ordonné
# Correspond aux différentes générations étudiées
d$generation <- factor(d$generation, levels = c("P15", "P25", "P27", "P30", "P50"))

# Transformation de la variable "lineage_group" en facteur
# Permet de regrouper les échantillons en deux catégories
d$lineage_group <- factor(d$lineage_group, levels = c("lineages_1_5", "lineages_6_10")) 

# Création du graphique :
# axe x correpsond aux générations (passages)
# axe y correspond au nombre de variants structuraux (SV)
# group permet de relier les points d'un même échantillon ( suivre leur evolution entre passage )
# color différencie les réplicats (échantillons) par couleur
p <- ggplot(d, aes(x = generation, y = n_SV, group = replicate, color = factor(replicate))) +
  
  # Ajout des lignes reliant les points
  geom_line(linewidth = 0.85, alpha = 0.85) +
  geom_point(size = 2.0) +
  
  # Création d'une grille de sous-graphes :
  # lignes = passages, colonnes = groupes de lignées
  # scales = "free_x" permet d'adapter l'axe x si besoin
  facet_grid(passage ~ lineage_group, scales = "free_x") +
  
  # Ajout des titres et labels
  labs(
    title = "Nombre de variants structuraux par echantillon et par generation",
    x = "Passage",
    y = "Nombre de SV",
    color = "Lignee"
  ) +
  theme_minimal(base_size = 12)

# Sauvegarde du graphique en fichier png dans le dossier de sortie
ggsave(outfile, p, width = 10, height = 6, dpi = 300)
