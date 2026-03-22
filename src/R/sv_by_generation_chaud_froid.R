library(ggplot2)
# Récupération du chemin du script 
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))  # Dossier du script

# Fichier d'entrée : table des variants SV
infile <- file.path(script_dir, "../../Files/plot_files/summary_table_sv_with_ORF.csv")
# Fichier de sortie csv
out_counts <- file.path(script_dir, "../../Files/plot_files/sv_per_generation_chaud_froid_summary.csv")
# Fichier de sortie : graphique enregistré dans le dossier Plots
out_plot <- file.path(script_dir, "../../Plots/sv_per_generation_chaud_froid.png")

# Ordre des générations
gen_order <- c("P15", "P25", "P27", "P30", "P50", "P65")

#lecture du fichier d'entrée
d <- read.csv(infile, stringsAsFactors = FALSE)
# Conversion de replicate en numérique 
d$replicate <- suppressWarnings(as.numeric(d$replicate))
# garder échantillons entre 1 et 10
d <- d[!is.na(d$replicate) & d$replicate >= 1 & d$replicate <= 10, ]
# garder uniquement les générations d'intérêt
d <- d[d$generation %in% gen_order, ]

# Comptage des SV par échantillon
sample_counts <- aggregate(
  position ~ generation + replicate,
  data = d,
  FUN = length  # compte le nombre de variants
)

# Renommer la colonne comptée
colnames(sample_counts)[colnames(sample_counts) == "position"] <- "n_sv"

sample_counts$condition <- ifelse(
  sample_counts$generation == "P15",
  "Avant choc",  # condition contrôle avant choc thermique
  ifelse(sample_counts$replicate <= 5,
  "Froid",       # echantillons 1 à 5
  "Chaud")       # echantillons 6 à 10
)

# Calcul des statistiques par génération et condition
summary_mean <- aggregate(n_sv ~ generation + condition, data = sample_counts, FUN = mean)
summary_sd <- aggregate(n_sv ~ generation + condition, data = sample_counts, FUN = sd)
summary_n <- aggregate(n_sv ~ generation + condition, data = sample_counts, FUN = length)
# Fusion des résultats (mean, sd, n)
summary_df <- Reduce(
  function(x, y) merge(x, y, by = c("generation", "condition")),
  list(summary_mean, summary_sd, summary_n)
)
# Renommage des colonnes
colnames(summary_df) <- c("generation", "condition", "mean_n_sv", "sd_n_sv", "n_samples")
# Calcul de l'erreur standard (SEM)
summary_df$sem_n_sv <- summary_df$sd_n_sv / sqrt(summary_df$n_samples)
# Ordre des facteurs pour un affichage cohérent
summary_df$generation <- factor(summary_df$generation, levels = gen_order)
summary_df$condition <- factor(summary_df$condition, levels = c("Avant choc", "Froid", "Chaud"))
# Tri des données
summary_df <- summary_df[order(summary_df$generation, summary_df$condition), ]
# Sauvegarde du tableau résumé
write.csv(summary_df, out_counts, row.names = FALSE)

# Création du graphique
p <- ggplot(summary_df, aes(x = generation, y = mean_n_sv, fill = condition)) +
  geom_col(
    position = position_dodge(width = 0.75),  # barres côte à côte
    width = 0.65,
    alpha = 0.85
  ) +
  geom_errorbar(
    aes(ymin = mean_n_sv - sem_n_sv, ymax = mean_n_sv + sem_n_sv),  
    position = position_dodge(width = 0.75),
    width = 0.15,
    linewidth = 0.8
  ) +
  scale_fill_manual(
    values = c(
      "Avant choc" = "#6B7280",
      "Froid" = "#1f78b4",  # bleu pour froid
      "Chaud" = "#e31a1c"   # rouge pour chaud
    )
  ) +
  labs(
    title = "Nombre moyen de variants structuraux entre P25 et P27",
    subtitle = "Comparaison des groupes froid et chaud",
    x = "Generation",
    y = "Nombre moyen de SV par echantillon",
    fill = "Condition"
  ) +
  theme_minimal(base_size = 12) +
  theme(panel.grid.minor = element_blank()) 

# Sauvegarde du graphique
ggsave(out_plot, p, width = 7.5, height = 5.5, dpi = 300)

# Messages de confirmation
cat("Sauvegardé:", out_counts, "\n")
cat("Sauvegardé:", out_plot, "\n")
