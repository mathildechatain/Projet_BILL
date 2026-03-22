library(ggplot2) 

# Récupération du chemin du script 
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))  # Dossier du script
# infile = fichier CSV récapitulatif des SV avec annotation ORF
infile <- file.path(script_dir, "../../Files/plot_files/summary_table_sv_with_ORF.csv")
# out_counts = fichier CSV de sortie
out_counts <- file.path(script_dir, "../../Files/plot_files/sv_per_generation_summary.csv")
# out_plot = chemin du fichier plot de sortie qui correspond au dossier Plots
out_plot <- file.path(script_dir, "../../Plots/sv_per_generation.png")

# Générations d'intérêt (ordre spécifique)
gen_order <- c("P15", "P25", "P27", "P30", "P50", "P65")

# Lecture du fichier CSV
d <- read.csv(infile, stringsAsFactors = FALSE)
# Conversion de la colonne replicate en numérique
d$replicate <- suppressWarnings(as.numeric(d$replicate)) 
# garder les replicates valides entre 1 et 10
d <- d[!is.na(d$replicate) & d$replicate >= 1 & d$replicate <= 10, ]
# garder uniquement les générations d'intérêt
d <- d[d$generation %in% gen_order, ]

# Comptage des SV par échantillon

# Calcul du nombre de SV (basé sur la colonne position) par génération et échantillon
sample_counts <- aggregate(
  position ~ generation + replicate,
  data = d,
  FUN = length 
)

# Renommage de la colonne "position" en "n_sv"
colnames(sample_counts)[colnames(sample_counts) == "position"] <- "n_sv"

# Moyenne du nombre de SV par génération
summary_mean <- aggregate(n_sv ~ generation, data = sample_counts, FUN = mean)
# Écart-type du nombre de SV par génération
summary_sd <- aggregate(n_sv ~ generation, data = sample_counts, FUN = sd)
# Nombre d'échantillons par génération
summary_n <- aggregate(n_sv ~ generation, data = sample_counts, FUN = length)
# Fusion des trois tableaux (mean, sd, n)
summary_df <- Reduce(
  function(x, y) merge(x, y, by = "generation"),
  list(summary_mean, summary_sd, summary_n)
)

# Renommage des colonnes
colnames(summary_df) <- c("generation", "mean_n_sv", "sd_n_sv", "n_samples")
# Calcul de l'erreur standard de la moyenne (SEM)
summary_df$sem_n_sv <- summary_df$sd_n_sv / sqrt(summary_df$n_samples)
# Transformation pour respecter l'ordre des générations
summary_df$generation <- factor(summary_df$generation, levels = gen_order)

# Tri du tableau selon cet ordre
summary_df <- summary_df[order(summary_df$generation), ]
# Sauvegarde du tableau résumé en CSV, qu'on peut retrouver dans le dossier plot_files
write.csv(summary_df, out_counts, row.names = FALSE)

# Construction du plot avec ggplot2
p <- ggplot(summary_df, aes(x = generation, y = mean_n_sv, group = 1)) +
  geom_line(linewidth = 1.0, color = "#2B6CB0") +  # ligne reliant les points
  geom_point(size = 2.5, color = "#2B6CB0") +      # points des moyennes
  geom_errorbar(
    aes(ymin = mean_n_sv - sem_n_sv, ymax = mean_n_sv + sem_n_sv),  # barres d'erreur (SEM)
    width = 0.15,
    linewidth = 0.8,
    color = "#2B6CB0"
  ) +
  labs(
    title = "Evolution du nombre moyen de variants structuraux par generation",
    x = "Generation",
    y = "Nombre moyen de SV par echantillon"
  ) +
  theme_minimal(base_size = 12) +  
  theme(panel.grid.minor = element_blank())  

# Sauvegarde du plot
ggsave(out_plot, p, width = 8, height = 5, dpi = 300)

# Messages de confirmation
cat("Sauvegardé", out_counts, "\n")
cat("Sauvegardé:", out_plot, "\n")