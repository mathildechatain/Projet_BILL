library(ggplot2) 

# Récupération du chemin du script
args <- commandArgs(trailingOnly = FALSE)
file_arg <- "--file="
script_path <- sub(file_arg, "", args[grep(file_arg, args)])
script_dir <- dirname(normalizePath(script_path))

# Fichier d'entrée : table des variants SV
infile <- file.path(script_dir, "../../Files/plot_files/summary_table_sv_with_ORF.csv")
# Fichier de sortie : graphique final
outfile <- file.path(script_dir, "../../Plots/orf_impact_choc_thermique_P25_P27.png")
# Lecture du CSV
d <- read.csv(infile, stringsAsFactors = FALSE)

# Filtrage des générations d'intérêt (avant/après choc)
d <- d[d$generation %in% c("P25", "P27"), ]
# Conversion des echantullons en caractères 
d$replicate <- as.character(d$replicate)

# Identification des echantillons présents dans les deux générations
common_replicates <- intersect(
  unique(d$replicate[d$generation == "P25"]),
  unique(d$replicate[d$generation == "P27"])
)
common_replicates <- sort(common_replicates)

# Filtrage pour ne garder que ces echantillons communs
d <- d[d$replicate %in% common_replicates, ]

# Création d'un identifiant unique pour chaque variant
d$variant_id <- paste(
  d$chrom, d$position, d$svtype, d$ref, d$alt, d$svlen,
  sep = "|"
)

# Fonction pour séparer les ORFs multiples
split_orfs <- function(x) {
  if (is.na(x) || x == "") {
    return("Intergenic")  # cas sans ORF
  }
  parts <- trimws(unlist(strsplit(x, ";", fixed = TRUE)))  # séparation
  parts[parts != ""]  # suppression des valeurs vides
}

# Expansion des ORFs (une ligne par ORF)
rows <- vector("list", nrow(d))
for (i in seq_len(nrow(d))) {
  orfs <- split_orfs(d$orf[i])
  rows[[i]] <- data.frame(
    generation = d$generation[i],
    replicate = d$replicate[i],
    variant_id = d$variant_id[i],
    orf = orfs,
    stringsAsFactors = FALSE
  )
}

# Fusion des lignes en un seul DataFrame
expanded <- do.call(rbind, rows)
# Liste pour stocker les résultats
result_list <- list()
idx <- 1

# Boucle sur chaque echantillon
for (rep in common_replicates) {
  rep_data <- expanded[expanded$replicate == rep, ]
  # Données pour chaque génération
  p25 <- rep_data[rep_data$generation == "P25", c("orf", "variant_id")]
  p27 <- rep_data[rep_data$generation == "P27", c("orf", "variant_id")]
  # Liste des ORFs présents
  all_orfs <- sort(unique(c(p25$orf, p27$orf)))
  for (orf_name in all_orfs) {
    # Variants associés à chaque ORF
    p25_ids <- unique(p25$variant_id[p25$orf == orf_name])
    p27_ids <- unique(p27$variant_id[p27$orf == orf_name])
    # Variants gagnés et perdus
    gained <- setdiff(p27_ids, p25_ids)
    lost <- setdiff(p25_ids, p27_ids)
    # Ajout des SV acquis
    if (length(gained) > 0) {
      result_list[[idx]] <- data.frame(
        replicate = rep,
        orf = orf_name,
        category = "SV acquis",
        n_sv = length(gained),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }

    # Ajout des SV perdus
    if (length(lost) > 0) {
      result_list[[idx]] <- data.frame(
        replicate = rep,
        orf = orf_name,
        category = "SV perdus",
        n_sv = length(lost),
        stringsAsFactors = FALSE
      )
      idx <- idx + 1
    }
  }
}

# Fusion des résultats
impact <- do.call(rbind, result_list)

# Exclusion des variants intergéniques
impact <- impact[impact$orf != "Intergenic", ]

# Somme des SV par ORF et catégorie
summary_orf <- aggregate(
  n_sv ~ orf + category,
  data = impact,
  FUN = sum
)

# Somme totale par ORF
total_orf <- aggregate(
  n_sv ~ orf,
  data = summary_orf,
  FUN = sum
)

# Sélection des ORFs les plus impactés
top_n <- 15
top_orfs <- head(total_orf[order(-total_orf$n_sv), "orf"], top_n)


plot_data <- summary_orf[summary_orf$orf %in% top_orfs, ]

# Ordre des ORFs pour affichage
orf_levels <- total_orf$orf[order(total_orf$n_sv)]
orf_levels <- orf_levels[orf_levels %in% top_orfs]

# Conversion en facteurs ordonnés
plot_data$orf <- factor(plot_data$orf, levels = orf_levels)
plot_data$category <- factor(plot_data$category, levels = c("SV perdus", "SV acquis"))



p <- ggplot(plot_data, aes(x = n_sv, y = orf, fill = category)) +
  geom_col(width = 0.75) +  
  labs(
    title = "ORFs les plus impactés par le choc thermique",
    subtitle = "Comparaison des SV acquis et perdus entre P25 et P27",
    x = "Nombre de SV acquis/perdus",
    y = "ORF",
    fill = "Catégorie"
  ) +
  scale_fill_manual(values = c("SV perdus" = "#d95f02", "SV acquis" = "#1b9e77")) +
  theme_minimal(base_size = 12)

# Sauvegarde du graphique
ggsave(outfile, p, width = 11, height = 8, dpi = 300)

# Message de confirmation
cat("Sauvegardé:", outfile, "\n")
