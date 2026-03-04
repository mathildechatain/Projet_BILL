# 1. Charger les données
variants <- read.table("../Files/plot_files/count_variants_sv.tsv", 
                       header = TRUE, sep = "\t", stringsAsFactors = FALSE)

# Vérifier les données
print(variants)

# 2. Séparer "File" en Generation et Condition
variants$Generation <- sub("_(chaud|froid)$", "", variants$File)
variants$Condition <- sub("^P[0-9]+_", "", variants$File)

# 3. Créer le graphique avec ggplot2
library(ggplot2)

p <- ggplot(variants, aes(x = Generation, y = N_variants, fill = Condition)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_minimal() +
  labs(title = "Nombre de variants par génération et condition",
       x = "Génération", y = "Nombre de variants") +
  scale_fill_manual(values = c("chaud" = "red", "froid" = "blue"))

# 4. Afficher le graphique
print(p)

# 5. Sauvegarder le graphique
ggsave("../Plots/count_variants_.png", plot = p, width = 6, height = 4, dpi = 300)
