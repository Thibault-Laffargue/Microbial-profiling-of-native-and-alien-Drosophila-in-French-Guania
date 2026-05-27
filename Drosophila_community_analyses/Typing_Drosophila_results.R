############## Typing Drosophila #############
##### Initializaion and import of data ######
setwd("~/Desktop/thèse/Microbiome Biodiversity/Bacteria and fungi/Typing Drosophila/All sequence")


install.packages(c("sf", "rnaturalearth", "rnaturalearthdata", "scatterpie"))


library(readxl)
library(writexl)
library(ggplot2)
library(plyr)
library(dplyr)
library(questionr)
library(lme4)
library(lmerTest)
library(lmtest)
library(car)
library(emmeans)
library(tidyverse)
library("gtsummary")
library(ggpubr)
library(rstatix)
library("entropy")
library(tidyr)
library(purrr)
library(Biostrings)
library(ape)
library(msa)
library(phangorn)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(scatterpie)
library(grid)
library(cowplot)
library(patchwork)



tab_brute <- read_excel("Blast_results_filtered.xlsx")
View(tab_brute)

Species <- unique(tab_brute$Species) #permet de créer un vecteur unique species qui contient chaque modalité
Species


tab_brute$Identification_results[tab_brute$Species == "Drosophila tropicalis"] <- "D. tropicalis" #Aucun doute possible avec autre willistoni
tab_brute$Identification_results[tab_brute$Species == "Drosophila nebulosa" & tab_brute$Percentage_Identity > 96] <- "D. nebulosa" #Aucun doute possible avec autre willistoni
tab_brute$Identification_results[tab_brute$Species == "Drosophila neocordata"] <- "D. neocordata" #Aucun doute possible avec autre saltans
tab_brute$Identification_results[tab_brute$Species == "Drosophila equinoxialis"] <- "D. equinoxialis" #Aucun doute possible avec autre willistoni

tab_brute <- tab_brute %>%
  mutate(
    Identification_results = as.character(Identification_results),
    Identification_results = coalesce(Identification_results, "Unidentified")
  )

##### Creat Group columns and Category columns ####

Species <- unique(tab_brute$Identification_results) #permet de créer un vecteur unique species qui contient chaque modalité
Species



# Group columns
tab_brute2 <- tab_brute %>%
  mutate(
    group_identification = case_when(
      Identification_results %in% c("Bipectinata complex sp.", "D. melanogaster", "D. ananassae") ~ "melanogaster",
      Identification_results %in% c("Saltans sp.", "D. neocordata", "D. sturtevanti") ~ "saltans",
      Identification_results %in% c("D. tropicalis", "Willistoni sp.", "D. paulistorum", "D. willistoni", "D. nebulosa","D. equinoxialis") ~ "willistoni",
      Identification_results == "Z. indianus" ~ "zaprionus",
      Identification_results == "Sc. latifasciaeformis" ~ "Scaptodrosophila",
      Identification_results == "Unidentified" ~ "Unidentified",
      TRUE ~ NA_character_
    )
  )
View(tab_brute2)



#Category

tab_final <- tab_brute2 %>%
  mutate(
    Category = case_when(
      group_identification %in% c("melanogaster", "zaprionus", "Scaptodrosophila") ~ "invasive",
      group_identification %in% c("saltans","willistoni") ~ "native",
      group_identification == "Unidentified" ~ "Unidentified",
      TRUE ~ NA_character_           # garde les NA s'il y en a
    )
  )
tab_final <- tab_final %>%
  mutate(Trap = str_extract(Sample_ID, "^[^S]+"))

info_trap <- read_excel("Info_trap.xlsx")
View(info_trap)

tab_final <- tab_final %>%
  select(-Locality) %>%  # Supprime l'ancienne colonne si elle existe déjà (facultatif mais propre)
  left_join(info_trap, by = "Trap")

tab_final <- tab_final %>%
  mutate(
    Identification_results = as.character(Identification_results)
  )

View(tab_final)

write_csv(tab_final, "Identification_results_vf.csv")

##################################### Graphical abstract: Category barplots ###############################################

library(dplyr)
library(ggplot2)

# Préparer les données : Cayenne vs autres localités
ga_category_data <- tab_final %>%
  filter(Category %in% c("invasive", "native")) %>%
  mutate(
    Locality_group = ifelse(Locality == "Cayenne", "Cayenne", "Other localities"),
    Category = factor(Category, levels = c("invasive", "native"))
  ) %>%
  count(Locality_group, Category, name = "Nbr_ind")

View(ga_category_data)

# Couleurs simples pour graphical abstract
category_colors <- c(
  invasive = "#C1440E",
  native = "#22427C"
)

# Barplot combiné avec deux panneaux
p_ga_category <- ggplot(
  ga_category_data,
  aes(x = Category, y = Nbr_ind, fill = Category)
) +
  geom_col(width = 0.65, color = "black") +
  geom_text(
    aes(label = Nbr_ind),
    vjust = -0.4,
    size = 6,
    fontface = "bold"
  ) +
  facet_wrap(~ Locality_group, nrow = 1) +
  scale_fill_manual(values = category_colors) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of individuals") +
  theme_minimal(base_size = 18) +
  theme(
    strip.text = element_text(size = 20, face = "bold"),
    axis.text.x = element_text(size = 16, face = "bold"),
    axis.text.y = element_text(size = 14, face = "bold"),
    axis.title.y = element_text(size = 18, face = "bold"),
    legend.position = "none",
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

p_ga_category

ggsave(
  "Graphical_abstract_category_barplots_Cayenne_vs_other_localities.pdf",
  p_ga_category,
  width = 8,
  height = 4
)


################## Graphiques préliminaire ###########


# Données de base
Species_by_Locality <- tab_final %>%
  dplyr::count(Identification_results, Locality, name = "Nbr_ind") %>%
  dplyr::rename(Sp = Identification_results, Loc = Locality)

# Définir l'ordre des espèces par catégorie
species_category <- tab_final %>%
  select(Identification_results, Category) %>%
  distinct()

invasive_species <- species_category %>%
  filter(Category == "invasive") %>%
  pull(Identification_results)

native_species <- species_category %>%
  filter(Category == "native") %>%
  pull(Identification_results)

unidentified_species <- species_category %>%
  filter(Category == "Unidentified") %>%
  pull(Identification_results)

Species_by_Locality$Sp <- factor(
  Species_by_Locality$Sp,
  levels = c(invasive_species, native_species, unidentified_species)
)

# Couleurs
colors <- setNames(
  c(
    colorRampPalette(c("#7B241C", "#C1440E", "#FFA500"))(length(invasive_species)),
    colorRampPalette(c("#22427C", "#1f9e89", "#228B22"))(length(native_species)),
    colorRampPalette("grey")(length(unidentified_species))
  ),
  c(invasive_species, native_species, unidentified_species)
)

# Proportions
Species_by_Locality <- Species_by_Locality %>%
  group_by(Loc) %>%
  mutate(Percent = Nbr_ind / sum(Nbr_ind)) %>%
  ungroup()

#passage en format large
pie_data <- Species_by_Locality %>%
  select(Loc, Sp, Percent) %>%
  tidyr::pivot_wider(names_from = Sp, values_from = Percent, values_fill = 0)

#Ajout des coordonnées GPS
coords <- tibble::tibble(
  Loc = c("Cayenne", "Kaw", "Bélizon", "Nouragues"),
  lon = c(-52.3333, -52.0344, -52.3817, -52.6872),
  lat = c(4.9333, 4.4776, 4.3083, 4.0788),
  x = c(-52.5, -51.5, -52.9, -54),  # positions décalées pour dessiner les camemberts
  y = c(5.5, 4.5, 2.5, 4.5)
)

pie_data <- left_join(pie_data, coords, by = "Loc")

# Calcul du total d'individus par localité et ajout
total_ind <- Species_by_Locality %>%
  group_by(Loc) %>%
  dplyr::summarise(Total = sum(Nbr_ind))
pie_data <- dplyr::left_join(pie_data, total_ind, by = "Loc")

# réordonne les colonnes selon l'ordre que tu veux
ordered_species <- c(invasive_species, native_species, unidentified_species)
pie_data <- pie_data[, c("Loc", ordered_species, "lon", "lat", "x", "y", "Total")]

# Crée des expressions de type expression(italic("...")) pour chaque espèce
species_labels <- setNames(
  lapply(ordered_species, function(sp) bquote(italic(.(sp)))),
  ordered_species
)

#Graph 
ggplot() +
  borders("world", fill = "gray90", colour = "black") +
  geom_segment( # flèches entre coord GPS (lon/lat) et camembert (x/y)
    data = pie_data,
    aes(x = lon, y = lat, xend = x, yend = y),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "gray30"
  ) +
  geom_point(#points au localité
    data = pie_data,
    aes(x = lon, y = lat),
    color = "red",
    size = 2
  ) +
  geom_text( # Nom des localités
    data = pie_data,
    aes(x = lon, y = lat, label = Loc),
    hjust = 1.1, vjust = 0.5,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  geom_scatterpie(
    data = pie_data,
    aes(x = x, y = y),
    cols = ordered_species,
    pie_scale = 10,
    legend_name = "Species"
  )+
  geom_label(
    data = pie_data,
    aes(x = x, y = y, label = Total),
    size = 4,
    fontface = "bold",
    fill = "gray90",    # couleur de fond
    color = "black"     # couleur du texte
  ) +
  scale_fill_manual(
    values = colors,
    labels = species_labels
  ) +
  coord_quickmap(xlim = c(-55, -51), ylim = c(2, 6), expand = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 15),       # taille du texte de la légende
    legend.title = element_text(size = 15, face = "bold")  # taille et style du titre
  ) +
  labs(title = "A - Distribution map of Drosophila species in French Guiana", fill = "Species", x = "Longitude", y = "Latitude")

ggsave(filename = "camemberts_par_localite.pdf", width = 10, height = 8)


# Fonds de carte seule
ggplot() +
  borders("world", fill = "gray90", colour = "black") +
  geom_point( # points aux localités
    data = pie_data,
    aes(x = lon, y = lat),
    color = "red",
    size = 2
  ) +
  geom_text( # noms des localités
    data = pie_data,
    aes(x = lon, y = lat, label = Loc),
    hjust = 1.1, vjust = 0.5,
    size = 6,
    fontface = "bold",
    color = "black"
  ) +
  coord_quickmap(xlim = c(-55, -51), ylim = c(2, 6), expand = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.position = "none"  # pas besoin de légende si pas de piecharts
  ) +
  labs(
    title = "A - Map of sampling localities in French Guiana",
    x = "Longitude",
    y = "Latitude"
  )

ggsave("fond_carte_sans_camemberts.pdf", width = 10, height = 8)


##### Version histo

# 1. Préparer les données
bar_data <- Species_by_Locality %>%
  left_join(coords, by = c("Loc")) %>%
  mutate(
    Sp = factor(Sp, levels = ordered_species),
    Loc = factor(Loc, levels = c("Cayenne", "Kaw", "Bélizon", "Nouragues"))
  )

# 2. Créer les mini-barplots
plots_list <- list()

for (loc in levels(bar_data$Loc)) {
  df_loc <- bar_data %>% filter(Loc == loc)
  total_ind <- sum(df_loc$Nbr_ind)
  
  p <- ggplot(df_loc, aes(x = Sp, y = Nbr_ind, fill = Sp)) +
    geom_col(width = 0.7, color = "black") +
    scale_fill_manual(values = colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
    labs(title = loc) +
    annotate(
      "label",
      x = length(unique(df_loc$Sp)) / 2 + 0.5,
      y = max(df_loc$Nbr_ind) * 0.9,
      label = paste0("Total: ", total_ind),
      size = 3.5,
      fontface = "bold",
      fill = "gray90",
      color = "black",
      label.size = 0.3
    ) +
    theme_minimal(base_size = 10) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 11, face = "bold"),
      legend.position = "none",
      plot.margin = margin(0, 0, 0, 0),
      # >>> Tick labels en gras sur l'axe Y
      axis.text.y = element_text(face = "bold")
    )
  
  plots_list[[loc]] <- ggplotGrob(p)
}

# 3. Carte de base
base_map <- ggplot() +
  borders("world", fill = "gray90", colour = "black") +
  geom_point(data = coords, aes(x = lon, y = lat), color = "red", size = 4) +
  geom_text(data = coords, aes(x = lon, y = lat, label = Loc),
            hjust = 1.1, vjust = 0.5, size = 6, fontface = "bold") +
  coord_quickmap(xlim = c(-55, -51), ylim = c(2, 6), expand = FALSE) +
  labs(x = "Longitude", y = "Latitude") +
  theme_minimal(base_size = 10) +  
  theme(
    axis.title = element_text(size = 12, face = "bold"),
    panel.grid = element_blank(),
    # >>> Tick labels en gras sur la carte (x et y)
    axis.text = element_text(size = 10, face = "bold")
  )

# 4. Positions verticales réparties dynamiquement
n_locs <- length(levels(bar_data$Loc))
heights <- rep(1 / n_locs, n_locs)
y_positions <- rev(cumsum(heights)) - heights
names(y_positions) <- levels(bar_data$Loc)

# 5. Légende indépendante
df_legend <- data.frame(Sp = factor(names(colors), levels = ordered_species))

legend_plot <- ggplot(df_legend, aes(x = Sp, fill = Sp)) +
  geom_bar() +
  scale_fill_manual(
    values = colors,
    labels = species_labels,
    name = "Species"
  ) +
  guides(fill = guide_legend(ncol = 3, title.position = "top")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.title = element_text(size = 13, face = "bold", hjust = 0.5),
    legend.text = element_text(size = 11)
  )

legend_grob <- cowplot::get_legend(legend_plot)

# 6. Colonne de gauche = carte + légende
left_column <- patchwork::wrap_elements(base_map) / patchwork::wrap_elements(legend_grob) +
  plot_layout(heights = c(0.7, 0.3))

# 7. Colonne de droite = barplots
barplots_stack <- cowplot::ggdraw()
for (loc in levels(bar_data$Loc)) {
  barplots_stack <- barplots_stack +
    draw_grob(plots_list[[loc]],
              x = 0, y = y_positions[loc],
              width = 1, height = 1 / n_locs)
}

# 8. Colonne flèche (grande flèche + texte)
arrow_plot <- ggplot() +
  geom_segment(aes(x = 0.5, xend = 0.5, y = 0.15, yend = 0.85),
               arrow = arrow(type = "closed", length = unit(0.4, "cm")),
               size = 3.5, color = "black") +
  annotate("text", x = 0.7, y = 0.5, label = "Anthropization gradient",
           angle = 90, size = 5, fontface = "bold", hjust = 0.5) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(0, 0, 0, 0),
    plot.background = element_rect(fill = "transparent", color = NA)
  )

# 9. Assembler les colonnes (carte+legende | barplots | flèche)
final_plot <- cowplot::plot_grid(
  cowplot::ggdraw(left_column),
  cowplot::ggdraw(barplots_stack),
  arrow_plot,
  ncol = 3,
  rel_widths = c(0.5, 0.43, 0.07),
  align = "h",
  axis = "tb"
)

# >>> 10. SUPPRIMÉ — pas de traits reliant carte et barplots

# 11. Titre global
final_with_title <- cowplot::ggdraw() +
  cowplot::draw_label(
    "B- Distribution map of Drosophila species in French Guiana",
    fontface = "bold",
    size = 16,
    x = 0.3,
    y = 0.98,
    hjust = 0.5
  ) +
  cowplot::draw_plot(final_plot, x = 0, y = 0, width = 1, height = 0.95)

# 12. Sauvegarde
ggsave("barplots_with_legend_below_map.pdf", plot = final_with_title, width = 12, height = 8)


##################################### Group percentage ###############################################

# Calcul du pourcentage d'individus par espèce au sein de chaque groupe taxonomique
species_percentage_group <- tab_final %>%
  group_by(group_identification, Identification_results) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(group_identification) %>%
  mutate(
    total = sum(n),
    percent = round((n / total) * 100, 2)
  ) %>%
  select(group_identification, Identification_results, percent)

View(species_percentage_group)

species_percentage_category <- tab_final %>%
  group_by(Category, Identification_results) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Category) %>%
  mutate(
    total = sum(n),
    percent = round((n / total) * 100, 2)
  ) %>%
  select(Category, Identification_results, percent)

View(species_percentage_category)

##################################### Statistics ######################################################


#### Calcul du nombre total d'individus et du nombre d'invasives par localité ###
prop_invasive <- tab_final %>%
  group_by(Locality) %>%
  summarise(
    total = n(),
    invasive = sum(Category == "invasive"),
    proportion_invasive = invasive / total
  )

View(prop_invasive)


### Global Chi² test : l'association entre Locality et Category est-elle significative ?

# Table de contingence

contingency_table <- tab_final %>%
  filter(Category %in% c("invasive", "native")) %>%  # on exclut les "Unidentified"
  dplyr::count(Locality, Category) %>%
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) %>%
  column_to_rownames("Locality") %>%
  as.matrix()

# Test du chi²
chisq_test <- chisq.test(contingency_table)

# Résultat du test
print(chisq_test)



#### exact fisher test pairwise ###

# Table résumée avec nb invasives et natives par localité
summary_tab <- tab_final %>%
  mutate(group = ifelse(Category == "invasive", "invasive", "native")) %>%
  dplyr::count(Locality, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

# Création de toutes les paires de localités
localities <- summary_tab$Locality
pairs <- combn(localities, 2, simplify = FALSE)

# Test de Fisher pour chaque paire
results <- map_df(pairs, function(pair) {
  # sous-table pour les 2 localités
  sub <- summary_tab %>% filter(Locality %in% pair)
  
  # table 2x2 pour fisher.test
  mat <- matrix(
    c(sub$invasive, sub$native),
    nrow = 2,
    byrow = FALSE,
    dimnames = list(Locality = sub$Locality, Group = c("invasive", "native"))
  )
  
  test <- fisher.test(mat)
  
  tibble(
    Locality1 = pair[1],
    Locality2 = pair[2],
    p_value = test$p.value
  )
})

# Optionnel : correction pour comparaisons multiples
results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # ou "bonferroni"
View(results)
