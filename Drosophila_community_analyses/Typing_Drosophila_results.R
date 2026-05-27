############## Typing Drosophila #############
##### Initialization and import of data ######
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

Species <- unique(tab_brute$Species) #create a vector of unique species with each variable
Species


tab_brute$Identification_results[tab_brute$Species == "Drosophila tropicalis"] <- "D. tropicalis" #can't be other willistoni
tab_brute$Identification_results[tab_brute$Species == "Drosophila nebulosa" & tab_brute$Percentage_Identity > 96] <- "D. nebulosa" #can't be other willistoni
tab_brute$Identification_results[tab_brute$Species == "Drosophila neocordata"] <- "D. neocordata" #can't be other saltans
tab_brute$Identification_results[tab_brute$Species == "Drosophila equinoxialis"] <- "D. equinoxialis" #can't be other willistoni

tab_brute <- tab_brute %>%
  mutate(
    Identification_results = as.character(Identification_results),
    Identification_results = coalesce(Identification_results, "Unidentified")
  )

##### Create Group columns and Category columns ####

Species <- unique(tab_brute$Identification_results) #create a vector of unique species with each variable
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
      TRUE ~ NA_character_           # Keep NA
    )
  )
tab_final <- tab_final %>%
  mutate(Trap = str_extract(Sample_ID, "^[^S]+"))

info_trap <- read_excel("Info_trap.xlsx")
View(info_trap)

tab_final <- tab_final %>%
  select(-Locality) %>%  # delete old col if already existing (optional but cleaner)
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

# prepare data: Cayenne vs others localities
ga_category_data <- tab_final %>%
  filter(Category %in% c("invasive", "native")) %>%
  mutate(
    Locality_group = ifelse(Locality == "Cayenne", "Cayenne", "Other localities"),
    Category = factor(Category, levels = c("invasive", "native"))
  ) %>%
  count(Locality_group, Category, name = "Nbr_ind")

View(ga_category_data)

# simple color for graphical abstract
category_colors <- c(
  invasive = "#C1440E",
  native = "#22427C"
)

# Barplot with combined panel
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


################## preliminary graphics ###########


# base data
Species_by_Locality <- tab_final %>%
  dplyr::count(Identification_results, Locality, name = "Nbr_ind") %>%
  dplyr::rename(Sp = Identification_results, Loc = Locality)

# define species order by category
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

# Colors
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

#large format
pie_data <- Species_by_Locality %>%
  select(Loc, Sp, Percent) %>%
  tidyr::pivot_wider(names_from = Sp, values_from = Percent, values_fill = 0)

#GPS data
coords <- tibble::tibble(
  Loc = c("Cayenne", "Kaw", "Bélizon", "Nouragues"),
  lon = c(-52.3333, -52.0344, -52.3817, -52.6872),
  lat = c(4.9333, 4.4776, 4.3083, 4.0788),
  x = c(-52.5, -51.5, -52.9, -54),  # positions of pie charts
  y = c(5.5, 4.5, 2.5, 4.5)
)

pie_data <- left_join(pie_data, coords, by = "Loc")

# number of individual by locality
total_ind <- Species_by_Locality %>%
  group_by(Loc) %>%
  dplyr::summarise(Total = sum(Nbr_ind))
pie_data <- dplyr::left_join(pie_data, total_ind, by = "Loc")

# Order col
ordered_species <- c(invasive_species, native_species, unidentified_species)
pie_data <- pie_data[, c("Loc", ordered_species, "lon", "lat", "x", "y", "Total")]

# Italic species name
species_labels <- setNames(
  lapply(ordered_species, function(sp) bquote(italic(.(sp)))),
  ordered_species
)

#Graph 
ggplot() +
  borders("world", fill = "gray90", colour = "black") +
  geom_segment( # arrow
    data = pie_data,
    aes(x = lon, y = lat, xend = x, yend = y),
    arrow = arrow(length = unit(0.15, "cm")),
    color = "gray30"
  ) +
  geom_point(# locality designed by point
    data = pie_data,
    aes(x = lon, y = lat),
    color = "red",
    size = 2
  ) +
  geom_text( # locality name
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
    fill = "gray90",    # background color
    color = "black"     # text color
  ) +
  scale_fill_manual(
    values = colors,
    labels = species_labels
  ) +
  coord_quickmap(xlim = c(-55, -51), ylim = c(2, 6), expand = FALSE) +
  theme_minimal() +
  theme(
    plot.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 15),       # legend text size
    legend.title = element_text(size = 15, face = "bold")  # Title size and style 
  ) +
  labs(title = "A - Distribution map of Drosophila species in French Guiana", fill = "Species", x = "Longitude", y = "Latitude")

ggsave(filename = "camemberts_par_localite.pdf", width = 10, height = 8)


# Map alone
ggplot() +
  borders("world", fill = "gray90", colour = "black") +
  geom_point( # Locality point
    data = pie_data,
    aes(x = lon, y = lat),
    color = "red",
    size = 2
  ) +
  geom_text( # Locality name
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
    legend.position = "none"  # no need for legend as we don't have the pie chart
  ) +
  labs(
    title = "A - Map of sampling localities in French Guiana",
    x = "Longitude",
    y = "Latitude"
  )

ggsave("fond_carte_sans_camemberts.pdf", width = 10, height = 8)


##### Version histogram

# 1. Prepare data
bar_data <- Species_by_Locality %>%
  left_join(coords, by = c("Loc")) %>%
  mutate(
    Sp = factor(Sp, levels = ordered_species),
    Loc = factor(Loc, levels = c("Cayenne", "Kaw", "Bélizon", "Nouragues"))
  )

# 2. create mini-barplot
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
      # >>> Tick labels in bold
      axis.text.y = element_text(face = "bold")
    )
  
  plots_list[[loc]] <- ggplotGrob(p)
}

# 3. MAP
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
    # >>> Tick labels in bold
    axis.text = element_text(size = 10, face = "bold")
  )

# 4. Dynamic vertical positions
n_locs <- length(levels(bar_data$Loc))
heights <- rep(1 / n_locs, n_locs)
y_positions <- rev(cumsum(heights)) - heights
names(y_positions) <- levels(bar_data$Loc)

# 5. legend
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

# 6. merging part
left_column <- patchwork::wrap_elements(base_map) / patchwork::wrap_elements(legend_grob) +
  plot_layout(heights = c(0.7, 0.3))

# 7. right col = barplots
barplots_stack <- cowplot::ggdraw()
for (loc in levels(bar_data$Loc)) {
  barplots_stack <- barplots_stack +
    draw_grob(plots_list[[loc]],
              x = 0, y = y_positions[loc],
              width = 1, height = 1 / n_locs)
}

# 8. arrow col
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

# 9. merge col again
final_plot <- cowplot::plot_grid(
  cowplot::ggdraw(left_column),
  cowplot::ggdraw(barplots_stack),
  arrow_plot,
  ncol = 3,
  rel_widths = c(0.5, 0.43, 0.07),
  align = "h",
  axis = "tb"
)

# >>> 10. No connecting lines between map and barplots

# 11. Global title
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

# 12. Save
ggsave("barplots_with_legend_below_map.pdf", plot = final_with_title, width = 12, height = 8)


##################################### Group percentage ###############################################

# Percentage of individuals of each species inside taxonomic species group
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


#### Total number of individuals and number of invasive individuals by locality ###
prop_invasive <- tab_final %>%
  group_by(Locality) %>%
  summarise(
    total = n(),
    invasive = sum(Category == "invasive"),
    proportion_invasive = invasive / total
  )

View(prop_invasive)


### Global Chi² test : does association between Locality and Category is significant?

# contingency table

contingency_table <- tab_final %>%
  filter(Category %in% c("invasive", "native")) %>%  # Unidentified excluded
  dplyr::count(Locality, Category) %>%
  pivot_wider(names_from = Category, values_from = n, values_fill = 0) %>%
  column_to_rownames("Locality") %>%
  as.matrix()

# Chi² test
chisq_test <- chisq.test(contingency_table)

# test Results
print(chisq_test)


#### exact fisher test pairwise ###

# resume table with number of invasive and native individuals by locality
summary_tab <- tab_final %>%
  mutate(group = ifelse(Category == "invasive", "invasive", "native")) %>%
  dplyr::count(Locality, group) %>%
  pivot_wider(names_from = group, values_from = n, values_fill = 0)

# create all locality pairwise
localities <- summary_tab$Locality
pairs <- combn(localities, 2, simplify = FALSE)

# Fisher test for each pair
results <- map_df(pairs, function(pair) {
  sub <- summary_tab %>% filter(Locality %in% pair)
  
  # 2x2 table for fisher.test
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

# Optional: multiple comparisons correction
results <- results %>%
  mutate(p_adj = p.adjust(p_value, method = "BH"))  # or "bonferroni"
View(results)
