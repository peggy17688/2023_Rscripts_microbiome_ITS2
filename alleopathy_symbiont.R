---
title: "Symbiont clade analysis from SymPortal"
author: "LeeLK"
output: html_document
---

library(tidyverse)
library(patchwork)
library(ggsci)

## Load data
raw <- read_csv("20210112_ITS2_relabun.csv")

## Clean `sample_id`
clean <- raw %>%
  mutate(sample_id = str_remove(sample_id, "_ITS2$")) %>%
  select(1:5)

## ITS2 type profiles
profile_levels <- names(raw)[6:ncol(raw)]

profile_long <- raw %>%
  mutate(sample_id = str_remove(sample_id, "_ITS2$")) %>%
  rename(macroalgal_treatment = `macroalgal_treatment type`) %>%
  pivot_longer(all_of(profile_levels),
               names_to = "profile", values_to = "abundance") %>%
  mutate(
    genus   = if_else(str_starts(profile, "C"), "Cladocopium", "Durusdinium"),
    profile = factor(profile, levels = profile_levels)
  )

## Profile colour palette

profiles_c <- profile_levels[str_starts(profile_levels, "C")]
profiles_d <- profile_levels[str_starts(profile_levels, "D")]

green_ramp  <- colorRampPalette(c("#1B5E20", "#A5D6A7"))  # dark → light green
orange_ramp <- colorRampPalette(c("#8C3A00", "#FFCC80"))  # dark → light orange

profile_cols <- c(
  setNames(green_ramp(length(profiles_c)),  profiles_c),
  setNames(orange_ramp(length(profiles_d)), profiles_d)
)

## Combined profile figure
treatment_order <- c("control", "5 cm", "direct")

make_profile_panel <- function(species_label) {
  profile_long %>%
    filter(sample_species == species_label,
           !macroalgal_treatment %in% c("pre_treatment", "field_sample")) %>%
    mutate(
      macroalgal_treatment = factor(macroalgal_treatment, levels = treatment_order),
      colony = factor(colony, levels = c(paste0("Satumu", 1:5), paste0("Kusu", 1:5)))
    ) %>%
    ggplot(aes(macroalgal_treatment, abundance, fill = profile)) +
    geom_col(position = "fill", width = 0.9) +
    facet_wrap(~ colony, nrow = 2) +
    scale_fill_manual(values = profile_cols, limits = profile_levels, drop = FALSE,
                      guide = guide_legend(ncol = 1)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Macroalgal treatment", y = "Proportional abundance",
         fill = "ITS2 type profile") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      strip.background = element_rect(fill = "grey90")
    )
}

p_acuta_prof    <- make_profile_panel("Pocillopora_acuta") +
  ggtitle(expression(italic("Pocillopora acuta")))
p_ampliata_prof <- make_profile_panel("Merulina_ampliata") +
  ggtitle(expression(italic("Merulina ampliata")))

combined_profile <- (p_acuta_prof / p_ampliata_prof) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  )

combined_profile

ggsave("figure_profile_composition.tiff", combined_profile,
       width = 230, height = 200, units = "mm")
ggsave("Fig6_llk_profile.pdf", combined_profile,
       width = 230, height = 200, units = "mm", dpi = 300, compression = "lzw")

## Load clade-level sequence abundances
clade_abund <- read_tsv(
  "158_280521_DBV_20210528T064539.seqs.absolute.abund_and_meta.txt"
) %>%
  select(sample_id = sample_name, C3:D2a) %>%
  filter(!is.na(sample_id))

## Tidy: join metadata to clade abundances
clade_levels_raw <- clade_abund %>% select(-sample_id) %>% names()
clade_levels     <- str_remove_all(clade_levels_raw, "_")

abund_long <- clade_abund %>%
  pivot_longer(-sample_id, names_to = "clade", values_to = "abundance") %>%
  mutate(
    genus = case_when(
      str_detect(clade, "^C[0-9]|_C$") ~ "Cladocopium",
      str_detect(clade, "^D[0-9]|_D$") ~ "Durusdinium"
    ),
    clade = str_remove_all(clade, "_"),                 # drop "_" from ITS2 name
    clade = factor(clade, levels = clade_levels)
  ) %>%
  left_join(clean, by = "sample_id") %>%
  rename(macroalgal_treatment = `macroalgal_treatment type`)
  
## Colour palette: top-10 clades per genus, rest as "Others"
top_n_per_genus <- 10

# Total sequence abundance per clade across the species/treatments plotted below,
# used to pick the most abundant clades to highlight.
clade_totals <- abund_long %>%
  filter(sample_species %in% c("Pocillopora_acuta", "Merulina_ampliata"),
         !macroalgal_treatment %in% c("pre_treatment", "field_sample")) %>%
  group_by(genus, clade) %>%
  summarise(tot = sum(abundance), .groups = "drop")

top_c <- clade_totals %>% filter(genus == "Cladocopium") %>%
  slice_max(tot, n = top_n_per_genus, with_ties = FALSE) %>% pull(clade) %>% as.character()
top_d <- clade_totals %>% filter(genus == "Durusdinium") %>%
  slice_max(tot, n = top_n_per_genus, with_ties = FALSE) %>% pull(clade) %>% as.character()

# Highlighted clades first (descending-abundance order), "Others" last.
display_levels <- c(top_c, top_d, "Others")

# Collapse every non-highlighted clade into a single "Others" bucket.
abund_long <- abund_long %>%
  mutate(clade_display = factor(
    if_else(as.character(clade) %in% c(top_c, top_d), as.character(clade), "Others"),
    levels = display_levels))

# ggsci qualitative palettes: NPG for Cladocopium, JCO for Durusdinium
# (each supplies up to 10 distinct colours), plus grey for "Others".
clade_cols <- c(
  setNames(pal_npg("nrc")(length(top_c)), top_c),
  setNames(pal_jco()(length(top_d)),      top_d),
  Others = "grey80"
)
## Combined publication figure
treatment_order <- c("control", "5 cm", "direct")

# One panel builder so A and B are guaranteed identical except for the data
make_panel <- function(species_label) {
  abund_long %>%
    filter(sample_species == species_label,
           !macroalgal_treatment %in% c("pre_treatment", "field_sample")) %>%
    mutate(
      macroalgal_treatment = factor(macroalgal_treatment, levels = treatment_order),
      colony = factor(colony, levels = c(paste0("Satumu", 1:5), paste0("Kusu", 1:5)))
    ) %>%
    ggplot(aes(macroalgal_treatment, abundance, fill = clade_display)) +
    geom_col(position = "fill", width = 0.9) +
    facet_wrap(~ colony, nrow = 2) +
    scale_fill_manual(values = clade_cols, limits = display_levels, drop = FALSE,
                      guide = guide_legend(ncol = 2, byrow = FALSE)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(x = "Macroalgal treatment", y = "Proportional abundance",
         fill = "ITS2 sequence") +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.key.size = unit(0.35, "cm"),
      legend.text = element_text(size = 7),
      strip.background = element_rect(fill = "grey90")
    )
}

p_acuta    <- make_panel("Pocillopora_acuta") +
  ggtitle(expression(italic("Pocillopora acuta")))
p_ampliata <- make_panel("Merulina_ampliata") +
  ggtitle(expression(italic("Merulina ampliata")))

combined <- (p_acuta / p_ampliata) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "A") &
  theme(
    legend.position = "right",
    plot.tag = element_text(face = "bold", size = 14)
  )

combined
