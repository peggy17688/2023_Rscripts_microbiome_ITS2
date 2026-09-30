library(dplyr)
library(phyloseq)
library(pairwiseAdonis)
merulina_initial <- subset_samples(merulina_rarefied_ps, macroalgal_treatment_type == "field_sample")
pocillopora_initial <- subset_samples(pocillopora_rarefied_ps, macroalgal_treatment_type == "field_sample")

merulina_initial_sqrt <- transform_sample_counts(merulina_initial, function(x) sqrt(x))
pocillopora_initial_sqrt <- transform_sample_counts(pocillopora_initial, function(x) sqrt(x))

#PERMANOVA/PERMDISP:
merulina_initial_bray_dist <- phyloseq::distance(merulina_initial_sqrt, method = "bray")
meta_merulina_initial_df <- as(sample_data(merulina_initial_sqrt), "data.frame")
merulina_initial_adonis <- adonis2(merulina_initial_bray_dist~location, data = meta_merulina_initial_df, permutations = 999)
merulina_initial_disp <- betadisper(merulina_initial_bray_dist, meta_merulina_initial_df$location)
anova(merulina_initial_disp)

pocillopora_initial_bray_dist <- phyloseq::distance(pocillopora_initial_sqrt, method = "bray")
meta_pocillopora_initial_df <- as(sample_data(pocillopora_initial_sqrt), "data.frame")
pociilopora_initial_adonis <- adonis2(pocillopora_initial_bray_dist~location, data = meta_pocillopora_initial_df, permutations = 999)
pocillopora_initial_disp <- betadisper(pocillopora_initial_bray_dist, meta_pocillopora_initial_df$location)
anova(pocillopora_initial_disp)

#NMDS plots (Bray-Curtis):
library(ggplot2)
set.seed(1)
merulina_initial_nmds <- metaMDS(merulina_initial_bray_dist, k = 2, trymax = 100)
merulina_initial_nmds$stress

nmds_merulina_initial_scores <- as.data.frame(scores(merulina_initial_nmds, display = "sites"))
nmds_merulina_initial_scores$location <- meta_merulina_initial_df$location
ggplot(nmds_merulina_initial_scores, aes(x = NMDS1, y = NMDS2, colour = location)) + geom_point(size = 3) + stat_ellipse(level = 0.95) + theme_minimal() + labs(title = "NMDS of Merulina ampliata initial timepoint microbiome composition by location") + annotate("text", x = min(nmds_merulina_initial_scores$NMDS1), y = max(nmds_merulina_initial_scores$NMDS2), label = paste0("Stress = ", merulina_initial_stress_val), hjust = -3.5, vjust = 1.5, size = 3.5)

pocillopora_initial_nmds <- metaMDS(pocillopora_initial_bray_dist, k = 2, trymax = 100)
pocillopora_initial_nmds$stress

nmds_pocillopora_initial_scores <- as.data.frame(scores(pocillopora_initial_nmds, display = "sites"))
nmds_pocillopora_initial_scores$location <- meta_pocillopora_initial_df$location
pocillopora_initial_stress_val <- round(pocillopora_initial_nmds$stress, 3)
ggplot(nmds_pocillopora_initial_scores, aes(x = NMDS1, y = NMDS2, colour = location)) + geom_point(size = 3) + stat_ellipse(level = 0.95) + theme_minimal() + labs(title = "NMDS of Pocillopora acuta initial timepoint microbiome composition by location") + annotate("text", x = min(nmds_pocillopora_initial_scores$NMDS1), y = max(nmds_pocillopora_initial_scores$NMDS2), label = paste0("Stress = ", pocillopora_initial_stress_val), hjust = -3.5, vjust = 1.5, size = 3.5)