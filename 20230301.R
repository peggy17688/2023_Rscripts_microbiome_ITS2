library(phyloseq)
library(vegan)
library(dplyr)
library(ggplot2)
library(DESeq2)
library(data.table)
library(forcats)
library(viridis)
library(microbiome)
library(tidyverse)
library(scales)
library(gridExtra)
library(ggpubr)
library(RColorBrewer)

#Plot relative abundances for P. acuta.
get_taxa_unique(raffles_pocillopora_rarefied_ps, "Phylum")
get_taxa_unique(kusu_pocillopora_rarefied_ps, "Phylum")
get_taxa_unique(raffles_pocillopora_rarefied_ps, "Genus")
get_taxa_unique(kusu_pocillopora_rarefied_ps, "Genus")

sample_data(raffles_pocillopora_rarefied_ps)$macroalgal_treatment_type <- factor(sample_data(raffles_pocillopora_rarefied_ps)$macroalgal_treatment_type)
sample_data(kusu_pocillopora_rarefied_ps)$macroalgal_treatment_type <- factor(sample_data(kusu_pocillopora_rarefied_ps)$macroalgal_treatment_type)

phylum_raffles_pocillopora <- raffles_pocillopora_rarefied_ps %>% aggregate_taxa(level = "Phylum") %>% microbiome::transform(transform = "compositional")
phylum_kusu_pocillopora <- kusu_pocillopora_rarefied_ps %>% aggregate_taxa(level = "Phylum") %>% microbiome::transform(transform = "compositional")
phylum_raffles_pocillopora
phylum_kusu_pocillopora
head(otu_table(phylum_raffles_pocillopora))
head(otu_table(phylum_kusu_pocillopora))

phylum_raffles_pocillopora_major <- raffles_pocillopora_rarefied_ps %>% aggregate_rare(level = "Phylum", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
phylum_kusu_pocillopora_major <- kusu_pocillopora_rarefied_ps %>% aggregate_rare(level = "Phylum", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
head(otu_table(phylum_raffles_pocillopora_major))
head(otu_table(phylum_kusu_pocillopora_major))
head(tax_table(phylum_raffles_pocillopora_major))
head(tax_table(phylum_kusu_pocillopora_major))

plot_raffles_pocillopora_phylum <- phylum_raffles_pocillopora_major %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_raffles_pocillopora_phylum <- plot_raffles_pocillopora_phylum + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_raffles_pocillopora_phylum <- plot_raffles_pocillopora_phylum + labs(title = "Pulau Satumu") + scale_fill_discrete(name = "Phylum")
plot_raffles_pocillopora_phylum <- plot_raffles_pocillopora_phylum + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_r_pocillopora <- plot_raffles_pocillopora_phylum + scale_fill_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Other" = "#CCCCCC", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))

plot_kusu_pocillopora_phylum <- phylum_kusu_pocillopora_major %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_kusu_pocillopora_phylum <- plot_kusu_pocillopora_phylum + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_kusu_pocillopora_phylum <- plot_kusu_pocillopora_phylum + labs(title = "Kusu Island") + scale_fill_discrete(name = "Phylum")
plot_kusu_pocillopora_phylum <- plot_kusu_pocillopora_phylum + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_k_pocillopora <- plot_kusu_pocillopora_phylum + scale_fill_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Other" = "#CCCCCC", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))

#Combine both plots.
legend_phylum <- cowplot::get_legend(plot_r_pocillopora)
grid.newpage()
grid.draw(legend_phylum)

plot_r_pocillopora_1 <- plot_r_pocillopora + theme(legend.position = "none")
plot_k_pocillopora_1 <- plot_k_pocillopora + theme(legend.position = "none")

relabun_pocillopora_20210112 <- ggarrange(plot_r_pocillopora_1, plot_k_pocillopora_1, legend_phylum, labels = c("A", "B"), nrow = 1)

#convert absolute to relative abundance.
#Pocillopora_Raffles:
pa_r <- tax_glom(raffles_pocillopora_rarefied_ps, "Phylum")
pa_r0 <- transform_sample_counts(pa_r, function(x)x/sum(x))
pa_r1 <- merge_samples(pa_r0, "macroalgal_treatment_type")
pa_r2 <- transform_sample_counts(pa_r1, function(x)x/sum(x))

#Pocillopora Kusu:
pa_k <- tax_glom(kusu_pocillopora_rarefied_ps, "Phylum")
pa_k0 <- transform_sample_counts(pa_k, function(x)x/sum(x))
pa_k1 <- merge_samples(pa_k0, "macroalgal_treatment_type")
pa_k2 <- transform_sample_counts(pa_k1, function(x)x/sum(x))

#Plot relative abundance based on Genus:
#Pocillopora:
genus_poc_raffles_genus_top <- raffles_pocillopora_rarefied_ps %>% aggregate_rare(level = "Genus", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
genus_poc_raffles_genus_top
head(otu_table(genus_poc_raffles_genus_top))
head(tax_table(genus_poc_raffles_genus_top))

genus_poc_kusu_genus_top <- kusu_pocillopora_rarefied_ps %>% aggregate_rare(level = "Genus", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
genus_poc_kusu_genus_top
head(otu_table(genus_poc_kusu_genus_top))
head(tax_table(genus_poc_kusu_genus_top))

nb.colors <- 55
mycolors <- colorRampPalette(brewer.pal(12, "Paired"))(nb.colors)

plot_poc_raffles_genus_top <- genus_poc_raffles_genus_top %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_poc_raffles_genus_top <- plot_poc_raffles_genus_top + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_poc_raffles_genus_top <- plot_poc_raffles_genus_top + labs(title = "Pulau Satumu") + scale_fill_discrete(name = "Genus") + scale_fill_manual(values = mycolors)
plot_poc_raffles_genus_top <- plot_poc_raffles_genus_top + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_poc_raffles_genus_top <- plot_poc_raffles_genus_top + theme(legend.position = "bottom")

plot_poc_kusu_genus_top <- genus_poc_kusu_genus_top %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_poc_kusu_genus_top <- plot_poc_kusu_genus_top + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_poc_kusu_genus_top <- plot_poc_kusu_genus_top + labs(title = "Kusu Island") + scale_fill_discrete(name = "Genus") + scale_fill_manual(values = mycolors)
plot_poc_kusu_genus_top <- plot_poc_kusu_genus_top + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_poc_kusu_genus_top <- plot_poc_kusu_genus_top + theme(legend.position = "bottom")

relabun_poc_genus <- ggarrange(plot_poc_raffles_genus_top, plot_poc_kusu_genus_top, labels = c("A", "B"), ncol = 1)

#Merulina:
genus_mer_raffles_genus_top <- raffles_merulina_rarefied_ps %>% aggregate_rare(level = "Genus", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
genus_mer_raffles_genus_top
head(otu_table(genus_mer_raffles_genus_top))
head(tax_table(genus_mer_raffles_genus_top))

genus_mer_kusu_genus_top <- kusu_merulina_rarefied_ps %>% aggregate_rare(level = "Genus", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
genus_mer_kusu_genus_top
head(otu_table(genus_mer_kusu_genus_top))
head(tax_table(genus_mer_kusu_genus_top))

plot_mer_raffles_genus_top <- genus_mer_raffles_genus_top %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_mer_raffles_genus_top <- plot_mer_raffles_genus_top + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_mer_raffles_genus_top <- plot_mer_raffles_genus_top + labs(title = "Pulau Satumu") + scale_fill_discrete(name = "Genus") + scale_fill_manual(values = mycolors)
plot_mer_raffles_genus_top <- plot_mer_raffles_genus_top + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_mer_raffles_genus_top <- plot_mer_raffles_genus_top + theme(legend.position = "bottom")

plot_mer_kusu_genus_top <- genus_mer_kusu_genus_top %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_mer_kusu_genus_top <- plot_mer_kusu_genus_top + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_mer_kusu_genus_top <- plot_mer_kusu_genus_top + labs(title = "Kusu Island") + scale_fill_discrete(name = "Genus") + scale_fill_manual(values = mycolors)
plot_mer_kusu_genus_top <- plot_mer_kusu_genus_top + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_mer_kusu_genus_top <- plot_mer_kusu_genus_top + theme(legend.position = "bottom")

relabun_mer_genus <- ggarrange(plot_mer_raffles_genus_top, plot_mer_kusu_genus_top, labels = c("A", "B"), ncol = 1)

poc_raffles_genus_top <- as.data.frame(otu_table(genus_poc_raffles_genus_top))
poc_kusu_genus_top <- as.data.frame(otu_table(genus_poc_kusu_genus_top))
mer_raffles_genus_top <- as.data.frame(otu_table(genus_mer_raffles_genus_top))
mer_kusu_genus_top <- as.data.frame(otu_table(genus_mer_kusu_genus_top))
write.csv(poc_raffles_genus_top, "20210112_poc_raffles_genus_top.csv", row.names = TRUE)
write.csv(poc_kusu_genus_top, "20210112_poc_kusu_genus_top.csv", row.names = TRUE)
write.csv(mer_raffles_genus_top, "20210112_mer_raffles_genus_top.csv", row.names = TRUE)
write.csv(mer_kusu_genus_top, "20210112_mer_kusu_genus_top.csv", row.names = TRUE)

##########

#Plot P. acuta richness plots
richness_pocillopora <- plot_richness(pocillopora_rarefied_ps, x = "location", color = "macroalgal_treatment_type", measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) + geom_boxplot()
richness_pocillopora <- richness_pocillopora + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
richness_pocillopora <- richness_pocillopora + theme(strip.background = element_blank())
richness_pocillopora <- richness_pocillopora + labs(x="", y="Alpha Diversity Index\n")
richness_pocillopora <- richness_pocillopora + theme(axis.text.x = element_text(angle = 45, hjust = 1))
richness_pocillopora <- richness_pocillopora + scale_x_discrete(limits=c("Pulau_Satumu", "Kusu_Island"), labels=c("Pulau Satumu", "Kusu Island"))
richness_pocillopora <- richness_pocillopora + theme(legend.title = element_blank())

p_poci <- richness_pocillopora
newTreatment <- c("field_sample", "pre_experiment", "control", "direct", "5cm")
p_poci$data$macroalgal_treatment_type <- as.character(p_poci$data$macroalgal_treatment_type)
p_poci$data$macroalgal_treatment_type <- factor(p_poci$data$macroalgal_treatment_type, levels = newTreatment)
print(p_poci)

p_poci <- p_poci + scale_color_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(legend.position = "bottom")

##########

#Plot relative abundances for M. ampliata.

get_taxa_unique(raffles_merulina_rarefied_ps, "Phylum")
get_taxa_unique(kusu_merulina_rarefied_ps, "Phylum")
get_taxa_unique(raffles_pocillopora_rarefied_ps, "Genus")
get_taxa_unique(kusu_merulina_rarefied_ps, "Genus")

sample_data(raffles_merulina_rarefied_ps)$macroalgal_treatment_type <- factor(sample_data(raffles_merulina_rarefied_ps)$macroalgal_treatment_type)
sample_data(kusu_merulina_rarefied_ps)$macroalgal_treatment_type <- factor(sample_data(kusu_merulina_rarefied_ps)$macroalgal_treatment_type)

phylum_raffles_merulina <- raffles_merulina_rarefied_ps %>% aggregate_taxa(level = "Phylum") %>% microbiome::transform(transform = "compositional")
phylum_kusu_merulina <- kusu_merulina_rarefied_ps %>% aggregate_taxa(level = "Phylum") %>% microbiome::transform(transform = "compositional")
head(otu_table(phylum_raffles_merulina))
head(otu_table(phylum_kusu_merulina))
head(tax_table(phylum_raffles_merulina))
head(tax_table(phylum_kusu_merulina))

phylum_raffles_merulina_major <- raffles_merulina_rarefied_ps %>% aggregate_rare(level = "Phylum", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
phylum_kusu_merulina_major <- kusu_merulina_rarefied_ps %>% aggregate_rare(level = "Phylum", detection = 1/100, prevalence = 50/100) %>% microbiome::transform(transform = "compositional")
head(otu_table(phylum_raffles_merulina_major))
head(otu_table(phylum_kusu_merulina_major))
head(tax_table(phylum_raffles_merulina_major))
head(tax_table(phylum_kusu_merulina_major))

plot_raffles_merulina_phylum <- phylum_raffles_merulina_major %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_raffles_merulina_phylum <- plot_raffles_merulina_phylum + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_raffles_merulina_phylum <- plot_raffles_merulina_phylum + labs(title = "Pulau Satumu") + scale_fill_discrete(name = "Phylum")
plot_raffles_merulina_phylum <- plot_raffles_merulina_phylum + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_r_merulina <- plot_raffles_merulina_phylum + scale_fill_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Other" = "#CCCCCC", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))

plot_kusu_merulina_phylum <- phylum_kusu_merulina_major %>% plot_composition(average_by = "macroalgal_treatment_type") + scale_y_continuous(labels = waiver())
plot_kusu_merulina_phylum <- plot_kusu_merulina_phylum + labs(x="", y = "Relative Abundance\n") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black"))
plot_kusu_merulina_phylum <- plot_kusu_merulina_phylum + labs(title = "Kusu Island") + scale_fill_discrete(name = "Phylum")
plot_kusu_merulina_phylum <- plot_kusu_merulina_phylum + scale_x_discrete(limits=c("field_sample", "pre_experiment", "control", "direct", "5cm"), labels=c("Field \n sample", "Before \n experiment", "Control", "Direct", "Water- \n mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_k_merulina <- plot_kusu_merulina_phylum + scale_fill_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Hydrogenedentes" = "#FF00FF", "Myxococcota" = "#CCCC33", "Other" = "#CCCCCC", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))

legend_phylum_merulina <- cowplot::get_legend(plot_k_merulina)
grid.newpage()
grid.draw(legend_phylum_merulina)

plot_r_merulina_1 <- plot_r_merulina + theme(legend.position = "none")
plot_k_merulina_1 <- plot_k_merulina + theme(legend.position = "none")

relabun_merulina_20210112 <- ggarrange(plot_r_merulina_1, plot_k_merulina_1, legend_phylum_merulina, labels = c("A", "B"), nrow = 1)

#convert absolute to relative abundance.
#Merulina_Raffles:
pm_r <- tax_glom(raffles_merulina_rarefied_ps, "Phylum")
pm_r0 <- transform_sample_counts(pm_r, function(x)x/sum(x))
pm_r1 <- merge_samples(pm_r0, "macroalgal_treatment_type")
pm_r2 <- transform_sample_counts(pm_r1, function(x)x/sum(x))

#Merulina_Kusu:
pm_k <- tax_glom(kusu_merulina_rarefied_ps, "Phylum")
pm_k0 <- transform_sample_counts(pm_k, function(x)x/sum(x))
pm_k1 <- merge_samples(pm_k0, "macroalgal_treatment_type")
pm_k2 <- transform_sample_counts(pm_k1, function(x)x/sum(x))

##########

#Plot richness of M. ampliata.
richness_merulina <- plot_richness(merulina_rarefied_ps, x = "location", color = "macroalgal_treatment_type", measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) + geom_boxplot()
richness_merulina <- richness_merulina + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
richness_merulina <- richness_merulina + theme(strip.background = element_blank())
richness_merulina <- richness_merulina + labs(x="", y="Alpha Diversity Index\n")
richness_merulina <- richness_merulina + theme(axis.text.x = element_text(angle = 45, hjust = 1))
richness_merulina <- richness_merulina + scale_x_discrete(limits=c("Pulau_Satumu", "Kusu_Island"), labels=c("Pulau Satumu", "Kusu Island"))
richness_merulina <- richness_merulina + theme(legend.title = element_blank())

p_meru <- richness_merulina
newTreatment <- c("field_sample", "pre_experiment", "control", "direct", "5cm")
p_meru$data$macroalgal_treatment_type <- as.character(p_meru$data$macroalgal_treatment_type)
p_meru$data$macroalgal_treatment_type <- factor(p_meru$data$macroalgal_treatment_type, levels = newTreatment)
print(p_meru)

p_meru <- p_meru + scale_color_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(legend.position = "bottom")

##########

#Rerun Kruskal-Wallis Tests based on location for P. acuta.
#Generate richness tables.
rich_raffles_pocillopora = estimate_richness(raffles_pocillopora_rarefied_ps)
rich_raffles_pocillopora = cbind(sample_data(raffles_pocillopora_rarefied_ps), rich_raffles_pocillopora)
rich_kusu_pocillopora = estimate_richness(kusu_pocillopora_rarefied_ps)
rich_kusu_pocillopora = cbind(sample_data(kusu_pocillopora_rarefied_ps), rich_kusu_pocillopora)

#Test normality.
hist(rich_raffles_pocillopora$Observed, main = "Observed", xlab = "", breaks = 10)
hist(rich_raffles_pocillopora$Chao1, main = "Chao1", xlab = "", breaks = 10)
hist(rich_raffles_pocillopora$Shannon, main = "Shannon", xlab = "", breaks = 10)
hist(rich_raffles_pocillopora$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)

shapiro.test(rich_raffles_pocillopora$Observed)
shapiro.test(rich_raffles_pocillopora$Chao1)
shapiro.test(rich_raffles_pocillopora$Shannon)
shapiro.test(rich_raffles_pocillopora$InvSimpson)

hist(rich_kusu_pocillopora$Observed, main = "Observed", xlab = "", breaks = 10)
hist(rich_kusu_pocillopora$Chao1, main = "Chao1", xlab = "", breaks = 10)
hist(rich_kusu_pocillopora$Shannon, main = "Shannon", xlab = "", breaks = 10)
hist(rich_kusu_pocillopora$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)

shapiro.test(rich_kusu_pocillopora$Observed)
shapiro.test(rich_kusu_pocillopora$Chao1)
shapiro.test(rich_kusu_pocillopora$Shannon)
shapiro.test(rich_kusu_pocillopora$InvSimpson)

#Run Kruskal-Wallis Tests and posthoc Dunn Test.
kruskal.test(Observed~sample_colony, data = rich_raffles_pocillopora)
kruskal.test(Observed~sample_type, data = rich_raffles_pocillopora)
kruskal.test(Observed~macroalgal_treatment_type, data = rich_raffles_pocillopora)
dunn.test::dunn.test(rich_raffles_pocillopora$Observed, rich_raffles_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Chao1~sample_colony, data = rich_raffles_pocillopora)
kruskal.test(Chao1~sample_type, data = rich_raffles_pocillopora)
kruskal.test(Chao1~macroalgal_treatment_type, data = rich_raffles_pocillopora)
dunn.test::dunn.test(rich_raffles_pocillopora$Chao1, rich_raffles_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Shannon~sample_colony, data = rich_raffles_pocillopora)
kruskal.test(Shannon~sample_type, data = rich_raffles_pocillopora)
kruskal.test(Shannon~macroalgal_treatment_type, data = rich_raffles_pocillopora)
dunn.test::dunn.test(rich_raffles_pocillopora$Shannon, rich_raffles_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(InvSimpson~sample_colony, data = rich_raffles_pocillopora)
kruskal.test(InvSimpson~sample_type, data = rich_raffles_pocillopora)
kruskal.test(InvSimpson~macroalgal_treatment_type, data = rich_raffles_pocillopora)
dunn.test::dunn.test(rich_raffles_pocillopora$InvSimpson, rich_raffles_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Observed~sample_colony, data = rich_kusu_pocillopora)
kruskal.test(Observed~sample_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Observed, rich_kusu_pocillopora$sample_type, method = "bonferroni")
kruskal.test(Observed~macroalgal_treatment_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Observed, rich_kusu_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Chao1~sample_colony, data = rich_kusu_pocillopora)
kruskal.test(Chao1~sample_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Chao1, rich_kusu_pocillopora$sample_type, method = "bonferroni")
kruskal.test(Chao1~macroalgal_treatment_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Chao1, rich_kusu_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Shannon~sample_colony, data = rich_kusu_pocillopora)
kruskal.test(Shannon~sample_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Shannon, rich_kusu_pocillopora$sample_type, method = "bonferroni")
kruskal.test(Shannon~macroalgal_treatment_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$Shannon, rich_kusu_pocillopora$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(InvSimpson~sample_colony, data = rich_kusu_pocillopora)
kruskal.test(InvSimpson~sample_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$InvSimpson, rich_kusu_pocillopora$sample_type, method = "bonferroni")
kruskal.test(InvSimpson~macroalgal_treatment_type, data = rich_kusu_pocillopora)
dunn.test::dunn.test(rich_kusu_pocillopora$InvSimpson, rich_kusu_pocillopora$macroalgal_treatment_type, method = "bonferroni")

#Generate median table.
test_raffles_pocillopora_means <- rich_raffles_pocillopora %>% group_by(macroalgal_treatment_type) %>% dplyr::summarise(median_Observed = median(Observed), mean_Observed = mean(Observed), median_Chao1 = median(Chao1), mean_Chao1 = mean(Chao1), median_Shannon = median(Shannon), mean_Shannon = mean(Shannon), median_InvSimpson = median(InvSimpson), mean_InvSimpson = mean(InvSimpson))
test_kusu_pocillopora_means <- rich_kusu_pocillopora %>% group_by(macroalgal_treatment_type) %>% dplyr::summarise(median_Observed = median(Observed), mean_Observed = mean(Observed), median_Chao1 = median(Chao1), mean_Chao1 = mean(Chao1), median_Shannon = median(Shannon), mean_Shannon = mean(Shannon), median_InvSimpson = median(InvSimpson), mean_InvSimpson = mean(InvSimpson))
write.csv(test_raffles_pocillopora_means, "20210112_raffles_pocillopora_medians.csv")
write.csv(test_kusu_pocillopora_means, "20210112_kusu_pocillopora_medians.csv")

##########

#Rerun Kruskal-Wallis Test based on location for M. ampliata.
#Generate richness tables.
rich_raffles_merulina = estimate_richness(raffles_merulina_rarefied_ps)
rich_raffles_merulina = cbind(sample_data(raffles_merulina_rarefied_ps), rich_raffles_merulina)
rich_kusu_merulin = estimate_richness(kusu_merulina_rarefied_ps)
rich_kusu_merulin = cbind(sample_data(kusu_merulina_rarefied_ps), rich_kusu_merulin)

#Test normality.
hist(rich_raffles_merulina$Observed, main = "Observed", xlab = "", breaks = 10)
hist(rich_raffles_merulina$Chao1, main = "Chao1", xlab = "", breaks = 10)
hist(rich_raffles_merulina$Shannon, main = "Shannon", xlab = "", breaks = 10)
hist(rich_raffles_merulina$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)

shapiro.test(rich_raffles_merulina$Observed)
shapiro.test(rich_raffles_merulina$Chao1)
shapiro.test(rich_raffles_merulina$Shannon)
shapiro.test(rich_raffles_merulina$InvSimpson)

hist(rich_kusu_merulin$Observed, main = "Observed", xlab = "", breaks = 10)
hist(rich_kusu_merulin$Chao1, main = "Chao1", xlab = "", breaks = 10)
hist(rich_kusu_merulin$Shannon, main = "Shannon", xlab = "", breaks = 10)
hist(rich_kusu_merulin$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)

shapiro.test(rich_kusu_merulin$Observed)
shapiro.test(rich_kusu_merulin$Chao1)
shapiro.test(rich_kusu_merulin$Shannon)
shapiro.test(rich_kusu_merulin$InvSimpson)

#Run Kruskal-Wallis Tests and posthoc Dunn's Test.
kruskal.test(Observed~sample_colony, data = rich_raffles_merulina)
kruskal.test(Observed~sample_type, data = rich_raffles_merulina)
kruskal.test(Observed~macroalgal_treatment_type, data = rich_raffles_merulina)

kruskal.test(Chao1~sample_colony, data = rich_raffles_merulina)
dunn.test::dunn.test(rich_raffles_merulina$Chao1, rich_raffles_merulina$sample_colony, method = "bonferroni")
kruskal.test(Chao1~sample_type, data = rich_raffles_merulina)
kruskal.test(Chao1~macroalgal_treatment_type, data = rich_raffles_merulina)

kruskal.test(Shannon~sample_colony, data = rich_raffles_merulina)
kruskal.test(Shannon~sample_type, data = rich_raffles_merulina)
kruskal.test(Shannon~macroalgal_treatment_type, data = rich_raffles_merulina)

kruskal.test(InvSimpson~sample_colony, data = rich_raffles_merulina)
kruskal.test(InvSimpson~sample_type, data = rich_raffles_merulina)
kruskal.test(InvSimpson~macroalgal_treatment_type, data = rich_raffles_merulina)

kruskal.test(Observed~sample_colony, data = rich_kusu_merulin)
kruskal.test(Observed~sample_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$Observed, rich_kusu_merulin$sample_type, method = "bonferroni")
kruskal.test(Observed~macroalgal_treatment_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$Observed, rich_kusu_merulin$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Chao1~sample_colony, data = rich_kusu_merulin)
kruskal.test(Chao1~sample_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$Chao1, rich_kusu_merulin$sample_type, method = "bonferroni")
kruskal.test(Chao1~macroalgal_treatment_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$Chao1, rich_kusu_merulin$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(Shannon~sample_colony, data = rich_kusu_merulin)
kruskal.test(Shannon~sample_type, data = rich_kusu_merulin)
kruskal.test(Shannon~macroalgal_treatment_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$Shannon, rich_kusu_merulin$macroalgal_treatment_type, method = "bonferroni")

kruskal.test(InvSimpson~sample_colony, data = rich_kusu_merulin)
kruskal.test(InvSimpson~sample_type, data = rich_kusu_merulin)
kruskal.test(InvSimpson~macroalgal_treatment_type, data = rich_kusu_merulin)
dunn.test::dunn.test(rich_kusu_merulin$InvSimpson, rich_kusu_merulin$macroalgal_treatment_type, method = "bonferroni")

#Generate sample means, median.
test_raffles_merulina_means <- rich_raffles_merulina %>% group_by(macroalgal_treatment_type) %>% dplyr::summarise(median_Observed = median(Observed), mean_Observed = mean(Observed), median_Chao1 = median(Chao1), mean_Chao1 = mean(Chao1), median_Shannon = median(Shannon), mean_Shannon = mean(Shannon), median_InvSimpson = median(InvSimpson), mean_InvSimpson = mean(InvSimpson))
test_kusu_merulina_means <- rich_kusu_merulin %>% group_by(macroalgal_treatment_type) %>% dplyr::summarise(median_Observed = median(Observed), mean_Observed = mean(Observed), median_Chao1 = median(Chao1), mean_Chao1 = mean(Chao1), median_Shannon = median(Shannon), mean_Shannon = mean(Shannon), median_InvSimpson = median(InvSimpson), mean_InvSimpson = mean(InvSimpson))
write.csv(test_raffles_merulina_means, "20210112_raffles_merulina_medians.csv")
write.csv(test_kusu_merulina_means, "20210112_kusu_merulina_medians.csv")

##########

#Plot NMDS plots.
#Generate data tables for plotting.
#For overall.
set.seed(1)
data.scores_overall<- vegan::scores(ordinate_bray_nmds_combined_overall)
data.scores_overall <- as.data.frame(scores(ordinate_bray_nmds_combined_overall$points))
data.scores_overall$sample_type <- df_combined_overall$sample_type
head(data.scores_overall)

#For combined P. acuta.
data.scores_pocillopora <- vegan::scores(ordinate_bray_nmds_pocillopora)
data.scores_pocillopora <- as.data.frame(scores(ordinate_bray_nmds_pocillopora$points))
data.scores_pocillopora$location <- df_pocillopora$location
data.scores_pocillopora$macroalgal_treatment_type <- df_pocillopora$macroalgal_treatment_type
head(data.scores_pocillopora)

#For combined M. ampliata.
data.scores_merulina <- vegan::scores(ordinate_bray_nmds_merulina)
data.scores_merulina <- as.data.frame(scores(ordinate_bray_nmds_merulina$points))
data.scores_merulina$location <- df_merulina$location
data.scores_merulina$macroalgal_treatment_type <- df_merulina$macroalgal_treatment_type
head(data.scores_merulina)

#Plot NMDS plot for overall.
plot_nmds_overall <- ggplot(data.scores_overall, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_overall$sample_type), shape=factor(data.scores_overall$sample_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_overall <- plot_nmds_overall + labs(shape = "Sample Type", color = "Sample Type")
plot_nmds_overall <- plot_nmds_overall + annotate("text", x = 1, y = 2, label = "Stress = 0.176")

p1 <- plot_nmds_overall + scale_color_manual(name = "Sample", labels = c("Field coral", "Tank coral", "Macroalgae", "Seawater", "Tank seawater"), limits = c("coral", "coral_tank", "macroalgae", "seawater", "seawater_tank"), values = c("red", "orange", "brown", "blue", "navy"))
p1 <- p1 + scale_shape_manual(name = "Sample", labels = c("Field coral", "Tank coral", "Macroalgae", "Seawater", "Tank seawater"), values = c(19, 17, 15, 3, 7))
p1 <- p1 + theme(legend.key = element_rect(fill = "white"))
p1 <- p1 + stat_ellipse(data = data.scores_overall , aes(x=MDS1, y=MDS2, color=sample_type))

#Plot NMDS plot for P. acuta.
r_k_treat_pocillopora <- subset(data.scores_pocillopora, macroalgal_treatment_type != "field_sample")
r_k_treat_pocillopora <- subset(r_k_treat_pocillopora, macroalgal_treatment_type != "pre_experiment")

plot_nmds_pocillopora <- ggplot(data.scores_pocillopora, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_pocillopora$location), shape=factor(data.scores_pocillopora$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_pocillopora <- plot_nmds_pocillopora + labs(shape = "Sample Type", color = "Location")
plot_nmds_pocillopora <- plot_nmds_pocillopora + annotate("text", x = 0.4, y = 0.5, label = "Stress = 0.133")

p2 <- plot_nmds_pocillopora + scale_shape_discrete(labels = c("Field Sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_experiment", "control", "direct", "5cm"))
p2 <- p2 + scale_color_discrete(labels = c("Pulau Satumu", "Kusu Island"), limits = c("Pulau_Satumu", "Kusu_Island"))
p2 <- p2 + theme(legend.key = element_rect(fill = "white"))

p2 <- p2 + stat_ellipse(data = r_k_treat_pocillopora, aes(x=MDS1, y=MDS2, color=location))

#Plot NMDS plot for M. ampliata.
r_k_treat_merulina <- subset(data.scores_merulina, macroalgal_treatment_type != "field_sample")
r_k_treat_merulina <- subset(r_k_treat_merulina, macroalgal_treatment_type != "pre_experiment")

plot_nmds_merulina <- ggplot(data.scores_merulina, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_merulina$location), shape=factor(data.scores_merulina$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_merulina <- plot_nmds_merulina + labs(shape = "Sample Type", color = "Location")
plot_nmds_merulina <- plot_nmds_merulina + annotate("text", x = 0.2, y = 0.55, label = "Stress = 0.146")

p3 <- plot_nmds_merulina + scale_shape_discrete(labels = c("Field Sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_experiment", "control", "direct", "5cm"))
p3 <- p3 + scale_color_discrete(labels = c("Pulau Satumu", "Kusu Island"), limits = c("Pulau_Satumu", "Kusu_Island"))
p3 <- p3 + theme(legend.key = element_rect(fill = "white"))

p3 <- p3 + stat_ellipse(data = r_k_treat_merulina, aes(x=MDS1, y=MDS2, color=location))

##########

#Try to merge the 2 P.acuta plots into 1 image.
#Extract legend.
p2_legend_mod <-p2 + theme(legend.direction = "horizontal", legend.position = "bottom", legend.title = element_blank())

legend_nmds_pocillopora <- cowplot::get_legend(p2_legend_mod)
grid.newpage()
grid.draw(legend_nmds_pocillopora)

library(magick)
library(patchwork)

plot_20200228_pocillopora <- image_read("C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112 mesocosm/20210112_new/20210112_new/beta-diversity/p1_mod.tiff") %>% image_ggplot()
plot_20210112_pocillopora <- image_read("C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112 mesocosm/20210112_new/20210112_new/beta-diversity/p2_mod.tiff") %>% image_ggplot()

plot_a <- ggarrange(plot_20200228_pocillopora, plot_20210112_pocillopora, labels = c("A", "B"), ncol = 2, nrow = 1)

plot_a_final <- cowplot::plot_grid(plot_a, legend_nmds_pocillopora, ncol = 1, rel_heights = c(0.95, 0.05))

##########
#Location-based PERMANOVA for treated microbiomes.
set.seed(1)

#Pocillopora:
pocillopora_raffles_tank_treated_ps <- subset_samples(pocillopora_raffles_ps, macroalgal_treatment_type != "field_sample")
pocillopora_raffles_tank_treated_ps <- subset_samples(pocillopora_raffles_tank_treated_ps, macroalgal_treatment_type != "pre_experiment")
pocillopora_kusu_tank_treated_ps <- subset_samples(pocillopora_kusu_ps, macroalgal_treatment_type != "field_sample")
pocillopora_kusu_tank_treated_ps <- subset_samples(pocillopora_kusu_tank_treated_ps, macroalgal_treatment_type != "pre_experiment")

pocillopora_raffles_tank_treated_ps_sqrt <- sqrt(otu_table(pocillopora_raffles_tank_treated_ps))
pocillopora_kusu_tank_treated_ps_sqrt <- sqrt(otu_table(pocillopora_kusu_tank_treated_ps))

bray_raffles_pocillopora_treated <- phyloseq::distance(pocillopora_raffles_tank_treated_ps_sqrt, method = "bray")
df_raffles_pocillopora_treated <-data.frame(sample_data(pocillopora_raffles_tank_treated_ps))
adonis2(bray_raffles_pocillopora_treated~macroalgal_treatment_type, strata = df_raffles_pocillopora_treated$sample_colony, data = df_raffles_pocillopora_treated)
beta_raffles_pocillopora_treated <- betadisper(bray_raffles_pocillopora_treated, df_raffles_pocillopora_treated$macroalgal_treatment_type)
permutest(beta_raffles_pocillopora_treated)

bray_kusu_pocillopora_treated <- phyloseq::distance(pocillopora_kusu_tank_treated_ps_sqrt, method = "bray")
df_kusu_pocillopora_treated <- data.frame(sample_data(pocillopora_kusu_tank_treated_ps))
adonis2(bray_kusu_pocillopora_treated~macroalgal_treatment_type, strata = df_kusu_pocillopora_treated$sample_colony, data = df_kusu_pocillopora_treated)
beta_kusu_pocillopora_treated <- betadisper(bray_kusu_pocillopora_treated, df_kusu_pocillopora_treated$macroalgal_treatment_type)
permutest(beta_kusu_pocillopora_treated)

#For Merulina.
merulina_raffles_tank_treated_ps <- subset_samples(merulina_raffles_ps, macroalgal_treatment_type != "field_sample")
merulina_raffles_tank_treated_ps <- subset_samples(merulina_raffles_tank_treated_ps, macroalgal_treatment_type != "pre_experiment")
merulina_kusu_tank_treated_ps <- subset_samples(merulina_kusu_ps, macroalgal_treatment_type != "field_sample")
merulina_kusu_tank_treated_ps <- subset_samples(merulina_kusu_tank_treated_ps, macroalgal_treatment_type != "pre_experiment")

merulina_raffles_tank_treated_ps_sqrt <- sqrt(otu_table(merulina_raffles_tank_treated_ps))
merulina_kusu_tank_treated_ps_sqrt <- sqrt(otu_table(merulina_kusu_tank_treated_ps))

bray_raffles_merulina_treated <- phyloseq::distance(merulina_raffles_tank_treated_ps_sqrt, method = "bray")
df_raffles_merulina_treated <- data.frame(sample_data(merulina_raffles_tank_treated_ps))
adonis2(bray_raffles_merulina_treated~macroalgal_treatment_type, strata = df_raffles_merulina_treated$sample_colony, data = df_raffles_merulina_treated)
beta_raffles_merulina_treated <- betadisper(bray_raffles_merulina_treated, df_raffles_merulina_treated$macroalgal_treatment_type)
permutest(beta_raffles_merulina_treated)

bray_kusu_merulina_treated <- phyloseq::distance(merulina_kusu_tank_treated_ps_sqrt, method = "bray")
df_kusu_merulina_treated <- data.frame(sample_data(merulina_kusu_tank_treated_ps))
adonis2(bray_kusu_merulina_treated~macroalgal_treatment_type, strata = df_kusu_merulina_treated$sample_colony, data = df_kusu_merulina_treated)
beta_kusu_merulina_treated <- betadisper(bray_kusu_merulina_treated, df_kusu_merulina_treated$macroalgal_treatment_type)
permutest(beta_kusu_merulina_treated)

#Raffles_pocillopora:
bray_raffles_poc <- phyloseq::distance(pocillopora_raffles_sqrt_ps, method = "bray")
df_raffles_poc <- data.frame(sample_data(pocillopora_raffles_sqrt_ps))
adonis2(bray_raffles_poc~macroalgal_treatment_type, strata = df_raffles_poc$sample_colony, data = df_raffles_poc)
pairwise.adonis2(bray_raffles_poc~macroalgal_treatment_type, strata = 'sample_colony', data = df_raffles_poc)
beta_raffles_poc <- betadisper(bray_raffles_poc, df_raffles_poc$macroalgal_treatment_type)
permutest(beta_raffles_poc)

#Raffles_merulina:
bray_raffles_mer <- phyloseq::distance(merulina_raffles_sqrt_ps, method = "bray")
df_raffles_mer <- data.frame(sample_data(merulina_raffles_sqrt_ps))
adonis2(bray_raffles_mer~macroalgal_treatment_type, strata = df_raffles_mer$sample_colony, data = df_raffles_mer)
pairwise.adonis2(bray_raffles_mer~macroalgal_treatment_type, strata = 'sample_colony', data = df_raffles_mer)
beta_raffles_mer <- betadisper(bray_raffles_mer, df_raffles_mer$macroalgal_treatment_type)
permutest(beta_raffles_mer)
TukeyHSD(beta_raffles_mer)

#Kusu_pocillopora:
bray_kusu_poc <- phyloseq::distance(pocillopora_kusu_sqrt_ps, method = "bray")
df_kusu_poc <- data.frame(sample_data(pocillopora_kusu_sqrt_ps))
adonis2(bray_kusu_poc~macroalgal_treatment_type, strata = df_kusu_poc$sample_colony, data = df_kusu_poc)
pairwise.adonis2(bray_kusu_poc~macroalgal_treatment_type, strata = 'sample_colony', data = df_kusu_poc)
beta_kusu_poc <- betadisper(bray_kusu_poc, df_kusu_poc$macroalgal_treatment_type)
permutest(beta_kusu_poc)

#Kusu_merulina:
bray_kusu_mer <- phyloseq::distance(merulina_kusu_sqrt_ps, method = "bray")
df_kusu_mer <- data.frame(sample_data(merulina_kusu_sqrt_ps))
adonis2(bray_kusu_mer~macroalgal_treatment_type, strata = df_kusu_mer$sample_colony, data = df_kusu_mer)
pairwise.adonis2(bray_kusu_mer~macroalgal_treatment_type, strata = 'sample_colony', data = df_kusu_mer)
beta_kusu_mer <- betadisper(bray_kusu_mer, df_kusu_mer$macroalgal_treatment_type)
permutest(beta_kusu_mer)
TukeyHSD(beta_kusu_mer)

##########
#Missing PERMANOVA for Merulina.
set.seed(1)
adonis2(bray_merulina~location + macroalgal_treatment_type, strata = df_merulina$sample_colony, data = df_merulina)

##########

#Plot location based NMDS plots for treated microbiomes.
set.seed(1)

#For Pocillopora.
ordinate_poc_raffles_tank_treated <- ordinate(pocillopora_raffles_tank_treated_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_poc_raffles_tank_treated <- vegan::scores(ordinate_poc_raffles_tank_treated)
data.scores_poc_raffles_tank_treated <- as.data.frame(scores(ordinate_poc_raffles_tank_treated$points))
data.scores_poc_raffles_tank_treated$macroalgal_treatment_type <- df_raffles_pocillopora_treated$macroalgal_treatment_type
head(data.scores_poc_raffles_tank_treated)

plot_nmds_poc_raffles_tank_treated <- ggplot(data.scores_poc_raffles_tank_treated, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_poc_raffles_tank_treated$macroalgal_treatment_type), shape=factor(data.scores_poc_raffles_tank_treated$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_poc_raffles_tank_treated <- plot_nmds_poc_raffles_tank_treated + labs(shape = "Contact Type", color = "Contact Type")
plot_nmds_poc_raffles_tank_treated <- plot_nmds_poc_raffles_tank_treated + annotate("text", x = 0.50, y = -0.2, label = "Stress = 0.124")

p4 <- plot_nmds_poc_raffles_tank_treated + scale_shape_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p4 <- p4 + scale_color_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p4 <- p4 + theme(legend.key = element_rect(fill = "white"))

ordinate_poc_kusu_tank_treated <- ordinate(pocillopora_kusu_tank_treated_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_poc_kusu_tank_treated <- vegan::scores(ordinate_poc_kusu_tank_treated)
data.scores_poc_kusu_tank_treated <- as.data.frame(scores(ordinate_poc_kusu_tank_treated$points))
data.scores_poc_kusu_tank_treated$macroalgal_treatment_type <- df_kusu_pocillopora_treated$macroalgal_treatment_type
head(data.scores_poc_kusu_tank_treated)

plot_nmds_poc_kusu_tank_treated <- ggplot(data.scores_poc_kusu_tank_treated, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_poc_kusu_tank_treated$macroalgal_treatment_type), shape=factor(data.scores_poc_kusu_tank_treated$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_poc_kusu_tank_treated <- plot_nmds_poc_kusu_tank_treated + labs(shape = "Contact Type", color = "Contact Type")
plot_nmds_poc_kusu_tank_treated <- plot_nmds_poc_kusu_tank_treated + annotate("text", x = 0.50, y = -0.2, label = "Stress = 0.143")

p5 <- plot_nmds_poc_kusu_tank_treated + scale_shape_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p5 <- p5 + scale_color_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p5 <- p5 + theme(legend.key = element_rect(fill = "white"))

#Combine plots.
plot_poc_tank_treated <- ggarrange(p4, p5, labels = c("A", "B"), ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

#For Merulina.
ordinate_mer_raffles_tank_treated <- ordinate(merulina_raffles_tank_treated_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_mer_raffles_tank_treated <- vegan::scores(ordinate_mer_raffles_tank_treated)
data.scores_mer_raffles_tank_treated <- as.data.frame(scores(ordinate_mer_raffles_tank_treated$points))
data.scores_mer_raffles_tank_treated$macroalgal_treatment_type <- df_raffles_merulina_treated$macroalgal_treatment_type
head(data.scores_mer_raffles_tank_treated)

plot_nmds_mer_raffles_tank_treated <- ggplot(data.scores_mer_raffles_tank_treated, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_mer_raffles_tank_treated$macroalgal_treatment_type), shape=factor(data.scores_mer_raffles_tank_treated$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_mer_raffles_tank_treated <- plot_nmds_mer_raffles_tank_treated + labs(shape = "Contact Type", color = "Contact Type")
plot_nmds_mer_raffles_tank_treated <- plot_nmds_mer_raffles_tank_treated + annotate("text", x = 0.25, y = -0.25, label = "Stress = 0.155")

p6 <- plot_nmds_mer_raffles_tank_treated + scale_shape_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p6 <- p6 + scale_color_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p6 <- p6 + theme(legend.key = element_rect(fill = "white"))

ordinate_mer_kusu_tank_treated <- ordinate(merulina_kusu_tank_treated_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_mer_kusu_tank_treated <- vegan::scores(ordinate_mer_kusu_tank_treated)
data.scores_mer_kusu_tank_treated <- as.data.frame(scores(ordinate_mer_kusu_tank_treated$points))
data.scores_mer_kusu_tank_treated$macroalgal_treatment_type <- df_kusu_merulina_treated$macroalgal_treatment_type
head(data.scores_mer_kusu_tank_treated)

plot_nmds_mer_kusu_tank_treated <- ggplot(data.scores_mer_kusu_tank_treated, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_mer_kusu_tank_treated$macroalgal_treatment_type), shape=factor(data.scores_mer_kusu_tank_treated$macroalgal_treatment_type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_mer_kusu_tank_treated <- plot_nmds_mer_kusu_tank_treated + labs(shape = "Contact Type", color = "Contact Type")
plot_nmds_mer_kusu_tank_treated <- plot_nmds_mer_kusu_tank_treated + annotate("text", x = 0.0, y = -0.2, label = "Stress = 0.116")

p7 <- plot_nmds_mer_kusu_tank_treated + scale_shape_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p7 <- p7 + scale_color_discrete(labels = c("Control", "Direct", "Water-mediated"), limits = c("control", "direct", "5cm"))
p7 <- p7 + theme(legend.key = element_rect(fill = "white"))

#Combine both plots together.
plot_mer_tank_treated <- ggarrange(p6, p7, labels = c("A", "B"), ncol = 1, nrow = 2, common.legend = TRUE, legend = "right")

##########

total_raw_reads = data.table(as(sample_data(combined_overall_ps), "data.frame"), TotalReads = sample_sums(combined_overall_ps), keep.rownames = TRUE)
write.csv(total_raw_reads, "20210112_overall_total_raw_reads.csv")
combined_coral_rarefied_reads = data.table(as(sample_data(combined_coral_rarefied_ps), "data.frame"), Reads = sample_sums(combined_coral_rarefied_ps), keep.rownames = TRUE)
write.csv(combined_coral_rarefied_reads, "20210112_combined_coral_rarefied_reads.csv")

##########
OTU_table_phylum_merulina_raffles <- as.data.frame(otu_table(phylum_raffles_merulina))
write.csv(OTU_table_phylum_merulina_raffles, "20210112_raffles_merulina_phylum.csv")

OTU_table_phylum_merulina_kusu <- as.data.frame(otu_table(phylum_kusu_merulina))
write.csv(OTU_table_phylum_merulina_kusu, "20210112_kusu_merulina_phylum.csv")

genus_raffles_merulina <- raffles_merulina_rarefied_ps %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_raffles_merulina
head(otu_table(genus_raffles_merulina))
OTU_table_genus_merulina_raffles <- as.data.frame(otu_table(genus_raffles_merulina))
write.csv(OTU_table_genus_merulina_raffles, "20210112_raffles_merulina_genus.csv")

genus_kusu_merulina <- kusu_merulina_rarefied_ps %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_kusu_merulina
head(otu_table(genus_kusu_merulina))
OTU_table_genus_merulina_kusu <- as.data.frame(otu_table(genus_kusu_merulina))
write.csv(OTU_table_genus_merulina_kusu, "20210112_kusu_merulina_genus.csv")

OTU_table_phylum_pocillopora_raffles <- as.data.frame(otu_table(phylum_raffles_pocillopora))
write.csv(OTU_table_phylum_pocillopora_raffles, "20210112_raffles_pocillopora_phylum.csv")

OTU_table_phylum_pocillopora_kusu <- as.data.frame(otu_table(phylum_kusu_pocillopora))
write.csv(OTU_table_phylum_pocillopora_kusu, "20210112_kusu_pocillopora_phylum.csv")

genus_raffles_pocillopora <- raffles_pocillopora_rarefied_ps %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_raffles_pocillopora
head(otu_table(genus_raffles_pocillopora))
OTU_table_genus_pocillopora_raffles <- as.data.frame(otu_table(genus_raffles_pocillopora))
write.csv(OTU_table_genus_pocillopora_raffles, "20210112_raffles_pocillopora_genus.csv")

genus_kusu_pocillopora <- kusu_pocillopora_rarefied_ps %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_kusu_pocillopora
head(otu_table(genus_kusu_pocillopora))
OTU_table_genus_pocillopora_kusu <- as.data.frame(otu_table(genus_kusu_pocillopora))
write.csv(OTU_table_genus_pocillopora_kusu, "20210112_kusu_pocillopora_genus.csv")


##########

#LME for alpha-diversity.
library(lme4)
library(lmerTest)
library(emmeans)
set.seed(1)

#Raffles_merulina:
mod1 <- lmer(log(Observed)~ macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina)
summary(mod1)
anova(mod1)
emmeans(mod1, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1))

mod1a <- lmer(log(Chao1)~ macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina)
summary(mod1a)
anova(mod1a)
plot(mod1a)
qqnorm(resid(mod1a))
qqline(resid(mod1a))

mod1b <- lmer(log(Shannon)~ macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina)
summary(mod1b)
anova(mod1b)
plot(mod1b)
qqnorm(resid(mod1b))
qqline(resid(mod1b))

mod1c <- lmer(log(InvSimpson)~ macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina)
summary(mod1c)
anova(mod1c)
plot(mod1c)
qqnorm(resid(mod1c))
qqline(resid(mod1c))

#Raffles_pocillopora:
mod2 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_pocillopora)
summary(mod2)
anova(mod2)
emmeans(mod2, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2))

mod2a <- lmer(log(Chao1)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_pocillopora)
summary(mod2a)
anova(mod2a)
emmeans(mod2a, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod2a)
qqnorm(resid(mod2a))
qqline(resid(mod2a))

mod2b <- lm(log(Shannon)~macroalgal_treatment_type + sample_colony, data = rich_raffles_pocillopora)
summary(mod2b)
anova(mod2b)
emmeans(mod2b, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod2b)

mod2c <- lm(log(InvSimpson)~macroalgal_treatment_type + sample_colony, data = rich_raffles_pocillopora)
summary(mod2c)
anova(mod2c)
emmeans(mod2c, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod2c)

#Kusu_merulina:
mod3 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulin)
summary(mod3)
anova(mod3)
emmeans(mod3, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod3)
qqnorm(resid(mod3))
qqline(resid(mod3))

mod3a <- lmer(log(Chao1)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulin)
summary(mod3a)
anova(mod3a)
emmeans(mod3a, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod3a)
qqnorm(resid(mod3a))
qqline(resid(mod3a))

mod3b <- lmer(log(Shannon)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulin)
summary(mod3b)
anova(mod3b)
emmeans(mod3b, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod3b)
qqnorm(resid(mod3b))
qqline(resid(mod3b))

mod3c <- lmer(log(InvSimpson)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulin)
summary(mod3c)
anova(mod3c)
emmeans(mod3c, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod3c)
qqnorm(resid(mod3c))
qqline(resid(mod3c))

#Kusu_pocillopora:
mod4 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora)
summary(mod4)
anova(mod4)
emmeans(mod4, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod4)
qqnorm(resid(mod4))
qqline(resid(mod4))

mod4a <- lm(log(Chao1)~macroalgal_treatment_type + sample_colony, data = rich_kusu_pocillopora)
summary(mod4a)
anova(mod4a)
emmeans(mod4a, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod4a)

mod4b <- lmer(log(Shannon)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora)
summary(mod4b)
anova(mod4b)
emmeans(mod4b, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod4b)
qqnorm(resid(mod4b))
qqline(resid(mod4b))

mod4c <- lmer(log(InvSimpson)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora)
summary(mod4c)
anova(mod4c)
emmeans(mod4c, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod4c)
qqnorm(resid(mod4c))
qqline(resid(mod4c))

#For treated samples:
merulina_raffles_tank_treated_ps_rarefied <- rarefy_even_depth(merulina_raffles_tank_treated_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(merulina_raffles_tank_treated_ps)), replace = FALSE)
merulina_kusu_tank_treated_ps_rarefied <- rarefy_even_depth(merulina_kusu_tank_treated_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(merulina_kusu_tank_treated_ps)), replace = FALSE)
pocillopora_raffles_tank_treated_ps_rarefied <- rarefy_even_depth(pocillopora_raffles_tank_treated_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(pocillopora_raffles_tank_treated_ps)), replace = FALSE)
pocillopora_kusu_tank_treated_ps_rarefied <- rarefy_even_depth(pocillopora_kusu_tank_treated_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(pocillopora_kusu_tank_treated_ps)), replace = FALSE)

rich_raffles_merulina_treated <- estimate_richness(merulina_raffles_tank_treated_ps_rarefied)
rich_raffles_merulina_treated <- cbind(rich_raffles_merulina_treated, sample_data(merulina_raffles_tank_treated_ps_rarefied))

rich_kusu_merulina_treated <- estimate_richness(merulina_kusu_tank_treated_ps_rarefied)
rich_kusu_merulina_treated <- cbind(rich_kusu_merulina_treated, sample_data(merulina_kusu_tank_treated_ps_rarefied))

rich_raffles_pocillopora_treated <- estimate_richness(pocillopora_raffles_tank_treated_ps_rarefied)
rich_raffles_pocillopora_treated <- cbind(rich_raffles_pocillopora_treated, sample_data(pocillopora_raffles_tank_treated_ps_rarefied))

rich_kusu_pocillopora_treated <- estimate_richness(pocillopora_kusu_tank_treated_ps_rarefied)
rich_kusu_pocillopora_treated <- cbind(rich_kusu_pocillopora_treated, sample_data(pocillopora_kusu_tank_treated_ps_rarefied))

#Raffles_merulina_treated:
mod5 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina_treated)
summary(mod5)
anova(mod5)
plot(mod5)
qqnorm(resid(mod5))
qqline(resid(mod5))

mod5a <- lmer(log(Chao1)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina_treated)
summary(mod5a)
anova(mod5a)
plot(mod5a)
qqnorm(resid(mod5a))
qqline(resid(mod5a))

mod5b <- lmer(log(Shannon)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina_treated)
summary(mod5b)
anova(mod5b)
plot(mod5b)
qqnorm(resid(mod5b))
qqline(resid(mod5b))

mod5c <- lmer(log(InvSimpson)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_merulina_treated)
summary(mod5c)
anova(mod5c)
plot(mod5c)
qqnorm(resid(mod5c))
qqline(resid(mod5c))

#Kusu_merulina_treated:
mod6 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulina_treated)
summary(mod6)
anova(mod6)
plot(mod6)
qqnorm(resid(mod6))
qqline(resid(mod6))

mod6a <- lmer(log(Chao1)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_merulina_treated)
summary(mod6a)
anova(mod6a)
plot(mod6a)
qqnorm(resid(mod6a))
qqline(resid(mod6a))

mod6b <- lm(log(Shannon)~macroalgal_treatment_type + sample_colony, data = rich_kusu_merulina_treated)
summary(mod6b)
anova(mod6b)
plot(mod6b)

mod6c <- lm(log(InvSimpson)~macroalgal_treatment_type + sample_colony, data = rich_kusu_merulina_treated)
summary(mod6c)
anova(mod6c)
emmeans(mod6c, list(pairwise~macroalgal_treatment_type), adjust = "tukey", type = "response")
plot(mod6c)

#Raffles_pocillopora_treated:
mod7 <- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_raffles_pocillopora_treated)
summary(mod7)
anova(mod7)
plot(mod7)
qqnorm(resid(mod7))
qqline(resid(mod7))

mod7a <- lm(log(Chao1)~macroalgal_treatment_type + sample_colony, data = rich_raffles_pocillopora_treated)
summary(mod7a)
anova(mod7a)
plot(mod7a)

mod7b <- lm(log(Shannon)~macroalgal_treatment_type + sample_colony, data = rich_raffles_pocillopora_treated)
summary(mod7b)
anova(mod7b)
plot(mod7b)

mod7c <- lm(log(InvSimpson)~macroalgal_treatment_type + sample_colony, data = rich_raffles_pocillopora_treated)
summary(mod7c)
anova(mod7c)
plot(mod7c)

#Kusu_pocillopora_treated:
mod8<- lmer(log(Observed)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora_treated)
summary(mod8)
anova(mod8)
plot(mod8)
qqnorm(resid(mod8))
qqline(resid(mod8))

mod8a<- lmer(log(Chao1)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora_treated)
summary(mod8a)
anova(mod8a)
plot(mod8a)
qqnorm(resid(mod8a))
qqline(resid(mod8a))

mod8b<- lmer(log(Shannon)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora_treated)
summary(mod8b)
anova(mod8b)
plot(mod8b)
qqnorm(resid(mod8b))
qqline(resid(mod8b))

mod8c<- lmer(log(InvSimpson)~macroalgal_treatment_type + (1|sample_colony), data = rich_kusu_pocillopora_treated)
summary(mod8c)
anova(mod8c)
plot(mod8b)
qqnorm(resid(mod8b))
qqline(resid(mod8b))

##########

#Plot combined plots for DESeq2 results.
#DESeq2 data already in database.

#For Merulina:
plot_deseq2_mr_c_d <- ggplot(res_mr1_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_mr_c_d <- plot_deseq2_mr_c_d + coord_flip()
plot_deseq2_mr_c_d <- plot_deseq2_mr_c_d + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_mr_c_d <- plot_deseq2_mr_c_d + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_mr_c_d <- plot_deseq2_mr_c_d + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_mr_c_5 <- ggplot(res_mr2_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_mr_c_5 <- plot_deseq2_mr_c_5 + coord_flip()
plot_deseq2_mr_c_5 <- plot_deseq2_mr_c_5 + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_mr_c_5 <- plot_deseq2_mr_c_5 + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_mr_c_5 <- plot_deseq2_mr_c_5 + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_mk_c_d <- ggplot(res_mk1_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_mk_c_d <- plot_deseq2_mk_c_d + coord_flip()
plot_deseq2_mk_c_d <- plot_deseq2_mk_c_d + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_mk_c_d <- plot_deseq2_mk_c_d + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_mk_c_d <- plot_deseq2_mk_c_d + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_mk_c_5 <- ggplot(res_mk2_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_mk_c_5 <- plot_deseq2_mk_c_5 + coord_flip()
plot_deseq2_mk_c_5 <- plot_deseq2_mk_c_5 + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_mk_c_5 <- plot_deseq2_mk_c_5 + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_mk_c_5 <- plot_deseq2_mk_c_5 + theme(legend.key = element_rect(fill = "white"))

ggarrange(plot_deseq2_mr_c_d, plot_deseq2_mr_c_5, plot_deseq2_mk_c_d, plot_deseq2_mk_c_5, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2)

#For Pocillopora:
plot_deseq2_pr_c_d <- ggplot(res_pr1_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_pr_c_d <- plot_deseq2_pr_c_d + coord_flip()
plot_deseq2_pr_c_d <- plot_deseq2_pr_c_d + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_pr_c_d <- plot_deseq2_pr_c_d + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_pr_c_d <- plot_deseq2_pr_c_d + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_pr_c_5 <- ggplot(res_pr2_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_pr_c_5 <- plot_deseq2_pr_c_5 + coord_flip()
plot_deseq2_pr_c_5 <- plot_deseq2_pr_c_5 + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_pr_c_5 <- plot_deseq2_pr_c_5 + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_pr_c_5 <- plot_deseq2_pr_c_5 + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_pk_c_d <- ggplot(res_pk1_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_pk_c_d <- plot_deseq2_pk_c_d + coord_flip()
plot_deseq2_pk_c_d <- plot_deseq2_pk_c_d + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_pk_c_d <- plot_deseq2_pk_c_d + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_pk_c_d <- plot_deseq2_pk_c_d + theme(legend.key = element_rect(fill = "white"))

plot_deseq2_pk_c_5 <- ggplot(res_pk2_sig, aes(x=Genus, y=log2FoldChange, color=Phylum)) + geom_point(size=2) + theme(axis.text.x = element_text(angle = 90, hjust = 0, vjust = 0))
plot_deseq2_pk_c_5 <- plot_deseq2_pk_c_5 + coord_flip()
plot_deseq2_pk_c_5 <- plot_deseq2_pk_c_5 + scale_color_manual(values = c("Acidobacteriota" = "#9999FF", "Actinobacteriota" = "#0066CC", "Bacteroidota" = "#66CCCC", "Bdellovibrionota" = "#009933", "Campilobacterota" = "#CC9900", "Chloroflexi" = "#FF9999", "Cyanobacteria" = "#FF3300", "Desulfobacterota" = "#FFCC33", "Firmicutes" = "#FFFF33", "Myxococcota" = "#CCCC33", "Nitrospirota" = "#66FF00", "Planctomycetota" = "#3333CC", "Proteobacteria" = "#330066", "Verrucomicrobiota" = "#CCC999"))
plot_deseq2_pk_c_5 <- plot_deseq2_pk_c_5 + theme(panel.background = element_blank(), panel.border = element_rect(fill = NA), axis.line = element_line(colour = "black"))
plot_deseq2_pk_c_5 <- plot_deseq2_pk_c_5 + theme(legend.key = element_rect(fill = "white"))

ggarrange(plot_deseq2_pr_c_d, plot_deseq2_pr_c_5, plot_deseq2_pk_c_d, plot_deseq2_pk_c_5, labels = c("A", "B", "C", "D"), nrow = 2, ncol = 2)
