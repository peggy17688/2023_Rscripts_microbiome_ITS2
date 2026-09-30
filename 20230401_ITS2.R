library(dada2)
library(phyloseq)
library(dplyr)
library(ggplot2)
library(vegan)
library(DECIPHER)
library(data.table)
library(ShortRead)
library(Biostrings)

#Process ITS2 sequences.
path <- "C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112_ITS2/raw_sequences_20210112_ITS2"
list.files(path)

fnFs <- sort(list.files(path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(path, pattern = "_R2_001.fastq.gz", full.names = TRUE))

FWD <- "GAATTGCAGAACTCCGTGAACC"
REV <- "CGGGTTCWCTTGTYTGACTTCATGC"

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients

fnFs.filtN <- file.path(path, "filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "filtN", basename(fnRs))
filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = FALSE)

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

#Trim reads with cutadapt.
cutadapt <- "C:/Users/PEIYIPEG001/Downloads/cutadapt.exe"
system2(cutadapt, args = "--version")

path.cut <- file.path(path, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC) 
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c(R1.flags, R2.flags, "-n", 2, # -n 2 required to remove FWD and REV from reads
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i], # output files
                             fnFs.filtN[i], fnRs.filtN[i])) # input files
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]]))

cutFs <- sort(list.files(path.cut, pattern = "_R1_001.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_R2_001.fastq.gz", full.names = TRUE))

get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)

#Inspect reads.
plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

filtFs <- file.path(path.cut, "filtered", basename(cutFs))
filtRs <- file.path(path.cut, "filtered", basename(cutRs))

out <- filterAndTrim(cutFs, filtFs, cutRs, filtRs, maxN = 0, maxEE = c(10, 10), 
                     truncQ = 2, minLen = 50, rm.phix = TRUE, compress = TRUE, multithread = FALSE)  # on windows, set multithread = FALSE

head(out)

#Learn error rates.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

plotErrors(errF, nominalQ = TRUE)

#Dereplicate identical reads.
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names

#Sample inference.
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)

#Merge reads.
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)

#Construct sequence table.
seqtab <- makeSequenceTable(mergers)
dim(seqtab)

#Remove chimeras.
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

table(nchar(getSequences(seqtab.nochim)))

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, 
                                                                       getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace
# sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", 
                     "nonchim")
rownames(track) <- sample.names
head(track)

#Assign taxonomy.
ITS2_taxa <- "C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112_ITS2/ITS2db_trimmed_derep_dada.fasta"
taxa <- assignTaxonomy(seqtab.nochim, ITS2_taxa, multithread = TRUE, tryRC = TRUE)

taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)

#Construction of phylogenetic tree.
library(phangorn)
sequences<-getSequences(seqtab.nochim)
names(sequences)<-sequences

alignment <- AlignSeqs(DNAStringSet(sequences), anchor=NA)
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) # Note, tip order != sequence order
fit = pml(treeNJ, data=phang.align)
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0))

write.csv(seqtab.nochim, "20210112_ITS2_seqtab.nochim.csv")
write.csv(taxa, "20210112_ITS2_taxa.csv")

#Create phyloseq object.
overallOTU_table <- read.table("C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112_ITS2/20210112_ITS2_seqtab.nochim.csv", header = T, row.names = 1, check.names = F, sep = ",")
overall_taxa <- as.matrix(read.table("C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112_ITS2/20210112_ITS2_taxa.csv", header = T, row.names = 1, check.names = F, sep = ","))
overall_sample <- read.table("C:/Users/PEIYIPEG001/Desktop/20230209_backup/20210112_ITS2/20210112_ITS2_metadata.csv", header = T, row.names = 1, check.names = F, sep = ",")

overallOTU_table_ps <- otu_table(overallOTU_table, taxa_are_rows = FALSE)
overall_sample_ps <- sample_data(overall_sample)
overall_taxa_ps <- tax_table(overall_taxa)

colnames(overallOTU_table_ps)
rownames(overall_sample_ps)
rownames(overallOTU_table_ps)
rownames(overall_taxa_ps)

overall_ITS2_ps <- phyloseq(overallOTU_table_ps, overall_taxa_ps, overall_sample_ps, phy_tree(fitGTR$tree))

set.seed(1)
phy_tree(overall_ITS2_ps) <- root(phy_tree(overall_ITS2_ps), sample(taxa_names(overall_ITS2_ps), 1), resolve.root = TRUE)
is.rooted(phy_tree(overall_ITS2_ps))

dna <- Biostrings::DNAStringSet(taxa_names(overall_ITS2_ps))
names(dna) <- taxa_names(overall_ITS2_ps)
overall_ITS2_ps <- merge_phyloseq(overall_ITS2_ps, dna)
taxa_names(overall_ITS2_ps) <- paste0("ASV", seq(ntaxa(overall_ITS2_ps)))

##########

library(dplyr)
library(vegan)
library(pairwiseAdonis)
library(ggplot2)
library(data.table)

merulina_ITS2_ps <- subset_samples(overall_ITS2_ps, sample_species == "Merulina_ampliata")
pocillopora_ITS2_ps <- subset_samples(overall_ITS2_ps, sample_species == "Pocillopora_acuta")
merulina_raffles_ps <- subset_samples(merulina_ITS2_ps, location == "Pulau_Satumu")
merulina_kusu_ps <- subset_samples(merulina_ITS2_ps, location == "Kusu_Island")
pocillopora_raffles_ps <- subset_samples(pocillopora_ITS2_ps, location == "Pulau_Satumu")
pocillopora_kusu_ps <- subset_samples(pocillopora_ITS2_ps, location == "Kusu_Island")

overall_ITS2_ps_rarefied = rarefy_even_depth(overall_ITS2_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(overall_ITS2_ps)), replace = F)
merulina_ITS2_ps_rarefied = rarefy_even_depth(merulina_ITS2_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(merulina_ITS2_ps)), replace = F)
pocillopora_ITS2_ps_rarefied = rarefy_even_depth(pocillopora_ITS2_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(pocillopora_ITS2_ps)), replace = F)
merulina_raffles_rarefied = rarefy_even_depth(merulina_raffles_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(merulina_raffles_ps)), replace = F)
merulina_kusu_rarefied = rarefy_even_depth(merulina_kusu_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(merulina_kusu_ps)), replace = F)
pocillopora_raffles_rarefied = rarefy_even_depth(pocillopora_raffles_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(pocillopora_raffles_ps)), replace = F)
pocillopora_kusu_rarefied = rarefy_even_depth(pocillopora_kusu_ps, rngseed = 1, sample.size = 0.95*min(sample_sums(pocillopora_kusu_ps)), replace = F)

overall_ITS2_ps_sqrt <- sqrt(otu_table(overall_ITS2_ps))
merulina_ITS2_ps_sqrt <- sqrt(otu_table(merulina_ITS2_ps))
pocillopora_ITS2_ps_sqrt <- sqrt(otu_table(pocillopora_ITS2_ps))
merulina_raffles_ps_sqrt <- sqrt(otu_table(merulina_raffles_ps))
merulina_kusu_ps_sqrt <- sqrt(otu_table(merulina_kusu_ps))
pocillopora_raffles_ps_sqrt <- sqrt(otu_table(pocillopora_raffles_ps))
pocillopora_kusu_ps_sqrt <- sqrt(otu_table(pocillopora_kusu_ps))

sdt1 = data.table(as(sample_data(overall_ITS2_ps), "data.frame"), TotalReads = sample_sums(overall_ITS2_ps), keep.rownames = T)
write.csv(sdt1, "20210112_ITS1_totalreads.csv")

sdt2 = data.table(as(sample_data(merulina_ITS2_ps_rarefied), "data.frame"), TotalReads = sample_sums(merulina_ITS2_ps_rarefied), keep.rownames = T)
write.csv(sdt2, "20210112_merulina_rarefied.csv")

sdt3 = data.table(as(sample_data(pocillopora_ITS2_ps_rarefied), "data.frame"), TotalReads = sample_sums(pocillopora_ITS2_ps_rarefied), keep.rownames = T)
write.csv(sdt3, "20210112_pocillopora_rarefied.csv")

sdt4 = data.table(as(sample_data(pocillopora_raffles_rarefied), "data.frame"), TotalReads = sample_sums(pocillopora_raffles_rarefied), keep.rownames = T)
sdt5 = data.table(as(sample_data(pocillopora_kusu_rarefied), "data.frame"), TotalReads = sample_sums(pocillopora_kusu_rarefied), keep.rownames = T)

sdt6 = data.table(as(sample_data(merulina_raffles_rarefied), "data.frame"), TotalReads = sample_sums(merulina_raffles_rarefied), keep.rownames = T)
sdt7 = data.table(as(sample_data(merulina_kusu_rarefied), "data.frame"), TotalReads = sample_sums(merulina_kusu_rarefied), keep.rownames = T)

##########
set.seed(1)

rich_merulina = estimate_richness(merulina_ITS2_ps_rarefied, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
rich_merulina <- cbind(sample_data(merulina_ITS2_ps_rarefied), rich_merulina)
rich_pocillopora = estimate_richness(pocillopora_ITS2_ps_rarefied, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
rich_pocillopora <- cbind(sample_data(pocillopora_ITS2_ps_rarefied), rich_pocillopora)
rich_overall = estimate_richness(overall_ITS2_ps_rarefied, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
rich_overall <- cbind(sample_data(overall_ITS2_ps_rarefied), rich_overall)

write.csv(rich_overall, "20210112_richness_overall.csv")

#Plot alpha-diversity.
plot_richness_merulina <- plot_richness(merulina_ITS2_ps_rarefied, x="location", color = "macroalgal_treatment.type", measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) + geom_boxplot()
plot_richness_merulina <- plot_richness_merulina + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot_richness_merulina <- plot_richness_merulina + theme(strip.background = element_blank())
plot_richness_merulina <- plot_richness_merulina + labs(x="", y="Alpha Diversity Index\n")
plot_richness_merulina <- plot_richness_merulina + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_richness_merulina <- plot_richness_merulina + scale_x_discrete(limits=c("Pulau_Satumu", "Kusu_Island"), labels=c("Pulau Satumu", "Kusu Island"))
plot_richness_merulina <- plot_richness_merulina + theme(legend.title = element_blank())
p1 <- plot_richness_merulina
newTreatment <- c("field_sample", "pre_treatment", "control", "direct", "5cm")
p1$data$macroalgal_treatment.type <- as.character(p1$data$macroalgal_treatment.type)
p1$data$macroalgal_treatment.type <- factor(p1$data$macroalgal_treatment.type, levels = newTreatment)
p1 <- p1 + scale_color_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(legend.position = "bottom")

plot_richness_pocillopora1 <- plot_richness(pocillopora_ITS2_ps_rarefied, x="location", color = "macroalgal_treatment.type", measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) + geom_boxplot()
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + theme(strip.background = element_blank())
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + labs(x="", y="Alpha Diversity Index\n")
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + scale_x_discrete(limits=c("Pulau_Satumu", "Kusu_Island"), labels=c("Pulau Satumu", "Kusu Island"))
plot_richness_pocillopora1 <- plot_richness_pocillopora1 + theme(legend.title = element_blank())
p2 <- plot_richness_pocillopora1
newTreatment <- c("field_sample", "pre_treatment", "control", "direct", "5cm")
p2$data$macroalgal_treatment.type <- as.character(p2$data$macroalgal_treatment.type)
p2$data$macroalgal_treatment.type <- factor(p2$data$macroalgal_treatment.type, levels = newTreatment)
p2 <- p2 + scale_color_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(legend.position = "bottom")

##########

rich_merulina_raffles = estimate_richness(merulina_raffles_rarefied)
rich_merulina_raffles = cbind(sample_data(merulina_raffles_rarefied), rich_merulina_raffles)
rich_merulina_kusu = estimate_richness(merulina_kusu_rarefied)
rich_merulina_kusu = cbind(sample_data(merulina_kusu_rarefied), rich_merulina_kusu)

rich_pocillopora_raffles = estimate_richness(pocillopora_raffles_rarefied, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
rich_pocillopora_raffles = cbind(sample_data(pocillopora_raffles_rarefied), rich_pocillopora_raffles)
rich_pocillopora_kusu = estimate_richness(pocillopora_kusu_rarefied, measures = c("Observed", "Chao1", "Shannon", "InvSimpson"))
rich_pocillopora_kusu = cbind(sample_data(pocillopora_kusu_rarefied), rich_pocillopora_kusu)

write.csv(rich_merulina_raffles, "20210112_merulina_raffles_ITS2_richness.csv")
write.csv(rich_merulina_kusu, "20210112_merulina_kusu_ITS2_richness.csv")
write.csv(rich_pocillopora_raffles, "20210112_pocillopora_raffles_ITS2_richness.csv")
write.csv(rich_pocillopora_kusu, "20210112_pocillopora_kusu_ITS2_richness.csv")

#For Raffles Merulina.
#Test normality.
hist(rich_merulina_raffles$Observed, main = "Observed", xlab = "", breaks = 10)
shapiro.test(rich_merulina_raffles$Observed)

hist(rich_merulina_raffles$Chao1, main = "Chao1", xlab = "", breaks = 10)
shapiro.test(rich_merulina_raffles$Chao1)

hist(rich_merulina_raffles$Shannon, main = "Shannon", xlab = "", breaks = 10)
shapiro.test(rich_merulina_raffles$Shannon)

hist(rich_merulina_raffles$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)
shapiro.test(rich_merulina_raffles$InvSimpson)

aov_test <- aov(Observed~macroalgal_treatment.type, rich_merulina_raffles, strata=rich_merulina_raffles$colony)
summary(aov_test)
TukeyHSD(aov_test)
aov_mer_ps2 <- aov(Observed~colony, rich_merulina_raffles)
summary(aov_mer_ps2)

aov_mer_ps3 <- aov(Chao1~macroalgal_treatment.type, rich_merulina_raffles)
summary(aov_mer_ps3)
TukeyHSD(aov_mer_ps3)
aov_mer_ps4 <- aov(Chao1~colony, rich_merulina_raffles)
summary(aov_mer_ps4)

aov_mer_ps5 <- aov(Shannon~macroalgal_treatment.type, rich_merulina_raffles)
summary(aov_mer_ps5)
aov_mer_ps6 <- aov(Shannon~colony, rich_merulina_raffles)
summary(aov_mer_ps6)

aov_mer_ps7 <- aov(InvSimpson~macroalgal_treatment.type, rich_merulina_raffles)
summary(aov_mer_ps7)
aov_mer_ps8 <- aov(InvSimpson~colony, rich_merulina_raffles)
summary(aov_mer_ps8)

kruskal.test(Observed~macroalgal_treatment.type, data = rich_merulina_raffles)
kruskal.test(Observed~colony, data = rich_merulina_raffles)

kruskal.test(Chao1~macroalgal_treatment.type, data = rich_merulina_raffles)
kruskal.test(Chao1~colony, data = rich_merulina_raffles)

kruskal.test(Shannon~macroalgal_treatment.type, data = rich_merulina_raffles)
kruskal.test(Shannon~colony, data = rich_merulina_raffles)

kruskal.test(InvSimpson~macroalgal_treatment.type, data = rich_merulina_raffles)
kruskal.test(InvSimpson~colony, data = rich_merulina_raffles)

#For Kusu Merulina.
hist(rich_merulina_kusu$Observed, main = "Observed", xlab = "", breaks = 10)
shapiro.test(rich_merulina_kusu$Observed)

hist(rich_merulina_kusu$Chao1, main = "Chao1", xlab = "", breaks = 10)
shapiro.test(rich_merulina_kusu$Chao1)

hist(rich_merulina_kusu$Shannon, main = "Shannon", xlab = "", breaks = 10)
shapiro.test(rich_merulina_kusu$Shannon)

hist(rich_merulina_kusu$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)
shapiro.test(rich_merulina_kusu$InvSimpson)

kruskal.test(Observed~macroalgal_treatment.type, data = rich_merulina_kusu)
kruskal.test(Observed~colony, data = rich_merulina_kusu)
dunn.test::dunn.test(rich_merulina_kusu$Observed, rich_merulina_kusu$colony, method = "bonferroni")

kruskal.test(Chao1~macroalgal_treatment.type, data = rich_merulina_kusu)
kruskal.test(Chao1~colony, data = rich_merulina_kusu)
dunn.test::dunn.test(rich_merulina_kusu$Chao1, rich_merulina_kusu$colony, method = "bonferroni")

kruskal.test(Shannon~macroalgal_treatment.type, data = rich_merulina_kusu)
kruskal.test(Shannon~colony, data = rich_merulina_kusu)
dunn.test::dunn.test(rich_merulina_kusu$Shannon, rich_merulina_kusu$colony, method = "bonferroni")

kruskal.test(InvSimpson~macroalgal_treatment.type, data = rich_merulina_kusu)
kruskal.test(InvSimpson~colony, data = rich_merulina_kusu)

#For Raffles Pocillopora.
hist(rich_pocillopora_raffles$Observed, main = "Observed", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_raffles$Observed)

hist(rich_pocillopora_raffles$Chao1, main = "Chao1", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_raffles$Chao1)

hist(rich_pocillopora_raffles$Shannon, main = "Shannon", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_raffles$Shannon)

hist(rich_pocillopora_raffles$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_raffles$InvSimpson)

kruskal.test(Observed~macroalgal_treatment.type, data = rich_pocillopora_raffles)
kruskal.test(Observed~colony, data = rich_pocillopora_raffles)
dunn.test::dunn.test(rich_pocillopora_raffles$Observed, rich_pocillopora_raffles$colony, method = "bonferroni")

kruskal.test(Chao1~macroalgal_treatment.type, data = rich_pocillopora_raffles)
kruskal.test(Chao1~colony, data = rich_pocillopora_raffles)
dunn.test::dunn.test(rich_pocillopora_raffles$Chao1, rich_pocillopora_raffles$colony, method = "bonferroni")

kruskal.test(Shannon~macroalgal_treatment.type, data = rich_pocillopora_raffles)
kruskal.test(Shannon~colony, data = rich_pocillopora_raffles)
dunn.test::dunn.test(rich_pocillopora_raffles$Shannon, rich_pocillopora_raffles$colony, method = "bonferroni")

kruskal.test(InvSimpson~macroalgal_treatment.type, data = rich_pocillopora_raffles)
kruskal.test(InvSimpson~colony, data = rich_pocillopora_raffles)
dunn.test::dunn.test(rich_pocillopora_raffles$InvSimpson, rich_pocillopora_raffles$colony, method = "bonferroni")

#For Kusu Pocillopora:
hist(rich_pocillopora_kusu$Observed, main = "Observed", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_kusu$Observed)

hist(rich_pocillopora_kusu$Chao1, main = "Chao1", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_kusu$Chao1)

hist(rich_pocillopora_kusu$Shannon, main = "Shannon", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_kusu$Shannon)

hist(rich_pocillopora_kusu$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)
shapiro.test(rich_pocillopora_kusu$InvSimpson)

kruskal.test(Observed~macroalgal_treatment.type, data = rich_pocillopora_kusu)
kruskal.test(Observed~colony, data = rich_pocillopora_kusu)

kruskal.test(Chao1~macroalgal_treatment.type, data = rich_pocillopora_kusu)
kruskal.test(Chao1~colony, data = rich_pocillopora_kusu)

kruskal.test(Shannon~macroalgal_treatment.type, data = rich_pocillopora_kusu)
kruskal.test(Shannon~colony, data = rich_pocillopora_kusu)

kruskal.test(InvSimpson~macroalgal_treatment.type, data = rich_pocillopora_kusu)
kruskal.test(InvSimpson~colony, data = rich_pocillopora_kusu)

#For overall.
hist(rich_overall$Observed, main = "Observed", xlab = "", breaks = 10)
shapiro.test(rich_overall$Observed)

hist(rich_overall$Chao1, main = "Chao1", xlab = "", breaks = 10)
shapiro.test(rich_overall$Chao1)

hist(rich_overall$Shannon, main = "Shannon", xlab = "", breaks = 10)
shapiro.test(rich_overall$Shannon)

hist(rich_overall$InvSimpson, main = "InvSimpson", xlab = "", breaks = 10)
shapiro.test(rich_overall$InvSimpson)

kruskal.test(Observed~macroalgal_treatment.type, data = rich_overall)
kruskal.test(Observed~colony, data = rich_overall)
kruskal.test(Observed~location, data = rich_overall)
kruskal.test(Observed~sample_species, data = rich_overall)
dunn.test::dunn.test(rich_overall$Observed, rich_overall$sample_species, method = "bonferroni")

kruskal.test(Chao1~macroalgal_treatment.type, data = rich_overall)
kruskal.test(Chao1~colony, data = rich_overall)
kruskal.test(Chao1~location, data = rich_overall)
kruskal.test(Chao1~sample_species, data = rich_overall)
dunn.test::dunn.test(rich_overall$Chao1, rich_overall$sample_species, method = "bonferroni")

kruskal.test(Shannon~macroalgal_treatment.type, data = rich_overall)
kruskal.test(Shannon~colony, data = rich_overall)
dunn.test::dunn.test(rich_overall$Shannon, rich_overall$colony, method = "bonferroni")
kruskal.test(Shannon~location, data = rich_overall)
kruskal.test(Shannon~sample_species, data = rich_overall)
dunn.test::dunn.test(rich_overall$Shannon, rich_overall$sample_species, method = "bonferroni")

kruskal.test(InvSimpson~macroalgal_treatment.type, data = rich_overall)
kruskal.test(InvSimpson~colony, data = rich_overall)
dunn.test::dunn.test(rich_overall$InvSimpson, rich_overall$colony, method = "bonferroni")
kruskal.test(InvSimpson~location, data = rich_overall)
kruskal.test(InvSimpson~sample_species, data = rich_overall)
dunn.test::dunn.test(rich_overall$InvSimpson, rich_overall$sample_species, method = "bonferroni")

##########

#Plot NMDS.
#For overall (Based on species).
df_overall <- data.frame(sample_data(overall_ITS2_ps))
ordinate_overall = ordinate(overall_ITS2_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_overall <- vegan::scores(ordinate_overall)
data.scores_overall <- as.data.frame(scores(ordinate_overall$points))
data.scores_overall$macroalgal_treatment.type <- df_overall$macroalgal_treatment.type
data.scores_overall$location <- df_overall$location
data.scores_overall$sample_species <- df_overall$sample_species
data.scores_overall$colony <- df_overall$colony
head(data.scores_overall)

plot_nmds_overall <- ggplot(data.scores_overall, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_overall$sample_species), shape=factor(data.scores_overall$location)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_overall <- plot_nmds_overall + labs(shape = "Location", color = "Species")
plot_nmds_overall <- plot_nmds_overall + annotate("text", x = 2, y = -1.5, label = "Stress = 0.103")
p3 <- plot_nmds_overall + scale_shape_discrete(labels = c("Pulau Satumu", "Kusu Island"), limits = c("Pulau_Satumu", "Kusu_Island"))
p3 <- p3 + scale_color_discrete(labels = c(expression(italic("Merulina ampliata")), expression(italic("Pocillopora acuta"))), limits = c("Merulina_ampliata", "Pocillopora_acuta"))
p3 <- p3 + theme(legend.key = element_rect(fill = "white"))
p3 <- p3 + stat_ellipse(data = data.scores_overall, aes(x=MDS1, y=MDS2, color=sample_species))

#For Merulina:
df_merulina <- data.frame(sample_data(merulina_ITS2_ps))
ordinate_merulina = ordinate(merulina_ITS2_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_merulina <- vegan::scores(ordinate_merulina)
data.scores_merulina <- as.data.frame(scores(ordinate_merulina$points))
data.scores_merulina$macroalgal_treatment.type <- df_merulina$macroalgal_treatment.type
data.scores_merulina$location <- df_merulina$location
data.scores_merulina$colony <- df_merulina$colony
head(data.scores_merulina)

plot_nmds_merulina <- ggplot(data.scores_merulina, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_merulina$location), shape=factor(data.scores_merulina$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_merulina <- plot_nmds_merulina + labs(shape = "Contact Type", color = "Location")
plot_nmds_merulina <- plot_nmds_merulina + annotate("text", x = -0.5, y = -0.5, label = "Stress = 0.182")
p4 <- plot_nmds_merulina + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p4 <- p4 + scale_color_discrete(labels = c("Pulau Satumu", "Kusu Island"), limits = c("Pulau_Satumu", "Kusu_Island"))
p4 <- p4 + theme(legend.key = element_rect(fill = "white"))
p4 <- p4 + stat_ellipse(data = data.scores_merulina, aes(x=MDS1, y=MDS2, color=location))

#For Pocillopora.
df_pocillopora <- data.frame(sample_data(pocillopora_ITS2_ps))
ordinate_pocillopora = ordinate(pocillopora_ITS2_ps_sqrt, method = "NMDS", distance = "bray")
data.scores_pocillopora <- vegan::scores(ordinate_pocillopora)
data.scores_pocillopora <- as.data.frame(scores(ordinate_pocillopora$points))
data.scores_pocillopora$macroalgal_treatment.type <- df_pocillopora$macroalgal_treatment.type
data.scores_pocillopora$location <- df_pocillopora$location
data.scores_pocillopora$colony <- df_pocillopora$colony
head(data.scores_pocillopora)

plot_nmds_pocillopora <- ggplot(data.scores_pocillopora, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_pocillopora$location), shape=factor(data.scores_pocillopora$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_pocillopora <- plot_nmds_pocillopora + labs(shape = "Contact Type", color = "Location")
plot_nmds_pocillopora <- plot_nmds_pocillopora + annotate("text", x = 1, y = -1.5, label = "Stress = 0.115")
p5 <- plot_nmds_pocillopora + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p5 <- p5 + scale_color_discrete(labels = c("Pulau Satumu", "Kusu Island"), limits = c("Pulau_Satumu", "Kusu_Island"))
p5 <- p5 + theme(legend.key = element_rect(fill = "white"))
p5 <- p5 + stat_ellipse(data = data.scores_pocillopora, aes(x=MDS1, y=MDS2, color=location))

#Location-based NMDS.
#For Raffles Merulina.
ordinate_merulina_raffles = ordinate(merulina_raffles_ps_sqrt, method = "NMDS", distance = "bray")
df_merulina_raffles <- data.frame(sample_data(merulina_raffles_ps))
data.scores_merulina_raffles <- vegan::scores(ordinate_merulina_raffles)
data.scores_merulina_raffles <- as.data.frame(scores(ordinate_merulina_raffles$points))
data.scores_merulina_raffles$macroalgal_treatment.type <- df_merulina_raffles$macroalgal_treatment.type
data.scores_merulina_raffles$colony <- df_merulina_raffles$colony
head(data.scores_merulina_raffles)

plot_nmds_merulina_raffles <- ggplot(data.scores_merulina_raffles, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_merulina_raffles$colony), shape=factor(data.scores_merulina_raffles$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_merulina_raffles <- plot_nmds_merulina_raffles + labs(shape = "Contact Type", color = "Colony")
plot_nmds_merulina_raffles <- plot_nmds_merulina_raffles + annotate("text", x = 1, y = -1.0, label = "Stress = 0.126")
p6 <- plot_nmds_merulina_raffles + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p6 <- p6 + scale_color_discrete(labels = c("R1", "R2", "R3", "R4", "R5"), limits = c("R1", "R2", "R3", "R4", "R5"))
p6 <- p6 + theme(legend.key = element_rect(fill = "white"))
p6 <- p6 + ggtitle("Pulau Satumu")
p6 <- p6 + stat_ellipse(data = data.scores_merulina_raffles, aes(x=MDS1, y=MDS2, color=colony))

#For Kusu Merulina.
ordinate_merulina_kusu = ordinate(merulina_kusu_ps_sqrt, method = "NMDS", distance = "bray")
df_merulina_kusu <- data.frame(sample_data(merulina_kusu_ps))
data.scores_merulina_kusu <- vegan::scores(ordinate_merulina_kusu)
data.scores_merulina_kusu <- as.data.frame(scores(ordinate_merulina_kusu$points))
data.scores_merulina_kusu$macroalgal_treatment.type <- df_merulina_kusu$macroalgal_treatment.type
data.scores_merulina_kusu$colony <- df_merulina_kusu$colony
head(data.scores_merulina_kusu)

plot_nmds_merulina_kusu <- ggplot(data.scores_merulina_kusu, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_merulina_kusu$colony), shape=factor(data.scores_merulina_kusu$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_merulina_kusu <- plot_nmds_merulina_kusu + labs(shape = "Contact Type", color = "Colony")
plot_nmds_merulina_kusu <- plot_nmds_merulina_kusu + annotate("text", x = 0.5, y = -0.5, label = "Stress = 0.150")
p7 <- plot_nmds_merulina_kusu + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p7 <- p7 + scale_color_discrete(labels = c("K1", "K2", "K3", "K4", "K5"), limits = c("K1", "K2", "K3", "K4", "K5"))
p7 <- p7 + theme(legend.key = element_rect(fill = "white"))
p7 <- p7 + ggtitle("Kusu Island")
p7 <- p7 + stat_ellipse(data = data.scores_merulina_kusu, aes(x=MDS1, y=MDS2, color=colony))

#For Raffles Pocillopora.
ordinate_pocillopora_raffles = ordinate(pocillopora_raffles_ps_sqrt, method = "NMDS", distance = "bray")
df_pocillopora_raffles <- data.frame(sample_data(pocillopora_raffles_ps))
data.scores_pocillopora_raffles <- vegan::scores(ordinate_pocillopora_raffles)
data.scores_pocillopora_raffles <- as.data.frame(scores(ordinate_pocillopora_raffles$points))
data.scores_pocillopora_raffles$macroalgal_treatment.type <- df_pocillopora_raffles$macroalgal_treatment.type
data.scores_pocillopora_raffles$colony <- df_pocillopora_raffles$colony
head(data.scores_pocillopora_raffles)

plot_nmds_pocillopora_raffles <- ggplot(data.scores_pocillopora_raffles, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_pocillopora_raffles$colony), shape=factor(data.scores_pocillopora_raffles$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_pocillopora_raffles <- plot_nmds_pocillopora_raffles + labs(shape = "Contact Type", color = "Colony")
plot_nmds_pocillopora_raffles <- plot_nmds_pocillopora_raffles + annotate("text", x = 1, y = -1.5, label = "Stress = 0.0652")
p8 <- plot_nmds_pocillopora_raffles + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p8 <- p8 + scale_color_discrete(labels = c("R1", "R2", "R3", "R4", "R5"), limits = c("R1", "R2", "R3", "R4", "R5"))
p8 <- p8 + theme(legend.key = element_rect(fill = "white"))
p8 <- p8 + ggtitle("Pulau Satumu")
p8 <- p8 + stat_ellipse(data = data.scores_pocillopora_raffles, aes(x=MDS1, y=MDS2, color=colony))

#For Kusu Pocillopora.
ordinate_pocillopora_kusu = ordinate(pocillopora_kusu_ps_sqrt, method = "NMDS", distance = "bray")
df_pocillopora_kusu <- data.frame(sample_data(pocillopora_kusu_ps))
data.scores_pocillopora_kusu <- vegan::scores(ordinate_pocillopora_kusu)
data.scores_pocillopora_kusu <- as.data.frame(scores(ordinate_pocillopora_kusu$points))
data.scores_pocillopora_kusu$macroalgal_treatment.type <- df_pocillopora_kusu$macroalgal_treatment.type
data.scores_pocillopora_kusu$colony <- df_pocillopora_kusu$colony
head(data.scores_pocillopora_kusu)

plot_nmds_pocillopora_kusu <- ggplot(data.scores_pocillopora_kusu, aes(x=NMDS1, y=NMDS2)) + geom_point(aes(MDS1, MDS2, color=factor(data.scores_pocillopora_kusu$colony), shape=factor(data.scores_pocillopora_kusu$macroalgal_treatment.type)), size=2) + theme(panel.background = element_rect(fill = NA,colour = "black"), axis.line = element_line(colour = "black"))
plot_nmds_pocillopora_kusu <- plot_nmds_pocillopora_kusu + labs(shape = "Contact Type", color = "Colony")
plot_nmds_pocillopora_kusu <- plot_nmds_pocillopora_kusu + annotate("text", x = -1, y = -1, label = "Stress = 0.0924")
p9 <- plot_nmds_pocillopora_kusu + scale_shape_discrete(labels = c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"), limits = c("field_sample", "pre_treatment", "control", "direct", "5cm"))
p9 <- p9 + scale_color_discrete(labels = c("K1", "K2", "K3", "K4", "K5"), limits = c("K1", "K2", "K3", "K4", "K5"))
p9 <- p9 + theme(legend.key = element_rect(fill = "white"))
p9 <- p9 + ggtitle("Kusu Island")
p9 <- p9 + stat_ellipse(data = data.scores_pocillopora_kusu, aes(x=MDS1, y=MDS2, color=colony))

##########

#PERMANOVA/PERMDISP.
set.seed(1)
bray_overall <- phyloseq::distance(overall_ITS2_ps_sqrt, method = "bray")
adonis2(bray_overall~macroalgal_treatment.type*sample_species*location, strata = df_overall$colony, data = df_overall)
adonis2(bray_overall~macroalgal_treatment.type + sample_species + location, strata = df_overall$colony, data = df_overall)
beta_overall_1 <- betadisper(bray_overall, df_overall$sample_species)
permutest(beta_overall_1)
TukeyHSD(beta_overall_1)
beta_overall_2 <- betadisper(bray_overall, df_overall$location)
permutest(beta_overall_2)
beta_overall_3 <- betadisper(bray_overall, df_overall$macroalgal_treatment.type)
permutest(beta_overall_3)

#Pocillopora only.
bray_pocillopora <- phyloseq::distance(pocillopora_ITS2_ps_sqrt, method = "bray")
adonis2(bray_pocillopora~macroalgal_treatment.type*location, strata = df_pocillopora$colony, data = df_pocillopora)
adonis2(bray_pocillopora~macroalgal_treatment.type + location, strata = df_pocillopora$colony, data = df_pocillopora)
pairwise.adonis2(bray_pocillopora~macroalgal_treatment.type, data = df_pocillopora, strata = 'colony')
pairwise.adonis2(bray_pocillopora~location, data = df_pocillopora, strata = 'colony')
beta_poc_1 <- betadisper(bray_pocillopora, df_pocillopora$macroalgal_treatment.type)
permutest(beta_poc_1)
TukeyHSD(beta_poc_1)
beta_poc_2 <- betadisper(bray_pocillopora, df_pocillopora$location)
permutest(beta_poc_2)

#Merulina only.
bray_merulina <- phyloseq::distance(merulina_ITS2_ps_sqrt, method = "bray")
adonis2(bray_merulina~macroalgal_treatment.type*location, strata = df_merulina$colony, data = df_merulina)
adonis2(bray_merulina~macroalgal_treatment.type + location, strata = df_merulina$colony, data = df_merulina)
beta_mer_1 <- betadisper(bray_merulina, df_merulina$macroalgal_treatment.type)
permutest(beta_mer_1)
beta_mer_2 <- betadisper(bray_merulina, df_merulina$location)
permutest(beta_mer_2)
TukeyHSD(beta_mer_2)

##########

#Location-based PERMANOVA.
#Raffles Merulina only.
bray_merulina_raffles <- phyloseq::distance(merulina_raffles_ps_sqrt, method = "bray")
adonis2(bray_merulina_raffles~macroalgal_treatment.type, strata = df_merulina_raffles$colony, data = df_merulina_raffles)
beta_raffles_mer <- betadisper(bray_merulina_raffles, df_merulina_raffles$macroalgal_treatment.type)
permutest(beta_raffles_mer)
adonis2(bray_merulina_raffles~colony, data = df_merulina_raffles)
beta_raffles_mer_2 <- betadisper(bray_merulina_raffles, df_merulina_raffles$colony)
permutest(beta_raffles_mer_2)
pairwise.adonis2(bray_merulina_raffles~colony, data = df_merulina_raffles)

#Kusu Merulina only.
bray_merulina_kusu <- phyloseq::distance(merulina_kusu_ps_sqrt, method = "bray")
adonis2(bray_merulina_kusu~macroalgal_treatment.type, strata = df_merulina_kusu$colony, data = df_merulina_kusu)
beta_kusu_mer <- betadisper(bray_merulina_kusu, df_merulina_kusu$macroalgal_treatment.type)
permutest(beta_kusu_mer)
adonis2(bray_merulina_kusu~colony, data = df_merulina_kusu)
pairwise.adonis2(bray_merulina_kusu~colony, data = df_merulina_kusu)
beta_kusu_mer_2 <- betadisper(bray_merulina_kusu, df_merulina_kusu$colony)
permutest(beta_kusu_mer_2)

#Raffles Pocillopora only.
bray_pocillopora_raffles <- phyloseq::distance(pocillopora_raffles_ps_sqrt, method = "bray")
adonis2(bray_pocillopora_raffles~macroalgal_treatment.type, strata = df_pocillopora_raffles$colony, data = df_pocillopora_raffles)
beta_raffles_poc <- betadisper(bray_pocillopora_raffles, df_pocillopora_raffles$macroalgal_treatment.type)
permutest(beta_raffles_poc)
adonis2(bray_pocillopora_raffles~colony, data = df_pocillopora_raffles)
pairwise.adonis2(bray_pocillopora_raffles~colony, data = df_pocillopora_raffles)
beta_raffles_poc_2 <- betadisper(bray_pocillopora_raffles, df_pocillopora_raffles$colony)
permutest(beta_raffles_poc_2)

#Kusu Pocillopora only.
bray_pocillopora_kusu <- phyloseq::distance(pocillopora_kusu_ps_sqrt, method = "bray")
adonis2(bray_pocillopora_kusu~macroalgal_treatment.type, strata = df_pocillopora_kusu$colony, data = df_pocillopora_kusu)
pairwise.adonis2(bray_pocillopora_kusu~macroalgal_treatment.type, strata = 'colony', data = df_pocillopora_kusu)
beta_kusu_poc <- betadisper(bray_pocillopora_kusu, df_pocillopora_kusu$macroalgal_treatment.type)
permutest(beta_kusu_poc)
adonis2(bray_pocillopora_kusu~colony, data = df_pocillopora_kusu)
pairwise.adonis2(bray_pocillopora_kusu~colony, data = df_pocillopora_kusu)
beta_kusu_poc_2 <- betadisper(bray_pocillopora_kusu, df_pocillopora_kusu$colony)
permutest(beta_kusu_poc_2)

##########

#Linear mixed models for alpha-diversity.
library(lme4)
library(lmerTest)
library(emmeans)

#Merulina:
#Overall:
mod1 <- lmer(log(Observed)~macroalgal_treatment.type + location + (1|colony), data = rich_merulina)
summary(mod1)
anova(mod1)
emmeans(mod1, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1))

mod1a <- lmer(log(Chao1)~macroalgal_treatment.type + location + (1|colony), data = rich_merulina)
summary(mod1a)
anova(mod1a)
emmeans(mod1a, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod1a)
qqnorm(resid(mod1a))
qqline(resid(mod1a))

mod1b <- lmer(log(Shannon)~macroalgal_treatment.type + location + (1|colony), data = rich_merulina)
summary(mod1b)
anova(mod1b)
emmeans(mod1b, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod1b)
qqnorm(resid(mod1b))
qqline(resid(mod1b))

mod1c <- lmer(log(InvSimpson)~macroalgal_treatment.type + location + (1|colony), data = rich_merulina)
summary(mod1c)
anova(mod1c)
plot(mod1c)
qqnorm(resid(mod1c))
qqline(resid(mod1c))

#Raffles:
mod2 <- lm(log(Observed)~macroalgal_treatment.type + colony , data = rich_merulina_raffles)
summary(mod2)
anova(mod2)
emmeans(mod2, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2))

mod2a <- lm(log(Chao1)~macroalgal_treatment.type + colony , data = rich_merulina_raffles)
summary(mod2a)
anova(mod2a)
emmeans(mod2a, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod2a)
qqnorm(resid(mod2a))
qqline(resid(mod2a))

mod2b <- lm(log(Shannon)~macroalgal_treatment.type + colony, data = rich_merulina_raffles)
summary(mod2b)
anova(mod2b)
plot(mod2b)
qqnorm(resid(mod2b))
qqline(resid(mod2b))

mod2c <- lmer(log(InvSimpson)~macroalgal_treatment.type + (1|colony), data = rich_merulina_raffles)
summary(mod2c)
anova(mod2c)
plot(mod2c)
qqnorm(resid(mod2c))
qqline(resid(mod2c))

#Kusu:
mod3 <- lmer(log(Observed)~macroalgal_treatment.type + (1|colony), data = rich_merulina_kusu)
summary(mod3)
anova(mod3)
plot(mod3)
qqnorm(resid(mod3))
qqline(resid(mod3))

mod3a <- lmer(log(Chao1)~macroalgal_treatment.type + (1|colony), data = rich_merulina_kusu)
summary(mod3a)
anova(mod3a)
plot(mod3a)
qqnorm(resid(mod3a))
qqline(resid(mod3a))

mod3b <- lmer(log(Shannon)~macroalgal_treatment.type + (1|colony), data = rich_merulina_kusu)
summary(mod3b)
anova(mod3b)
plot(mod3b)
qqnorm(resid(mod3b))
qqline(resid(mod3b))

mod3c <- lmer(log(InvSimpson)~macroalgal_treatment.type + (1|colony), data = rich_merulina_kusu)
summary(mod3c)
anova(mod3c)
plot(mod3c)
qqnorm(resid(mod3c))
qqline(resid(mod3c))

#Pocillopora:
#Overall:
mod4 <- lmer(log(Observed)~macroalgal_treatment.type + location + (1|colony), data = rich_pocillopora)
summary(mod4)
anova(mod4)
emmeans(mod4, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod4)
qqnorm(resid(mod4))
qqline(resid(mod4))

mod4a <- lmer(log(Chao1)~macroalgal_treatment.type + location + (1|colony), data = rich_pocillopora)
summary(mod4a)
anova(mod4a)
emmeans(mod4a, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod4a)
qqnorm(resid(mod4a))
qqline(resid(mod4a))

mod4b <- lmer(log(Shannon)~macroalgal_treatment.type + location + (1|colony), data = rich_pocillopora)
summary(mod4b)
anova(mod4b)
plot(mod4b)
qqnorm(resid(mod4b))
qqline(resid(mod4b))

mod4c <- lmer(log(InvSimpson)~macroalgal_treatment.type + location + (1|colony), data = rich_pocillopora)
summary(mod4c)
anova(mod4c)
emmeans(mod4c, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod4c)
qqnorm(resid(mod4c))
qqline(resid(mod4c))

#Raffles:
mod5 <- lmer(log(Observed)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_raffles)
summary(mod5)
anova(mod5)
emmeans(mod5, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod5)
qqnorm(resid(mod5))
qqline(resid(mod5))

mod5a <- lmer(log(Chao1)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_raffles)
summary(mod5a)
anova(mod5a)
emmeans(mod5a, list(pairwise~macroalgal_treatment.type), adjust = "tukey", type = "response")
plot(mod5a)
qqnorm(resid(mod5a))
qqline(resid(mod5a))

mod5b <- lmer(log(Shannon)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_raffles)
summary(mod5b)
anova(mod5b)
plot(mod5b)
qqnorm(resid(mod5b))
qqline(resid(mod5b))

mod5c <- lmer(log(InvSimpson)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_raffles)
summary(mod5c)
anova(mod5c)
plot(mod5c)
qqnorm(resid(mod5c))
qqline(resid(mod5c))

#Kusu:
mod6 <- lm(log(Observed)~macroalgal_treatment.type + colony, data = rich_pocillopora_kusu)
summary(mod6)
anova(mod6)
plot(mod6)
qqnorm(resid(mod6))
qqline(resid(mod6))

mod6a <- lm(log(Chao1)~macroalgal_treatment.type + colony, data = rich_pocillopora_kusu)
summary(mod6a)
anova(mod6a)
plot(mod6a)
qqnorm(resid(mod6a))
qqline(resid(mod6a))

mod6b <- lmer(log(Shannon)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_kusu)
summary(mod6b)
anova(mod6b)
plot(mod6b)
qqnorm(resid(mod6b))
qqline(resid(mod6b))

mod6c <- lmer(log(InvSimpson)~macroalgal_treatment.type + (1|colony), data = rich_pocillopora_kusu)
summary(mod6c)
anova(mod6c)
plot(mod6c)
qqnorm(resid(mod6c))
qqline(resid(mod6c))

##########

#Overall relative abundance:
get_taxa_unique(overall_ITS2_ps_rarefied, "Genus")
sample_data(overall_ITS2_ps_rarefied)$sample_species <- factor(sample_data(overall_ITS2_ps_rarefied)$sample_species)
sample_data(overall_ITS2_ps_rarefied)$macroalgal_treatment.type <- factor(sample_data(overall_ITS2_ps_rarefied)$macroalgal_treatment.type)
sample_data(overall_ITS2_ps_rarefied)$colony <- factor(sample_data(overall_ITS2_ps_rarefied)$colony)
sample_data(overall_ITS2_ps_rarefied)$location <- factor(sample_data(overall_ITS2_ps_rarefied)$location)

genus_ITS2_combined <- overall_ITS2_ps_rarefied %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_ITS2_combined
head(otu_table(genus_ITS2_combined))
head(sample_data(genus_ITS2_combined))

genus_ITS2_abund <- genus_ITS2_combined %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
head(genus_ITS2_abund)

all_ITS2_combined <- genus_ITS2_abund %>% 
  select(Genus, Abundance, sample_species, location, macroalgal_treatment.type, colony) %>%
  filter(Abundance != 0) %>%
  mutate(Genus = as.character(Genus))
head(all_ITS2_combined)

genus_ITS2_combined_plot <- all_ITS2_combined %>% 
  select(Genus, Abundance, sample_species, location, macroalgal_treatment.type, colony) %>%
  group_by(sample_species, location, macroalgal_treatment.type) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(sample_species, location, macroalgal_treatment.type, Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    totalSum,
    Real_ab = Abundance/totalSum) %>%
  unique()
head(genus_ITS2_combined_plot)
max(genus_ITS2_combined_plot$Real_ab)
length(unique(genus_ITS2_combined_plot$Genus))

plot_genus_ITS2_all <- ggplot(genus_ITS2_combined_plot) + geom_col(mapping = aes(x = macroalgal_treatment.type, y = Real_ab, fill = Genus), position = "stack", show.legend = TRUE) + 
  facet_grid(vars(location), vars(sample_species)) + 
  ylab("Relative Abundance") +
  xlab(NULL) + 
  theme_bw()

genus_ITS2_plot_extract_mod <- genus_ITS2_combined_plot
genus_ITS2_plot_extract_mod$location <- factor(genus_ITS2_plot_extract_mod$location, labels = c(Pulau_Satumu = "Pulau Satumu", Kusu_Island = "Kusu Island"))
genus_ITS2_plot_extract_mod$sample_species <- factor(genus_ITS2_plot_extract_mod$sample_species, labels = c(Merulina_ampliata = "Merulina ampliata", Pocillopora_acuta = "Pocillopora acuta"))

plot_genus_ITS2_all <- ggplot(genus_ITS2_plot_extract_mod) + geom_col(mapping = aes(x = macroalgal_treatment.type, y = Real_ab, fill = Genus), position = "stack", show.legend = TRUE) + 
  facet_grid(vars(location), vars(sample_species)) + 
  ylab("Relative Abundance") +
  xlab(NULL) +
  theme_bw()

p_genus_ITS2_all <- plot_genus_ITS2_all
p_genus_ITS2_all <- p_genus_ITS2_all + theme(strip.text.x.top = element_text(face = "italic"))
p_genus_ITS2_all <- p_genus_ITS2_all + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
p_genus_ITS2_all <- p_genus_ITS2_all + theme(strip.background = element_blank())

newTreatment <- c("field_sample", "pre_treatment", "control", "direct", "5cm")
p_genus_ITS2_all$data$macroalgal_treatment.type <- as.character(p_genus_ITS2_all$data$macroalgal_treatment.type)
p_genus_ITS2_all$data$macroalgal_treatment.type <- factor(p_genus_ITS2_all$data$macroalgal_treatment.type, levels = newTreatment)
print(p_genus_ITS2_all)
p_genus_ITS2_all <- p_genus_ITS2_all + scale_x_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))

#Split by colony:
#Pocillopora:
get_taxa_unique(merulina_ITS2_ps_rarefied, "Genus")
sample_data(merulina_ITS2_ps_rarefied)$macroalgal_treatment.type <- factor(sample_data(merulina_ITS2_ps_rarefied)$macroalgal_treatment.type)
sample_data(merulina_ITS2_ps_rarefied)$colony <- factor(sample_data(merulina_ITS2_ps_rarefied)$colony)
sample_data(merulina_ITS2_ps_rarefied)$location <- factor(sample_data(merulina_ITS2_ps_rarefied)$location)

genus_mer_ITS2_combined <- merulina_ITS2_ps_rarefied %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_mer_ITS2_combined
head(otu_table(genus_mer_ITS2_combined))
head(sample_data(genus_mer_ITS2_combined))

genus_mer_ITS2_abund <- genus_mer_ITS2_combined %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
head(genus_mer_ITS2_abund)

all_mer_ITS2 <- genus_mer_ITS2_abund %>% 
  select(Genus, Abundance, sample_species, location, macroalgal_treatment.type, colony) %>%
  filter(Abundance != 0) %>%
  mutate(Genus = as.character(Genus))
head(all_mer_ITS2)

genus_ITS2_mer_plot <- all_mer_ITS2 %>% 
  select(Genus, Abundance, location, macroalgal_treatment.type, colony) %>%
  group_by(colony, location, macroalgal_treatment.type) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(colony, location, macroalgal_treatment.type, Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    totalSum,
    Real_ab = Abundance/totalSum) %>%
  unique()
head(genus_ITS2_mer_plot)
max(genus_ITS2_mer_plot$Real_ab)
length(unique(genus_ITS2_mer_plot$Genus))

plot_genus_mer_ITS2_all <- ggplot(genus_ITS2_mer_plot) + geom_col(mapping = aes(x = macroalgal_treatment.type, y = Real_ab, fill = Genus), position = "stack", show.legend = TRUE) + 
  facet_wrap(vars(colony), ncol = 5) + 
  ylab("Relative Abundance") +
  xlab(NULL) + 
  theme_bw()

newTreatment <- c("field_sample", "pre_treatment", "control", "direct", "5cm")
plot_genus_mer_ITS2_all$data$macroalgal_treatment.type <- as.character(plot_genus_mer_ITS2_all$data$macroalgal_treatment.type)
plot_genus_mer_ITS2_all$data$macroalgal_treatment.type <- factor(plot_genus_mer_ITS2_all$data$macroalgal_treatment.type, levels = newTreatment)
print(plot_genus_mer_ITS2_all)
plot_genus_mer_ITS2_all <- plot_genus_mer_ITS2_all + scale_x_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_genus_mer_ITS2_all <- plot_genus_mer_ITS2_all + labs(title = "Merulina ampliata")
plot_genus_mer_ITS2_all <- plot_genus_mer_ITS2_all + theme(plot.title = element_text(face = "italic"))
plot_genus_mer_ITS2_all <- plot_genus_mer_ITS2_all + theme(legend.position = "bottom")

#Pocillopora:
get_taxa_unique(pocillopora_ITS2_ps_rarefied, "Genus")
sample_data(pocillopora_ITS2_ps_rarefied)$macroalgal_treatment.type <- factor(sample_data(pocillopora_ITS2_ps_rarefied)$macroalgal_treatment.type)
sample_data(pocillopora_ITS2_ps_rarefied)$colony <- factor(sample_data(pocillopora_ITS2_ps_rarefied)$colony)
sample_data(pocillopora_ITS2_ps_rarefied)$location <- factor(sample_data(pocillopora_ITS2_ps_rarefied)$location)

genus_poc_ITS2_combined <- pocillopora_ITS2_ps_rarefied %>% aggregate_taxa(level = "Genus") %>% microbiome::transform(transform = "compositional")
genus_poc_ITS2_combined
head(otu_table(genus_poc_ITS2_combined))
head(sample_data(genus_poc_ITS2_combined))

genus_poc_ITS2_abund <- genus_poc_ITS2_combined %>% 
  tax_glom(taxrank = "Genus") %>%
  transform_sample_counts(function(x) {x/sum(x)}) %>%
  psmelt()
head(genus_poc_ITS2_abund)

all_poc_ITS2 <- genus_poc_ITS2_abund %>% 
  select(Genus, Abundance, sample_species, location, macroalgal_treatment.type, colony) %>%
  filter(Abundance != 0) %>%
  mutate(Genus = as.character(Genus))
head(all_poc_ITS2)

genus_ITS2_poc_plot <- all_poc_ITS2 %>% 
  select(Genus, Abundance, location, macroalgal_treatment.type, colony) %>%
  group_by(colony, location, macroalgal_treatment.type) %>%
  mutate(totalSum = sum(Abundance)) %>%
  ungroup() %>%
  group_by(colony, location, macroalgal_treatment.type, Genus) %>%
  summarise(
    Abundance = sum(Abundance),
    totalSum,
    Real_ab = Abundance/totalSum) %>%
  unique()
head(genus_ITS2_poc_plot)
max(genus_ITS2_poc_plot$Real_ab)
length(unique(genus_ITS2_poc_plot$Genus))

plot_genus_poc_ITS2_all <- ggplot(genus_ITS2_poc_plot) + geom_col(mapping = aes(x = macroalgal_treatment.type, y = Real_ab, fill = Genus), position = "stack", show.legend = TRUE) + 
  facet_wrap(vars(colony), ncol = 5) + 
  ylab("Relative Abundance") +
  xlab(NULL) + 
  theme_bw()

newTreatment <- c("field_sample", "pre_treatment", "control", "direct", "5cm")
plot_genus_poc_ITS2_all$data$macroalgal_treatment.type <- as.character(plot_genus_poc_ITS2_all$data$macroalgal_treatment.type)
plot_genus_poc_ITS2_all$data$macroalgal_treatment.type <- factor(plot_genus_poc_ITS2_all$data$macroalgal_treatment.type, levels = newTreatment)
print(plot_genus_poc_ITS2_all)
plot_genus_poc_ITS2_all <- plot_genus_poc_ITS2_all + scale_x_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated")) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
plot_genus_poc_ITS2_all <- plot_genus_poc_ITS2_all + labs(title = "Pocillopora acuta")
plot_genus_poc_ITS2_all <- plot_genus_poc_ITS2_all + theme(plot.title = element_text(face = "italic"))
plot_genus_poc_ITS2_all <- plot_genus_poc_ITS2_all + theme(legend.position = "bottom")

##########

#Combined alpha-diversity figure:
fig_ITS2_rich <- plot_richness(overall_ITS2_ps_rarefied, x = "location", color = "macroalgal_treatment.type", measures = c("Observed", "Chao1", "Shannon", "InvSimpson")) + geom_boxplot()
fig_ITS2_rich <- fig_ITS2_rich + facet_grid(variable ~ sample_species, scales = "free_y")
fig_ITS2_rich <- fig_ITS2_rich + theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
fig_ITS2_rich <- fig_ITS2_rich + theme(strip.background = element_blank())
fig_ITS2_rich <- fig_ITS2_rich + labs(x="", y="Alpha Diversity Index\n")
fig_ITS2_rich <- fig_ITS2_rich + theme(axis.text.x = element_text(angle = 45, hjust = 1))
fig_ITS2_rich <- fig_ITS2_rich + scale_x_discrete(limits=c("Pulau_Satumu", "Kusu_Island"), labels=c("Pulau Satumu", "Kusu Island"))

newTreatment <- c("field_sample", "pre_experiment", "control", "direct", "5cm")
fig_ITS2_rich <- fig_ITS2_rich
fig_ITS2_rich$data$macroalgal_treatment.type <- as.character(fig_ITS2_rich$data$macroalgal_treatment.type)
fig_ITS2_rich$data$macroalgal_treatment.type <- factor(fig_ITS2_rich$data$macroalgal_treatment.type, levels = newTreatment)
print(fig_ITS2_rich)

fig_ITS2_rich <- fig_ITS2_rich + scale_color_discrete(labels= c("Field sample", "Before experiment", "Control", "Direct", "Water-mediated"))
fig_ITS2_rich <- fig_ITS2_rich + theme(legend.position = "bottom", legend.direction = "horizontal")
fig_ITS2_rich <- fig_ITS2_rich + theme(legend.title = element_blank())

labels <- c(Pocillopora_acuta = "Pocillopora acuta", Merulina_ampliata = "Merulina ampliata")
fig_ITS2_rich <- fig_ITS2_rich + facet_grid(variable ~ sample_species, labeller = labeller(sample_species = labels), scales = "free_y")
fig_ITS2_rich <- fig_ITS2_rich + theme(strip.text.x.top = element_text(face = "italic"))
