# set work directory
setwd('e:/permafrost')
wd_its <- file.path(getwd(),"data/ITS")
if (!dir.exists(wd_its)) {
  dir.create(wd_its)
}
# set saving directory
save.dir <- file.path(getwd(),"result/ITS")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}

# loading packages
library(phyloseq)
library(ape)
library(Biostrings)
# read data
## metadata
metadata <- read.delim(file.path(wd_its, "./metadata.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$sample_id)
## dna sequences
dna_seqs <- readDNAStringSet(file.path(wd_its, "./otus.fa"), format = "fasta", nrec = -1L, 
                             skip = 0L, seek.first.rec = FALSE,
                             use.names = TRUE)

## tree
tree <- read_tree(file.path(wd_its, "./otus.nwk"))
## otu table
otu <- read.delim(file.path(wd_its, "./otutab.txt"), header = T, row.names = 1, sep = "\t")
otu = otu[, metadata$sample_id[metadata$sample_id %in% colnames(otu)]]
## rarefy out table
otu_rare <- read.delim(file.path(wd_its, "./otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_rare = otu_rare[, metadata$sample_id[metadata$sample_id %in% colnames(otu_rare)]]
## tax
tax <- read.delim(file.path(wd_its, "./taxonomy.txt"), header = T, row.names = 1, sep = "\t")
tax <- as.matrix(tax)
# phyloseq object
otu <- otu_table(otu, taxa_are_rows = TRUE)
otu_rare <- otu_table(otu_rare, taxa_are_rows = TRUE)
tax <- tax_table(tax)
meta_dat <- sample_data(metadata)
phylo <- phyloseq(otu, tax, meta_dat, tree, dna_seqs)
phylo_rare <- phyloseq(otu_rare, tax, meta_dat, tree, dna_seqs)
# phylo_SUR <- subset_samples(phylo, layer == 'SUR') 
# phylo_SUB <- subset_samples(phylo, layer == 'SUB') 
# phylo_PL <- subset_samples(phylo, layer == 'PL') 

# alpha diversity
library(phyloseq)
library(picante)
library(microbiome)
library(tidyverse)
library(ggpubr)
library(ggplot2)
alpha_div <- estimate_richness(phylo_rare, measures = c("Observed", "Chao1", 'Shannon', 'Simpson'))
pd <- pd(t(otu), tree, include.root = F)
alpha_div <- cbind(Layers =metadata$layer, alpha_div, Faith = pd$PD,
                   Evenness = alpha_div$Shannon/log(alpha_div$Observed)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
# mpd <- picante::mpd(t(otu), cophenetic(tree), abundance.weighted = T)
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
p <- ggplot(alpha_div, aes(x = Layers, y = Evenness)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Evenness', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())
ggsave(file.path(save.dir, "./figs/alpha/alpha_boxplot_evenness.pdf"), 
       p, width = 89, height = 89, units = "mm")


p <- ggplot(alpha_div, aes(x = Layers, y = Faith)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Faith index', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())
ggsave(file.path(save.dir, "./figs/alpha/alpha_boxplot_Faith.pdf"),
       p, width = 89, height = 89, units = "mm")

# beta diversity
library(vegan)
library(pairwiseAdonis)
## taxonomy
### PERMANOVA test
tax_dist <- as.matrix(vegdist(t(otu), "bray" ))
adonis2(tax_dist~layer, data = metadata)
set.seed(123)
pairwise.adonis(tax_dist, metadata$layer)
### PCoA plot with bray-curties as distance
ord.tax <-  cmdscale(tax_dist,  k = 2, eig = T, add = T)
round(ord.tax$eig*100/sum(ord.tax$eig),1)[c(1,2)]
pcoa_tax_plot <- data.frame(Layers = metadata$layer, scores(ord.tax)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Dim1, y = Dim2, shape = Layers, color = Layers)) + 
  geom_point(size = 1, alpha = 0.8) + 
  stat_ellipse(geom = "polygon", aes(fill = Layers), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  labs(x=paste("PCoA1 (", format(100 * ord.tax$eig[1] / sum(ord.tax$eig), digits = 3), "%)", sep = ""),
       y=paste("PCoA2 (", format(100 * ord.tax$eig[2] / sum(ord.tax$eig), digits = 3), "%)", sep = "")) +
  theme(axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

ggsave(file.path(save.dir, "./figs/beta/PCoA_tax_bray.pdf"),
       pcoa_tax_plot, width = 89, height = 59, units = "mm")
### difference in taxonomic variance among layers
beta_tax_plot <- sapply(unique(metadata$layer), function(x) usedist::dist_subset(tax_dist, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Layers, y = distance)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Taxonomic variance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

ggsave(file.path(save.dir, "./figs/beta/tax_variance.pdf"),
       beta_tax_plot, width = 89, height = 89, units = "mm")

## phylogeny
### PERMANOVA test
phy_dist <- as.matrix(UniFrac(phylo, weighted = TRUE, normalized = TRUE, parallel = T, fast = TRUE))
adonis2(phy_dist~layer, data = metadata)
set.seed(123)
pairwise.adonis(phy_dist, metadata$layer)
### PCoA plot with unweighted UniFrac as distance
ord.phy <- ordinate(phylo, method = "PCoA", distance = "unifrac", weighted = TRUE)
PCoA_phy_plot <- plot_ordination(phylo, ord.phy, color = "layer") + theme(aspect.ratio=1)
PCoA_unifrac_plot <- data.frame(Layers = metadata$layer, ord.phy$vectors[, 1:2]) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Axis.1, y = Axis.2, shape = Layers, color = Layers)) + 
  geom_point(size = 1, alpha = 0.8) + 
  stat_ellipse(geom = "polygon", aes(fill = Layers), alpha = 0.2, show.legend = FALSE, level = 0.95) +
  labs(x=paste("PCoA1 (", format(100 * ord.phy$values[1, 2], digits = 3), "%)", sep = ""),
       y=paste("PCoA2 (", format(100 * ord.phy$values[2, 2], digits = 3), "%)", sep = "")) +
  theme(axis.title = element_text(size = 8, colour = "black"),
        axis.text = element_text(size = 6, colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

ggsave(file.path(save.dir, "./figs/beta/PCoA_unifrac_plot.pdf"),
       PCoA_unifrac_plot, width = 89, height = 59, units = "mm")
### determine the betaMNTD
#### Method 1
require(picante)
beta.mntd.weighted <- as.matrix(comdistnt(otu), cophenetic(tree), abundance.weighted = T)

#### Method 2 (time efficiency)
p.dist.mat <- cophenetic(tree)
get.presents <- function(x) {
  names(x[x > 0])
}
list.of.names <- apply(t(otu), 1, get.presents)
Dnn.apply.funtion <- function(x) {
  tmp.function <- function(z) {
    mean(c(apply(p.dist.mat[x, z], MARGIN = 1, min, na.rm = T), 
           apply(p.dist.mat[x, z], MARGIN = 2, min, na.rm = T)))
  }
  lapply(list.of.names, FUN = tmp.function)
}
beta.mntd.weighted <- do.call(cbind, lapply(list.of.names, Dnn.apply.funtion))
write.csv(beta.mntd.weighted, file.path(save.dir, './tables/beta.mntd.weighted.csv'))

### difference in MNTD among layers
beta_MNTD_plot <- sapply(unique(metadata$layer), function(x) usedist::dist_subset(beta.mntd.weighted, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Layers, y = distance)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'beta-MNTD', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

ggsave(file.path(save.dir, "./figs/beta/beta_mntd.pdf"),
       beta_MNTD_plot, width = 89, height = 89, units = "mm")

# composition
## determine the Order compositions within top 10 phylums
subphylo <- tax_glom(phylo, 'Phylum')
subphylo.rel <- microbiome::transform(subphylo, "compositional")
ra.tab <- otu_table(subphylo.rel)
subtaxa_tab <- tax_table(subphylo.rel)[, 2]
subtaxa_tab <- data.frame(subtaxa_tab, ra.tab) %>% group_by(Phylum) %>% 
  summarise(across(everything(), sum)) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(15, MRA) %>%
  select(., -c('MRA')) %>% 
  mutate(Phylum = factor(Phylum, levels = Phylum)) %>%
  tidyr::pivot_longer(cols = -c(Phylum), names_to = "sample_id", values_to = 'rel_abun') %>%
  mutate(layer = sapply(stringr::str_split(sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL')))

## boxplot shows the community composition
p <- ggplot(subtaxa_tab, aes(Phylum, rel_abun)) + 
  geom_boxplot(width = 0.5, aes(fill = layer)) +
  stat_compare_means(aes(group = layer), label = "p.signif") +
  #scale_fill_manual(values= cols) +
  labs(x = 'Phylum', y = 'Relative abundance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 8, colour = "black"),
        axis.text.x = element_text(size = 7, colour = "black", 
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 8, colour = "black"),
        legend.title = element_text(size = 8),
        legend.text = element_text(size = 7),
        panel.grid = element_blank())
ggsave(file.path(save.dir, "./figs/composition/tax_composition.pdf"), p, width = 189, height = 89, units = "mm")  
