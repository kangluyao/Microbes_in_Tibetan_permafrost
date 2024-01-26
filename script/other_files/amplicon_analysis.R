#################################################################
# alpha diversity
library(phyloseq)
library(picante)
library(microbiome)
library(tidyverse, quietly = TRUE, warn.conflicts = F)
library(ggpubr)
library(ggplot2)
alpha_div <- estimate_richness(phylo_rare, measures = c("Observed", "Chao1", 'Shannon', 'Simpson'))
pd <- pd(t(otu), tree, include.root = F)
alpha_div <- cbind(Layers = metadata$Layer, alpha_div, Faith = pd$PD,
                   Evenness = alpha_div$Shannon/log(alpha_div$Observed)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
# mpd <- picante::mpd(t(otu), cophenetic(tree), abundance.weighted = T)
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
p <- ggplot(alpha_div, aes(x = Layers, y = Evenness)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons, paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Evenness', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())
# ggsave(file.path(save.dir, "./figs/alpha/alpha_boxplot_evenness.pdf"), 
#        p, width = 89, height = 89, units = "mm")


p <- ggplot(alpha_div, aes(x = Layers, y = Faith)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons,  paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Faith index', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())
# ggsave(file.path(save.dir, "./figs/alpha/alpha_boxplot_Faith.pdf"),
#        p, width = 89, height = 89, units = "mm")
##########################################################################
# beta diversity
library(vegan)
library(pairwiseAdonis)
## taxonomy
### PERMANOVA test
tax_dist <- as.matrix(vegdist(t(otu), "bray" ))
adonis2(tax_dist ~ Layer, data = metadata)
set.seed(123)
pairwise.adonis(tax_dist, metadata$Layer)
### PCoA plot with bray-curties as distance
ord.tax <-  cmdscale(tax_dist,  k = 2, eig = T, add = T)
round(ord.tax$eig*100/sum(ord.tax$eig),1)[c(1,2)]
pcoa_tax_plot <- data.frame(Layers = metadata$Layer, scores(ord.tax)) %>%
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

# ggsave(file.path(save.dir, "./figs/beta/PCoA_tax_bray.pdf"),
#        pcoa_tax_plot, width = 89, height = 59, units = "mm")
### difference in taxonomic variance among layers
beta_tax_plot <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(tax_dist, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Layers, y = distance)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons,  paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Taxonomic variance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

# ggsave(file.path(save.dir, "./figs/beta/tax_variance.pdf"),
#        beta_tax_plot, width = 89, height = 89, units = "mm")

## phylogeny
### PERMANOVA test
phy_dist <- as.matrix(UniFrac(phylo, weighted = TRUE, normalized = TRUE, parallel = T, fast = TRUE))
adonis2(phy_dist ~ Layer, data = metadata)
set.seed(123)
pairwise.adonis(phy_dist, metadata$Layer)
### PCoA plot with unweighted UniFrac as distance
ord.phy <- ordinate(phylo, method = "PCoA", distance = "unifrac", weighted = TRUE)
PCoA_phy_plot <- plot_ordination(phylo, ord.phy, color = "layer") + theme(aspect.ratio = 1)
PCoA_unifrac_plot <- data.frame(Layers = metadata$Layer, ord.phy$vectors[, 1:2]) %>%
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

# ggsave(file.path(save.dir, "./figs/beta/PCoA_unifrac_plot.pdf"),
#        PCoA_unifrac_plot, width = 89, height = 59, units = "mm")
### determine the betaMNTD
#### Method 1
require(picante)
beta.mntd.weighted <- as.matrix(comdistnt(t(data.frame(otu)), cophenetic(tree), abundance.weighted = T))

#### Method 2 (time efficiency)
library(doParallel)  
cores <- detectCores() - 2  

p.dist.mat <- cophenetic(tree)
get.presents <- function(x) {
  names(x[x > 0])
}
list.of.names <- apply(otu, 2, get.presents)
Dnn.apply.funtion <- function(x) {
  tmp.function <- function(z) {
    mean(c(apply(p.dist.mat[x, z], MARGIN = 1, min, na.rm = T), 
         apply(p.dist.mat[x, z], MARGIN = 2, min, na.rm = T)))
  }
  lapply(list.of.names, FUN = tmp.function)
}
beta.mntd.weighted <- do.call(cbind, lapply(list.of.names, Dnn.apply.funtion))

beta.mntd.weighted <- read.table(file.path(save.dir, './tables/beta.mntd.weighted.txt'),
                                 header = T, row.names = 1)
### difference in MNTD among layers
beta_MNTD_plot <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(beta.mntd.weighted, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Layers, y = distance)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons,  paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'beta-MNTD', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

# ggsave(file.path(save.dir, "./figs/beta/beta_mntd.pdf"),
#        beta_MNTD_plot, width = 89, height = 89, units = "mm")

##########################################################################
## functional analysis
### PERMANOVA test
fun_dist <- as.matrix(vegdist(t(ko_tpm_table), "bray" ))
adonis2(fun_dist ~ Layer, data = metadata)
set.seed(123)
pairwise.adonis(fun_dist, metadata$Layer)
### PCoA plot
ord.fun <-  cmdscale(fun_dist,  k = 2, eig = T, add = T)
round(ord.fun$eig*100/sum(ord.fun$eig),1)[c(1,2)]
library(ggplot2)
pcoa_fun_plot <- data.frame(Layers = metadata$Layer, scores(ord.fun)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Dim1, y = Dim2, shape = Layers, color = Layers)) + 
  geom_point(size = 1.5, alpha = 0.8) + 
  labs(x=paste("PCoA1 (", format(100 * ord.fun$eig[1] / sum(ord.fun$eig), digits = 3), "%)", sep = ""),
       y=paste("PCoA2 (", format(100 * ord.fun$eig[2] / sum(ord.fun$eig), digits = 3), "%)", sep = "")) +
  theme(axis.title = element_text(size = 6, colour = "black"),
        axis.text = element_text(size = 5, colour = "black"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 5),
        legend.key = element_blank(),
        panel.grid = element_blank(),
        panel.background = element_blank(), 
        panel.border = element_rect(fill = NA, colour = "black"))

# ggsave(file.path(save.dir, "./figs/beta/PCoA_fun_bray.pdf"),
#        pcoa_fun_plot, width = 89, height = 59, units = "mm")
### difference in functional variance among layers
beta_fun_plot <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(fun_dist, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = Layers, y = distance)) + 
  geom_boxplot(width = 0.5, aes(fill = Layers))+
  stat_compare_means(comparisons = my_comparisons,  paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Functional variance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())

# ggsave(file.path(save.dir, "./figs/beta/fun_variance.pdf"),
#        beta_fun_plot, width = 89, height = 89, units = "mm")

# couple or decouple effect between taxonomic, phylogenetic and functional profile
pcoa1_tax <- data.frame(Layers = metadata$Layer, scores(ord.tax)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
pcoa1_phy <- data.frame(Layers = metadata$Layer, ord.phy$vectors[, 1:2]) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
pcoa1_fun <- data.frame(Layers = metadata$Layer, scores(ord.fun)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))

pcoa1_dat <- data.frame(pcoa1_tax[, 1:2], pcoa1_phy[,2], pcoa1_fun[,2])
colnames(pcoa1_dat) <- c("Layers", "Taxnomic_PCoA1", "Phylogenetic_PCoA1", "Functional_PCoA1")

install.packages('ggpmisc')
library(ggpmisc)

p_linear <- ggplot(pcoa1_dat, aes(x = Taxnomic_PCoA1, y = Functional_PCoA1, fill = Layers)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = Layers)) +
  # geom_smooth(method="lm", size=1, se=T, colour='black') +
  stat_poly_line(colour='black') +
  stat_poly_eq(aes(label = paste(after_stat(eq.label),
                                 after_stat(rr.label), sep = "*\", \"*"), colour = Layers)) +
  # scale_color_manual(values = c('#1b9e77', '#d95f02', "red")) +
  ylab("Taxnomic PCoA1")+xlab("Functional PCoA1") +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color = 'black', size = 6),
        axis.ticks.length = unit(0.2, "lines"), axis.ticks = element_line(color = 'black'),
        axis.line = element_blank(), 
        axis.text = element_text(colour = 'black',size = 5),
        strip.text = element_text(size = 6))
p_linear

dist_tax <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(tax_dist, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
dist_phy <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(beta.mntd.weighted, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
dist_fun <- sapply(unique(metadata$Layer), function(x) usedist::dist_subset(fun_dist, grep(x, metadata$sample_id, value = T))) %>%
  data.frame() %>% gather("Layers", "distance") %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))

dist_dat <- data.frame(dist_tax[, 1:2], dist_phy[,2], dist_fun[,2])
colnames(dist_dat) <- c("Layers", "Taxnomic_distance", "Phylogenetic_distance", "Functional_distance")
p_linear <- ggplot(dist_dat, aes(x = Phylogenetic_distance, y = Functional_distance, fill = Layers)) +
  geom_point(size=3.5, alpha=0.8, aes(colour = Layers)) +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  scale_color_manual(values = c('#1b9e77', '#d95f02', "red")) +
  # facet_wrap( .~ Layers, scales = "free_x", ncol = 2) +
  ylab("Taxnomic PCoA1")+xlab("Phylogenetic PCoA1") +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')

##################################################
# composition
## determine the compositions within top 10 phylums
subphylo <- tax_glom(phylo, 'Phylum')
subphylo.rel  = transform_sample_counts(subphylo, function(x) x / sum(x) )
ntaxa(subphylo.rel)
ra.tab <- otu_table(subphylo.rel)
sum(ra.tab[, 1])
subtaxa_tab <- tax_table(subphylo.rel)[, 2]
## boxplot shows the community composition
box_plot <- data.frame(subtaxa_tab, ra.tab) %>% group_by(Phylum) %>% 
  summarise(across(everything(), sum)) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(10, MRA) %>%
  select(., -c('MRA')) %>% 
  mutate(Phylum = factor(Phylum, levels = Phylum)) %>%
  tidyr::pivot_longer(cols = -c(Phylum), names_to = "sample_id", values_to = 'rel_abun') %>%
  mutate(layer = sapply(stringr::str_split(sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(Phylum, rel_abun*100)) + 
  geom_boxplot(width = 0.5, aes(fill = layer)) +
  stat_compare_means(aes(group = layer),  paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols) +
  labs(x = 'Phylum', y = 'Relative abundance (%)', fill='Layers') +
  theme_bw() +
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text.x = element_text(size = 10, colour = "black", 
                                   angle = 45, vjust = 1, hjust = 1),
        axis.text.y = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 10),
        legend.text = element_text(size = 8),
        panel.grid = element_blank())

## barplot shows the community composition
bar_plot <- data.frame(subtaxa_tab, ra.tab) %>% 
  mutate(MRA = rowMeans(select(., colnames(ra.tab)))) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(10, MRA) %>%
  select(., -c('MRA')) %>% 
  bind_rows(summarise_all(., ~if(is.numeric(.)) 1-sum(.) else "Other")) %>%
  mutate(Phylum = factor(Phylum, levels = Phylum)) %>%
  pivot_longer(cols = -c(Phylum), names_to = "sample_id", values_to = 'rel_abun') %>%
  mutate(layer = sapply(stringr::str_split(sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = sample_id, y = 100*rel_abun, fill = Phylum))+
  geom_bar(stat = "identity") +
  scale_fill_brewer(palette = "Paired")+
  scale_y_continuous(expand = c(0,0))+
  labs(x = 'Sample', y = 'Mean relative abundance (%)') +
  theme_linedraw() + 
  guides(fill = guide_legend(keywidth = 0.5, keyheight = 0.5)) +
  facet_grid(~layer, scales = "free_x", space = "free_x") + 
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size = 6, colour = "black"),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 4), 
        axis.text.y = element_text(size = 5), 
        panel.spacing = unit(0, "lines"),
        strip.text = element_blank(),
        strip.background = element_blank(),
        legend.title = element_text(size = 6), 
        legend.text = element_text(size = 5))

ggsave(file.path(save.dir, "./figs/composition/tax_composition_10.pdf"), box_plot, width = 189, height = 95, units = "mm")  
ggsave(file.path(save.dir, "./figs/composition/tax_composition_bar.pdf"), bar_plot, width = 120, height = 65, units = "mm")  

##################################################################
# unique otus profile among three layers
# Load the library
library(limma)

# Generate example data
sur_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(select(., grep('SUR', metadata$sample_id, value = T)))) %>%
  filter(rowsum > 0) %>%
  rownames()
sub_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(select(., grep('SUB', metadata$sample_id, value = T)))) %>%
  filter(rowsum > 0) %>%
  rownames()
pl_venn <- otu %>% data.frame() %>%
  mutate(rowsum = rowSums(select(., grep('PL', metadata$sample_id, value = T)))) %>%
  filter(rowsum > 0) %>%
  rownames()

# What are the possible letters in the universe?
universe <- sort(unique(c(sur_venn, sub_venn, pl_venn)))

# Generate a matrix, with the sets in columns and possible letters on rows
Counts <- matrix(0, nrow=length(universe), ncol=3)
# Populate the said matrix
for (i in 1:length(universe)) {
  Counts[i,1] <- universe[i] %in% sur_venn
  Counts[i,2] <- universe[i] %in% sub_venn
  Counts[i,3] <- universe[i] %in% pl_venn
}

# Specify the colors for the sets
# Prepare a palette of 3 colors:
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
n = 3
cols = gg_color_hue(n)

# Name the columns with the sample names
colnames(Counts) <- c("SUR","SUB","PL")
vennDiagram(vennCounts(Counts), circle.col = cols)

# library
library(VennDiagram)

#Make the plot
venn_plot <- venn.diagram(
  x = list(
    otu %>% data.frame() %>%
      mutate(rowsum = rowSums(select(., grep('SUR', metadata$sample_id, value = T)))) %>%
      filter(rowsum > 0) %>%
      rownames(), 
    otu %>% data.frame() %>%
      mutate(rowsum = rowSums(select(., grep('SUB', metadata$sample_id, value = T)))) %>%
      filter(rowsum > 0) %>%
      rownames(), 
    otu %>% data.frame() %>%
      mutate(rowsum = rowSums(select(., grep('PL', metadata$sample_id, value = T)))) %>%
      filter(rowsum > 0) %>%
      rownames()
  ),
  category.names = c("SUR" , "SUB" , "PL"),
  filename = NULL,
  # output = F,
  # imagetype="png" ,
  # height = 480 ,
  # width = 480 ,
  # resolution = 300,
  # compression = "lzw",
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = cols,
  # Numbers
  cex = 1,
  fontfamily = "sans",
  # Set names
  cat.cex = 1,
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  cat.col = 'black',
  rotation = 1
)
grid::grid.draw(venn_plot)
# grid.newpage()
# pdf(file= file.path(save.dir, "OTU_Group_venn.pdf"), height = 3, width = 3)
# grid.draw(venn_plot)
# dev.off()

plt <- venn.diagram(
  filename = NULL,
  cex = 1,
  cat.cex = 1,
  lwd = 2,
)

####################################################################
# null model
phylo_sur <- subset_samples(phylo, layer == 'SUR')
phylo_sur <- prune_taxa(taxa_sums(phylo_sur) >= 1, phylo_sur)
phylo_sub <- subset_samples(phylo, layer == 'SUB')
phylo_sub <- prune_taxa(taxa_sums(phylo_sub) >= 1, phylo_sub)
phylo_pl <- subset_samples(phylo, layer == 'PL')
phylo_pl <- prune_taxa(taxa_sums(phylo_pl) >= 1, phylo_pl)

## extract the community matrix (row: sample; column: species) and phylogenetic tree
comm_sur <- t(otu_table(phylo_sur))
tree_sur <- phy_tree(phylo_sur)

comm_sub <- t(otu_table(phylo_sub))
tree_sub <- phy_tree(phylo_sub)

comm_pl <- t(otu_table(phylo_pl))
tree_pl <- phy_tree(phylo_pl)

#基于 Ning et al (2020) 的方法的群落构建分析
#ses.cut = 1.96，以 βNRI=1.96 作为同质和异质选择的阈值；rc.cut=0.95，以 RC=0.95 作为扩散和漂变的阈值
#bin.size.limit 用来指定 bin 中的分类单元数,这个数据集太小因此使用使用 5，实际情况中可根据系统发育信号测试或尝试一些设置选择合理的 bin，作者的经验是 12、24、48 等
#rand 指定随机化次数构建零分布，nworker 用来多线程运行以提升计算速率
#更多详情 ?icamp.big
set.seed(123)
icamp.sur.out <- icamp.big(comm = comm_sur, tree = tree_sur, pd.wd = paste0(save.dir,"/tables/null_model/sur"), 
                           ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24, 
                           rand = 1000, nworker = 8)
## Heterogeneous.Selection, Homogeneous.Selection, Dispersal.Limitation, Homogenizing.Dispersal, Drift.and.Others.
head(icamp.sur.out$CbMPDiCBraya)
write.csv(icamp.sur.out$CbMPDiCBraya, 
          file.path(save.dir, './tables/null_model/sur/iCAMP.process.CbMPDiCBraya.csv'))

set.seed(123)
icamp.sub.out <- icamp.big(comm = comm_sub, tree = tree_sub, pd.wd = paste0(save.dir,"/tables/null_model/sub"), 
                           ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24, 
                           rand = 1000, nworker = 8)
head(icamp.sub.out$CbMPDiCBraya)
write.csv(icamp.sub.out$CbMPDiCBraya, 
          file.path(save.dir, './tables/null_model/sub/iCAMP.process.CbMPDiCBraya.csv'))


set.seed(123)
icamp.pl.out <- icamp.big(comm = comm_pl, tree = tree_pl, 
                          pd.wd = paste0(save.dir,"/tables/null_model/pl"),
                          ses.cut = 1.96, rc.cut = 0.95, bin.size.limit = 24,
                          rand = 1000, nworker = 8)
head(icamp.pl.out$CbMPDiCBraya)
write.csv(icamp.pl.out$CbMPDiCBraya, 
          file.path(save.dir, './tables/null_model/pl/iCAMP.process.CbMPDiCBraya.csv'))

## arrange the data
null_sur <- read.csv(file.path(save.dir, './tables/null_model/sur/iCAMP.process.CbMPDiCBraya.csv'),
                     header = T, row.names = 1, stringsAsFactors = F)
null_sub <- read.csv(file.path(save.dir, './tables/null_model/sub/iCAMP.process.CbMPDiCBraya.csv'),
                     header = T, row.names = 1, stringsAsFactors = F)
null_pl <- read.csv(file.path(save.dir, './tables/null_model/pl/iCAMP.process.CbMPDiCBraya.csv'),
                    header = T, row.names = 1, stringsAsFactors = F)

null_df <- rbind(cbind(layer = rep('SUR', nrow(null_sur)), null_sur[, 3:7]),
                 cbind(layer = rep('SUB', nrow(null_sub)), null_sub[, 3:7]),
                 cbind(layer = rep('PL', nrow(null_pl)), null_pl[, 3:7]))

## alluvial diagram
library(ggalluvial)
null_plot <- null_df %>%
  group_by(layer) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))) %>%
  ggplot(aes(y = value, x = layer)) + 
  geom_flow(aes(alluvium = process), alpha = 0.9, lty = 2, fill = "white",
            color = "black", curve_type = "linear", width = 0.5) + 
  geom_col(aes(fill = process), width = 0.5, color = "black") + 
  labs(x = 'Layers', y = 'Relative importance', fill = 'Processes') +
  # scale_fill_manual(values = c("#000000", "#294e63", "#496a80", "#7c98ac", "#b3c4d2")) + 
  scale_y_continuous(expand = c(0, 0)) + theme_classic() +
  theme(axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"))
ggsave(file.path(save.dir, './figs/null_model/null_stacked_plot.pdf'), null_plot, width = 6, height = 4)

## boxplot
library(ggpubr)
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
p <- null_df %>%
  select(c('layer', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))) %>%
  ggplot(aes(layer, value))+
  geom_boxplot(width = 0.5, aes(fill = layer))+
  facet_grid(. ~ process, scales = 'free_x', space = 'free_x') +
  stat_compare_means(comparisons = my_comparisons, paired = TRUE, 
                     p.adjust.method = "BH", label = "p.signif") +
  #scale_fill_manual(values= cols) +
  labs(x = 'Layers', y = 'Relative importance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(file.path(save.dir, './figs/null_model/null_box_plot.pdf'), p, width = 8, height = 4)

############################################################
# network analysis
library(phyloseq)
## pruned the OTUs with the mean relative abundance less than 0.1%
phylo_rel <- transform_sample_counts(phylo, function(x) x / sum(x) )
phylo_rel_core <- phyloseq::filter_taxa(phylo_rel, function(x) mean(x) > .0001, TRUE)
ntaxa(phylo_rel_core)
phylo_rel <- microbiome::transform(phylo, "compositional")
phylo_rel_core <- microbiome::core(phylo_rel, detection = 0.1/100, prevalence = 1/66)
ntaxa(phylo_rel_core)
physeqaF <- phylo_rel_core
physeqaF

## divided the filtered OTU table into three groups (i.e., SUR, SUB, and PL)
phylo_sur_net <- subset_samples(physeqaF, Layer == "SUR")
phylo_sur_net <- prune_taxa(taxa_sums(phylo_sur_net) > 0, phylo_sur_net) 

phylo_sub_net <- subset_samples(physeqaF, Layer == "SUB")
phylo_sub_net <- prune_taxa(taxa_sums(phylo_sub_net) > 0, phylo_sub_net) 

phylo_pl_net <- subset_samples(physeqaF, Layer == "PL")
phylo_pl_net <- prune_taxa(taxa_sums(phylo_pl_net) > 0, phylo_pl_net) 

## extract the otu talbe
comm.sur.net <- otu_table(phylo_sur_net)
comm.sub.net <- otu_table(phylo_sub_net)
comm.pl.net <- otu_table(phylo_pl_net)

# #prepare the format of input data for Molecular Ecological Network Analyses (MENA) pipeline (http://ieg4.rccc.ou.edu/mena)
# comm.sur.net[comm.sur.net == 0] <- ""
# comm.sub.net[comm.sub.net == 0] <- ""
# comm.pl.net[comm.pl.net == 0] <- ""
# 
# write.table(comm.sur.net, file.path(save.dir, 'tables/network/16S/comm_sur_net.txt'))
# write.table(comm.sub.net, file.path(save.dir, 'tables/network/16S/comm_sub_net.txt'))
# write.table(comm.pl.net, file.path(save.dir, 'tables/network/16S/comm_pl_net.txt'))
# ## then import the otu table into MENA pipeline

##To compare the difference in empirical network 
##indices between high-altitude and high-latitude 
##thermokarst lakes, the student t-test was employed 
##using the standard deviations derived from corresponding random networks
library(BSDA)
#Average clustering coefficient
tsum.test(mean.x = 0.643, s.x = 0.006, n.x = 118,
          mean.y = 0.653, s.y = 0.003, n.y = 118)

#Average path distance
tsum.test(mean.x = 2.692, s.x = 0.018, n.x = 118,
          mean.y = 5.080, s.y = 0.015, n.y = 118)

tsum.test(mean.x = 0.544, s.x = 0.003, n.x = 118,
          mean.y = 0.868, s.y = 0.004, n.y = 118)

net.cal<-function(x){  # x: community matrix with species as rows, sample as cols
  library(WGCNA)
  library(impute)
  library(preprocessCore)
  library(igraph)
  x <- t(x)
  cp <- corAndPvalue(x,
                   use = "pairwise.complete.obs",
                   alternative = c("two.sided")
  )
  # To get a matrix of the corresponding FDR, use:
  cp.r <- cp$cor  # extract the r matrix
  cp.p <- apply(cp$p, 2, p.adjust, method = "fdr")  # extract the p value matrix
  
  # define the threshold of correlation strength
  cp.r[cp.p > 0.01 | abs(cp.r) < 0.8] = 0
  igraph <- graph_from_adjacency_matrix(cp.r, mode = "undirected", weighted = TRUE, diag = FALSE)
  igraph <- simplify(igraph)  # remove the self loops
  bad.vs <- V(igraph)[degree(igraph) == 0] 
  igraph <- delete.vertices(igraph, bad.vs) # remove the nodes with 0 degree
  
  #该模式下，边权重代表了相关系数
  #由于权重通常为正值，因此最好取个绝对值，相关系数重新复制一列
  E(igraph)$correlation <- E(igraph)$weight
  E(igraph)$weight <- abs(E(igraph)$weight)
  
  results<-list(cp.r, igraph)
  return(results)
}
## Construct igraph project
net_sur <- net.cal(comm.sur.net)
net_sub <- net.cal(comm.sub.net)
net_pl <- net.cal(comm.pl.net)

## obtain the network
igraph_sur <- net_sur[[2]]
igraph_sur   #10 19

igraph_sub <- net_sub[[2]]
igraph_sub   #204 499

igraph_pl <- net_pl[[2]]
igraph_pl   #204 499

# calculate vulnerability of each node
source('https://raw.githubusercontent.com/Mengting-Maggie-Yuan/warming-network-complexity-stability/master/Fig3_and_S7.stability/info.centrality.R')
sur.node.vul<-info.centrality.vertex(igraph_sur)
sub.node.vul<-info.centrality.vertex(igraph_sub)
pl.node.vul<-info.centrality.vertex(igraph_pl)
max(sur.node.vul)
max(sub.node.vul)
max(pl.node.vul)


# Network attack analysis
# Libraries
library(SpiecEasi)
require(reshape2)
library(phyloseq)
library(igraph)
# Calculating the natural connectivity from adjacency matrix
ncc <- function(ig) {
  evals <- eigen(ig)$value
  nc <- log(mean(exp(evals)))
}

# Calculating the natural connectivity from adjacency matrix of a graph
natcon <- function(ig) {
  adj <- get.adjacency(ig)
  evals <- eigen(adj)$value
  nc <- log(mean(exp(evals)))
}

# Targeted attack ordered by betweenness
nc.attackbetweenness <- function(ig) {
  hubord <- order(rank(betweenness(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}

# Targeted attack ordered by node degree.
nc.attackdegree <- function(ig) {
  hubord <- order(rank(degree(ig)), decreasing=TRUE)
  sapply(1:round(vcount(ig)*.8), function(i) {
    ind <- hubord[1:i]
    tmp <- delete_vertices(ig, V(ig)$name[ind])
    natcon(tmp)
  })
}


# Node removals
attack<-function (adj.mat, node.sup)
{
  n.nodes <- dim(adj.mat)[1]
  adj.mat[node.sup, ] <- rep(0, n.nodes)
  adj.mat[, node.sup] <- rep(0, n.nodes)
  nc<-ncc(adj.mat)
  list(new.mat = adj.mat, nc=nc)
}

nc.deg.rmt_sur <- nc.attackdegree(igraph_sur)
nc.bet.rmt_sur <- nc.attackbetweenness(igraph_sur)

nc.deg.rmt_sub <- nc.attackdegree(igraph_sub)
nc.bet.rmt_sub <- nc.attackbetweenness(igraph_sub)

nc.deg.rmt_pl <- nc.attackdegree(igraph_pl)
nc.bet.rmt_pl <- nc.attackbetweenness(igraph_pl)

nc.deg.rmt_sur <- cbind(rm_pro = c(1:length(nc.deg.rmt_sur))/vcount(igraph_sur), nat.connet = nc.deg.rmt_sur)
nc.bet.rmt_sur <- cbind(rm_pro = c(1:length(nc.bet.rmt_sur))/vcount(igraph_sur), nat.connet = nc.bet.rmt_sur)

nc.deg.rmt_sub <- cbind(rm_pro = c(1:length(nc.deg.rmt_sub))/vcount(igraph_sub), nat.connet = nc.deg.rmt_sub)
nc.bet.rmt_sub <- cbind(rm_pro = c(1:length(nc.bet.rmt_sub))/vcount(igraph_sub), nat.connet = nc.bet.rmt_sub)

nc.deg.rmt_pl <- cbind(rm_pro = c(1:length(nc.deg.rmt_pl))/vcount(igraph_pl), nat.connet = nc.deg.rmt_pl)
nc.bet.rmt_pl <- cbind(rm_pro = c(1:length(nc.bet.rmt_pl))/vcount(igraph_pl), nat.connet = nc.bet.rmt_pl)


setwd(save.wd)
write.table(pa.nc.deg.rmt, file="pa-deg-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(pa.nc.bet.rmt, file="pa-bet-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(tp.nc.deg.rmt, file="tp-deg-att-rmt.txt", 
            sep="\t", quote=F, row.names = T, col.names = T)
write.table(tp.nc.bet.rmt, file="tp-bet-att-rmt.txt",
            sep="\t", quote=F, row.names = T, col.names = T)



robut.df.rmt <- data.frame(rbind(cbind(Layer = rep('SUR', sum(nrow(nc.deg.rmt_sur), nrow(nc.bet.rmt_sur))),
                                       type = c(rep('Proportion of removed nodes', nrow(nc.deg.rmt_sur)), 
                                                rep('Proportion of removed betweenness', nrow(nc.bet.rmt_sur))),
                                       nat.connet = rbind(nc.deg.rmt_sur, nc.bet.rmt_sur)),
                                 cbind(Layer = rep('SUB', sum(nrow(nc.deg.rmt_sub), nrow(nc.bet.rmt_sub))),
                                       type = c(rep('Proportion of removed nodes', nrow(nc.deg.rmt_sub)), 
                                                rep('Proportion of removed betweenness', nrow(nc.bet.rmt_sub))),
                                       nat.connet = rbind(nc.deg.rmt_sub, nc.bet.rmt_sub)),
                                 cbind(Layer = rep('PL', sum(nrow(nc.deg.rmt_pl), nrow(nc.bet.rmt_pl))),
                                       type = c(rep('Proportion of removed nodes', nrow(nc.deg.rmt_pl)), 
                                                rep('Proportion of removed betweenness', nrow(nc.bet.rmt_pl))),
                                       nat.connet = rbind(nc.deg.rmt_pl, nc.bet.rmt_pl))))
robut.df.rmt$Layer <- factor(robut.df.rmt$Layer, ordered = T, 
                              levels = c('SUR', 'SUB', 'PL'))
robut.df.rmt$type <- factor(robut.df.rmt$type, ordered = T, 
                            levels = c('Proportion of removed nodes', 'Proportion of removed betweenness'))

robut.df.rmt$rm_pro <- as.numeric(robut.df.rmt$rm_pro)
robut.df.rmt$nat.connet <- as.numeric(robut.df.rmt$nat.connet)

robutness_network_plot <- ggplot(robut.df.rmt, aes(x = rm_pro, y = nat.connet, colour = Layer))+
  geom_point() +
  scale_color_manual(values= c('#0073c2', '#efc000', '#868686')) +
  labs(x = NULL, y = 'Natural connectivity') +
  facet_wrap(~type, scales = 'free') +
  theme_bw() +
  theme(strip.text = element_text(size = 14),
        axis.line = element_line(colour = "black"),
        axis.title = element_text(size = 14, colour = 'black'),
        axis.text = element_text(size = 12, colour = 'black'),
        axis.text.x = element_text(colour='black', size = 12, 
                                   angle = 45, hjust = 1),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        legend.position = c(0.9, 0.8),
        panel.grid = element_blank())
robutness_network_plot

#community.plot
cairo_pdf(file.path(save.dir, "/figs/network/robutness_network.pdf"), width = 8, height = 4, onefile = TRUE)
robutness_network_plot
dev.off()

## modularity (cluster_fast_greedy)
fc_sur <- cluster_fast_greedy(igraph_sur)
modul_sur <- modularity(igraph_sur, membership(fc_sur))
modul_sur  # 0.7165607

fc_sub <- cluster_fast_greedy(igraph_sub)
modul_sub <- modularity(igraph_sub, membership(fc_sub))
modul_sub  # 0.6715852

fc_pl <- cluster_fast_greedy(igraph_pl)
modul_pl <- modularity(igraph_pl, membership(fc_pl))
modul_pl  # 0.7311247

#igraph提供了可以被 gephi 或 cytoscape 等直接识别的格式
#graphml 格式，可使用 gephi 软件打开并进行可视化编辑
write.graph(igraph_sur, file.path(save.dir, 'tables/network/16S/network_sur.graphml'), format = 'graphml')
write.graph(igraph_sub, file.path(save.dir, 'tables/network/16S/network_sub.graphml'), format = 'graphml')
write.graph(igraph_pl, file.path(save.dir, 'tables/network/16S/network_pl.graphml'), format = 'graphml')

library(RCy3)  # ensure that cytoscape app is running
createNetworkFromIgraph(igraph_sur,"Igraph_sur")
createNetworkFromIgraph(igraph_sub,"Igraph_sub")
createNetworkFromIgraph(igraph_pl,"Igraph_pl")

