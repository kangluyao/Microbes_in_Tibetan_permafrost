#extract the top 500 ASVs in each layer
top_n_sel <- function(i, N) {
  as.data.frame(otu[, grep(i, colnames(otu))]) %>% 
    mutate(mean_abundance = rowMeans(.)) %>% 
    arrange(desc(mean_abundance)) %>% top_n(500, mean_abundance) %>% 
    select(grep('SUR', colnames(otu))) %>% rownames()
}
otu_top500_sur <- top_n_sel('SUR', 500)
otu_top500_sub <- top_n_sel('SUB', 500)
otu_top500_pl <- top_n_sel('PL', 500)
otu_top500 <- full_join(otu_top500_sur, otu_top500_sub, otu_top500_pl)

otu_top500_sur <- as.data.frame(otu[, grep('SUR', colnames(otu))]) %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% top_n(500, mean_abundance) %>% 
  select(grep('SUR', colnames(otu))) %>% rownames()

phylo_rel <- transform_sample_counts(phylo, function(x) x / sum(x) )
phylo_rel_core <- phyloseq::filter_taxa(phylo_rel, function(x) mean(x) > .0001, TRUE)
ntaxa(phylo_rel_core)


#读取 OTU 丰度表和样本分组
otu_specif <- data.frame(otu_table(phylo_rel_core), Taxonomy = data.frame(tax_table(phylo_rel_core))$Phylum) 
##计算各组样本中，OTU 的特异性(specificity)和占有率(occupancy)
Nindividuals_S <- rep(0, nrow(otu_specif))
for (i in unique(metadata$Layer)) {
  otu_group_i <- otu_specif[ ,rownames(subset(metadata, Layer == i))]
  Nindividuals_S <- Nindividuals_S + rowMeans(otu_group_i)  #计算 Nindividuals S
}

spec_occu <- NULL
for (i in unique(metadata$Layer)) {
  otu_group_i <- otu_specif[ ,rownames(subset(metadata, Layer == i))]
  Nindividuals_SH <- apply(otu_group_i, 1, mean)  #计算 Nindividuals SH
  Specificity <- Nindividuals_SH / Nindividuals_S  #计算 Specificity
  Nsites_H <- ncol(otu_group_i)  #计算 Nsites H
  Nsites_SH <- apply(otu_group_i, 1, function(x) sum(x>0))  #计算 Nsites SH
  Occupancy <- Nsites_SH / Nsites_H  #计算 Occupancy
  spec_occu_group_i <- data.frame(layer = i, OTU = rownames(otu_group_i), 
                                  Specificity = Specificity, Occupancy = Occupancy, 
                                  Abundance_mean = rowMeans(otu_group_i), Taxonomy = otu_specif$Taxonomy)
  spec_occu <- rbind(spec_occu, spec_occu_group_i)  #合并各组统计
}
spec_occu <- spec_occu %>% mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL')))
head(spec_occu)  #该数据框包含各组中各 OTU 的名称、特异性、占有率、平均丰度、类群信息等

#输出统计表格
#write.table(spec_occu, 'spec_occu.txt', row.names = FALSE, sep = '\t', quote = FALSE)

##识别各组样本中的特化种（specialist species）
#在上述统计表格中直接根据特异性和占有率 ≥0.7 做筛选，保留下的 OTU 即为特化种
spec_occu_specialist <- subset(spec_occu, Specificity >= 0.7 & Occupancy >= 0.7) %>% 
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL')))
head(spec_occu_specialist)

#输出统计表格
#write.table(spec_occu_specialist, 'spec_occu_specialist.txt', row.names = FALSE, sep = '\t', quote = FALSE)

#绘制 SPEC-OCCU 图, 添加阈值线
library(ggplot2)

p <- spec_occu %>% mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(Occupancy, Specificity)) +
  # geom_point(aes(size = log10(Abundance_mean), color = 'Taxonomy'), alpha = 0.5) +
  geom_jitter(aes(size = log10(Abundance_mean), colour = Specificity >= 0.7 & Occupancy >= 0.7), alpha = 0.8) +  #如果觉得普通散点图中的点重叠严重不好看，可以仿照 Gweon et al (2021) 使用抖动点图来展示
  scale_colour_manual(values = c("grey", "#5f9ea0")) +
  scale_size(breaks = c(-1, -2, -3, -4), labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4})), range = c(0, 1)) +
  #scale_color_manual(values = c('#E7272E', '#F59D1F', '#768CC5', '#9BC648', '#794779', '#A19E9D'), limits = c('Proteobacteria', 'Actinobacteria', 'Acidobacteria', 'Bacteroidetes', 'Firmicutes', 'Others')) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', color = 'gray30'), 
        legend.key = element_blank(), legend.position = "bottom") +
  facet_wrap(~layer, ncol = 3) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  labs(x = 'Occupancy', y = 'Specificity', size = 'Average relative abundance') + # add color option, color = 'Taxonomy'
  coord_cartesian(clip = 'off') +
  geom_segment(aes(x = 0.7, xend = 1, y = 0.7, yend = 0.7), linetype = 2)+
  geom_segment(aes(x = 0.7, xend = 0.7, y = 0.7, yend = 1), linetype = 2)

p <- ggplot(spec_occu_specialist, aes(x = Occupancy, y = Specificity, colour = layer, size = log10(Abundance_mean))) +
  geom_point(data = spec_occu, colour = "#0087A8", alpha = .2) +
  geom_point() + 
  scale_size(breaks = c(-1, -2, -3, -4), labels = c(expression(10^{-1}), expression(10^{-2}), expression(10^{-3}), expression(10^{-4})), range = c(0, 1)) +
  theme(axis.title = element_text(colour = 'black'),
        axis.text = element_text(colour = 'black'),
        panel.grid = element_blank(), 
        panel.background = element_rect(fill = 'white', color = 'gray30'), 
        legend.key = element_blank(), legend.position = "bottom") +
  facet_wrap(~layer, ncol = 3) +
  scale_x_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  scale_y_continuous(breaks = c(0, 0.5, 1), expand = c(0, 0), limit = c(0, 1)) +
  labs(x = 'Occupancy', y = 'Specificity', size = 'Average relative abundance') + # add color option, color = 'Taxonomy'
  coord_cartesian(clip = 'off') +
  geom_segment(aes(x = 0.7, xend = 1, y = 0.7, yend = 0.7), color = "black", linetype = 2)+
  geom_segment(aes(x = 0.7, xend = 0.7, y = 0.7, yend = 1), color = "black", linetype = 2)

cairo_pdf(filename = file.path(save.dir, "./figs/diff_abun/spec_occu_plot.pdf"), 
          width = 8.5, height = 4.1, onefile = TRUE)
p
dev.off()

ggsave(file.path(save.dir, "./figs/diff_abun/spec_occu_plot.pdf"), p, width = 180, height = 90, units = "mm")

# tree图显示特化种的具体类群
# devtools::install_github("ropensci/taxa", force = TRUE)
library(metacoder)
metacode_tax <- read.table(file.path(wd_16s, 'metacode_tax.txt'), header=T, sep="\t", comment.char="")
obj <- parse_tax_data(metacode_tax,
                      class_cols = "lineage", # the column that contains taxonomic information
                      class_sep = ";", # The character used to separate taxa in the classification
                      class_regex = "^(.+)__(.+)$", # Regex identifying where the data for each taxon is
                      class_key = c(tax_rank = "info", # A key describing each regex capture group
                                    tax_name = "taxon_name"))
#comparison between PL and SUR
#transform to proportion data
obj$data$tax_data <- calc_obs_props(obj, "tax_data")
spec_sur <- rownames(otu) %in% subset(spec_occu_specialist, layer == 'SUR')$OTU
spec_sub <- rownames(otu) %in% subset(spec_occu_specialist, layer == 'SUB')$OTU
spec_pl <- rownames(otu) %in% subset(spec_occu_specialist, layer == 'PL')$OTU
obj_sur <- filter_obs(obj, "tax_data", spec_sur, drop_taxa = T)
print(obj_sur)
obj_sub <- filter_obs(obj, "tax_data", spec_sub, drop_taxa = T)
print(obj_sub)
obj_pl <- filter_obs(obj, "tax_data", spec_pl, drop_taxa = T)
print(obj_pl)
#Getting per-taxon information
obj_sur$data$tax_abund <- calc_taxon_abund(obj_sur, "tax_data", cols = metadata$sample_id)
obj_sub$data$tax_abund <- calc_taxon_abund(obj_sub, "tax_data", cols = metadata$sample_id)
obj_pl$data$tax_abund <- calc_taxon_abund(obj_pl, "tax_data", cols = metadata$sample_id)
#calculate the number of samples that have reads for each taxon
obj_sur$data$tax_occ <- calc_n_samples(obj_sur, "tax_abund", groups = metadata$layer, cols = metadata$sample_id)
print(obj_sur)
obj_sub$data$tax_occ <- calc_n_samples(obj_sub, "tax_abund", groups = metadata$layer, cols = metadata$sample_id)
print(obj_sub)
obj_pl$data$tax_occ <- calc_n_samples(obj_pl, "tax_abund", groups = metadata$layer, cols = metadata$sample_id)
print(obj_pl)

#Plotting taxonomic data
set.seed(1) # This makes the plot appear the same each time it is run 
p_sur <- obj_sur %>% 
  metacoder::filter_taxa(obj_sur$data$class_data$taxon_id[obj_sur$data$class_data$tax_rank == 'o'], 
                         supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs) # node_color_range = c("grey", "#F8766D")
ggsave(file.path(save.dir, "./figs/diff_abun/spec_sur_tree.pdf"), p_sur, width = 89, height = 59, units = "mm")
set.seed(1) 
p_sub = obj_sub %>% 
  metacoder::filter_taxa(obj_sub$data$class_data$taxon_id[obj_sub$data$class_data$tax_rank == 'o'], 
                         supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs)
ggsave(file.path(save.dir, "./figs/diff_abun/spec_sub_tree.pdf"), p_sub, width = 89, height = 59, units = "mm")
set.seed(1) 
p_pl = obj_pl %>% 
  metacoder::filter_taxa(obj_pl$data$class_data$taxon_id[obj_pl$data$class_data$tax_rank == 'o'], 
                         supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_size = n_obs,
            node_color = n_obs)
ggsave(file.path(save.dir, "./figs/diff_abun/spec_pl_tree.pdf"), p_pl, width = 89, height = 59, units = "mm")

spec_occu_tree <- cowplot::plot_grid(p_sur, p_sub, p_pl, ncol = 3)
ggsave(file.path(save.dir, "./figs/diff_abun/spec_occu_tree.pdf"), spec_occu_tree, width = 210, height = 120, units = "mm")










# the relation ship between the envs and specificity OTUs
speci_otus <- c(subset(spec_occu_specialist, layer == 'SUR')$OTU, 
                subset(spec_occu_specialist, layer == 'SUB')$OTU,
                subset(spec_occu_specialist, layer == 'PL')$OTU)
phylo_speci <- prune_taxa(speci_otus, phylo_rel)
ps2 <- tax_glom(phylo_speci, "Order", NArm = TRUE)


#Use code snippet straight provided by authors of PhyloSeq to export phyloseq object to an EdgeR object
phyloseq_to_edgeR = function(physeq, method="RLE", ...){
  require("edgeR")
  require("phyloseq")
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
  } 
  # Now turn into a DGEList
  y = DGEList(counts=x, genes=taxonomy, remove.zeros = TRUE, ...)
  # Calculate the normalization factors
  z = calcNormFactors(y, method=method)
  # Check for division by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    stop("Something wrong with edgeR::calcNormFactors on this data,
         non-finite $norm.factors, consider changing `method` argument")
  }
  # Estimate dispersions
  return(estimateTagwiseDisp(estimateCommonDisp(z)))
}

# Make normalized phyloseq object (ps2) into an edgeR object. It needs a grouping factor. We use location.
dge <- phyloseq_to_edgeR(ps2)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a <- sample_data(ps2)

#So many things don't understand what kind of variables they need to be. Make sure they understand
layer <- as.numeric(a$layer)
moisture <- as.numeric(a$moisture)
pH <- as.numeric(a$pH)
Clay <- as.numeric(a$Clay)
Silt <- as.numeric(a$Silt) 
Sand <- as.numeric(a$Sand)
TN <- as.numeric(a$TN)
SOC <- as.numeric(a$SOC)
NH4_N <- as.numeric(a$NH4) 
NO3_N <- as.numeric(a$NO3)
DIN <- as.numeric(a$DIN)
DON <- as.numeric(a$DON)

# Design for my linear model, must avoid the linear dependencies between your variables
design <-model.matrix(~ moisture + pH + Clay + TN + SOC + NH4_N + NO3_N + DIN + DON) 

# EdgeR needs to calculate dispersion again after you've fed it the design. Why? I don't know.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x <- calcNormFactors(dge, method="RLE")
x <- estimateGLMCommonDisp(dge, design)
x <- estimateGLMTrendedDisp(dge, design)
x <- estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef = 2:10)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-table[,1:9]
table<-as.data.frame(table)

q<-lrt$genes
table<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
             OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2<-table[apply(table[1:9], 1, function(x) any(abs(x)>1.96)), ]

write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted<-reshape::melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter<-unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(phylo_speci)) %in% filter), phylo_speci)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=tax$Genus, label3=tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette

library(RColorBrewer)

colourCount = length(unique(tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)

p<-ggtree(f) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))
p + xlim (0, 3) 

test<-as.matrix(table2[1:15])
tested<-apply(test, 2, function(x){cut(x, br=c(-8, -6, -4, -1.96, 1.96, 4, 6, 8))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(7, "Spectral"))(7)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))

# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#F2F2F2" # -1.96 to 1.96
heatmapcols[2] = "#AABBDD" #"#D7E3F4" # -1.96 to -4
heatmapcols[3] = "#112288" # -4 to -6 - dk blue
heatmapcols[4] = "#FFCCCF"
heatmapcols[5] = "#D9444D"
heatmapcols[6] = "#881111" 
heatmapcols[7] = 
  
  
  p1<-gheatmap(p, tested, offset = 1.5, width = 3.5, font.size=2.5, 
               colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1

