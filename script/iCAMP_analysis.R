save.dir = "E:/permafrost/result"
if(!dir.exists(save.dir)){dir.create(save.dir)}
setwd(save.dir)
library(iCAMP)
#extract the top 500 ASVs in each layer
library(tidyverse)
top_n_sel <- function(i, N) {
  as.data.frame(otu[, grep(i, colnames(otu))]) %>% 
    mutate(mean_abundance = rowMeans(.)) %>% 
    arrange(desc(mean_abundance)) %>% top_n(500, mean_abundance) %>% 
    select(grep(i, colnames(otu))) %>% rownames()
}
otu_top500_sur <- top_n_sel('SUR', 500)
otu_top500_sub <- top_n_sel('SUB', 500)
otu_top500_pl <- top_n_sel('PL', 500)
otu_top500 <- full_join(otu_top500_sur, otu_top500_sub, otu_top500_pl)

otu_top500_sur <- as.data.frame(otu[, grep('SUR', colnames(otu))]) %>% 
  mutate(mean_abundance = rowMeans(.)) %>% 
  arrange(desc(mean_abundance)) %>% top_n(500, mean_abundance) %>% 
  select(grep('SUR', colnames(otu))) %>% rownames()

phylo_sur <- subset_samples(phylo, layer == 'SUR')
phylo_sur <- prune_taxa(taxa_sums(phylo_sur)>=1, phylo_sur)
phylo_sub <- subset_samples(phylo, layer == 'SUB')
phylo_sub <- prune_taxa(taxa_sums(phylo_sub)>=1, phylo_sub)
phylo_pl <- subset_samples(phylo, layer == 'PL')
phylo_pl <- prune_taxa(taxa_sums(phylo_pl)>=1, phylo_pl)

#以 iCAMP 包的内置数据集 example.data 为例(行样本, 列 OTU)，包含一个小型微生物群落数据及它们的系统发育树
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
#异质选择（Heterogeneous.Selection）、同质选择（Homogeneous.Selection）、扩散限制（Dispersal.Limitation）、同质扩散（Homogenizing.Dispersal）以及漂变（Drift.and.Others）的样本对得分
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


null_sur <- read.csv(file.path(save.dir, './tables/null_model/sur/iCAMP.process.CbMPDiCBraya.csv'),
                     header = T, row.names = 1, stringsAsFactors = F)
null_sub <- read.csv(file.path(save.dir, './tables/null_model/sub/iCAMP.process.CbMPDiCBraya.csv'),
                     header = T, row.names = 1, stringsAsFactors = F)
null_pl <- read.csv(file.path(save.dir, './tables/null_model/pl/iCAMP.process.CbMPDiCBraya.csv'),
                    header = T, row.names = 1, stringsAsFactors = F)

null_df <- rbind(cbind(layer = rep('SUR', nrow(null_sur)), null_sur[, 3:7]),
                 cbind(layer = rep('SUB', nrow(null_sub)), null_sub[, 3:7]),
                 cbind(layer = rep('PL', nrow(null_pl)), null_pl[, 3:7]))


library(ggalluvial)

#alluvial diagram
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

#boxplot
library(ggpubr)
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
p <- null_df %>%
  select(c('layer', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))) %>%
  ggplot(aes(layer, value))+
  geom_boxplot(width = 0.5, aes(fill = layer))+
  facet_grid(. ~ process, scales = 'free_x', space = 'free_x') +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols) +
  labs(x = 'Layers', y = 'Relative importance', fill='Layers') +
  theme_bw() +
  theme(axis.title = element_text(colour = "black"),
        axis.text = element_text(colour = "black"),
        panel.grid = element_blank(),
        panel.spacing = unit(0, "lines"))
ggsave(file.path(save.dir, './figs/null_model/null_box_plot.pdf'), p, width = 8, height = 4)




# Dissimilarity-Overlap analysis
# Dissimilarity
d <- phyloseq::distance(microbiome::transform(phylo_pl, "compositional"), "jsd", parallel=TRUE)

# Overlap
o <- microbiome::overlap(phylo_pl, detection = 0.2/100)

# Compare
dvec <- d[lower.tri(d)]
ovec <- o[lower.tri(o)]

# Assess rough correlation
cc <- cor(dvec, ovec, method = "spearman", use = "pairwise.complete.obs")

# Scatterplot
plot(dvec, ovec, pch = 2, main = paste("Spearman rho", round(cc, 2)), las = 1, xlab = "Dissimilarity (Jensen-Shannon)", ylab = "Overlap")










#stegen
library(iCAMP)

#以 iCAMP 包的内置数据集 example.data 为例
data(example.data)
comm <- example.data$comm  #包含 20 行样本的 30 列 OTU
tree <- example.data$tree  #30 个 OTU 序列的系统发育树

#根据进化树获取系统发育距离
pd <- cophenetic(tree) 

#基于 Stegen et al (2013 and 2015) 的方法的群落构建分析
#更多详情 ?qpen，本示例默认以 βNTI=2 和 RC=0.95 作为不同过程的阈值，随机化 1000 次构建零分布，8 线程执行
set.seed(123)
qpen.out <- qpen(comm = comm, pd = pd, sig.bNTI = 2, sig.rc = 0.95, rand.time = 1000, nworker = 8)

qpen.out$ratio  #异质选择（Heterogeneous.Selection）、同质选择（Homogeneous.Selection）、扩散限制（Dispersal.Limitation）、同质扩散（Homogenizing.Dispersal）以及漂变（Undominated）的样本对占比
head(qpen.out$result)  #各样本对的 bMNTD、BC、bNTI、RC 值以及最重要的生态过程

#输出
#write.csv(qpen.out, 'qpen.out.csv', row.names = FALSE)