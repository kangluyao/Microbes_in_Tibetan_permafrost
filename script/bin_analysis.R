# set work directory
setwd('e:/permafrost')
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")
# load data from TDbook, including tree_hmptree, 
# df_tippoint (the abundance and types of microbes),
# df_ring_heatmap (the abundance of microbes at different body sites),
# and df_barplot_attr (the abundance of microbes of greatest prevalence)

#loading packages
library(phyloseq)
library(picante)
library(microbiome)
library(tidyverse)
library(readxl)
library(ggpubr)
library(ggplot2)
wd_fun <- file.path(getwd(),"data/metagenome/")

tree.file <- file.path(wd_fun, 'binning_70/result/gtdb_tree/bacteria/tax.unrooted.tree')
abundance_tab.file <- file.path(wd_fun, "binning_70/result/bin_abundance_table.tab")
tax.file <- file.path(wd_fun, "binning_70/result/tax.txt")
bin_fun_file <- file.path(wd_fun, "metabolic/METABOLIC_70/metabolic_heatmap.csv")
metabolic_output_file <- file.path(wd_fun, "metabolic/METABOLIC_70/METABOLIC_result.xlsx")

tree <- read.tree(tree.file)
abundance_tab <- read_delim(abundance_tab.file, col_names = T)
colnames(abundance_tab)[1] <- 'ID'
tax_tab <- read_delim(tax.file, col_names = T)
bin_fun_tab <- read_delim(bin_fun_file, col_names = T)
metabolic_tab <- read_excel(metabolic_output_file)


dat1 <- tax_tab
# count the number of MAGs of each Class
dat1 %>% select(c("Class")) %>%
  mutate(number = c(rep(1, nrow(dat1)))) %>%
  group_by(Class) %>%
  summarise(across(everything(), sum)) %>%
  arrange(desc(number)) %>%
  mutate(Class = gsub("c__", "", Class)) %>%
  mutate(Class = factor(Class, levels = Class)) %>%
  ggplot(aes(x = Class, y = number)) +
  geom_bar(stat = "identity") +
  scale_y_continuous(limits = c(0, 38), expand = c(0, 0)) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"))
  
  
  coord_polar()

dat1 %>% select(c("Class")) %>%
  mutate(number = c(rep(1, nrow(dat1)))) %>%
  group_by(Class) %>%
  summarise(across(everything(), sum)) %>%
  arrange(desc(number)) %>%
  mutate(Class = gsub("c__", "", Class)) %>%
  mutate(Class = factor(Class, levels = Class)) %>%
  ggplot(aes(x = Class, y = number)) +
  geom_bar(stat = "identity", fill = "#0073c2", width = 0.65) +
  scale_y_continuous(limits = c(0, 38), expand = c(0, 0)) +
  labs(x = "Class", y = "Number of MAGs") +
  theme_bw() +
  theme(axis.title = element_text(color = "black", size = 12),
        axis.text.x = element_text(colour = "black", angle = 90, vjust = 0.5, hjust=1),
        axis.text.y = element_text(colour = "black"))



  geom_bar(stat = 'identity', width = 1) +
  coord_polar(theta = 'x') +
  #scale_fill_manual(values = c('#8DD3C7', '#FFFFB3', '#BEBADA', '#FB8072', '#80B1D3', '#FDB462', '#B3DE69', '#FCCDE5', '#BC80BD', '#CCEBC5', 'gray')) +
  theme(panel.grid = element_blank(), 
        panel.background = element_blank(), 
        axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        plot.title = element_text(hjust = 0.5))

  

dat2 <- abundance_tab %>% 
  mutate(SUR = rowMeans(select(., grep('SUR', colnames(abundance_tab), value = T)))) %>%
  mutate(SUB = rowMeans(select(., grep('SUB', colnames(abundance_tab), value = T)))) %>%
  mutate(PL = rowMeans(select(., grep('PL', colnames(abundance_tab), value = T)))) %>%
  select(c(ID, SUR, SUB, PL)) %>%
  pivot_longer(cols = -c(ID), names_to = "Layers", values_to = 'rel_abun') %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
# write the table for itol annotation
write.csv(dat2, "E:/permafrost/data/metagenome/binning_70/result/itol_70/DIY_manual/abundance_annotation.csv")
##metabolic fun heatmap
dat3 <- bin_fun_tab %>% select(-c(1, 3, 4)) %>%
  group_by(Function_final) %>% 
  summarise(across(everything(), sum)) %>%
  filter(Function_final %in% c("Aromatics degradation", "Cellulose degrading", "Hemicullulose debranching",
                               "Endohemicellulases", "Other oligosaccharide degrading", "Amylolytic enzymes",
                               "Chitin degrading", "Methane oxidation - Partculate methane monooxygenase",
                               "Methane oxidation - Soluble methane monoxygenase", "Methane production")) %>%
  pivot_longer(cols = -c(Function_final), names_to = "ID", values_to = 'presence_or_absent') %>%
  mutate(presence_or_absent = ifelse(presence_or_absent >= 1, 1, 0))

dat4 <- metabolic_tab %>% 
  select(c('Category', 'Function', grep('Hmm.presence', colnames(metabolic_tab)))) %>%
  filter(Category %in% c('Ethanol fermentation', 'Fatty acid degradation', 'Aromatics degradation',
                         'Complex carbon degradation', 'Fermentation', 'Methane metabolism', 
                         'Nitrogen cycling', 'Sulfur cycling', 'Hydrogenases', 'As cycling',
                         'Selenate reduction', 'Iron cycling', 'Manganese cycling')) %>% 
  mutate(across(starts_with("s_"), ~ifelse(. == "Present", 1, 0))) %>%
  group_by(Category, Function) %>% 
  summarise(across(everything(), sum)) %>%
  filter(rowSums(across(where(is.numeric)))!=0) %>% 
  pivot_longer(cols = -c(Category, Function), names_to = "ID", values_to = 'presence_or_absent') %>%
  mutate(ID = gsub(".Hmm.presence", "", ID)) %>%
  mutate(presence_or_absent = ifelse(presence_or_absent >= 1, 1, 0))

# write the table for itol annotation
write.csv(dat4, "E:/permafrost/data/metagenome/binning_70/result/itol_70/DIY_manual/annotation.csv")

# count the number of MAGs of each metabolic parhway (Carbon & Nitrigen)
C_N_metapathway <- c("Amylolytic enzymes", "Cellulose degrading", "Endohemicellulases",  
                     "Chitin degrading", "Phenol => Benzoyl-CoA", "Fatty acid degradation",
                     "Methane oxidation - Partculate methane monooxygenase",
                     "Methane oxidation - Soluble methane monoxygenase", 
                     "Ammonia oxidation", "Nitrate reduction", "Nitric oxide reduction",
                     "Nitrite oxidation", "Nitrite reduction to ammonia", "Nitrite reduction",
                     "Nitrous oxide reduction", "N2 fixation")

other_metapathway <- c("Acetogenesis", "Lactate utilization", "Pyruvate oxidation",  
                       "Pyruvate <=> acetyl-CoA + formate", "Acetate to acetyl-CoA", 
                       "Sulfate reduction", "Sulfide oxidation", "Sulfite reduction",
                       "Sulfur oxidation", "Thiosulfate disproportionation", 
                       "Thiosulfate oxidation", "Iron oxidation", "Iron reduction",
                       "Metal (Iron/Manganese) reduction", "Arsenate reduction for detoxification",
                       "Arsenite methylation", "Arsenite oxidation for detoxification and energy metabolism",
                       "Dissimilatory arsenate reduction", "FeFe hydrogenase", "Ni-Fe Hydrogenase",
                       "Manganese oxidation", "Selenate reduction")

dat4 %>% select(c(1, 2, 4)) %>%
  group_by(Category, Function) %>% 
  summarise(across(everything(), sum)) %>%
  filter(Function %in% C_N_metapathway) %>% 
  mutate(Function = factor(Function, levels = C_N_metapathway)) %>%
  ggplot(aes(x = Function, y = presence_or_absent, fill = Function)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FFFFB3", "#FB8072", "#FDB462", "#80B1D3", "#BEBADA", "#B3DE69", "#FCCDE5", "#A6CEE3",
                               "#8DD3C7", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C", "#FDBF6F", "#FF7F00", "#1F78B4")) +
  scale_y_continuous(limits = c(0, 200), expand = c(0, 0)) +
  theme_classic()

# count the number of MAGs of each metabolic parhway (others)
dat4 %>% select(c(1, 2, 4)) %>%
  group_by(Category, Function) %>% 
  summarise(across(everything(), sum)) %>%
  filter(Function %in% other_metapathway) %>% 
  mutate(Function = factor(Function, levels = other_metapathway)) %>%
  ggplot(aes(x = Function, y = presence_or_absent, fill = Function)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#6A3D9A", "#FCCDE5", "#1F78B4", "#A6CEE3", "#B3DE69", "#E69F00", "#56B4E9", "#009E73",
                               "#F0E442", "#0072B2", "#D55E00", "#FB9A99", "#E31A1C", "#FDBF6F", "#8DD3C7", "#FFFFB3", 
                               "#BEBADA", "#FB8072", "#B2DF8A", "#33A02C", "#FF7F00", "#000000")) +
  scale_y_continuous(limits = c(0, 250), expand = c(0, 0)) +
  theme_classic() +
  theme(legend.position = "none")
  
dat1 <- df_tippoint
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
library(ggtree)
library(ggstar)
library(ggnewscale)
library(ggtreeExtra)

nodedf <- data.frame(id = c(269, 261, 26, 357, 356, 354, 120, 119, 380, 414, 415, 196, 430, 199, 242, 235, 254, 251), 
                    Phylum = c("Actinobacteriota", "Firmicutes", "Atribacterota", "Chloroflexota",
                               "Patescibacteria", "Caldisericota", "Deinococcota", "Thermodesulfobiota",
                               "Acidobacteriota", "Desulfobacterota", "Methylomirabilota", "Nitrospirota",
                               "Proteobacteria", "Campylobacterota", "Bacteroidota", "Gemmatimonadota",
                               "Verrucomicrobiota", "Spirochaetota"))


nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.0, 0.997, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
ggtree(tree, branch.length='none', layout='circular', open.angle = 30) +
  geom_text(aes(label=node), size = 2, hjust=-.3) +
  geom_tiplab(size = 2, align = TRUE, linesize = 0.5, offset = 0.5, hjust = 0)

p <- ggtree(tree, branch.length='none', layout='circular', open.angle = 20, size = 0.15) +
  geom_hilight(data = nodedf, mapping=aes(node = id, fill = Phylum),
               alpha = 0.3, color = "grey50", size = 0.05)

p <- ggtree(tree, layout="fan", size = 0.15, open.angle = 5) +
  geom_hilight(data = nodedf, mapping=aes(node = id, fill = Phylum),
               alpha = 0.3, color = "grey50", size = 0.05)
  # geom_cladelab(data=labdf, 
  #               mapping=aes(node=node, 
  #                           label=label,
  #                           offset.text=pos),
  #               hjust=0.5,
  #               angle="auto",
  #               barsize=NA,
  #               horizontal=FALSE, 
  #               fontsize=1.4,
  #               fontface="italic"
  # )

# p <- p %<+% dat1 + geom_star(
#   mapping=aes(fill=Phylum),
#   position="identity",starstroke=0.1) +
#   scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
#                              "#800080", "#9ACD32","#D15FEE","#FFC0CB",
#                              "#EE6A50","#8DEEEE", "#006400","#800000",
#                              "#B0171F","#191970", "red"),
#                     guide=guide_legend(keywidth = 0.5,
#                                        keyheight = 0.5, order=1,
#                                        override.aes=list(starshape=15)),
#                     na.translate=FALSE)+
#   scale_starshape_manual(values=c(15, 1),
#                          guide=guide_legend(keywidth = 0.5,
#                                             keyheight = 0.5, order=2),
#                          na.translate=FALSE)+
#   scale_size_continuous(range = c(1, 2.5),
#                         guide = guide_legend(keywidth = 0.5,
#                                              keyheight = 0.5, order=3,
#                                              override.aes=list(starshape=15)))

p <- p + new_scale_fill() + 
  geom_fruit(data = dat2, geom = geom_tile, 
             mapping = aes(y = ID, x = Layers, alpha = rel_abun, fill = "red"), 
             color = "grey50", offset = 0.01, size = 0.02) + 
  scale_alpha_continuous(range = c(0, 1), 
                         guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 5)) + 
  # scale_fill_manual(values = c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
  #                              "#ceaa3b", "#7a3260", "#432d7c", "#c6d179"), 
  #                   guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)) + 
  geom_treescale(fontsize = 3, linesize = 0.3, x = 0.9, y = 0.1) + 
  theme(legend.position = c(1.2, 0.5), 
        legend.background = element_rect(fill = NA), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        legend.spacing.y = unit(0.04, "cm"))

p <- p + new_scale_fill() + 
  geom_fruit(data = dat4, geom = geom_tile, 
             mapping = aes(y = ID, x = Function, alpha = presence_or_absent, fill = Function), 
             color = "grey50", offset = 0.01, size = 0.02) + 
  scale_alpha_continuous(range = c(0, 1), 
                         guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 5)) + 
  # scale_fill_manual(values = c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
  #                              "#ceaa3b", "#7a3260", "#432d7c", "#c6d179"), 
  #                   guide = guide_legend(keywidth = 0.3, keyheight = 0.3, order = 4)) + 
  geom_treescale(fontsize = 3, linesize = 0.3, x = 0.9, y = 0.1) + 
  theme(legend.position = c(1.2, 0.5), 
        legend.background = element_rect(fill = NA), 
        legend.title = element_text(size = 8), 
        legend.text = element_text(size = 6), 
        legend.spacing.y = unit(0.04, "cm"))
p

p + layout_rectangular() + 
  theme(legend.position=c(.05, .7))


plot.name <- paste0(save.dir, "/figs/bin/bin_70_fun_1.pdf")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8, height = 5.5, onefile = TRUE)
p
dev.off()

# relationship of bin with carbon flux
bin_flux_file = file.path(wd_fun, '/bin_flux.csv')
bin_flux <- read_delim(bin_flux_file, col_names = T)
p <- ggplot(bin_flux, aes(x = s_SUB_15_15, y = Cumulative_5)) +
  geom_point(size=3.5, alpha=0.8, colour= '#0000ff') +
  geom_smooth(method="lm", size=1, se=T, colour='black') +
  xlab("Bin abundance (TPM)") + ylab("CO2") +
  theme_bw() +
  theme(panel.grid=element_blank(), 
        axis.title = element_text(color='black',size=14),
        axis.ticks.length = unit(0.2,"lines"), axis.ticks = element_line(color='black'),
        axis.line = element_blank(), 
        axis.text.y = element_text(colour='black',size=12),
        axis.text.x = element_text(colour='black', size = 12),
        strip.text = element_text(size = 14),
        legend.position='none')

ggsave(file.path(save.dir, "./figs/bin/bin_flux.pdf"),
       p, width = 110, height = 120, units = "mm")


#  Generate a network plot with reactions as node, and edges are MAGs, and are colored by taxonomic group.
sur_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_sur/METABOLIC_Figures_Input/functional_network_input.txt")
sub_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_sub/METABOLIC_Figures_Input/functional_network_input.txt")
pl_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_pl/METABOLIC_Figures_Input/functional_network_input.txt")

network.plots.folder <- file.path(wd_fun, "metabolic/results/network")
table_sur <- read.csv(sur_input_table, header=T, sep="\t")
table_sub <- read.csv(sub_input_table, header=T, sep="\t")
table_pl <- read.csv(pl_input_table, header=T, sep="\t")

#Change the column names

#install.packages("ggraph")
library(ggraph)
library(igraph)
library(tidyverse)
library(tidygraph)

my_graph_sur <- table_sur[, c(2,3,4,5)] %>% 
  graph_from_data_frame()
my_graph_sub <- table_sub[, c(2,3,4,5)] %>% 
  graph_from_data_frame()
my_graph_pl <- table_pl[, c(2,3,4,5)] %>% 
  graph_from_data_frame()

deg_sur <- degree(my_graph_sur, mode="all")
deg_sub <- degree(my_graph_sub, mode="all")
deg_pl <- degree(my_graph_pl, mode="all")

# color
sur_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#432d7c", "#c6d179", "#8f379a", "#c988d1")
sub_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#a152dd", "#5139c2", "#ceaa3b", "#432d7c",
               "#c6d179", "#8f379a", "#70d68c", "#d9432f", "#6ad5be", "#d5416a",
               "#6a75d5", "#c988d1", "#598939")
pl_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
              "#ceaa3b", "#7a3260", "#432d7c", "#c6d179", "#8f379a", "#70d68c",
              "#d9432f", "#6ad5be", "#76c2d7", "#d87a71", "#6a75d5", "#836834",
              "#c988d1", "#598939")
# this generates the whole community plot
community.plot_sur <- table_sur[, c(2,3,4,5)] %>% 
  graph_from_data_frame() %>% 
  ggraph(layout = "linear", circular = TRUE) +
  geom_edge_arc(alpha = .25, 
                aes(width = Coverage.value.average., color=as.factor(Taxonomic.Group))) +
  scale_edge_colour_manual(values = sur_color) +
  geom_node_point(aes(size = 0.02*deg_sur), color = "black", alpha = 0.75) +
  geom_node_text(aes(label = name),  color="black", repel = TRUE)+
  #theme_graph()+
  labs(title = 'Metabolic connections within dataset', 
       subtitle = 'No scaling')
community.plot_sub <- table_sub[, c(2,3,4,5)] %>% 
  graph_from_data_frame() %>% 
  ggraph(layout = "linear", circular = TRUE) +
  geom_edge_arc(alpha = .25, 
                aes(width = Coverage.value.average., color=as.factor(Taxonomic.Group))) +
  scale_edge_colour_manual(values = sub_color) +
  geom_node_point(aes(size = 0.02*deg_sub), color = "black", alpha = 0.75) +
  geom_node_text(aes(label = name),  color="black", repel = TRUE)+
  #theme_graph()+
  labs(title = 'Metabolic connections within dataset', 
       subtitle = 'No scaling')
community.plot_pl <- table_pl[, c(2,3,4,5)] %>% 
  graph_from_data_frame() %>% 
  ggraph(layout = "linear", circular = TRUE) +
  geom_edge_arc(alpha = .25, 
                aes(width = Coverage.value.average., color=as.factor(Taxonomic.Group))) +
  scale_edge_colour_manual(values = pl_color) +
  geom_node_point(aes(size = 0.02*deg_pl), color = "black", alpha = 0.75) +
  geom_node_text(aes(label = name),  color="black", repel = TRUE)+
  #theme_graph()+
  labs(title = 'Metabolic connections within dataset', 
       subtitle = 'No scaling')

#community.plot
plot.name <- paste0(network.plots.folder, "/SurCommunityPlot.PDF")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8, height = 5.5, onefile = TRUE)
community.plot_sur
dev.off()

plot.name <- paste0(network.plots.folder, "/SubCommunityPlot.PDF")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8, height = 5.5, onefile = TRUE)
community.plot_sub
dev.off()

plot.name <- paste0(network.plots.folder, "/PlCommunityPlot.PDF")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8, height = 5.5, onefile = TRUE)
community.plot_pl
dev.off()

#sankey plot





#LEGEND_COLORS	
#(#d64e9e	#6cd54c	#dd49d1	#c8dd41	#a152dd	#4d4040	#5139c2	#ceaa3b	
#432d7c	#c6d179	#8f379a	#70d68c	#d9432f	#6ad5be	#d5416a	
#76c2d7	#d87a71	#6a75d5	#836834	#c988d1	#598939	#7a3260	
#bed3b3	#8f372e	#6082b3	#d47c35	#312749	#d4ac8b	#314825
#LEGEND_LABELS	
# Acidobacteriota(#d64e9e)	 Actinobacteriota(#6cd54c)	
# Alphaproteobacteria(#dd49d1) Atribacterota(#c8dd41)	
# Bacteroidota(#a152dd)	 Caldisericota(#4d4040)	
# Campylobacterota(#5139c2) Chloroflexota(#ceaa3b)
# Deinococcota (#7a3260)
# Desulfobacterota(#432d7c)	 Firmicutes(#c6d179)	
# Gammaproteobacteria(#8f379a)  Gemmatimonadota(#70d68c)	
# Halobacteriota(#d9432f) Methylomirabilota(#6ad5be)
# Nitrospirota(#d5416a) Patescibacteria(#76c2d7)
# Spirochaetota(#d87a71) Thermodesulfobiota(#6a75d5)
# Thermoplasmatota(#836834) Thermoproteota(#c988d1)
# Verrucomicrobiota(#598939)
sur_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#432d7c", "#c6d179", "#8f379a", "#c988d1")
sub_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#a152dd", "#5139c2", "#ceaa3b", "#432d7c",
               "#c6d179", "#8f379a", "#70d68c", "#d9432f", "#6ad5be", "#d5416a",
               "#6a75d5", "#c988d1", "#598939")
pl_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
              "#ceaa3b", "#7a3260", "#432d7c", "#c6d179", "#8f379a", "#70d68c",
              "#d9432f", "#6ad5be", "#76c2d7", "#d87a71", "#6a75d5", "#836834",
              "#c988d1", "#598939")

all_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
               "#5139c2", "#ceaa3b", "#7a3260", "#432d7c", "#c6d179", "#8f379a",
               "#70d68c", "#d9432f", "#6ad5be", "#d5416a", "#76c2d7", "#d87a71",
               "#6a75d5", "#836834", "#c988d1", "#598939")
  

