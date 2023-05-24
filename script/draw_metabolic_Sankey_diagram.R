# Create directory to hold plots ---------------------------
sankey.plots.folder <- file.path(wd_fun, "metabolic/results/sankey") # Name of new output folder
if (!dir.exists(sankey.plots.folder)) {
  dir.create(sankey.plots.folder)
}


#  Generate a network plot with reactions as node, and edges are MAGs, and are colored by taxonomic group.
sur_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_sur/METABOLIC_Figures_Input/Metabolic_Sankey_diagram_input.txt")
sub_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_sub/METABOLIC_Figures_Input/Metabolic_Sankey_diagram_input.txt")
pl_input_table <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_pl/METABOLIC_Figures_Input/Metabolic_Sankey_diagram_input.txt")


# Load packages ---------------------------
library(ggalluvial)
library(ggthemes)

# Load the energy flow diagram input ---------------------------
sankey.plots.folder <- file.path(wd_fun, "metabolic/results/sankey")
table_sur <- read.table(sur_input_table, header = F, sep = "\t")
colnames(table_sur) <- c("Taxa", "Reaction", "Freq")
table_sub <- read.table(sub_input_table, header = F, sep = "\t")
colnames(table_sub) <- c("Taxa", "Reaction", "Freq")
table_pl <- read.table(pl_input_table, header = F, sep = "\t")
colnames(table_pl) <- c("Taxa", "Reaction", "Freq")

# Rename the categories ---------------------------
table_sur$Category <- ifelse(grepl("C-S", table_sur$Reaction), "Carbon", 
                               ifelse(grepl("N-S", table_sur$Reaction), "Nitrogen",
                                      ifelse(grepl("S-S", table_sur$Reaction), "Sulfur",
                                             ifelse(grepl("O-S", table_sur$Reaction), "Others",""))))
table_sur$Category <- as.factor(table_sur$Category)

table_sub$Category <- ifelse(grepl("C-S", table_sub$Reaction), "Carbon", 
                               ifelse(grepl("N-S", table_sub$Reaction), "Nitrogen",
                                      ifelse(grepl("S-S", table_sub$Reaction), "Sulfur",
                                             ifelse(grepl("O-S", table_sub$Reaction), "Others",""))))
table_sub$Category <- as.factor(table_sub$Category)

table_pl$Category <- ifelse(grepl("C-S", table_pl$Reaction), "Carbon", 
                               ifelse(grepl("N-S", table_pl$Reaction), "Nitrogen",
                                      ifelse(grepl("S-S", table_pl$Reaction), "Sulfur",
                                             ifelse(grepl("O-S", table_pl$Reaction), "Others",""))))
table_pl$Category <- as.factor(table_pl$Category)

table_all <- rbind(cbind(Layer = rep("SUR", nrow(table_sur)), table_sur),
                   cbind(Layer = rep("SUB", nrow(table_sub)), table_sub),
                   cbind(Layer = rep("PL", nrow(table_pl)), table_pl)) %>%
  mutate(Layer = factor(Layer, levels = c("SUR", "SUB", "PL")))

taxa_names <- table_all %>% select(Taxa, Freq) %>%
  group_by(Taxa) %>%
  summarise(sum = sum(Freq)) %>%
  arrange(desc(sum)) %>%
  pull(Taxa)
reaction_orders <- c(sort(grep('C-', unique(table_all$Reaction), value = T), decreasing = F),
                     sort(grep('N-', unique(table_all$Reaction), value = T), decreasing = F),
                     sort(grep('S-S', unique(table_all$Reaction), value = T), decreasing = F),
                     sort(grep('O-', unique(table_all$Reaction), value = T), decreasing = F))

table_all <- table_all %>%
  mutate(Taxa = factor(Taxa, levels = taxa_names)) %>%
  mutate(Category = factor(Category, levels = c('Carbon', 'Nitrogen', 
                                                'Sulfur', 'Others'))) %>%
  mutate(Reaction = factor(Reaction, levels = reaction_orders))
# Load packages and confirm format of table ---------------------------
is_alluvia_form(as.data.frame(table_sur), axes = 1:3, silent = TRUE)
is_alluvia_form(as.data.frame(table_sub), axes = 1:3, silent = TRUE)
is_alluvia_form(as.data.frame(table_pl), axes = 1:3, silent = TRUE)
is_alluvia_form(as.data.frame(table_all), axes = 1:4, silent = TRUE)

# Create plot ---------------------------
all_color <- c("#d64e9e", "#6cd54c", "#dd49d1", "#c8dd41", "#a152dd", "#4d4040",
               "#5139c2", "#ceaa3b", "#7a3260", "#432d7c", "#c6d179", "#8f379a",
               "#70d68c", "#d9432f", "#6ad5be", "#d5416a", "#76c2d7", "#d87a71",
               "#6a75d5", "#836834", "#c988d1", "#598939")

all_color <- c( "#FF7F00", "#B3DE69", "#1F78B4", "#FB8072", "#FDB462", "#80B1D3", "#BEBADA", "#c988d1", "#FCCDE5", "#A6CEE3",
"#8DD3C7", "#B2DF8A", "#33A02C", "#FB9A99", "#FFFFB3", "#FDBF6F", "#E31A1C", "#d5416a", "#76c2d7", "#d87a71",
"#6a75d5", "#836834", "#598939")
alluvial.plot.all <- ggplot(as.data.frame(table_all),
                            aes(y = Freq, axis1= Layer, axis2= Taxa, 
                                axis3 = Reaction, axis4 = Category)) +
  geom_alluvium(aes(fill = Taxa), width = 1/12)+
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_text(stat = "stratum", infer.label = TRUE) +
  # geom_label(stat = "stratum", infer.label = TRUE) +
  # label.strata was deprecated so I changed it to infer.label
  # scale_x_discrete(limits = c("Layer", "Taxa", "Reaction", "Category"), 
  #                  expand = c(.05, .05)) +
  scale_fill_manual(values = all_color) +
  ggtitle("Reactions and taxa") +
  theme_bw() +
  guides(fill = guide_legend(ncol = 1))

alluvial.plot.sur <- ggplot(as.data.frame(table_sur),
       aes(y = Freq, axis1= Taxa, axis2 = Reaction, axis3 = Category)) +
  geom_alluvium(aes(fill = Taxa), width = 1/12)+
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
# label.strata was deprecated so I changed it to infer.label
  scale_x_discrete(limits = c("Reaction", "Taxa","Category"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Reactions and taxa")+
  theme_bw()

alluvial.plot.sub <- ggplot(as.data.frame(table_sub),
                            aes(y = Freq, axis1= Taxa, axis2 = Reaction, axis3 = Category)) +
  geom_alluvium(aes(fill = Taxa), width = 1/12)+
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
  # label.strata was deprecated so I changed it to infer.label
  scale_x_discrete(limits = c("Reaction", "Taxa","Category"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Reactions and taxa")+
  theme_bw()

alluvial.plot.pl <- ggplot(as.data.frame(table_pl),
                            aes(y = Freq, axis1= Taxa, axis2 = Reaction, axis3 = Category)) +
  geom_alluvium(aes(fill = Taxa), width = 1/12)+
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  geom_label(stat = "stratum", infer.label = TRUE) +
  # label.strata was deprecated so I changed it to infer.label
  scale_x_discrete(limits = c("Reaction", "Taxa","Category"), expand = c(.05, .05)) +
  #scale_fill_brewer(type = "qual", palette = "Set1") +
  ggtitle("Reactions and taxa")+
  theme_bw()

# Save plot ---------------------------
plot.name <- paste(sankey.plots.folder ,"/","sankey.all_2.plot.pdf", sep="")
pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
alluvial.plot.all
dev.off()

plot.name <- paste(sankey.plots.folder ,"/","sankey.sur.plot.pdf", sep="")
pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
alluvial.plot.sur
dev.off()

plot.name <- paste(sankey.plots.folder ,"/","sankey.sub.plot.pdf", sep="")
pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
alluvial.plot.sub
dev.off()

plot.name <- paste(sankey.plots.folder ,"/","sankey.pl.plot.pdf", sep="")
pdf(file = plot.name, width = 11, height = 8.5, onefile=FALSE)
alluvial.plot.pl
dev.off()


#
remotes::install_github("davidsjoberg/ggsankey")
library(ggsankey)
df <- table_all %>%
  +     make_long(Layer, Taxa, Reaction, Category)

library(ggplot2)
library(dplyr)
ggplot(df, aes(x = x,
               next_x = next_x,
               node = node,
               next_node = next_node,
               fill = factor(node))) +
  geom_sankey() +
  theme_sankey(base_size = 16) +
  theme(legend.position = "none")
