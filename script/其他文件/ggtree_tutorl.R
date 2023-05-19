library(ggtree)
set.seed(2015-12-21)
tree <- rtree(30)
p <- ggtree(tree) + xlim(NA, 8)

p + geom_cladelab(node=45, label="test label") +
  geom_cladelab(node=34, label="another clade")
p + geom_cladelab(node=45, label="test label", align=TRUE,  
                  offset = .2, textcolor='red', barcolor='red') +
  geom_cladelab(node=34, label="another clade", align=TRUE, 
                offset = .2, textcolor='blue', barcolor='blue')
p + geom_cladelab(node=45, label="test label", align=TRUE, angle=270, 
                  hjust='center', offset.text=.5, barsize=1.5, fontsize=8) +
  geom_cladelab(node=34, label="another clade", align=TRUE, angle=45)
p + geom_cladelab(node=34, label="another clade", align=TRUE, 
                  geom='label', fill='lightblue')

dat <- data.frame(node = c(45, 34), 
                  name = c("test label", "another clade"))
# The node and label is required when geom="text" 
## or geom="label" or geom="shadowtext".
p1 <- p + geom_cladelab(data = dat, 
                        mapping = aes(node = node, label = name, color = name), 
                        fontsize = 3)

dt <- data.frame(node = c(45, 34), 
                 image = c("7fb9bea8-e758-4986-afb2-95a2c3bf983d", 
                           "0174801d-15a6-4668-bfe0-4c421fbe51e8"), 
                 name = c("specie A", "specie B"))

# when geom="phylopic" or geom="image", the image of aes is required.
p2 <- p + geom_cladelab(data = dt, 
                        mapping = aes(node = node, label = name, image = image), 
                        geom = "phylopic", imagecolor = "black", 
                        offset=1, offset.text=0.5)

# The color or size of image also can be mapped.
p3 <- p + geom_cladelab(data = dt, 
                        mapping = aes(node = node, label = name, 
                                      image = image, color = name), 
                        geom = "phylopic", offset = 1, offset.text=0.5)

ggtree(tree, layout="daylight") + 
  geom_text(aes(label=node), hjust=-.3) +
  geom_cladelab(node=35, label="test label", angle=0, 
                fontsize=8, offset=.5, vjust=.5)  + 
  geom_cladelab(node=55, label='another clade', 
                angle=-95, hjust=.5, fontsize=8)

p + geom_tiplab() +
  geom_strip('t10', 't30', barsize=2, color='red', 
             label="associated taxa", offset.text=.1) + 
  geom_strip('t1', 't18', barsize=2, color='blue', 
             label = "another label", offset.text=.1)

nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
ggtree(tree) + 
  geom_hilight(node=21, fill="steelblue", alpha=.6) +
  geom_hilight(node=17, fill="darkgreen", alpha=.6) 

ggtree(tree, layout="circular") + 
  geom_hilight(node=21, fill="steelblue", alpha=.6) +
  geom_hilight(node=23, fill="darkgreen", alpha=.6)

## type can be 'encircle' or 'rect'
p + geom_hilight(node=55, linetype = 3) + 
  geom_hilight(node=35, fill='darkgreen', type="rect")

ggtree(tree) +
  geom_balance(node=16, fill='steelblue', color='white', alpha=0.6, extend=1) +
  geom_balance(node=19, fill='darkgreen', color='white', alpha=0.6, extend=1) 


## using external data
d <- data.frame(node=c(17, 21), type=c("A", "B"))
ggtree(tree) + geom_hilight(data=d, aes(node=node, fill=type),
                            type = "roundrect")

## using data stored in the tree object
library(treeio)
x <- read.nhx(system.file("extdata/NHX/ADH.nhx", package="treeio"))
ggtree(x) + geom_hilight(mapping=aes(subset = node %in% c(10, 12), 
                                     fill = S),
                         type = "gradient", gradient.direction = 'rt',
                         alpha = .8) +
  scale_fill_manual(values=c("steelblue", "darkgreen"))

# 5.2.3 Taxa connection 
p1 <- ggtree(tree) + geom_tiplab() + geom_taxalink(taxa1='A', taxa2='E') + 
  geom_taxalink(taxa1='F', taxa2='K', color='red', linetype = 'dashed',
                arrow=arrow(length=unit(0.02, "npc")))

p2 <- ggtree(tree, layout="circular") + 
  geom_taxalink(taxa1='A', taxa2='E', color="grey", alpha=0.5, 
                offset=0.05, arrow=arrow(length=unit(0.01, "npc"))) + 
  geom_taxalink(taxa1='F', taxa2='K', color='red', 
                linetype = 'dashed', alpha=0.5, offset=0.05,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1="L", taxa2="M", color="blue", alpha=0.5, 
                offset=0.05, hratio=0.8, 
                arrow=arrow(length=unit(0.01, "npc"))) + 
  geom_tiplab()
# when the tree was created using reverse x, 
# we can set outward to FALSE, which will generate the inward curve lines.
p3 <- ggtree(tree, layout="inward_circular", xlim=150) +
  geom_taxalink(taxa1='A', taxa2='E', color="grey", alpha=0.5, 
                offset=-0.2, outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1='F', taxa2='K', color='red', linetype = 'dashed', 
                alpha=0.5, offset=-0.2, outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_taxalink(taxa1="L", taxa2="M", color="blue", alpha=0.5, 
                offset=-0.2, outward=FALSE,
                arrow=arrow(length=unit(0.01, "npc"))) +
  geom_tiplab(hjust=1) 

dat <- data.frame(from=c("A", "F", "L"), 
                  to=c("E", "K", "M"), 
                  h=c(1, 1, 0.1), 
                  type=c("t1", "t2", "t3"), 
                  s=c(2, 1, 2))
p4 <- ggtree(tree, layout="inward_circular", xlim=c(150, 0)) +
  geom_taxalink(data=dat, 
                mapping=aes(taxa1=from, 
                            taxa2=to, 
                            color=type, 
                            size=s), 
                ncp=10,
                offset=0.15) + 
  geom_tiplab(hjust=1) +
  scale_size_continuous(range=c(1,3))
plot_list(p1, p2, p3, p4, ncol=2, tag_levels='A')

# 5.2.4 Uncertainty of evolutionary inference
file <- system.file("extdata/MEGA7", "mtCDNA_timetree.nex", package = "treeio")
x <- read.mega(file)
p1 <- ggtree(x) + geom_range('reltime_0.95_CI', color='red', size=3, alpha=.3)
p2 <- ggtree(x) + geom_range('reltime_0.95_CI', color='red', size=3, 
                             alpha=.3, center='reltime')  
p3 <- p2 + scale_x_range() + theme_tree2()

#6.1 Viewing Selected Clade
library(ggtree)
nwk <- system.file("extdata", "sample.nwk", package="treeio")
tree <- read.tree(nwk)
p <- ggtree(tree) + geom_tiplab()
viewClade(p, MRCA(p, "I", "L"))

# 6.2 Scaling Selected Clade
tree2 <- groupClade(tree, c(17, 21))
p <- ggtree(tree2, aes(color=group)) + theme(legend.position='none') +
  scale_color_manual(values=c("black", "firebrick", "steelblue"))
scaleClade(p, node=17, scale=.1) 

# 6.3 Collapsing and Expanding Clade
p2 <- p %>% collapse(node=21) + 
  geom_point2(aes(subset=(node==21)), shape=21, size=5, fill='green')
p2 <- collapse(p2, node=23) + 
  geom_point2(aes(subset=(node==23)), shape=23, size=5, fill='red')
print(p2)
expand(p2, node=23) %>% expand(node=21)

# Triangles are often used to represent the collapsed clade and ggtree also supports it.
p2 <- p + geom_tiplab()
node <- 21
collapse(p2, node, 'max') %>% expand(node)
collapse(p2, node, 'min') %>% expand(node)
collapse(p2, node, 'mixed') %>% expand(node)

collapse(p, 21, 'mixed', fill='steelblue', alpha=.4) %>% 
  collapse(23, 'mixed', fill='firebrick', color='blue')
scaleClade(p, 23, .2) %>% collapse(23, 'min', fill="darkgreen")  

# 6.4 Grouping Taxa
data(iris)
rn <- paste0(iris[,5], "_", 1:150)
rownames(iris) <- rn
d_iris <- dist(iris[,-5], method="man")

tree_iris <- ape::bionj(d_iris)
grp <- list(setosa     = rn[1:50],
            versicolor = rn[51:100],
            virginica  = rn[101:150])

p_iris <- ggtree(tree_iris, layout = 'circular', branch.length='none')
groupOTU(p_iris, grp, 'Species') + aes(color=Species) +
  theme(legend.position="right")

tree_iris <- groupOTU(tree_iris, grp, "Species")
ggtree(tree_iris, aes(color=Species), layout = 'circular', 
       branch.length = 'none') + 
  theme(legend.position="right")

# 7.1 Mapping Data to The tree Structure
library(ggimage)
library(ggtree)
library(TDbook)

# load `tree_boots`, `df_tip_data`, and `df_inode_data` from 'TDbook'
p <- ggtree(tree_boots) %<+% df_tip_data + xlim(-.1, 4)
p2 <- p + geom_tiplab(offset = .6, hjust = .5) +
  geom_tippoint(aes(shape = trophic_habit, color = trophic_habit, 
                    size = mass_in_kg)) + 
  theme(legend.position = "right") + 
  scale_size_continuous(range = c(3, 10))

p2 %<+% df_inode_data + 
  geom_label(aes(label = vernacularName.y, fill = posterior)) + 
  scale_fill_gradientn(colors = RColorBrewer::brewer.pal(3, "YlGnBu"))

# 7.2 Aligning Graph to the Tree Based on the Tree Structure
library(ggtree)
library(TDbook)

## load `tree_nwk`, `df_info`, `df_alleles`, and `df_bar_data` from 'TDbook'
tree <- tree_nwk
snps <- df_alleles
snps_strainCols <- snps[1,] 
snps<-snps[-1,] # drop strain names
colnames(snps) <- snps_strainCols

gapChar <- "?"
snp <- t(snps)
lsnp <- apply(snp, 1, function(x) {
  x != snp[1,] & x != gapChar & snp[1,] != gapChar
})
lsnp <- as.data.frame(lsnp)
lsnp$pos <- as.numeric(rownames(lsnp))
lsnp <- tidyr::gather(lsnp, name, value, -pos)
snp_data <- lsnp[lsnp$value, c("name", "pos")]

## visualize the tree 
p <- ggtree(tree) 

## attach the sampling information data set 
## and add symbols colored by location
p <- p %<+% df_info + geom_tippoint(aes(color=location))

## visualize SNP and Trait data using dot and bar charts,
## and align them based on tree structure
p + geom_facet(panel = "SNP", data = snp_data, geom = geom_point, 
               mapping=aes(x = pos, color = location), shape = '|') +
  geom_facet(panel = "Trait", data = df_bar_data, geom = geom_col, 
             aes(x = dummy_bar_value, color = location, 
                 fill = location), orientation = 'y', width = .6) +
  theme_tree2(legend.position=c(.05, .85))


# 7.3 Visualize a Tree with an Associated Matrix
beast_file <- system.file("examples/MCC_FluA_H3.tree", package="ggtree")
beast_tree <- read.beast(beast_file)

genotype_file <- system.file("examples/Genotype.txt", package="ggtree")
genotype <- read.table(genotype_file, sep="\t", stringsAsFactor=F)
colnames(genotype) <- sub("\\.$", "", colnames(genotype))
p <- ggtree(beast_tree, mrsd="2013-01-01") + 
  geom_treescale(x=2008, y=1, offset=2) + 
  geom_tiplab(size=2)
gheatmap(p, genotype, offset=5, width=0.5, font.size=3, 
         colnames_angle=-45, hjust=0) +
  scale_fill_manual(breaks=c("HuH3N2", "pdm", "trig"), 
                    values=c("steelblue", "firebrick", "darkgreen"), name="genotype")

p <- ggtree(beast_tree, mrsd="2013-01-01") + 
  geom_tiplab(size=2, align=TRUE, linesize=.5) + 
  theme_tree2()
gheatmap(p, genotype, offset=8, width=0.6, 
         colnames=FALSE, legend_title="genotype") +
  scale_x_ggtree() + 
  scale_y_continuous(expand=c(0, 0.3))
# 7.3.1 Visualize a tree with multiple associated matrices
nwk <- system.file("extdata", "sample.nwk", package="treeio")

tree <- read.tree(nwk)
circ <- ggtree(tree, layout = "circular")

df <- data.frame(first=c("a", "b", "a", "c", "d", "d", "a", 
                         "b", "e", "e", "f", "c", "f"),
                 second= c("z", "z", "z", "z", "y", "y", 
                           "y", "y", "x", "x", "x", "a", "a"))
rownames(df) <- tree$tip.label

df2 <- as.data.frame(matrix(rnorm(39), ncol=3))
rownames(df2) <- tree$tip.label
colnames(df2) <- LETTERS[1:3]


p1 <- gheatmap(circ, df, offset=.8, width=.2,
               colnames_angle=95, colnames_offset_y = .25) +
  scale_fill_viridis_d(option="D", name="discrete\nvalue")


library(ggnewscale)
p2 <- p1 + new_scale_fill()
gheatmap(p2, df2, offset=15, width=.3,
         colnames_angle=90, colnames_offset_y = .25) +
  scale_fill_viridis_c(option="A", name="continuous\nvalue")

# 7.4 Visualize a Tree with Multiple Sequence Alignments
library(TDbook)

# load `tree_seq_nwk` and `AA_sequence` from 'TDbook'
p <- ggtree(tree_seq_nwk) + geom_tiplab(size=3)
msaplot(p, AA_sequence, offset=3, width=2)

p <- ggtree(tree_seq_nwk, layout='circular') + 
  geom_tiplab(offset=4, align=TRUE) + xlim(NA, 12)
msaplot(p, AA_sequence, window=c(120, 200))

# 7.5 Composite Plots
library(ggplot2)
library(ggtree)

set.seed(2019-10-31)
tr <- rtree(10)

d1 <- data.frame(
  # only some labels match
  label = c(tr$tip.label[sample(5, 5)], "A"),
  value = sample(1:6, 6))

d2 <- data.frame(
  label = rep(tr$tip.label, 5),
  category = rep(LETTERS[1:5], each=10),
  value = rnorm(50, 0, 3)) 

g <- ggtree(tr) + geom_tiplab(align=TRUE) + hexpand(.01)

p1 <- ggplot(d1, aes(label, value)) + geom_col(aes(fill=label)) + 
  geom_text(aes(label=label, y= value+.1)) +
  coord_flip() + theme_tree2() + theme(legend.position='none')

p2 <- ggplot(d2, aes(x=category, y=label)) + 
  geom_tile(aes(fill=value)) + scale_fill_viridis_c() + 
  theme_minimal() + xlab(NULL) + ylab(NULL)

#If we align them using cowplot, the composite plots are not aligned properly as we anticipated (Figure 7.6A).

cowplot::plot_grid(g, p2, p1, ncol=3) 

#Using aplot, it will do all the dirty work for us and all the subplots are aligned properly as demonstrated in Figure 7.6B.

library(aplot)
p2 %>% insert_left(g) %>% insert_right(p1, width=.5) 

# ggtreeExtra for Presenting Data on a Circular Layout
library(ggtreeExtra)
library(ggtree)
library(phyloseq)
library(dplyr)

data("GlobalPatterns")
GP <- GlobalPatterns
GP <- prune_taxa(taxa_sums(GP) > 600, GP)
sample_data(GP)$human <- get_variable(GP, "SampleType") %in%
  c("Feces", "Skin")
mergedGP <- merge_samples(GP, "SampleType")
mergedGP <- rarefy_even_depth(mergedGP,rngseed=394582)
mergedGP <- tax_glom(mergedGP,"Order")

melt_simple <- psmelt(mergedGP) %>%
  filter(Abundance < 120) %>%
  select(OTU, val=Abundance)

p <- ggtree(mergedGP, layout="fan", open.angle=10) + 
  geom_tippoint(mapping=aes(color=Phylum), 
                size=1.5,
                show.legend=FALSE)
p <- rotate_tree(p, -90)

p <- p +
  geom_fruit(
    data=melt_simple,
    geom=geom_boxplot,
    mapping = aes(
      y=OTU,
      x=val,
      group=label,
      fill=Phylum,
    ),
    size=.2,
    outlier.size=0.5,
    outlier.stroke=0.08,
    outlier.shape=21,
    axis.params=list(
      axis       = "x",
      text.size  = 1.8,
      hjust      = 1,
      vjust      = 0.5,
      nbreak     = 3,
    ),
    grid.params=list()
  ) 

p <- p +
  scale_fill_discrete(
    name="Phyla",
    guide=guide_legend(keywidth=0.8, keyheight=0.8, ncol=1)
  ) +
  theme(
    legend.title=element_text(size=9), 
    legend.text=element_text(size=7) 
  )
p

#  10.3 Aligning Multiple Graphs to the Tree for Multi-dimensional Data 
library(ggtreeExtra)
library(ggtree)
library(treeio)
library(tidytree)
library(ggstar)
library(ggplot2)
library(ggnewscale)
library(TDbook)

# load data from TDbook, including tree_hmptree, 
# df_tippoint (the abundance and types of microbes),
# df_ring_heatmap (the abundance of microbes at different body sites),
# and df_barplot_attr (the abundance of microbes of greatest prevalence)
tree <- tree_hmptree
dat1 <- df_tippoint
dat2 <- df_ring_heatmap
dat3 <- df_barplot_attr

# adjust the order
dat2$Sites <- factor(dat2$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)","Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
dat3$Sites <- factor(dat3$Sites, 
                     levels=c("Stool (prevalence)", "Cheek (prevalence)",
                              "Plaque (prevalence)", "Tongue (prevalence)",
                              "Nose (prevalence)", "Vagina (prevalence)",
                              "Skin (prevalence)"))
# extract the clade label information. Because some nodes of tree are
# annotated to genera, which can be displayed with high light using ggtree.
nodeids <- nodeid(tree, tree$node.label[nchar(tree$node.label)>4])
nodedf <- data.frame(node=nodeids)
nodelab <- gsub("[\\.0-9]", "", tree$node.label[nchar(tree$node.label)>4])
# The layers of clade and hightlight
poslist <- c(1.6, 1.4, 1.6, 0.8, 0.1, 0.25, 1.6, 1.6, 1.2, 0.4,
             1.2, 1.8, 0.3, 0.8, 0.4, 0.3, 0.4, 0.4, 0.4, 0.6,
             0.3, 0.4, 0.3)
labdf <- data.frame(node=nodeids, label=nodelab, pos=poslist)

# The circular layout tree.
p <- ggtree(tree, layout="fan", size=0.15, open.angle=5) +
  geom_hilight(data=nodedf, mapping=aes(node=node),
               extendto=6.8, alpha=0.3, fill="grey", color="grey50",
               size=0.05) +
  geom_cladelab(data=labdf, 
                mapping=aes(node=node, 
                            label=label,
                            offset.text=pos),
                hjust=0.5,
                angle="auto",
                barsize=NA,
                horizontal=FALSE, 
                fontsize=1.4,
                fontface="italic"
  )

p <- p %<+% dat1 + geom_star(
  mapping=aes(fill=Phylum, starshape=Type, size=Size),
  position="identity",starstroke=0.1) +
  scale_fill_manual(values=c("#FFC125","#87CEFA","#7B68EE","#808080",
                             "#800080", "#9ACD32","#D15FEE","#FFC0CB",
                             "#EE6A50","#8DEEEE", "#006400","#800000",
                             "#B0171F","#191970"),
                    guide=guide_legend(keywidth = 0.5, 
                                       keyheight = 0.5, order=1,
                                       override.aes=list(starshape=15)),
                    na.translate=FALSE)+
  scale_starshape_manual(values=c(15, 1),
                         guide=guide_legend(keywidth = 0.5, 
                                            keyheight = 0.5, order=2),
                         na.translate=FALSE)+
  scale_size_continuous(range = c(1, 2.5),
                        guide = guide_legend(keywidth = 0.5, 
                                             keyheight = 0.5, order=3,
                                             override.aes=list(starshape=15)))

p <- p + new_scale_fill() +
  geom_fruit(data=dat2, geom=geom_tile,
             mapping=aes(y=ID, x=Sites, alpha=Abundance, fill=Sites),
             color = "grey50", offset = 0.04,size = 0.02)+
  scale_alpha_continuous(range=c(0, 1),
                         guide=guide_legend(keywidth = 0.3, 
                                            keyheight = 0.3, order=5)) +
  geom_fruit(data=dat3, geom=geom_bar,
             mapping=aes(y=ID, x=HigherAbundance, fill=Sites),
             pwidth=0.38, 
             orientation="y", 
             stat="identity",
  ) +
  scale_fill_manual(values=c("#0000FF","#FFA500","#FF0000",
                             "#800000", "#006400","#800080","#696969"),
                    guide=guide_legend(keywidth = 0.3, 
                                       keyheight = 0.3, order=4))+
  geom_treescale(fontsize=2, linesize=0.3, x=4.9, y=0.1) +
  theme(legend.position=c(0.93, 0.5),
        legend.background=element_rect(fill=NA),
        legend.title=element_text(size=6.5),
        legend.text=element_text(size=4.5),
        legend.spacing.y = unit(0.02, "cm"),
  )
p

p + layout_rectangular() + 
  theme(legend.position=c(.05, .7))
