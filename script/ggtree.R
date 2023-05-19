library(ape)
library(ggplot2)
library("phytools") # for sims and ASRs
library("ggtree")

#extract the tree at Class rank 
phylo_rel <- microbiome::transform(phylo, "compositional")
phylo_rel_core <- microbiome::core(phylo_rel, detection = 0.1/100, prevalence = 7/66)
ntaxa(phylo_rel_core)
core_tree <- phy_tree(phylo_rel_core)
# x1 <- tax_glom(phylo_rare, taxrank = "Genus")
# labels <- as.vector(tax_table(phylo_rel_core)[,6])
labels <- as.vector(rownames(tax_table(phylo_rel_core)))
# add tip labels (making room with xlim first), node labels, background color, 
# branch colors (based on branch legths), and a legend for the branch colors
ggtree(core_tree, branch.length='none', layout='circular', open.angle = 20)






#phylogenetic tree
library(ggtreeExtra)
ps1 <- prune_taxa(taxa_names(phylo_rel_core), phylo)
taxa_sums(ps1)
melt_simple <- data.frame(Phylum = data.frame(tax_table(ps1))[2],
                          Genus = data.frame(tax_table(ps1))[6],
                          ASV = taxa_names(ps1),
                          abundance = taxa_sums(ps1))
f<-phy_tree(ps1)
core.tax = as.data.frame(tax_table(ps1))
# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label = mytree, label2=core.tax$Genus, label3=core.tax$Phylum, stringsAsFactors = F)
#foo$label2[foo$label2 == NA] <- 'uncultured'
# Make a less putrid color palette
library(RColorBrewer)
colourCount = length(unique(core.tax$Phylum))
Mycolors <- colorRampPalette(brewer.pal(11, "Spectral"))(colourCount)
names(Mycolors) <- levels(as.factor(core.tax$Phylum))
colScale <- scale_colour_manual(name = "grp",values = Mycolors)
p <- ggtree(f, branch.length='none', layout='circular', open.angle = 35) %<+% foo + 
  geom_tiplab(size=2, align=TRUE, linesize=0.5, aes(label=label2), offset=5.5, hjust =0) +
  geom_tippoint(size=2, aes(label=label2, group=label3, color=label3) ) + 
  theme_tree2() + 
  scale_color_manual(values=colorRampPalette(brewer.pal(11, "Spectral"))(colourCount))

p_tree <- p + geom_fruit(
  data = melt_simple,
  geom = geom_bar,
  mapping = aes(x = abundance, y = ASV, fill = 'red4'),
  orientation="y",
  stat="identity",
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
p_tree
#edgeR analysis
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
dge = phyloseq_to_edgeR(ps1, group = a$layer)
#The crunching to follow is much easier if metadata is pulled out this way into object "a"
a = sample_data(ps1)
a$layer = as.factor(a$layer)
a$layer = relevel(a$layer, ref="SUR")

#So many things don't understand what kind of variables they need to be. Make sure they understand
layer <- as.numeric(a$layer)
# Design for my linear model
design <-model.matrix(~ layer) 

# EdgeR needs to calculate dispersion again after you've fed it the design.
#This step projectile vomits if chloride or conductivity aren't sqrt transformed.
x = calcNormFactors(dge, method="RLE")
x = estimateGLMCommonDisp(dge, design)
x = estimateGLMTrendedDisp(dge, design)
x = estimateGLMTagwiseDisp(dge, design)

fit <-glmFit(x, design)

# grab the coefficients I care about
lrt <- glmLRT(fit, coef=2)

# lrt to z-scores
table<-lrt$table
table<-apply(table, 2, function(x) scale(x, center = TRUE, scale = TRUE))
rownames(table) = rownames(lrt$table)
table<-as.data.frame(table)

q<-lrt$genes
table1<-cbind(table, Kingdom = q$Kingdom, Phylum = q$Phylum, Class = q$Class, Order = q$Order,
              Family = q$Family, Genus = q$Genus, OTU = rownames(table))

# Sometimes want to do this to check what has at least one z score beyond threshold
table2 <- table1[abs(table1$logFC) > 1.96, ]

#write.csv(table2, "./tables/EdgeR-Zscores-Genera.csv")
melted <- reshape2::melt(table2)
melt2<-subset(melted, abs(melted$value) > 1.96 )

filter <- unique(melt2$OTU)
mini<-prune_taxa((rownames(otu_table(water_physeq)) %in% filter), water_physeq)
f<-phy_tree(mini)
tax = as.data.frame(tax_table(mini))

tax$Genus = gsub("D_5__", "", tax$Genus)
tax$Phylum = gsub("D_1__", "", tax$Phylum)

# change tip labels to genera / need to gsub
library(ggtree)
mytree = phy_tree(f)$tip.label
foo = data.frame(label=mytree, label2=paste(tax$Genus, tax$ASV, sep='_'), label3=tax$Phylum, stringsAsFactors = F)
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

test<-as.matrix(table2[1:20])
tested<-apply(test, 2, function(x){cut(x, br=c(-14, -4, -2, 2, 4, 14))})
rownames(tested) = rownames(table2)
colnames(tested) = gsub("logFC.", "", colnames(tested))
colnames(tested) = gsub("logCPM", "Rotorua", colnames(tested))

heatmapcols <-colorRampPalette(brewer.pal(6, "Spectral"))(6)
names(heatmapcols) <- levels(as.factor(tested[1:10000]))
heatmapcols <- c(heatmapcols[5],heatmapcols[4],heatmapcols[2],heatmapcols[3],heatmapcols[1],heatmapcols[6])
# Make sure to check that the order of these is correct if this is rerunhe
heatmapcols[1] = "#881111"
heatmapcols[2] = "#FFCCCF"
heatmapcols[3] = "#F2F2F2"
heatmapcols[4] = "#AABBDD"
heatmapcols[5] = "#112288"
heatmapcols[6] = "#112288"

mycols <- c("#881111","#D9444D","#FFCCCF","#F2F2F2","#AABBDD","#112288")
tested <- tested[,c('chitinolysis', 'cellulolysis', 'fermentation', 'methanogenesis', 
                    'methanotrophy', 'methylotrophy',  'pH', 'MAP', 'MAT',
                    'Comp4', 'Comp3', 'S275_295', 'SUVA254', 'DOC', 'NH4_N',
                    'Na', 'K', 'Ca', 'Conductivity', 'Depth')] 
p1<-gheatmap(p, tested, offset = 0.25, width = 3, font.size=2.5, 
             colnames_angle=-90, colnames_position = "top", hjust=1) + 
  scale_fill_manual(values=heatmapcols)
p1 <-p1 + theme(legend.position="right")
p1
#ggsave(filename="output/edgeR-Genus.pdf", plot=p, width=8, height=10, units="in")
