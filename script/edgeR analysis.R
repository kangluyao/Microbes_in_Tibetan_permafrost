save.dir <- file.path(getwd(),"result")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}
# 1. alternately edgeR analysi.
library("edgeR")

#The field in the class definition file that defines the classes of the data.
data_classes <- "layer"

# Load Expression Data
RNASeq <- ko_count_table

# Load subtype information
RNASeq[1:5, 1:5]
classDefinitions_RNASeq <- meta_dat
classDefinitions_RNASeq[1:5, 1:3]

# Filter Data (RNA-seq read counts are converted to CPM values 
# and genes with CPM > 1 in at least 50 of the samples are 
# retained for further study, a gene mush have at least 50 measurements 
# with more than 1 CPM in one of the classes to be included in the analysis)
cpms <- cpm(RNASeq)
keep <- rowSums(cpms > 1) >= 1
counts <- RNASeq[keep,]

# Normalization and Dispersion
# create data structure to hold counts and subtype information for each sample.
d <- DGEList(counts = counts, group = classDefinitions_RNASeq$SUBTYPE)
# Normalize the data
d <- calcNormFactors(d)

# create multidimensional scaling(MDS) plot. The command below will automatically
# generate the plot containing all samples where each subtype is a different color.
# Ideally there should be a good separation between the different classes.
mds_filename <- file.path(save.dir, "figs/metagenome/mdsplot_allsamples.png")
png(filename = mds_filename)
mds_output <- plotMDS(d, labels = NULL, pch = 1,
                      col= c("darkgreen", "red", "orange")[factor(classDefinitions_RNASeq$layer)],
                      xlim = c(-2.5,4), ylim = c(-2.5,4))
legend("topright",
       legend = levels(factor(classDefinitions_RNASeq$layer)),
       pch = c(1), col = c("darkgreen","red", "orange"), title = "Class",
       bty = 'n', cex = 0.75)
dev.off()

#calculate dispersion
d <- estimateCommonDisp(d)
d <- estimateTagwiseDisp(d)

#compare SUR layer to the remaining two layers.
classes <- factor(classDefinitions_RNASeq$layer)
modelDesign <- model.matrix(~ 0 + classes)
contrast_sur <- makeContrasts(
  survsrest = "classesSUR-(classesSUB + classesPL)/2", levels = modelDesign)
fit_glm <- glmFit(d, modelDesign)
survsrest <- glmLRT(fit_glm , contrast = contrast_sur)
tt_survsrest <- topTags(survsrest, n = nrow(d))
# Create g:Profiler input list
# tt <- tt_exact_test
#get the indices of scored dataset that have FDR < 0.05
select_sur_genes = rownames(tt_survsrest)[which(tt_survsrest$table$FDR < 0.01 & tt_survsrest$table$logFC > 1)]
#output how many genes there are in the set that have FDR < 0.05
length(select_sur_genes)

#compare SUB layer to the remaining two layers.
contrast_sub <- makeContrasts(
  subvsrest = "classesSUB-(classesSUR + classesPL)/2", levels = modelDesign)
fit_glm <- glmFit(d, modelDesign)
subvsrest <- glmLRT(fit_glm , contrast = contrast_sub)
tt_subvsrest <- topTags(subvsrest, n = nrow(d))
#get the indices of scored dataset that have FDR < 0.05
select_sub_genes = rownames(tt_subvsrest)[which(tt_subvsrest$table$FDR < 0.01 & tt_subvsrest$table$logFC > 1)]
#output how many genes there are in the set that have FDR < 0.05
length(select_sub_genes)

#compare PL layer to the remaining two layers.
contrast_pl <- makeContrasts(
  plvsrest = "classesPL-(classesSUR + classesSUB)/2", levels = modelDesign)
fit_glm <- glmFit(d, modelDesign)
plvsrest <- glmLRT(fit_glm , contrast = contrast_pl)
tt_plvsrest <- topTags(plvsrest, n = nrow(d))
#get the indices of scored dataset that have FDR < 0.05
select_pl_genes = rownames(tt_plvsrest)[which(tt_plvsrest$table$FDR < 0.01 & tt_plvsrest$table$logFC > 1)]
#output how many genes there are in the set that have FDR < 0.05
length(select_pl_genes)

length(intersect(select_sur_genes, select_sub_genes))
length(intersect(select_sur_genes, select_pl_genes))
length(intersect(select_sub_genes, select_pl_genes))
############################################################################################################
# 2. Volcano plot

############################################################################################################
# 3. KEGG decoder
kegg_decoder_fun <- function(layer, kegg_sig) {
  samples <- grep(layer, meta_dat$sample_id, value = T)
  kegg_decoder_df = data.frame(Sample = c(), KEGG = c())
  for (kegg_id in kegg_sig) {
    for (sample in samples) {
      if (ko_count_table[rownames(ko_count_table) == kegg_id, sample] > 0) {
        kegg_decoder_df = rbind(kegg_decoder_df, data.frame(Sample = sample, KEGG = kegg_id))
      }
    }
  }
  return(kegg_decoder_df)
}
kegg_decoder_sur_df <- kegg_decoder_fun('SUR', select_sur_genes)
kegg_decoder_sub_df <- kegg_decoder_fun('SUB', select_sub_genes)
kegg_decoder_pl_df <- kegg_decoder_fun('PL', select_pl_genes)

kegg_decoder_sur_df$Sample <- sapply(kegg_decoder_sur_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_sub_df$Sample <- sapply(kegg_decoder_sub_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_pl_df$Sample <- sapply(kegg_decoder_pl_df$Sample, FUN = function(x){sub('_', '', x)})

# write.table(kegg_decoder_sur_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_sur_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_sub_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_sub_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_pl_df,
#             file.path(save.dir, 'tables/kegg/edgeR/KEGG_decoder_pl_df_in.tsv'), sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)

# in the Data folder using the CryoBiome-kegg-decoder env: 
# python KEGG-decoder.py -i KEGG_decoder_sur_df_in.tsv -o KEGG_decoder_sur_df_output.tsv
# python KEGG-decoder.py -i KEGG_decoder_sub_df_in.tsv -o KEGG_decoder_sub_df_output.tsv
# python KEGG-decoder.py -i KEGG_decoder_pl_df_in.tsv -o KEGG_decoder_pl_df_output.tsv
# python -m pip install xxxx
KEGG_decoder_sur <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_sur_df_output.tsv"), 
                               sep = "\t", header = TRUE)
KEGG_decoder_sub <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_sub_df_output.tsv"),
                               sep = "\t", header = TRUE)
KEGG_decoder_pl <- read.table(file.path(save.dir, "tables/kegg/edgeR/KEGG_decoder_pl_df_output.tsv"),
                              sep = "\t", header = TRUE)
KEGG_decoder <- rbind(KEGG_decoder_sur, KEGG_decoder_sub, KEGG_decoder_pl)

colsums <- colSums(KEGG_decoder[, !colnames(KEGG_decoder) %in% c("Function")])
KEGG_decoder <- KEGG_decoder[, c("Function", names(colsums[colsums > 0]))]
names(KEGG_decoder)[names(KEGG_decoder) == "Function"] <- "Sample"
KEGG_decoder <- reshape2::melt(KEGG_decoder, variable.name = "Function", value.name = "Completion")

# Renaming pathways for plotting
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "TCA.Cycle"] <- "TCA Cycle"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "NAD.P.H.quinone.oxidoreductase"] <- "NAD(P)H-quinone oxidoreductase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "X3.Hydroxypropionate.Bicycle"] <- "3-Hydroxypropionate Bicycle"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "basic.endochitinase.B"] <- "Basic endochitinase B"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "DMS.dehydrogenase"] <- "DMS dehydrogenase"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "thiamin.biosynthesis"] <- "Thiamin biosynthesis"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "transporter..vitamin.B12"] <- "Transporter: vitamin B12"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Flagellum"] <- "Flagellum"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Photosystem.II"] <- "Photosystem II"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Photosystem.I"] <- "Photosystem I"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Cytochrome.b6.f.complex"] <- "Cytochrome b6/f complex"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Entner.Doudoroff.Pathway"] <- "Entner-Doudoroff Pathway"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Mixed.acid..Ethanol..Acetate.to.Acetylaldehyde"] <- "Mixed acid: Ethanol, Acetate to Acetylaldehyde"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Mixed.acid..Ethanol..Acetyl.CoA.to.Acetylaldehyde..reversible."] <- "Mixed acid: Ethanol, Acetyl-CoA to Acetylaldehyde (reversible)"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Naphthalene.degradation.to.salicylate"] <- "Naphthalene degradation to salicylate"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Curli.fimbriae.biosynthesis"] <- "Curli fimbriae biosynthesis"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Competence.related.core.components"] <- "Competence-related core components"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "CP.lyase.operon"] <- "CP-lyase operon"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.III.Secretion"] <- "Type III Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.IV.Secretion"] <- "Type IV Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "Type.Vabc.Secretion"] <- "Type Vabc Secretion"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "serine"] <- "Serine"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "tyrosine"] <- "Tyrosine"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "starch.glycogen.degradation"] <- "Starch/glycogen degradation"
levels(KEGG_decoder$Function)[levels(KEGG_decoder$Function) == "end.product.astaxanthin"] <- "End-product astaxanthin"

KEGG_decoder$mean_comp = vapply(1:nrow(KEGG_decoder), function(x) mean(KEGG_decoder$Completion[KEGG_decoder$Function == KEGG_decoder$Function[x]]), FUN.VALUE = numeric(1))
KEGG_decoder$Function <- factor(KEGG_decoder$Function, levels=unique((KEGG_decoder$Function)[order(KEGG_decoder$mean_comp)]))
meta_dat$Sample <- sapply(meta_dat$sample_id, FUN = function(x){sub('_', '', x)})
KEGG_decoder$layer = sapply(1:nrow(KEGG_decoder), function(x) meta_dat$layer[meta_dat$Sample == KEGG_decoder$Sample[x]])
KEGG_decoder$Sample <- factor(KEGG_decoder$Sample, levels = meta_dat$Sample, ordered = T)
KEGG_decoder$layer <- factor(KEGG_decoder$layer, levels = c('SUR', 'SUB', 'PL'), ordered = T)

# plot
library(ggplot2)
p <- ggplot(KEGG_decoder) + 
  geom_tile(aes(x = as.factor(Sample), y = as.factor(Function), fill = Completion)) +
  xlab("") + ylab("") + 
  theme_linedraw() + 
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 10)) +
  scale_fill_gradientn(colours = c("white", "#00A087FF"), values = c(0, 1)) + 
  facet_grid(~layer, scales = "free_x", space = "free_x") + 
  theme(axis.title.x = element_blank(), panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5), 
        axis.text.y = element_text(size = 6), 
        panel.spacing = unit(0, "lines"),
        legend.title = element_text(size = 6),
        legend.text = element_text(size = 6),
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave(file.path(save.dir, '/figs/metagenome/kegg_enrichment_edgeR.pdf'), 
       p, width = 10, height = 5)

############################################################################################################
# 4. enrichment analysis with clusterProfiler package
library(clusterProfiler)
KO_pathway <- read.csv("e:/permafrost/data/metagenome/fun/keggKO_ko_metabolism.csv", header = T)
termgene <- KO_pathway[, c(2, 1)]
termname <- KO_pathway[, c(2, 3)]
enrichpathway_sur <- enricher(select_sur_genes, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_sub <- enricher(select_sub_genes, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_pl <- enricher(select_pl_genes, TERM2GENE = termgene, TERM2NAME = termname)

write.csv(enrichpathway_sur, file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_sur_edgeR.csv"))
write.csv(enrichpathway_sub,  file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_sub_edgeR.csv"))
write.csv(enrichpathway_pl,  file.path(save.dir, "tables/kegg/clusterProfiler/enrichpathway_pl_edgeR.csv"))
