
#set work directory
save.dir <- file.path(getwd(),"result")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}
setwd(save.dir)
#loading packages
library(tidyverse)
library(EnhancedVolcano)
library(DESeq2)

############################################################################################################
# 1. DESeq2 analysis
meta_dat <- metadata
meta_dat$layer <- as.factor(meta_dat$layer)
KEGG_dds <- DESeqDataSetFromMatrix(countData=round(ko_count_table + 1), colData=meta_dat, design=~layer)
KEGG_deseq <- DESeq(KEGG_dds)
#t treatment, c control
kegg_sig_fun <- function(t, c) {
  KEGG_res <- results(KEGG_deseq, contrast = c("layer", t, c))
  KEGG_res$padj[is.na(KEGG_res$padj)] = 1
  KEGG_sig_in_t = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]
  KEGG_sig_in_c = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange < -1)]
  KEGG_sig <- list(KEGG_sig_in_t, KEGG_sig_in_c)
  names(KEGG_sig) <- c(t, c)
  return(KEGG_sig)
}
kegg_sig_pl_sur <- kegg_sig_fun('PL', 'SUR')
kegg_sig_pl_sub <- kegg_sig_fun('PL', 'SUB')
kegg_sig_sub_sur <- kegg_sig_fun('SUB', 'SUR')

kegg_sig_in_sur <- intersect(kegg_sig_pl_sur[['SUR']], kegg_sig_sub_sur[['SUR']])
length(kegg_sig_in_sur)
write.table(kegg_sig_in_sur, file.path(save.dir, './tables/kegg/DESeq/KOs_sig_in_sur.txt'),
            row.names = F, col.names = F, quote = F)
kegg_sig_in_sub <- intersect(kegg_sig_sub_sur[['SUB']], kegg_sig_pl_sub[['SUB']])
length(kegg_sig_in_sub)
write.table(kegg_sig_in_sub, file.path(save.dir, './tables/kegg/DESeq/KOs_sig_in_sub.txt'),
            row.names = F, col.names = F, quote = F)
kegg_sig_in_pl <- intersect(kegg_sig_pl_sur[['PL']], kegg_sig_pl_sub[['PL']])
length(kegg_sig_in_pl)
write.table(kegg_sig_in_pl, file.path(save.dir, './tables/kegg/DESeq/KOs_sig_in_pl.txt'),
            row.names = F, col.names = F, quote = F)
# KEGG_significant <- c(kegg_sig_in_sur, kegg_sig_in_sub, kegg_sig_in_pl)
# KEGG_significant <- kegg_sig_in_sub
# KEGG_res <- results(KEGG_deseq, contrast = c("layer", "SUB", "SUR"))
# KEGG_res$padj[is.na(KEGG_res$padj)] = 1
# KEGG_significant = rownames(KEGG_res)[(KEGG_res$padj < 0.05) & (KEGG_res$log2FoldChange > 1)]
# length(KEGG_significant)
# write.csv(as.data.frame(KEGG_res), file = 'Data/KEGG_deseq_results.csv')

############################################################################################################
# 2. DESeq2 volcano plot
# EnhancedVolcano(KEGG_res,
#                 lab = rownames(KEGG_res),
#                 pCutoff = 0.05,
#                 FCcutoff = 1,
#                 col = c("grey", "grey30", "grey30", "#23A671"),
#                 x = 'log2FoldChange',
#                 title = NULL,
#                 subtitle = NULL,
#                 caption = NULL,
#                 y = 'padj')
# ggsave('3_Functional_analysis/3_1_Enrichment/KEGG_enriched_layer.pdf', width = 7, height = 7)

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
kegg_decoder_sur_df <- kegg_decoder_fun('SUR', kegg_sig_in_sur)
kegg_decoder_sub_df <- kegg_decoder_fun('SUB', kegg_sig_in_sub)
kegg_decoder_pl_df <- kegg_decoder_fun('PL', kegg_sig_in_pl)

kegg_decoder_sur_df$Sample <- sapply(kegg_decoder_sur_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_sub_df$Sample <- sapply(kegg_decoder_sub_df$Sample, FUN = function(x){sub('_', '', x)})
kegg_decoder_pl_df$Sample <- sapply(kegg_decoder_pl_df$Sample, FUN = function(x){sub('_', '', x)})

# write.table(kegg_decoder_sur_df, file='tables/kegg/KEGG_decoder_sur_df_in.tsv', sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_sub_df, file='tables/kegg/KEGG_decoder_sub_df_in.tsv', sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)
# write.table(kegg_decoder_pl_df, file='tables/kegg/KEGG_decoder_pl_df_in.tsv', sep='\t',  col.names=FALSE, row.names=FALSE , quote=FALSE)

# in the Data folder using the CryoBiome-kegg-decoder env: KEGG-decoder -i KEGG_decoder_df_in.tsv -o KEGG_decoder_output
KEGG_decoder_sur <- read.table("tables/kegg/KEGG_decoder_sur_output.tsv", sep = "\t", header = TRUE)
KEGG_decoder_sub <- read.table("tables/kegg/KEGG_decoder_sub_output.tsv", sep = "\t", header = TRUE)
KEGG_decoder_pl <- read.table("tables/kegg/KEGG_decoder_pl_output.tsv", sep = "\t", header = TRUE)
KEGG_decoder <- rbind(KEGG_decoder_sur, KEGG_decoder_sub, KEGG_decoder_pl)
# write.table(KEGG_decoder, file='tables/kegg/KEGG_decode_output.tsv', sep='\t',  col.names=T, row.names=F , quote=FALSE)

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
        strip.text = element_blank(),
        strip.background = element_blank())
ggsave('./figs/metagenome/kegg_enrichment.pdf', p, width = 10, height = 4)

############################################################################################################
# 4. enrichment analysis with clusterProfiler package
library(clusterProfiler)
KO_pathway <- read.csv("e:/permafrost/data/metagenome/fun/keggKO_ko_metabolism.csv", header = T)
termgene <- KO_pathway[, c(2, 1)]
termname <- KO_pathway[, c(2, 3)]
enrichpathway_sur <- enricher(kegg_sig_in_sur, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_sub <- enricher(kegg_sig_in_sub, TERM2GENE = termgene, TERM2NAME = termname)
enrichpathway_pl <- enricher(kegg_sig_in_pl, TERM2GENE = termgene, TERM2NAME = termname)

write.csv(enrichpathway_sur, "enrichpathway_sur.csv")
write.csv(enrichpathway_sub, "enrichpathway_sub.csv")
write.csv(enrichpathway_pl, "enrichpathway_pl.csv")
