# set work directory
setwd('e:/permafrost')
wd_16s <- file.path(getwd(),"data/16S")
# if (!dir.exists(wd_16s)) {
#   dir.create(wd_16s)
# }
wd_fun <- file.path(getwd(),"data/metagenome")
save.dir <- file.path(getwd(),"result")
# loading packages
library(phyloseq)
library(ape)
library(Biostrings)
library(tidyverse)
# read data
## metadata
metadata <- read.delim(file.path(wd_16s, "metadata_final.txt"), header = T, sep = "\t")
rownames(metadata) <- (metadata$sample_id)
## dna sequences
dna_seqs <- readDNAStringSet(file.path(wd_16s, "otus.fa"), format = "fasta", nrec = -1L, 
                             skip = 0L, seek.first.rec = FALSE,
                             use.names = TRUE)

## tree
tree <- read_tree(file.path(wd_16s, "otus.nwk"))
## otu table
otu <- read.delim(file.path(wd_16s, "otutab.txt"), header = T, row.names = 1, sep = "\t")
otu <- otu[, metadata$sample_id[metadata$sample_id %in% colnames(otu)]]
## rarefy out table
otu_rare <- read.delim(file.path(wd_16s, "otutab_rare.txt"), header = T, row.names = 1, sep = "\t")
otu_rare <- otu_rare[, metadata$sample_id[metadata$sample_id %in% colnames(otu_rare)]]
## tax
tax <- read.delim(file.path(wd_16s, "taxonomy.txt"), header = T, row.names = 1, sep = "\t")
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
# read functional table
wd_fun <- file.path(getwd(),"data/metagenome")
# if (!dir.exists(wd_fun)) {
#   dir.create(wd_fun)
# }
ko_tpm_table <- read.delim(file.path(wd_fun, "fun/eggnog.KEGG_ko.raw.tpm.txt"), 
                           header = T, row.names = 1, sep = "\t")
ko_count_table <- read.delim(file.path(wd_fun, "fun/eggnog.KEGG_ko.raw.counts.txt"), 
                             header = T, row.names = 1, sep = "\t")
ko_tpm_table <- ko_tpm_table[, metadata$sample_id[metadata$sample_id %in% colnames(ko_tpm_table)]]
ko_count_table <- ko_count_table[, metadata$sample_id[metadata$sample_id %in% colnames(ko_count_table)]]
