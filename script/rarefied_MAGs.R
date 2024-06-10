source("script/read_data.R")
# tree.file <- file.path(wd_fun, 'binning_70/result/gtdb_tree/bacteria/tax.unrooted.tree')
abundance_tab.file <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_all/bin_abundance_table.tab")
tax.file <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_all/tax.txt")
# bin genome information
bin_infor_file <- file.path(wd_fun, "metabolic/METABOLIC_70_layer/METABOLIC_70_all/Widb.csv")
metabolic_output_file <- file.path(wd_fun, "/metabolic/METABOLIC_70_layer/METABOLIC_70_all/metabolic.csv")
# tree <- read.tree(tree.file)
abundance_tab <- read_delim(abundance_tab.file, col_names = T)
colnames(abundance_tab)[1] <- 'ID'
tax_tab <- read_delim(tax.file, col_names = T)
bin_infor <- read.csv(bin_infor_file) %>%
  select(c("genome", "completeness", "contamination", "size")) %>%
  mutate(Layers = sapply(stringr::str_split(genome, "_",  n = 4), `[`, 2)) %>%
  mutate(Layers = factor(Layers, levels = c('SUR', 'SUB', 'PL')))
metabolic_tab <- read.csv(metabolic_output_file) %>%
  select(c('Category', 'Function', grep('Function.presence', colnames(.))))

for(i in 1:99) {
  set.seed(i)
  bin_infor %>%
    filter(Layers == "SUB") %>%
    pull(genome) %>%
    sample(., 34, replace = FALSE) %>%
    write.table(file = file.path(paste("E:/permafrost/revision/resample/bins_names_sub", i, ".txt", sep = "")), quote = F, row.names = F)
}


for(i in 1:99) {
  set.seed(i)
  bin_infor %>%
    filter(Layers == "PL") %>%
    pull(genome) %>%
    sample(., 34, replace = FALSE) %>%
    write.table(file = file.path(paste("E:/permafrost/revision/resample/bins_names_pl", i, ".txt", sep = "")), quote = F, row.names = F)
}
