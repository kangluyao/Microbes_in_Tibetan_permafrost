# set work directory
wd_fun = "E:/permafrost/data/metagenome"
setwd(wd_fun)
library(tidyverse)
# read data
ko_tax_file <- "./fun/parse_dat.txt"
ko_description <- './fun/KO_description.txt'
ko_9_levels <- './fun/ko_9_levels.txt'
ko_tax_table <- read_delim(ko_tax_file, col_names = T)
ko_description_table <- read_delim(ko_description, col_names = T)
KO <- read_delim(ko_9_levels, col_names = T)

#
ko_enrich_sur <- KO %>%
  filter(KO %in% kegg_sig_in_sur)
unique(ko_enrich_sur$L2)

tax_enrich_sur <- ko_tax_table %>%
  filter(KO %in% kegg_sig_in_sur)


ko_enrich_sub <- KO %>%
  filter(KO %in% kegg_sig_in_sub)
unique(ko_enrich_sub$L2)
ko_enrich_pl <- KO %>%
  filter(KO %in% kegg_sig_in_pl)
unique(ko_enrich_pl$L2)
