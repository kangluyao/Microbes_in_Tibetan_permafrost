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

set.seed(123)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub.txt", quote = F, row.names = F)

set.seed(1)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub1.txt", quote = F, row.names = F)

set.seed(2)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub2.txt", quote = F, row.names = F)

set.seed(3)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub3.txt", quote = F, row.names = F)

set.seed(4)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub4.txt", quote = F, row.names = F)

set.seed(5)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub5.txt", quote = F, row.names = F)

set.seed(6)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub6.txt", quote = F, row.names = F)

set.seed(7)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub7.txt", quote = F, row.names = F)

set.seed(8)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub8.txt", quote = F, row.names = F)

set.seed(9)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub9.txt", quote = F, row.names = F)

set.seed(10)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub10.txt", quote = F, row.names = F)

set.seed(11)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub11.txt", quote = F, row.names = F)

set.seed(12)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub12.txt", quote = F, row.names = F)

set.seed(13)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub13.txt", quote = F, row.names = F)

set.seed(14)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub14.txt", quote = F, row.names = F)

set.seed(15)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub15.txt", quote = F, row.names = F)

set.seed(16)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub16.txt", quote = F, row.names = F)


set.seed(17)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub17.txt", quote = F, row.names = F)

set.seed(18)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub18.txt", quote = F, row.names = F)

set.seed(19)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub19.txt", quote = F, row.names = F)

set.seed(20)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub20.txt", quote = F, row.names = F)

set.seed(21)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub21.txt", quote = F, row.names = F)

set.seed(22)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub22.txt", quote = F, row.names = F)

set.seed(23)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub23.txt", quote = F, row.names = F)


set.seed(24)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub24.txt", quote = F, row.names = F)


set.seed(25)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub25.txt", quote = F, row.names = F)


set.seed(26)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub26.txt", quote = F, row.names = F)


set.seed(27)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub27.txt", quote = F, row.names = F)


set.seed(28)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub28.txt", quote = F, row.names = F)


set.seed(29)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub29.txt", quote = F, row.names = F)


set.seed(30)
bin_infor %>%
  filter(Layers == "SUB") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_sub30.txt", quote = F, row.names = F)

for(i in 31:99) {
  set.seed(i)
  bin_infor %>%
    filter(Layers == "SUB") %>%
    pull(genome) %>%
    sample(., 34, replace = FALSE) %>%
    write.table(file = file.path(paste("E:/permafrost/revision/resample/bins_names_sub", i, ".txt", sep = "")), quote = F, row.names = F)
}



for(i in 31:99) {
  set.seed(i)
  bin_infor %>%
    filter(Layers == "PL") %>%
    pull(genome) %>%
    sample(., 34, replace = FALSE) %>%
    write.table(file = file.path(paste("E:/permafrost/revision/resample/bins_names_pl", i, ".txt", sep = "")), quote = F, row.names = F)
}




set.seed(123)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl.txt", quote = F, row.names = F)

set.seed(1)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl1.txt", quote = F, row.names = F)


set.seed(2)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl2.txt", quote = F, row.names = F)

set.seed(3)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl3.txt", quote = F, row.names = F)

set.seed(4)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl4.txt", quote = F, row.names = F)

set.seed(5)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl5.txt", quote = F, row.names = F)

set.seed(6)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl6.txt", quote = F, row.names = F)

set.seed(7)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl7.txt", quote = F, row.names = F)

set.seed(8)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl8.txt", quote = F, row.names = F)

set.seed(9)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl9.txt", quote = F, row.names = F)

set.seed(10)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl10.txt", quote = F, row.names = F)

set.seed(11)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl11.txt", quote = F, row.names = F)

set.seed(12)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl12.txt", quote = F, row.names = F)

set.seed(13)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl13.txt", quote = F, row.names = F)

set.seed(14)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl14.txt", quote = F, row.names = F)

set.seed(15)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl15.txt", quote = F, row.names = F)

set.seed(16)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl16.txt", quote = F, row.names = F)

set.seed(17)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl17.txt", quote = F, row.names = F)

set.seed(18)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl18.txt", quote = F, row.names = F)

set.seed(19)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl19.txt", quote = F, row.names = F)

set.seed(20)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl20.txt", quote = F, row.names = F)

set.seed(21)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl21.txt", quote = F, row.names = F)

set.seed(22)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl22.txt", quote = F, row.names = F)

set.seed(23)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl23.txt", quote = F, row.names = F)

set.seed(24)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl24.txt", quote = F, row.names = F)

set.seed(25)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl25.txt", quote = F, row.names = F)

set.seed(26)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl26.txt", quote = F, row.names = F)

set.seed(27)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl27.txt", quote = F, row.names = F)

set.seed(28)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl28.txt", quote = F, row.names = F)

set.seed(29)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl29.txt", quote = F, row.names = F)


set.seed(30)
bin_infor %>%
  filter(Layers == "PL") %>%
  pull(genome) %>%
  sample(., 34, replace = FALSE) %>%
  write.table(file = "E:/permafrost/revision/bins_names_pl30.txt", quote = F, row.names = F)
