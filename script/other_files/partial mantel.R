########################################################################################
## partial mantel test function for the relationship between the envs and community
partial.mantel.fun <- function(phylo) {
  env.table <- data.frame(sample_data(phylo))
  #env.table <- env.table[complete.cases(env.table), ]
  otu_table <- as.matrix(t(otu_table(phylo)))
  otu_table_hel <- decostand(otu_table, 'hellinger')
  otu_table_hel_dist <- vegdist(otu_table_hel, 'bray',upper=F)
  df <- NULL
  vars <- c("MAT", "AI", "Moisture",	"pH",	"Clay",	"Silt",	"Sand",	"TC", 
            "SOC", "DOC", "TN",	"NH4_N",	"NO3_N",	"DON",	
            "LCP1",	"LCP2",	"RCP")
  for (x in vars) {
    y.dist <- vegdist(scale(env.table[,x]), 'euclidean', na.rm = T)
    z.dist <- vegdist(scale(env.table[ , setdiff(vars, x)]), 'euclidean', na.rm = T)
    mode <- mantel.partial(otu_table_hel_dist, y.dist, z.dist, 
                           method = "pearson", permutations = 999, na.rm = T)
    r <- mode$statistic
    p <- mode$signif
    tmp <- data.frame(env = x, r = r, p.value = p)
    if(is.null(df))
      df <- tmp
    else
      df <- rbind(df ,tmp)
  }
  return(df)
}

phylo_sur <- subset_samples(phylo, Layer == 'SUR')
phylo_sur <- prune_taxa(taxa_sums(phylo_sur)>=1, phylo_sur)
phylo_sub <- subset_samples(phylo, Layer == 'SUB')
phylo_sub <- prune_taxa(taxa_sums(phylo_sub)>=1, phylo_sub)
phylo_pl <- subset_samples(phylo, Layer == 'PL')
phylo_pl <- prune_taxa(taxa_sums(phylo_pl)>=1, phylo_pl)

set.seed(123)
par.mant.sur <- partial.mantel.fun(phylo_sur)
set.seed(123)
par.mant.sub <- partial.mantel.fun(phylo_sub)
set.seed(123)
par.mant.pl <- partial.mantel.fun(phylo_pl)
## PLOT
## devtools::install_github('hannet91/ggcor')
library(ggcor)
par.man.tibble <- tibble(spec = c(rep('SUR', nrow(par.mant.sur)), rep('SUB', nrow(par.mant.sub)), rep('PL', nrow(par.mant.pl))), 
                         rbind(par.mant.sur, par.mant.sub, par.mant.pl))

vars <- c("MAT", "AI", "Moisture",	"pH",	"Clay",	"Silt",	"Sand",	"TC", 
          "SOC", "DOC", "TN",	"NH4_N",	"NO3_N",	"DON",	
          "LCP1",	"LCP2",	"RCP")
env.table <- sample_data(phylo)[ , vars]
mantel02 <- par.man.tibble %>% 
  mutate(r = cut(r, breaks = c(-Inf, 0.2, 0.5, Inf), 
                 labels = c("<0.20", "0.20-0.5", ">0.50"),
                 right = FALSE),
         p.value = cut(p.value, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
                       labels = c("<0.001", "0.001-0.01", "0.01-0.05", ">0.05"),
                       right = T))
p <- quickcor(env.table, type = "lower") + geom_square() + 
  add_link(mantel02, mapping = aes(colour = p.value, size = r),
           diag.label = TRUE) +
  scale_color_manual(values = c('#d95f02', '#1b9e77', '#3C5488FF', 'grey')) +
  scale_size_manual(values = c(0.5, 2, 3)) +
  geom_diag_label() + remove_axis("x")

ggsave(file.path(save.dir, './figs/env_effect/partail_matel_plot.pdf'), p, width = 9, height = 6, units = "in")
