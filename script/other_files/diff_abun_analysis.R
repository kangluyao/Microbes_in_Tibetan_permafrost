# 1.Kruskal-Wallis test
groups <- meta_dat$layer
kw.cal<-function(abund_table){
  #Apply normalisation (either use relative or log-relative transformation)
  #data<-log(abund_table/rowSums(abund_table))
  #data<-log((abund_table+1)/(rowSums(abund_table)+dim(abund_table)[2]))
  data<-as.data.frame(abund_table)
  
  #Reference: http://www.bigre.ulb.ac.be/courses/statistics_bioinformatics/practicals/microarrays_berry_2010/berry_feature_selection.html
  kruskal.wallis.table <- data.frame()
  for (i in 1:dim(data)[2]) {
    ks.test <- kruskal.test(data[,i], g=groups)
    # Store the result in the data frame
    kruskal.wallis.table <- rbind(kruskal.wallis.table,
                                  data.frame(id=names(data)[i],
                                             p.value=ks.test$p.value
                                  ))
    # Report number of values tested
    #cat(paste("Kruskal-Wallis test for ",names(data)[i]," ", i, "/", 
    #          dim(data)[2], "; p-value=", ks.test$p.value,"\n", sep=""))
  }
  
  
  kruskal.wallis.table$E.value <- kruskal.wallis.table$p.value * dim(kruskal.wallis.table)[1]
  
  kruskal.wallis.table$FWER <- pbinom(q=0, p=kruskal.wallis.table$p.value, 
                                      size=dim(kruskal.wallis.table)[1], lower.tail=FALSE)
  
  kruskal.wallis.table <- kruskal.wallis.table[order(kruskal.wallis.table$p.value,
                                                     decreasing=FALSE), ]
  kruskal.wallis.table$q.value.factor <- dim(kruskal.wallis.table)[1] / 1:dim(kruskal.wallis.table)[1]
  kruskal.wallis.table$q.value <- kruskal.wallis.table$p.value * kruskal.wallis.table$q.value.factor
  return(kruskal.wallis.table)
}

kw.plot<-function(table){
  plot(table$p.value,
       table$E.value,
       main='Multitesting corrections',
       xlab='Nominal p-value',
       ylab='Multitesting-corrected statistics',
       log='xy',
       col='blue',
       panel.first=grid(col='#BBBBBB',lty='solid'))
  lines(table$p.value,
        table$FWER,
        pch=20,col='darkgreen', type='p'
  )
  lines(table$p.value,
        table$q.value,
        pch='+',col='darkred', type='p'
  )
  kruskal.wallis.alpha=0.01
  abline(h=kruskal.wallis.alpha, col='red', lwd=2)
  legend('topleft', legend=c('E-value', 'p-value', 'q-value'), col=c('blue', 'darkgreen','darkred'), lwd=2,bg='white',bty='o')
}

diff.cal<-function(table){
  kruskal.wallis.alpha=0.01
  last.significant.element <- max(which(table$q.value <= kruskal.wallis.alpha))
  selected <- 1:last.significant.element
  diff.cat.factor <- table$id[selected]
  diff.cat <- as.vector(diff.cat.factor)#获得差异显著的gene的name
  return(diff.cat)
}

kw_C_N <- kw.cal(C_N_path2)
kw.plot(kw_C_N)
kw_C_N.diff <- diff.cal(kw_C_N)
kw_C_N.diff <- kw_C_N.diff[order(match(kw_C_N.diff, sel_ko$Enzyme_protein_encoded))]
length(kw_C_N.diff)

library(reshape2)
kw_dat <- data.frame(layer = meta_dat$layer, C_N_path2[ , c(kw_C_N.diff)])
mydf <- melt(kw_dat, measure.vars=names(kw_dat)[-1])
mydf$layer <- factor(mydf$layer, levels = c('SUR', 'SUB', 'PL'))
mydf$variable <- factor(mydf$variable, levels = kw_C_N.diff)

hacky_df <- sel_ko %>% filter(Enzyme_protein_encoded %in% kw_C_N.diff) %>%
  select(Subcategory, Enzyme_protein_encoded)
library(gtable)
library(grid)

d <- data.frame(fruit = rep(c("apple", "orange", "plum", "banana", "pear", "grape")), 
                farm = rep(c(0,1,3,6,9,12), each=6), 
                weight = rnorm(36, 10000, 2500), 
                size=rep(c("small", "large")))

dummy <- ggplot(data = d, aes(x = farm, y = weight))+ facet_wrap(~fruit) + 
  geom_rect(aes(fill=size), xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf) +
  theme_minimal()



library(ggplot2)
p1 <- ggplot(mydf,aes(x = layer, y = value)) +
  geom_boxplot(aes(fill = layer)) +
  facet_wrap(~variable, nrow = 4, scales = "free")

# plot code
p1 <-
  ggplot(mydf) +    # don't specify x and y here.  Otherwise geom_rect will complain.
  geom_rect(
    data = hacky_df,
    aes(xmin = -Inf, xmax = Inf,
        ymin = -Inf, ymax = Inf,     # totally defined by trial-and-error
        fill = Subcategory, alpha = 0.4)) +
  geom_boxplot(aes(x = layer, y = value, fill = layer)) +     
  # coord_cartesian(clip = "off", ylim = c(10, 35)) +
  facet_wrap(~variable, nrow = 4, scales = "free")
  # scale_fill_manual(values = c("area" = "green", "bat" = "red", "vege" = "blue", "indus" = "black")) +

# Extract only the legend from "p1" plot
g_legend <- function(p1){ 
  tmp <- ggplot_gtable(ggplot_build(p1)) 
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box") 
  legend <- tmp$grobs[[leg]] 
  return(legend)} 
# Assign the legend to a separate object
facet.legend <- g_legend(p1)

plot_new


# dataset for plotting
df <- mtcars %>% gather(-mpg, key = "var", value = "value")

# dataset for facet label colors
hacky_df <- data.frame(
  var = c("am", "carb", "cyl", "disp", "drat", "gear", "hp", "qsec", "vs", "wt"),
  var_color = c("area", "indus", "indus", "bat", "bat", "bat", "area", "indus", "vege", "vege")
)






































# kw without e and q value
library(dplyr)
library(ggplot2)
kw_dat <- data.frame(layer = meta_dat$layer, C_N_path2)
mydf <- kw_dat %>% tidyr::pivot_longer(cols = -c(layer), 
                                       names_to = "pathway", 
                                       values_to = 'abundance') 
pv <- mydf %>% group_by(pathway) %>% 
  summarize(p.value = kruskal.test(abundance ~ layer)$p.value)
# library(rcompanion)
# mod1 <- kruskal.test(narB ~ layer, data = kw_dat)
# mod1$p.value
# mod2 <- rcompanion::epsilonSquared(x = kw_dat$narB, g = kw_dat$layer)
# 
#   
#   melt(kw_dat, measure.vars=names(kw_dat)[-1])
# pv <- mydf %>% group_by(pathway) %>%
#   summarize(p.value = kruskal.test(abundance ~ layer)$p.value)

pv <- pv[pv$p.value < 0.05, ]
kw_C_N.diff <- pv$pathway[order(match(pv$pathway, sel_ko$Enzyme_protein_encoded))]
mydf$layer <- factor(mydf$layer, levels = c('SUR', 'SUB', 'PL'))
mydf$pathway <- factor(mydf$pathway, levels = kw_C_N.diff)

ggplot(mydf,aes(x = layer, y = abundance)) +
  geom_boxplot(aes(fill = layer)) +
  facet_wrap(~pathway, nrow = 6, scales = "free")


# ANCOMBC
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")
source('E:/R/permafrost/ancom.R')



# zero_cut = 0.99 means that only genera present in at least 21 samples are taken into account
feature_table = merged_table; sample_var = "Sample"; group_var = 'Cryosphere'
out_cut = 0.05; zero_cut = 0.995; lib_cut = 1000; neg_lb = FALSE
prepro = feature_table_pre_process(feature_table, merged_metadata, sample_var, group_var, out_cut, zero_cut, lib_cut, neg_lb)
feature_table = prepro$feature_table # Preprocessed feature table
meta_data = prepro$meta_data # Preprocessed metadata
struc_zero = prepro$structure_zeros # Structural zero info

# Step 2: ANCOM
main_var = 'Cryosphere'; p_adj_method = "BH"; alpha = 0.05
adj_formula = NULL; rand_formula = "~ 1 | Dataset"
control = lmeControl(maxIter = 100, msMaxIter = 100, opt = "optim")
res = ANCOM(feature_table, meta_data, struc_zero, main_var, p_adj_method, alpha, adj_formula, rand_formula, control=control)

# Step 3: Volcano Plot
# Number of taxa except structural zeros
n_taxa = ifelse(is.null(struc_zero), nrow(feature_table), sum(apply(struc_zero, 1, sum) == 0))
# Cutoff values for declaring differentially abundant taxa
cut_off = c(0.9 * (n_taxa -1), 0.8 * (n_taxa -1), 0.7 * (n_taxa -1), 0.6 * (n_taxa -1))
names(cut_off) = c("detected_0.9", "detected_0.8", "detected_0.7", "detected_0.6")

# Annotation data W=1277
colnames(res$fig$data)[colnames(res$fig$data) == 'x'] = 'CLR_mean_diff'
colnames(res$fig$data)[colnames(res$fig$data) == 'y'] = 'W'
dat_ann = data.frame(x = min(res$fig$data$CLR_mean_diff), y = cut_off["detected_0.7"], label = "W[0.7]")

res$fig$data$DA = 'Others'
res$fig$data$DA[(res$fig$data$W>dat_ann$y)&(res$fig$data$CLR_mean_diff>0)] = 'Overrepresented'

ggplot(res$fig$data) + geom_point(aes(x=CLR_mean_diff,y=W,color=DA)) + scale_color_manual(values=c('dimgrey','#1A7A7F')) + theme_linedraw() +
  geom_hline(yintercept = dat_ann$y, linetype = "dashed") + ylab('W statistic') + xlab('CLR mean difference') + theme(legend.position = 'none') +
  geom_text(data = dat_ann, aes(x = x, y = y, label = label), 
            size = 4, vjust = -0.5, hjust = 0, color = "orange", parse = TRUE) +
  geom_vline(xintercept = 0,size=0.2)
ggsave('2_Genus_analysis/2_1_Diff_abund/2_1_Amplicon_ancom_res.pdf',width = 4,height = 4)

over_genera = as.character(res$fig$data$taxa_id[(res$fig$data$CLR_mean_diff>0) & (res$fig$data$W > dat_ann$y)])
over_data = res$fig$data[(res$fig$data$CLR_mean_diff>0)&(res$fig$data$W>dat_ann$y),]

write.csv(res$fig$data, file='Data/Amplicon_ancom_res.csv')

# Most overrepresented genera barplot w=1251 for 0.7
over_data = read.csv('Data/Amplicon_ancom_res.csv')
over_data$Genus = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; g__')[[1]][2], FUN.VALUE = character(1))
over_data$Family = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][5], FUN.VALUE = character(1))
over_data$Order = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][4], FUN.VALUE = character(1))
over_data$Class = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][3], FUN.VALUE = character(1))
over_data$Phylum = vapply(as.character(over_data$taxa_id), function(x) strsplit(x, split='; ')[[1]][2], FUN.VALUE = character(1))
over_data$Phylum = gsub('p__','',over_data$Phylum)

table(over_data$Phylum[over_data$DA == 'Overrepresented'])

# phylum: c('Acidobacteriota','Actinobacteriota','Bacteroidota','Chloroflexi','Cyanobacteria','Firmicutes','Others,'Patescibacteria','Planctomycetota','Proteobacteria','Verrucomicrobiota')
# colors: c('#B09C85FF','#F39B7FFF','#DC0000FF','#91D1C2FF','#00A087FF','#7E6148FF','grey','dimgrey','#4DBBD5FF','#3C5488FF','#8491B4FF')
most_over = over_data[(over_data$CLR_mean_diff > 0.35) & (over_data$W >= 1269),]
ggplot(most_over) + 
  geom_bar(aes(x=CLR_mean_diff,y=reorder(as.factor(Genus), CLR_mean_diff),fill=Phylum), stat='identity') + ylab('') + xlab('CLR mean difference') +
  theme_linedraw() + theme(panel.grid = element_line(colour = 'darkgrey'), axis.text.y = element_text(face = "italic")) + 
  scale_fill_manual(values=c('#F39B7FFF','#DC0000FF','#00A087FF','grey','dimgrey','#3C5488FF','#8491B4FF'))
ggsave('2_Genus_analysis/2_1_Diff_abund/2_1_Ancom_most_over.pdf', width=6.5, height = 6)


# Create Taxonomic tree at the family level with the edge size as the number of overrepresented genera
library(metacoder)
tree_data = parse_tax_data(data.frame(Taxonomy= over_data$taxa_id[(over_data$CLR_mean_diff > 0) & (over_data$W >= 1269)]),
                           class_cols = "Taxonomy", # the column that contains taxonomic information
                           class_sep = "; ", # The character used to separate taxa in the classification
                           class_regex = "^(.*)__(.*)$", # Regex identifying where the data for each taxon is
                           class_key = c(tax_rank = "info", # A key describing each regex capture group
                                         tax_name = "taxon_name"))
tree_data %>% filter_taxa(tree_data$data$class_data$taxon_id[tree_data$data$class_data$tax_rank == 'c'], supertaxa = TRUE) %>%
  filter_taxa(n_obs > 1, supertaxa = TRUE) %>%
  heat_tree(node_label = taxon_names,
            node_color_range = c('grey','#B09C85FF','#00A087FF'),
            node_size=n_obs,
            node_color=n_obs,
            node_color_trans='log10',
            node_size_trans='log10',
            node_color_axis_label = "Cryo. genera #",
            margin_size = c(0.1,0.1,0.1,0.1),
            initial_layout = "re",
            overlap_avoidance = 10,
            output_file = "2_Genus_analysis/2_1_Diff_abund/2_1_Ancom_over_genera.pdf")

