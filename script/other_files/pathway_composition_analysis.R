#set work directory
wd_fun = 'E:/permafrost/data/metagenome' #work directory
save.dir = "E:/permafrost/result" #Output directory
setwd(wd_fun)

#loading packages
library(tidyverse)
#reading data
funcprof <- read.csv('fun/funcprof_abund.csv', header = T)


##ggplot2 堆叠面积图
p <- funcprof %>% 
  pivot_longer(-Pathway, names_to = 'sample_id', values_to = 'Abundance') %>%
  mutate(layer = sapply(stringr::str_split(sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL'), ordered = T)) %>%
  mutate(sample_id = rep(1:66, 11)) %>% 
  group_by(layer, sample_id) %>%
  mutate(prop = (Abundance / sum(Abundance)*100)) %>%
  ggplot(aes(x = sample_id, y = prop, fill = Pathway)) +
  geom_area() +
  labs(x = 'Sample', y = 'Proportion (%)') +
  scale_fill_brewer(palette = "Paired") +
  scale_y_continuous(expand = c(0,0)) +
  facet_grid(~ layer, scales = 'free_x', space = 'free_x') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 11, colour = "black"),
        axis.text.x = element_blank(),
        strip.text = element_text(size = 12, colour = "black"),
        axis.ticks.x = element_blank(),
        legend.title = element_text(size = 12), 
        legend.text = element_text(size = 11),
        legend.position = "left",
        legend.key.size = unit(0.7, "cm"),
        panel.spacing = unit(0, "lines"))

#save the picture
plot.name <- paste0(save.dir, "/figs/metagenome/pathway_composition.pdf")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8, height = 5.5, onefile = TRUE)
p
dev.off()
