#environmental factors test among layers
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
env_vars <- c("Clay", "Silt", "Sand", "pH", "Moisture", "TC", "SOC", "DOC", "LCP1", "LCP2", "RCP",
              "TN", "NH4_N", "NO3_N", "DON")

meta_dat %>% as.tibble() %>% select(env_vars) %>%
  group_by(Layer) %>%
  summarise_all(list(mean = mean, se = plotrix::std.error)) %>%
  write.csv(., file = file.path(save.dir, "tables/env_comparison/env_compare.csv"))

p <- meta_dat %>% as.tibble() %>% select("Layer", env_vars) %>%
  pivot_longer(cols = -c(Layer), names_to = "env", values_to = "value") %>% 
  mutate(Layer = factor(Layer, levels = c("SUR", "SUB", "PL"))) %>% 
  mutate(env = factor(env, levels = env_vars)) %>%
  ggbarplot(x = "Layer", y = "value", fill = "Layer", 
          facet.by = "env", palette = "jco",
          add = "mean_se", position = position_dodge(0.8)) + 
  facet_wrap(vars(env), scales = "free", ncol = 5) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  stat_compare_means(comparisons=my_comparisons, paired = TRUE, label = "p.signif",
                     test = "wilcox.test", size = 4, p.adjust.method = "BH",
                     step.increase = 0.22) +
  theme_linedraw() +
  theme(strip.text = element_text(color = 'black'),
        strip.background = element_rect(fill = "#aaa9a9"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        legend.position = "none")

plot.name <- paste0(save.dir, "/figs/env_effect/env_paired_wilcox_test_abc.pdf")
print(plot.name)
cairo_pdf(filename = plot.name, width = 8.6, height = 7, onefile = TRUE)
p
dev.off()

library(tidyverse)
library(rstatix)
meta_dat %>% as.tibble() %>% select(-c(1, 2, 4:6)) %>%
  pivot_longer(cols = -c(Layer), names_to = "env", values_to = "value") %>% 
  group_by(env) %>% 
  pairwise_wilcox_test(value ~ Layer, p.adjust.method = "BH")
  summarize(out = pairwise.wilcox.test(value, as.factor(Layer), p.adjust.method="BH")$p.value)
