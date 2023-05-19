########################################################################################
#heatmap for alpha diversity and environmental factors
library(ggplot2)
library(tidyr)
library(dplyr)
library(tibble)
cor_div_env_tab <- data.frame(metadata[-c(1, 2, 4:6)], Evenness = alpha_div$Evenness)
#You can use sel_env to specify the variables you want to use and sel_env_label to specify the labes for the pannel
sel_env <- c("MAT", "AI", "Moisture",	"pH",	"Clay",	"Silt",	"Sand",	"TC", 
           "SOC", "DOC", "TN",	"NH4_N",	"NO3_N",	"DON",	
           "LCP1",	"LCP2",	"RCP")

#
groups <- cor_div_env_tab$Layer
x <- data.frame(cor_div_env_tab[, sel_env])
y <- data.frame(Evenness = cor_div_env_tab$Evenness)
#You can use kendall, spearman, or pearson below:
method<-"pearson"
#Now calculate the correlation between individual Taxa and the environmental data
df <- NULL
for (i in colnames(x)) {
  for (j in colnames(y)) {
    for (k in unique(groups)) {
      a <- x[groups == k, i, drop = F]
      b <- y[groups == k, j, drop = F]
      tmp <- c(i, j, cor(a[complete.cases(b), ], b[complete.cases(b), ], use = "everything", method = method), 
               cor.test(a[complete.cases(b), ], b[complete.cases(b), ], method = method)$p.value, k)
      if (is.null(df)) {
        df <- tmp
      } else {
        df <- rbind(df, tmp)
      }
    }
  }
}

df <- data.frame(row.names = NULL, df)
colnames(df) <- c("Env", "div", "Rho", "Pvalue", "Layer")
df$Pvalue <- as.numeric(as.character(df$Pvalue))
df$AdjPvalue <- rep(0, dim(df)[1])
df$Rho <- as.numeric(as.character(df$Rho))

# Now we generate the labels for signifant values
df$Significance <- cut(df$Pvalue, breaks = c(-Inf, 0.001, 0.01, 0.05, Inf), label = c("***", "**", "*", ""))

# We ignore NAs
df <- df[complete.cases(df), ]
df <- df[, c("Env", "Layer", "Rho", "Significance")]
# We want to reorganize the Env data based on they appear
df$Env <- factor(df$Env, levels = sel_env)
df$Layer <- factor(df$Layer, levels = rev(c("SUR", "SUB", "PL")))

# plot
row_num = length(levels(df$Layer))
p <- ggplot(aes(x = Env, y = as.numeric(Layer), fill = Rho), data = df) + 
  xlim(c("", sel_env)) +
  ylim(c(-row_num/1.5, row_num + 1)) +
  geom_tile() + ylab("") + 
  scale_fill_gradient2(low = "#2C7BB6", mid = "white", high = "#D7191C", limits = c(-1, 1)) + 
  annotate(x = "", y = 1:row_num, label = levels(df$Layer), size = 2.5, geom = "text") +
  geom_text(aes(label = Significance), color = "black", size = 4)

p <- p + coord_polar(start = -0.15, theta = "x") + 
  theme_bw() + 
  theme(axis.text = element_text(color = "black"),
        legend.position = c(0.5, 0.5), 
        legend.key.size = unit(0.5, "cm"))

ggsave(file.path(save.dir, './figs/env_effect/env_alpha_plot.pdf'), p, width = 9, height = 6, units = "in")
