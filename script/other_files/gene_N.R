#set work directory
setwd('e:/permafrost')
wd_fun <- file.path(getwd(),"data/metagenome")
if (!dir.exists(wd_fun)) {
  dir.create(wd_fun)
}
save.dir <- file.path(getwd(),"result")
if (!dir.exists(save.dir)) {
  dir.create(save.dir)
}

library(tidyverse)
N_genes <- read_csv(file.path(wd_fun, 'flux/N_cycle.csv'), col_names = T, skip = 1)

#You can use sel_env to specify the variables you want to use and sel_env_label to specify the labes for the pannel
names <- c("narB", "nasA",	"nasB",	"nirA",	"napA",	"napB",	"narG",	"narH",	"narI",	"nirK",
           "nirS",	"norB",	"norC",	"nosZ",	"nirB",	"nirD",	"nrfA",	"nrfH",	"nrtA",	"nrtB",
           "nrtC",	"nrtD",	"nifD",	"nifH",	"nifK",	"amoA",	"amoB",	"amoC",	"hao",	"nxrA",	"nxrB")

###N过程组间差异
library(ggplot2)
library(ggpubr)
my_comparisons <- list( c("SUR", "SUB"), c("SUB", "PL"), c("SUR", "PL") )
p <- N_genes[ ,c(2, 4:8)] %>% pivot_longer(cols = -layer, 
                                           names_to = 'N_process', values_to = 'value') %>%
  mutate(layer = factor(layer, levels = c('SUR', 'SUB', 'PL'))) %>%
  ggplot(aes(x = layer, y = value)) + 
  geom_boxplot(width = 0.5, aes(fill = layer))+
  stat_compare_means(comparisons = my_comparisons, label = "p.signif") +
  #scale_fill_manual(values= cols)+
  labs(x = 'Layers', y = 'Value', fill='Layers') +
  facet_wrap(.~N_process, scales = 'free') +
  theme_bw() +
  theme(axis.title = element_text(size = 12, colour = "black"),
        axis.text = element_text(size = 10, colour = "black"),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        panel.grid = element_blank())
ggsave(file.path(save.dir, "./figs/metagenome/N_process.pdf"), p, width = 89, height = 89, units = "mm")


library(vegan)
x<-data.frame(N_genes[, 4:8])
y<-data.frame(N_genes[, -c(1:8)])
colnames(y) = names
groups<-c(rep('All',22))

#You can use kendall, spearman, or pearson below:
method<-"pearson"


#Now calculate the correlation between individual Taxa and the environmental data
df<-NULL
for(i in colnames(x)){
  for(j in colnames(y)){
    for(k in unique(groups)){
      a<-x[groups==k,i,drop=F]
      b<-y[groups==k,j,drop=F]
      tmp<-c(i,j,cor(a[complete.cases(b),],b[complete.cases(b),],
                     use="everything",method=method),
             cor.test(a[complete.cases(b),],b[complete.cases(b),],method=method)$p.value,k)
      if(is.null(df)){
        df<-tmp  
      }
      else{
        df<-rbind(df,tmp)
      }    
    }
  }
}

df<-data.frame(row.names=NULL,df)
colnames(df)<-c("Index","Env","Correlation","Pvalue","Type")
df$Pvalue<-as.numeric(as.character(df$Pvalue))
df$AdjPvalue<-rep(0,dim(df)[1])
df$Correlation<-as.numeric(as.character(df$Correlation))
#You can adjust the p-values for multiple comparison using Benjamini & Hochberg (1995):
# 1 -> donot adjust
# 2 -> adjust Env + Type (column on the correlation plot)
# 3 -> adjust Taxa + Type (row on the correlation plot for each type)
# 4 -> adjust Taxa (row on the correlation plot)
# 5 -> adjust Env (panel on the correlation plot)
adjustment_label<-c("NoAdj","AdjEnvAndType","AdjTaxaAndType","AdjTaxa","AdjEnv")
adjustment<-5

if(adjustment==1){
  df$AdjPvalue<-df$Pvalue
} else if (adjustment==2){
  for(i in unique(df$Env)){
    for(j in unique(df$Type)){
      sel<-df$Env==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==3){
  for(i in unique(df$Taxa)){
    for(j in unique(df$Type)){
      sel<-df$Taxa==i & df$Type==j
      df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
    }
  }
} else if (adjustment==4){
  for(i in unique(df$Taxa)){
    sel<-df$Taxa==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
} else if (adjustment==5){
  for(i in unique(df$Env)){
    sel<-df$Env==i
    df$AdjPvalue[sel]<-p.adjust(df$Pvalue[sel],method="BH")
  }
}

#Now we generate the labels for significant values
df$Significance<-cut(df$Pvalue, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), label=c("***", "**", "*", ""))

#We ignore NAs
df<-df[complete.cases(df),]

#We want to reorganize the Env data based on they appear
df$Env<-factor(df$Env, levels= names)
# Renaming pathways for plotting
df$Index[df$Index == "Net.N.mineralization"] <- "Net N mineralization"
df$Index[df$Index == "Net.nitrification"] <- "Net nitrification"
df$Index[df$Index == "Gross.N.mineralization"] <- "Gross N mineralization"
df$Index[df$Index == "Microbial.immobilization"] <- "Microbial immobilization"
df$Index[df$Index == "Gross.nitrification"] <- "Gross nitrification"

df$Index<-factor(df$Index,levels = c("Net N mineralization",	"Gross N mineralization", 
                                     "Net nitrification",	"Gross nitrification", 
                                     "Microbial immobilization"))

df$r_squar <- (df$Correlation)^2
df

#heatmap
p_N_env <- ggplot(aes(x = Env, y = Index, fill = Correlation), data = df) + 
  geom_tile() + 
  scale_fill_gradient2(low = "#009E73", mid = "white", high = "#E69F00") + 
  geom_text(aes(label = paste(round(Correlation, 2), Significance)), 
            color = "black", size = 3) + 
  labs(y = "Carbon properties", x = "pathway", fill = "Pearson R^2") + 
  theme(axis.title = element_blank(), 
        axis.text = element_text(size = 12, colour = "black"), 
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
print(p_N_env)

ggplot(melted,aes(x=value,y=Cumulative_11))+
  geom_point(shape=19, size=2,colour='tomato',alpha=0.8)+
  geom_smooth(method="lm", size=1.5, se=T,colour='black') +
  ylab('Cumulative')+
  scale_y_continuous(limits = c(-1,6.5))+
  facet_wrap(.~ variable , scales="free",ncol = 4)+
  mytheme

