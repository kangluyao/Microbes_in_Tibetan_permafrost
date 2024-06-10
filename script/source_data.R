# Load necessary library
library(dplyr)

# Function to perform pairwise Wilcoxon tests and adjust p-values
pairwise_wilcox_test <- function(data, response_vars, group_var) {
  results <- list()
  
  for (response in response_vars) {
    test_result <- pairwise.wilcox.test(data[[response]], data[[group_var]], p.adjust.method = "none", paired = TRUE)
    results[[response]] <- test_result$p.value
  }
  
  # Flatten the list of p-values and adjust using BH method
  all_p_values <- unlist(results)
  adjusted_p_values <- p.adjust(all_p_values, method = "BH")
  
  # Reshape adjusted p-values back to original structure
  index <- 1
  adjusted_results <- list()
  
  for (response in response_vars) {
    p_matrix <- results[[response]]
    adjusted_p_matrix <- matrix(adjusted_p_values[index:(index + length(p_matrix) - 1)], nrow = nrow(p_matrix))
    rownames(adjusted_p_matrix) <- rownames(p_matrix)
    colnames(adjusted_p_matrix) <- colnames(p_matrix)
    adjusted_results[[response]] <- adjusted_p_matrix
    index <- index + length(p_matrix)
  }
  
  return(adjusted_results)
}


wilcox_test <- function(data, response_vars, group_var) {
  results <- list()
  
  for (response in response_vars) {
    test_result <- pairwise.wilcox.test(data[[response]], data[[group_var]], p.adjust.method = "none", paired = F)
    results[[response]] <- test_result$p.value
  }
  
  # Flatten the list of p-values and adjust using BH method
  all_p_values <- unlist(results)
  adjusted_p_values <- p.adjust(all_p_values, method = "BH")
  
  # Reshape adjusted p-values back to original structure
  index <- 1
  adjusted_results <- list()
  
  for (response in response_vars) {
    p_matrix <- results[[response]]
    adjusted_p_matrix <- matrix(adjusted_p_values[index:(index + length(p_matrix) - 1)], nrow = nrow(p_matrix))
    rownames(adjusted_p_matrix) <- rownames(p_matrix)
    colnames(adjusted_p_matrix) <- colnames(p_matrix)
    adjusted_results[[response]] <- adjusted_p_matrix
    index <- index + length(p_matrix)
  }
  
  return(adjusted_results)
}

# Alpha diversity
# Function to calculate standard error
se <- function(x) {
  sd(x) / sqrt(length(x))
}
# Calculate summary statistics using dplyr
summary_stats_alpha <- alpha_div %>%
  group_by(Layers) %>%
  summarise(across(c("Shannon", "Faith"), list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))


# Print the summary statistics
print(summary_stats_alpha)

# Specify response variables and group variable
response_vars <- c("Shannon", "Faith")
group_var <- "Layers"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_alpha <- pairwise_wilcox_test(alpha_div, response_vars, group_var)

# Print the adjusted results
adjusted_results_alpha


# Beta diversity
# prepare the beta data
beta_data <- data.frame(sapply(unique(metadata$Layer), function(x) usedist::dist_subset(tax_dist, grep(x, metadata$sample_id, value = T))) %>%
                          data.frame() %>% gather("Layers", "distance"),
                        sapply(unique(metadata$Layer), function(x) usedist::dist_subset(beta.mntd.weighted, grep(x, metadata$sample_id, value = T))) %>% 
                          data.frame() %>% gather("Layers", "distance"),
                        sapply(unique(metadata$Layer), function(x) usedist::dist_subset(fun_dist, grep(x, metadata$sample_id, value = T))) %>% 
                          data.frame() %>% gather("Layers", "distance"))[, c(1, 2, 4, 6)] %>%
  mutate(Layers = factor(Layers, levels = c("SUR", "SUB", "PL")))
colnames(beta_data) <- c("Layers", "Tax_distance", "Phylogenetic distance", "Functional distance")

# Calculate summary statistics using dplyr
summary_stats_beta <- beta_data %>%
  group_by(Layers) %>%
  summarise(across(c("Tax_distance", "Phylogenetic distance", "Functional distance"), list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# Print the summary statistics
print(summary_stats_beta)


# Specify response variables and group variable
response_vars <- c("Tax_distance", "Phylogenetic distance", "Functional distance")
group_var <- "Layers"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_beta <- pairwise_wilcox_test(beta_data, response_vars, group_var)

# Print the adjusted results
adjusted_results_beta




# taxa composition
compo_table <- data.frame(subtaxa_tab, ra.tab) %>% group_by(Phylum) %>% 
  summarise(across(everything(), sum)) %>% 
  mutate(MRA = rowMeans(.[, colnames(ra.tab)])) %>%
  arrange(desc(MRA)) %>% dplyr::top_n(10, MRA) %>%
  dplyr::select(-MRA) %>%
  column_to_rownames("Phylum") %>% t() %>% data.frame() %>%
  mutate(Layers = sapply(stringr::str_split(metadata$sample_id, "_",  n = 2), `[`, 1)) %>%
  mutate(Layers = factor(Layers, levels = c("SUR", "SUB", "PL")))
  
# Calculate summary statistics using dplyr
top_phlya <- c("Proteobacteria", "Actinobacteriota", "Chloroflexi", "Acidobacteriota",
               "Bacteroidota", "Firmicutes", "Desulfobacterota", "Planctomycetota",	
               "Verrucomicrobiota", "Gemmatimonadota")
summary_stats_compo <- compo_table %>%
  group_by(Layers) %>%
  summarise(across(top_phlya, list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# Print the summary statistics
print(summary_stats_compo)
# write.csv(summary_stats_compo, file = "C:/Users/dell/Desktop/summary_stats_compo.csv")

# Specify response variables and group variable
response_vars <- top_phlya
group_var <- "Layers"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_compo <- pairwise_wilcox_test(compo_table, response_vars, group_var)

# Print the adjusted results
adjusted_results_compo

# null model
process <- c('Heterogeneous.Selection', 'Homogeneous.Selection', 
             'Dispersal.Limitation', 'Homogenizing.Dispersal', 'Drift.and.Others')

summary_stats_null_process <- null_df %>%
  mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))) %>%
  group_by(layer) %>%
  summarise(across(process, list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# Print the summary statistics
print(summary_stats_null_process)
write.csv(summary_stats_null_process, file = "C:/Users/dell/Desktop/summary_stats_null_process.csv")

# Specify response variables and group variable
response_vars <- process
group_var <- "layer"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_null_process <- wilcox_test(null_df %>%
                                               mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))), 
                                             response_vars, group_var)

# Print the adjusted results
adjusted_results_null_process


# each layer
summary_stats_null_eachlayer <- null_df %>%
  dplyr::select(c('layer', 'Heterogeneous.Selection', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Homogenizing.Dispersal', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  mutate(layer = factor(layer, levels = c("SUR", "SUB", "PL"))) %>%
  mutate(process = factor(process, levels = c('Dispersal.Limitation', 'Drift.and.Others',
                                              'Homogeneous.Selection', 'Homogenizing.Dispersal', 'Heterogeneous.Selection'))) %>%
  group_by(layer, process) %>%
  summarise(across("value", list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# Print the summary statistics
print(summary_stats_null_eachlayer)
write.csv(summary_stats_null_eachlayer, file = "C:/Users/dell/Desktop/summary_stats_null_eachlayer.csv")


# prepare the data
null_df_sur <- null_df %>%
  dplyr::select(c('layer', 'Heterogeneous.Selection', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Homogenizing.Dispersal', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  filter(layer == "SUR") %>%
  mutate(process = factor(process, levels = c('Dispersal.Limitation', 'Drift.and.Others',
                                              'Homogeneous.Selection', 'Homogenizing.Dispersal', 'Heterogeneous.Selection')))

null_df_sub <- null_df %>%
  dplyr::select(c('layer', 'Heterogeneous.Selection', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Homogenizing.Dispersal', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  filter(layer == "SUB") %>%
  mutate(process = factor(process, levels = c('Dispersal.Limitation', 'Drift.and.Others',
                                              'Homogeneous.Selection', 'Homogenizing.Dispersal', 'Heterogeneous.Selection')))

null_df_pl <- null_df %>%
  dplyr::select(c('layer', 'Heterogeneous.Selection', 'Homogeneous.Selection', 'Dispersal.Limitation', 'Homogenizing.Dispersal', 'Drift.and.Others')) %>%
  pivot_longer(cols = -c(layer), names_to = "process", values_to = "value") %>%
  filter(layer == "PL") %>%
  mutate(process = factor(process, levels = c('Dispersal.Limitation', 'Drift.and.Others',
                                              'Homogeneous.Selection', 'Homogenizing.Dispersal', 'Heterogeneous.Selection')))
  
# Specify response variables and group variable
response_vars <- "value"
group_var <- "process"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_null_sur <- wilcox_test(null_df_sur, response_vars, group_var)
adjusted_results_null_sub <- wilcox_test(null_df_sub, response_vars, group_var)
adjusted_results_null_pl <- wilcox_test(null_df_pl, response_vars, group_var)

# Print the adjusted results
adjusted_results_null_sur
adjusted_results_null_sub
adjusted_results_null_pl



# bin size
bin_size_data <- bin_infor %>% data.frame() %>%
  dplyr::select(Layers, size) %>%
  mutate(Layers = factor(Layers, levels = c("SUR", "SUB", "PL")))

summary_stats_bin_size <- bin_size_data %>%
  group_by(Layers) %>%
  summarise(across("size", list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))

# Print the summary statistics
print(summary_stats_bin_size)
# write.csv(summary_stats_bin_size, file = "C:/Users/dell/Desktop/summary_stats_bin_size.csv")

# Specify response variables and group variable
response_vars <- "size"
group_var <- "Layers"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_bin_size <- wilcox_test(bin_size_data, response_vars, group_var)

# Print the adjusted results
adjusted_results_bin_size



# human footprint
summary_stats_HF <- HF_df %>%
  mutate(Type = factor(Type, levels = c("Sampling_sites", "Cities"))) %>%
  group_by(Type ) %>%
  summarise(across("HF_2000_2015", list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))
# Print the summary statistics
print(summary_stats_HF)


# Environmental variables variance
env_vars <- c("Clay", "Silt", "pH",  "Moisture", "SOC", "DOC", "LCP1", "LCP2", "RCP", "NH4_N", "NO3_N", "DON")

envdata <- meta_dat %>% as.tibble() %>% dplyr::select(c("Layer", env_vars)) %>%
  mutate(Layer = factor(Layer, levels = c("SUR", "SUB", "PL")))

summary_stats_env <- envdata %>%
  group_by(Layer) %>%
  summarise(across(env_vars, list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x)
  ), .names = "{.col}_{.fn}"))
# Print the summary statistics
print(summary_stats_env)

# Specify response variables and group variable
response_vars <- env_vars
group_var <- "Layer"

# Perform the pairwise Wilcoxon tests and get adjusted p-values
adjusted_results_env <- pairwise_wilcox_test(envdata, response_vars, group_var)

# Print the adjusted results
adjusted_results_env


# MAGs cheracteristics
MAGs_traits <- c("Clean_data", "completeness", "contamination", "Recovery_rate_MAGs")
summary_stats_MAGs <- MAG_tab %>% group_by(Layer, Index) %>%
  summarise(across("Value", list(
    mean = ~ mean(.x, na.rm = TRUE),
    sd = ~ sd(.x, na.rm = TRUE),
    se = ~ se(.x),
    min = ~ min(.x, na.rm = TRUE),
    q1 = ~ quantile(.x, 0.25, na.rm = TRUE),
    median = ~ median(.x, na.rm = TRUE),
    q3 = ~ quantile(.x, 0.75, na.rm = TRUE),
    max = ~ max(.x, na.rm = TRUE)
  ), .names = "{.col}_{.fn}"))
# Print the summary statistics
print(summary_stats_MAGs)

tt_survsrest_source <- tt_survsrest$table %>%
  rownames_to_column(var = "gene")

tt_subvsrest_source <- tt_subvsrest$table %>%
  rownames_to_column(var = "gene")

tt_plvsrest_source <- tt_plvsrest$table %>%
  rownames_to_column(var = "gene")

fig4data <- tt_survsrest_source %>%
  right_join(tt_subvsrest_source, by = "gene") %>%
  right_join(tt_plvsrest_source, by = "gene") %>%
  # filter(gene %in% c(C_names, fermentation_names, N_names, S_names, Other_names)) %>%
  arrange(factor(gene, levels = c(C_names, fermentation_names, N_names, S_names, Other_names))) %>%
  write.table(., file = "C:/Users/dell/Desktop/tt_all.csv", row.names = T, sep = ',', quote = FALSE)

tt_survsrest_source <- tt_survsrest$table[c(C_names, fermentation_names, N_names, S_names, Other_names), ]

write.table(table_all, file = "C:/Users/dell/Desktop/table_all.csv", row.names = T, sep = ',', quote = FALSE)

