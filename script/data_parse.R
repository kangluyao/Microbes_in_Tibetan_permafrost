
ko_tax_TPM.file <- 'E:/permafrost/bioinfo_results/result/eggnog/ko_tax_TPM.txt'
ko_tax_TPM_table <- read.table(ko_tax_TPM.file, header = TRUE, sep = "\t", as.is = TRUE, 
                   stringsAsFactors = FALSE, comment.char = "", check.names = FALSE)
nrow(ko_tax_TPM_table)

ko_tax_TPM_table[1:250, 1]
parse_data <- function(data) {
  require(dplyr)
  data = as.data.frame(data)
  df = NULL
  for (i in 1:nrow(data)) {
    ks = data[i, 1]
    ks = unlist(strsplit(ks, ','))
    if (length(ks) > 1) {
      temp_rows = data.frame(ko = ks, data[, -1] %>% slice(rep(i, length(ks))))
    } else {
      temp_rows = data[i, ]
    } 
    if (is.null(df)) {
      df <- temp_rows
    } else {
      df <- rbind(df, temp_rows)
    }
  }
  return(df)
}

library(foreach)
library(doParallel)
# 查看自己电脑的物理核心数量
detectCores(logical = F)
# 创建一个集群并注册
cl <- makeCluster(8)
registerDoParallel(cl)

# 启动并行计算
x1 <- foreach(i = 1:100000, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
  ks = ko_tax_TPM_table[i, 1]
  ks = unlist(strsplit(ks, ','))
  if (length(ks) > 1) {
    temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
  } else {
    temp_rows = ko_tax_TPM_table[i, ]
  } 
  return(temp_rows)
}
write.table(x1, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_1.txt', sep = '\t', quote = F)
rm(x1)
x2 <- foreach(i = 100001:200000, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
    ks = ko_tax_TPM_table[i, 1]
    ks = unlist(strsplit(ks, ','))
    if (length(ks) > 1) {
      temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
    } else {
      temp_rows = ko_tax_TPM_table[i, ]
    } 
    return(temp_rows)
}
write.table(x2, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_2.txt', sep = '\t', quote = F)
rm(x2)
x3 <- foreach(i = 200001:300000, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
  ks = ko_tax_TPM_table[i, 1]
  ks = unlist(strsplit(ks, ','))
  if (length(ks) > 1) {
    temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
  } else {
    temp_rows = ko_tax_TPM_table[i, ]
  } 
  return(temp_rows)
}
write.table(x3, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_3.txt', sep = '\t', quote = F)
rm(x3)
x4 <- foreach(i = 300001:400000, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
  ks = ko_tax_TPM_table[i, 1]
  ks = unlist(strsplit(ks, ','))
  if (length(ks) > 1) {
    temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
  } else {
    temp_rows = ko_tax_TPM_table[i, ]
  } 
  return(temp_rows)
}
write.table(x4, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_4.txt', sep = '\t', quote = F)
rm(x4)
x5 <- foreach(i = 400001:500000, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
  ks = ko_tax_TPM_table[i, 1]
  ks = unlist(strsplit(ks, ','))
  if (length(ks) > 1) {
    temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
  } else {
    temp_rows = ko_tax_TPM_table[i, ]
  } 
  return(temp_rows)
}
write.table(x5, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_5.txt', sep = '\t', quote = F)
rm(x5)
x6 <- foreach(i = 500001:630161, .combine = rbind, .packages = 'dplyr', .errorhandling = "pass") %dopar% {
  ks = ko_tax_TPM_table[i, 1]
  ks = unlist(strsplit(ks, ','))
  if (length(ks) > 1) {
    temp_rows = data.frame(ko = ks, ko_tax_TPM_table[, -1] %>% slice(rep(i, length(ks))))
  } else {
    temp_rows = ko_tax_TPM_table[i, ]
  } 
  return(temp_rows)
}
write.table(x6, file = 'E:/permafrost/bioinfo_results/result/eggnog/parse_dat_6.txt', sep = '\t', quote = F)
rm(x6)

# 在计算结束后别忘记关闭集群
stopImplicitCluster()
stopCluster(cl)
