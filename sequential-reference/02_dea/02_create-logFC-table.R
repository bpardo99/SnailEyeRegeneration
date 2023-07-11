# Create logFC tables from DEA results
# Brenda Pardo
# 2023-05-15

#Libraries
library(dplyr)
library(purrr)
library(readr)
library(here)


#Load logFC tables for each comparison
files <- paste0("eTest_s_", c(paste0(c(3, 6, 9, 12, 15, 21, 28), "dpa"), "intact"), "-s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa-seq-ref.csv.gz")
et <- map(files, function(x){
  f<- read.csv(here(paste0("sequential-reference/02_dea/tables/", x)))
})
names(et) <- files

#Get logFC for each gene in each comparison
logfc<- map(et, function(x){
  col <- data.frame(x) %>%
    dplyr::select(logFC, gene_id) 
})

logfc.t<- logfc[[1]] %>%
  dplyr::rename(s3dpa_s1dpa = logFC)  %>%
  left_join(y= logfc[[2]], by= "gene_id") %>%
  dplyr::rename(s6dpa_s3dpa = logFC) %>%
  left_join(y= logfc[[3]], by= "gene_id") %>%
  dplyr::rename(s9dpa_s6dpa = logFC)  %>%
  left_join(y= logfc[[4]], by= "gene_id") %>%
  dplyr::rename(s12dpa_s9dpa = logFC)  %>%
  left_join(y= logfc[[5]], by= "gene_id") %>%
  dplyr::rename(s15dpa_s12dpa = logFC)  %>%
  left_join(y= logfc[[6]], by= "gene_id") %>%
  dplyr::rename(s21dpa_s15dpa = logFC)  %>%
  left_join(y= logfc[[7]], by= "gene_id") %>%
  dplyr::rename(s28dpa_s21dpa = logFC)  %>%
  left_join(y= logfc[[8]], by= "gene_id") %>%
  dplyr::rename(sintact_s28dpa = logFC) 

st <- paste0(c(paste0("s", c(3, 6, 9, 12, 15, 21, 28), "dpa"), "sintact"), "_s", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa")
logfc.t <- logfc.t %>% dplyr::select(gene_id, all_of(st))

#Save table with logFC
write_csv(logfc.t, here("sequential-reference/02_dea/tables/logFC-seq-ref.csv.gz"))
