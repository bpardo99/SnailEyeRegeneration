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
st <- paste0("s_", c(paste0(c(3, 6, 9, 12, 15, 21, 28), "dpa"), "intact"), "-s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa")
et <- pmap(list(files, st), function(x, st){
  f<- read.csv(here(paste0("sequential-reference/02_dea/tables/", x)))
  colnames(f)[2:5] <- paste0(st, ".", colnames(f)[2:5])
  return(f)
  })
names(et) <- files


et.t<-et[[1]] %>%
  left_join(y= et[[2]], by= "gene_id") %>%
  left_join(y= et[[3]], by= "gene_id") %>%
  left_join(y= et[[4]], by= "gene_id") %>%
  left_join(y= et[[5]], by= "gene_id") %>%
  left_join(y= et[[6]], by= "gene_id") %>%
  left_join(y= et[[7]], by= "gene_id") %>%
  left_join(y= et[[8]], by= "gene_id")


#Save table with logFC
write_csv(et.t, here("sequential-reference/02_dea/tables/logFC-logCPM-Pval-seq-ref.csv.gz"))
