# Create logFC tables from DEA results
# Brenda Pardo
# 2023-05-15

#Libraries
library(dplyr)
library(purrr)
library(readr)
library(here)


#Load logFC tables for each comparison
files <- c(paste0("eTest_s_", c(3, 6, 9, 12, 15, 21, 28), "dpa-s_1dpa-1dpa-ref.csv.gz"), "eTest_s_intact-s_1dpa-1dpa-ref.csv.gz")
st<- c(paste0("s_", c(3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact")
et <- pmap(list(files, st), function(x, st){
  f<- read_csv(here(paste0("1dpa-reference/02_dea/tables/", x)))
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
write_csv(et.t, here("1dpa-reference/02_dea/tables/logFC-logCPM-Pval-1dpa-ref.csv.gz"))
