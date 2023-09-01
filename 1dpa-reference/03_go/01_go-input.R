# GO enrichment analysis using 1dpa time-point as the reference.
# Brenda Pardo
# 2023-05-15
# lmd: 2023-06-09

#Libraries
library(purrr)
library(stringr)
library(readr)
library(here)

##Create input tables

comp<- paste0(c(paste0(c(3, 6, 9, 12, 15, 21, 28), "dpa"), "intact"), "-s_1dpa")
path<- here("1dpa-reference/02_dea/tables/")
files<- paste0("de-up_s_", comp,"-1dpa-ref_lfc")
files<- c(paste0(files, 0, ".csv.gz"), 
          paste0(files, 1, ".csv.gz"), 
          paste0(files, 2, ".csv.gz"))

go.input<- map(files, function(x){
  f<- read.csv(paste0(path, x))
  de<- f %>%
    dplyr::select(gene_id)
  
  nam.fil<- x %>%
    str_replace(pattern = "de-up", replacement = "de-up-go")%>%
    str_replace(pattern = ".csv.gz", replacement = ".txt.gz")
  
  write_csv(de, 
            file=here(paste0("1dpa-reference/03_go/input/", nam.fil)))

})

