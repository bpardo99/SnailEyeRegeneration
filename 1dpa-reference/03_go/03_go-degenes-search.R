# Serch list of go terms and their DE genes
# Brenda Pardo
# 2023-05-15
# lmd: 2023-06-09: ejr

#Libraries
library(dplyr)
library(stringr)
library(purrr)
library(here)
library(tidyverse)

#Comparisons
st<- c(paste0(c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "intact")
comp.set<- paste0("de-up-go_s_", st[2:9], "-s_1dpa-1dpa-ref_lfc0.txt.gz")

#List of go term-genes that will be the base of the table 
go.tab<- read.table(here("01_raw-data/genes2go-ref.txt.gz"), sep="\t",
                    header = TRUE, quote = "")

#GO breakdown
go.br<- read_tsv(here("1dpa-reference/03_go/output/go-breakdown_lfc0.txt.gz"), header = TRUE, quote= "", comment.char = "")

  all.comp <- map(comp.set, function(comp){
    
    #Get go terms just for one comparison
    go.sub <- go.br %>%
      filter(cluster== comp) 
    go.sub2<- unique(go.sub$ID)
    #306 go terms for example comparison
    
    #filter base list, by -my go list-
    go.t.sub<- go.tab %>%
      filter(go_term %in% go.sub2)
    
    #Confirm right amount of go terms
    # length(unique(go.t.sub$go_term))
    # > length(unique(go.t.sub$go_term))
    # [1] 305
    
    ## Filter base list by upregulated genes in the comp
    #Get upregulated genes
    comp.de<- comp %>%
      str_replace("-go", "") %>%
      str_replace("txt", "csv.gz")
    
    de.sub<- read_csv(here(paste0("1dpa-reference/02_dea/tables/", comp.de)), header = TRUE)
    de.sub<- unique(de.sub$gene_id)
    
    #Get number of genes
    # length(de.sub)
    # > length(de.sub)
    # [1] 2178
    
    #filter base list, by -my de genes upregulated list-
    go.t.sub2<- go.t.sub %>%
      filter(gene_id %in% de.sub)
    
    #Check number of genes
    # length(unique(go.t.sub$go_term))
    # > length(unique(go.t.sub$go_term))
    # [1] 305
    #Check number of go terms
    # length(unique(go.t.sub$gene_id))
    # > length(unique(go.t.sub$gene_id))
    # [1] 929
    d <- dim(go.t.sub2)[1]
    comparison <- rep(comp.de, times=d)
    
    go.t.sub3<- cbind(go.t.sub2, comparison)
    
  })
  
  df <- do.call("rbind", all.comp)
  
  write_tav(df, file= here("1dpa-reference/03_go/go-degenes-search.txt.gz"), quote=FALSE, row.names = FALSE)
  