# Create input fot GO figure
# Brenda Pardo
# 2023-05-15
# lmd: 2023-06-09: ejr

#Libraries
library(purrr)
library(dplyr)
library(tidyr)
#library(xlsx)
library(stringr)
library(here)
library(tidyverse)


##Divide GO terms by comparison and create input for GO-figure algoritm

#GO breakdown
go<- read_tsv(here("1dpa-reference/03_go/output/go-breakdown_lfc0.txt.gz"), header = TRUE, quote= "", comment.char = "")
##Stages
st<- c(paste0(c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "intact")
#Comparisons 
comp.set<- paste0("de-up-go_s_", st[2:9], "-s_1dpa-1dpa-ref_lfc0.txt.gz")
  
all.comp <- map(comp.set, function(comp){
    
    #Get go terms just for one comparison
    go.sub <- go %>%
      filter(cluster== comp) %>%
      separate(GeneRatio, c("no.gen", "tot.gen"), sep="/") %>%
      mutate(rat.gen= (as.numeric(no.gen)*100) /as.numeric(tot.gen)) %>%
      arrange(desc(rat.gen))
    
  write_tsv(go.sub, file= here(paste0("1dpa-reference/04_go-figure/input/go-breakdown_", comp)), row.names = FALSE, quote = FALSE)
  
  #Make table to give as an input to go-figure
  #Standard input consists of a tabular format file with two columns: GO term, P-value
  go.fig<- go.sub %>%
    # select(ID, p.adjust) 
    select_if(names(.) %in% c('ID', 'p.adjust'))

  comp2<- comp %>%
    str_replace("de-up-go_s_", "") %>%
    str_replace("-1dpa-ref_lfc0.txt.gz", "")%>%
    str_replace("s_", "") 
  
  write_tsv(go.fig, file=here(paste0("1dpa-reference/04_go-figure/input/input_file_", comp2, ".tsv.gz")), row.names = FALSE, col.names= FALSE, 
              quote = FALSE)
  
})



  
