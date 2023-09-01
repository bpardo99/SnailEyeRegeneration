# Create input fot GO figure
# Brenda Pardo
# 2023-05-15

#Libraries
library(purrr)
library(dplyr)
library(tidyr)
#library(xlsx)
library(readr)
library(stringr)
library(here)


##Divide GO terms by comparison and create input for GO-figure algoritm

#GO breakdown
go<- read.table(here("sequential-reference/03_go/output/go-breakdown_lfc0.txt.gz"), sep="\t", header = TRUE, quote= "", comment.char = "")
##Stages
st <- paste0(c(paste0("s_", c(3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact"), "-s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa")
#Comparisons 
comp.set<- paste0("de-up-go_", st, "-seq-ref_lfc0.txt.gz")


all.comp <- map(comp.set, function(comp){
    
    #Get go terms just for one comparison
    go.sub <- go %>%
      filter(cluster== comp) %>%
      separate(GeneRatio, c("no.gen", "tot.gen"), sep="/") %>%
      mutate(rat.gen= (as.numeric(no.gen)*100) /as.numeric(tot.gen)) %>%
      arrange(desc(rat.gen))
    
  write.table(go.sub, file= here(paste0("sequential-reference/04_go-figure/input/go-breakdown_", comp)), row.names = FALSE, quote = FALSE, sep = '\t')
  
  #Make table to give as an input to go-figure
  #Standard input consists of a tabular format file with two columns: GO term, P-value
  go.fig<- go.sub %>%
    # select(ID, p.adjust) 
    select_if(names(.) %in% c('ID', 'p.adjust'))

  comp2<- comp %>%
    str_replace("de-up-go_s_", "") %>%
    str_replace("-seq-ref_lfc0.txt", "")%>%
    str_replace("s_", "") 
  
  write.table(go.fig, file=here(paste0("sequential-reference/04_go-figure/input/input_file_", comp2, ".tsv.gz")), row.names = FALSE, col.names= FALSE, 
              quote = FALSE, sep = '\t')
  
})



  
