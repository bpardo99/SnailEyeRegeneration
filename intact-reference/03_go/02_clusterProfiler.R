# Use clusterprofiler to run GO enrichment analysis in upregulated DE genes
# Brenda Pardo
# 2023-05-15
# Script written by ejr: 2021-11-25

#Libraries
library(tidyverse)
library(clusterProfiler)
library(GOSemSim)
library(cowplot)
library(enrichplot)
library(org.Pcanaliculata.eg.db)
library(here)

#to install Pcan ORGDB
#install.packages("/home/ejr/sciproj/SCI-003759-SBPCAN/genomes/GCF_003073045.1_ASM307304v1/analysis/orgpkg/org.Pcanaliculata.eg.db", repos=NULL)
#pcGO <- godata('org.Pcanaliculata.eg.db', ont="MF", keytype = "GID")



pth <- here("intact-reference/03_go/input/")
setwd(pth)

comp<- paste0(paste0(c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "-s_intact")
files<- paste0("de-up-go_s_", comp,"-intact-ref_lfc")

#lfc thresholfs
lfc_list <- c(0:2)

go.lfc <-map(lfc_list, function(lfc){
  file_list<- paste0(files, lfc, ".txt")
                
  go <- map(file_list, function(x) {
    
    gene_list <- read_tsv(x)   
    
    enrich_out <- map(c("MF", "CC", "BP"), function(y) {
      ego <- enrichGO(gene      = gene_list$gene_id,
                      OrgDb         = org.Pcanaliculata.eg.db,
                      keyType       = "GID", 
                      ont           = y,
                      pAdjustMethod = "BH",
                      pvalueCutoff  = 0.01,
                      qvalueCutoff  = 0.05,
                      readable      = TRUE)
      print(paste0(x, "_", y))
      if (!is.null(ego)) {
        if((length(ego@termsim) > 0)) {
        p <- barplot(ego, showCategory=20)
        filename = here(paste0("intact-reference/03_go/output/", x, y, "_bar.png"))
        save_plot(filename, p, base_width=12, base_height=6)
      }    
      as.data.frame(ego) %>% mutate(ont = y) %>% mutate(cluster = x)
      }
    })
    
    do.call("rbind", enrich_out)
  })

  df <- do.call("rbind", go)

  write_tsv(df, 
          file= here(paste0("intact-reference/03_go/output/go-breakdown_lfc", lfc, ".txt.gz")))
})

