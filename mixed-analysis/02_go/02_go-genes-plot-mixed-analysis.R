# GO plot using subset of GO terms selected from go-figure analysis and manual 
# selection (By Alice and Brenda). Manually selected genes that belong to each 
# go term are shown
# 05-25-2023

#Libraries
library(openxlsx)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(ggpubr)
library(cowplot)
library(readr)
library(here)


#Load selected go terms
bp.df<- read.xlsx(here("mixed-analysis/02_go/20230401-terms-from-gofigure-selection-BP.xlsx"))
bp<- bp.df$go_id #IDs
bp.tit<- bp.df$term #Terms

mf.df<- read.xlsx(here("mixed-analysis/02_go/20230401-terms-from-gofigure-selection-MF.xlsx"))
mf<- mf.df$go_id #IDs
mf.tit<- mf.df$term #Terms

cc.df<- read.xlsx(here("mixed-analysis/02_go/20230401-terms-from-gofigure-selection-CC.xlsx"))
cc<- cc.df$go_id #IDs
cc.tit<- cc.df$term #Terms

term.cat<- list(bp, mf, cc)
term.catn<- c("BP", "MF", "CC")
tit.cat<- list(bp.tit, mf.tit, cc.tit)

#All go enriched
#references<- c("1dpa", "intact", "sequential")
references<- c("1dpa", "sequential")
breakdown<- map(references, function(x) {
  # read in GO enrichment
  go <- read_tsv(here(paste0(x, "-reference/03_go/output/go-breakdown_lfc0.txt")))
})
breakdown2<- do.call("rbind", breakdown)

enriched<- map(term.cat, function(t){
  enr<- breakdown2 %>%
    filter(ID %in% t)
})
enriched2<- do.call("rbind", enriched)
enriched3<- enriched2$ID %>%
  unique() %>%
  length() #all the list of go terms form part of the enriched in 1dpa and sequential analysis
l.bp.cc.mf<- c(bp.tit, mf.tit, cc.tit)


#Load selected list of genes
#Wound and division
wound.div<- read.xlsx(here("1dpa-reference/02_dea/selection-woundhealing-celldivision-aai.xlsx"), colNames = FALSE)
wound.div<- wound.div$X1
#TFs
tf<- read.xlsx(here("1dpa-reference/02_dea/selection-markers-tfs-aai.xlsx"), colNames = FALSE)
tf<- tf$X1
genl<- unique(c(wound.div, tf))
#Filter gene selection by if they are DE
references<- c("1dpa", "intact", "sequential")
de.l<- map(references, function(reference){
  file<-read.table(here(paste0(reference, "-reference/02_dea/tables/de-up_lfc2.csv")), sep="\t", quote = "", header= TRUE)
})
de.l<- do.call("rbind", de.l) %>%
  unique()

genl.d<- data.frame(genl) %>%
  rename(gene_id=genl)
genl.d.f<- genl.d %>% 
  filter(gene_id %in% de.l$gene_id)
genl<- genl.d.f$gene_id


#Load de go-gene data base
go.db<- read.table(here("mixed-analysis/02_go/pcan_go_20210924_with_ancestors.txt"), sep="\t", quote = "", header= TRUE)



#Read references
ref<- read.table(here("01_raw-data/gene-ref.txt"), sep="\t", quote = "", header= TRUE)
ref$id.desc= paste0(ref$gene_id,"_", ref$description)
ref <- ref %>%
  select(gene_id, id.desc)

#Load zscore table with reference
tpm.z<- read.csv(here("02_processed-data/tpm-zscores.csv"), header = TRUE)
zsc <- tpm.z %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()

#Filter DE gene list from go.db by selected 
a<- go.db %>%
  filter(gene_id %in% genl) 
#Basically we are asking if a gene is DE and what DE enriched go-terms include it





# #Test to make sure all the genes in the list Alice gave me are DE in any comparison
# de.1dpa<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/1dpa-reference/02_dea/tables/de-up_lfc0.csv", sep="\t", quote = "", header= TRUE)
# de.int<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/intact-reference/02_dea/tables/de-up_lfc0.csv", sep="\t", quote = "", header= TRUE)
# de.seq<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/sequential-reference/02_dea/tables/de-up_lfc0.csv", sep="\t", quote = "", header= TRUE)
# de.conf.1dpa<- de.1dpa %>%
#   filter(gene_id %in% genl)
# de.conf.int<- de.int %>%
#   filter(gene_id %in% genl)
# de.conf.seq<- de.seq %>%
#   filter(gene_id %in% genl)
# de.conf.all<- rbind(de.conf.1dpa, de.conf.int, de.conf.seq)
# de.conf.all.uniq<- unique(de.conf.all) 
# #Just 90 genes of Alice's list are de in comparisons with lfc>=0 and thus I will filter the list



g<-  pmap(list(term.cat, tit.cat, term.catn), function(term.c, tit.c, term.cn){
  # term.c<- term.cat[[1]]
  # tit.c<- tit.cat[[1]]
  # term.cn<- term.catn[[1]]
  
  d<- pmap(list(term.c, tit.c), function(term, tit){
    # term<- term.c[[1]]
    # tit<- tit.c[[1]]
  
    #Filter list of selected and de genes by go-term
    b<- a %>%
      filter(go_term %in% term)
    
    #Get gene list to get its expression
    c<- b$gene_id %>%
      unique()
    
    #Get its expression
    d.p<- zsc %>%
      filter(gene_id %in% c)
    
    if(dim(d.p)[1] == 0){
      
      print(paste0("dimension is 0-", tit))
      
      return(NULL)
      
      
    } else if(dim(d.p)[1] < 2){
      
      print(paste0("dimension is < 2-", tit))
      
      d.p2 <- d.p %>% 
        mutate(go_term= rep(tit, dim(d.p)[1])) %>%
        pivot_longer(starts_with("s"), names_to = "st", values_to = "zscore")
      
      ##Arrange data frames
      d.p3 <- d.p2 %>% 
        mutate(order = factor(id.desc)) 
      
      return(d.p3)
      
  
    } else{
      print(paste0("we have a plot-", tit))
      
      m<- d.p
      rownames(m)<- m$id.desc
      m$gene_id<- m$id.desc<- NULL
      m <- as.matrix(m)
      ##Cluster
      h <- hclust(dist(m), method = 'ward.D2')
      row.order <- rownames(m[h$order,])
      
      ##Pivot dataframes
      d.p2 <- d.p %>% 
        mutate(go_term= rep(tit, dim(d.p)[1])) %>%
      pivot_longer(starts_with("s"), names_to = "st", values_to = "zscore")
      
      ##Arrange data frames
      d.p3 <- d.p2 %>% 
        mutate(order = factor(id.desc, levels = row.order)) %>%
        arrange(order) 
      
      return(d.p3)
    }
    
  })  
    
  e <- d[!sapply(d,is.null)] 
  f<- do.call("rbind", e)
    
  
  
  ##Plots 
  mypalette<-rev(brewer.pal(11,"RdBu"))
  lev<- rev(c(paste0("s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact"))
  
  plot <- ggplot(f, aes(y=factor(st, levels = lev), x=order)) +
    geom_tile(aes(fill= zscore)) +
    scale_fill_gradientn(colors=mypalette) +
    theme_cowplot()+
    theme(axis.text.x = element_text(angle= 90, hjust= 0.99, vjust = 0.95), panel.spacing.x = unit(0.5, "pt")) +
    facet_grid(~go_term, scales = "free_x", space = "free_x") +
    ylab("Stages") +
    xlab("Genes") #+
    #theme(
      #strip.text.x = element_text(
        #size = 3, color = "black"))
  
  ggsave(here(paste0("mixed-analysis/02_go/figures/go-genes-plot-", term.cn, "-mixed-reference-lfc2.pdf")), plot, height=10, width=20, limitsize=FALSE)
    
})
  
