# Heatmaps for comparing expression between regeneration and development
# Plot of upregulated DE genes at timeppoints of regeneration. I plot expression 
# in regeneration as well as in development. Genes are clustered based on their 
# expression across development and regeneration
# Brenda Pardo
# 2023-05-15

#Libraries
library(dplyr)
library(tidyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(purrr)
library(here)


#Read gene reference
ref=read_tsv(here("01_raw-data/gene-ref.txt.gz"), header=TRUE, quote = "")
ref$id.desc= paste0(ref$gene_id,"_", ref$description)
ref <- ref %>%
  select(gene_id, id.desc) 

#Read logFC regeneration 
reg.lfc<- read_tsv(here("1dpa-reference/02_dea/tables/logFC-1dpa-ref.csv.gz"), sep= ",", header=TRUE)
reg.lfc <- reg.lfc %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()

#Read zscores regeneration
reg.z<- read.table(here("02_processed-data/tpm-zscores.csv.gz"), sep= ",", header=TRUE)
colnames(reg.z) = gsub("_", "", colnames(reg.z))
reg.z <- reg.z %>%
  dplyr::rename(gene_id=geneid) %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()
reg.exp<- list(reg.lfc, reg.z)

#Read logFC development
d<-read_tsv(here("01_raw-data/dev_expr_mixed_ref_20220310.txt.gz"), header=TRUE)
dev.lfc<- d %>%
  select(target_id, contrast, logFC) %>%
  mutate(contrast = str_replace(contrast, "Stage ", "st")) %>% 
  pivot_wider(names_from = contrast, values_from = logFC) %>%
  dplyr::rename(h= Hatchling) %>%
  dplyr::rename(gene_id= target_id) %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()

#Read zscores development
d<-read.table(here("02_processed-data/dev-cpm-zscores.csv.gz"), sep= ",", header=TRUE)
dev.z<- d %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()
dev.exp<- list(dev.lfc, dev.z)

exp.lbl<- c("logFC", "zscores")

#Set levels
lv.lfc<- c(paste0("s", c(3, 6, 9, 12, 15, 21, 28), "dpa"), "sintact", paste0("st", 2:10), "h")
lv.z<- c(paste0("s", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "sintact", paste0("st", 2:10), "h")
lv<- list(lv.lfc, lv.z)

#1. Extract DE genes at all time points during regeneration
path<- here("1dpa-reference/02_dea/tables/")
st<- c(paste0(c(3, 6, 9, 12, 15, 21, 28), "dpa"), "intact")
list.de<- paste0(path, "de-up_s_", st, "-s_1dpa-1dpa-ref_lfc2.csv.gz")


#Plor sizes
p.height<- c(40, rep(170, times=3), 220, rep(250, times=3))
#Loop for different expression matrixes
exp<- pmap(list(reg.exp, dev.exp, exp.lbl, lv), function(reg.exp, dev.exp, exp.lbl, lv){
  #Loop for different comparisons
  filt.dev<- pmap(list(list.de, st, p.height), function(comp, st, p.height){
    comp.r<-read.table(comp,  sep= ",", header=TRUE)
    genes<- comp.r$gene_id
    
    #Regeneration, subset data for plotting 
    d.p.reg<- reg.exp %>%
      filter(gene_id %in% genes)
    
    #Development, subset data for plotting 
    d.p.dev<- dev.exp %>%
      filter(gene_id %in% genes)
    
    #Join tables
    d.p <- full_join(d.p.reg, d.p.dev, by= c("gene_id", "id.desc")) %>%
      na.omit() #I HAD TO RESTRICT TO PLOT GENES THAT ARE PRESENT IN BOTH LISTS
    
    #Matrix for clustering based on development and regeneration
    m<- as.data.frame(d.p)
    rownames(m)<- m$id.desc
    m$gene_id<- m$id.desc<- NULL
    m <- as.matrix(m)
    #Cluster
    h <- hclust(dist(m), method = 'ward.D2')
    row.order <- rownames(m[h$order,])
    
    #Pivot both reg and dev dataframes
    d.p2 <- d.p %>% pivot_longer(starts_with("s") | starts_with("h"), names_to = "compa", values_to = "logFC")
    
    #Arrange data frames
    d.p3 <- d.p2 %>% 
      mutate(order = factor(id.desc, levels = row.order)) %>%
      arrange(order) %>%
      #Set maximum values for color range
      mutate(logFC = replace(logFC, logFC >= 4, 4))  %>%
      mutate(logFC = replace(logFC, logFC <= -4, -4))
    
    mypalette<-rev(brewer.pal(11,"RdBu"))
    #Plot
    av <- max(abs(d.p3$logFC))
    
    plot <- ggplot(d.p3, aes(y=order, x=factor(compa, levels = lv))) +
      geom_tile(aes(fill= logFC)) +
      scale_fill_gradientn(colors=mypalette, limits= c(-av, av)) +
      theme(axis.text.x = element_text(angle=90)) +
      ylab("Genes") +
      xlab("Stages")+ 
      labs(fill = exp.lbl)
    
    ggsave(plot, filename=here(paste0("1dpa-reference/05_overlap-reg-dev/figures/reg-deup_", st,"-1dpa-ref-logFC2_dev-mixref_clusregdev-", exp.lbl, ".pdf")), width=10, height= p.height, limitsize = FALSE)
  })
})

