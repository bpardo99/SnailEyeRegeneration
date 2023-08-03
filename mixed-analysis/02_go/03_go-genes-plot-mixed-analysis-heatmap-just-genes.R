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

#Load selected list of genes
#Wound and division
# wound.div<- read.xlsx("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/1dpa-reference/02_dea/selection-woundhealing-celldivision-aai.xlsx", colNames = FALSE)
# wound.div<- wound.div$X1
# #TFs
# tf<- read.xlsx("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/1dpa-reference/02_dea/selection-markers-tfs-aai.xlsx", colNames = FALSE)
# tf<- tf$X1

aai<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/mixed-analysis/02_go/go-genes-plot-mixed-reference-just-selected-genes-de-lfc0_AAI_selection.txt", sep="\t", quote = "", header= TRUE)
aai<- aai$gene_id

#genl<- unique(c(wound.div, tf))
genl<- unique(aai)
#Filter gene selection by if they are DE
references<- c("1dpa", "intact", "sequential")
de.l<- map(references, function(reference){
  file<-read.table(paste0("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/", reference, "-reference/02_dea/tables/de-up_lfc0.csv"), sep="\t", quote = "", header= TRUE)
})
de.l<- do.call("rbind", de.l) %>%
  unique()

genl.d<- data.frame(genl) %>%
  dplyr::rename(gene_id=genl)
genl.d.f<- genl.d %>% 
  filter(gene_id %in% de.l$gene_id)
genl<- genl.d.f$gene_id


#Read references
ref<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/01_raw-data/gene-ref.txt", sep="\t", quote = "", header= TRUE)
ref$id.desc= paste0(ref$gene_id,"_", ref$description)
ref <- ref %>%
  select(gene_id, id.desc) %>%
  mutate(id.desc = str_replace(id.desc, "-like", ""))

#Load zscore table with reference
tpm.z<- read.csv("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/02_processed-data/tpm-zscores.csv", header = TRUE)
tpm.lp<- read.csv("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/02_processed-data/tpm-mean-log1p.csv", header = TRUE)

exp.list<- list(tpm.lp, tpm.z)
exp.nam<- list("log1p.tpm", "zscore")

l<- pmap(list(exp.list, exp.nam), function(exp, nam){
  
  zsc <- exp %>%
    left_join(y=ref, by= "gene_id") %>%
    na.omit()
d.p<- zsc %>%
  filter(gene_id %in% genl)

m<- d.p
rownames(m)<- m$id.desc
m$gene_id<- m$id.desc<- NULL
m <- as.matrix(m)
##Cluster
h <- hclust(dist(m), method = 'ward.D2')
row.order <- rownames(m[h$order,])

##Pivot dataframes
d.p2 <- d.p %>% 
  pivot_longer(starts_with("s"), names_to = "st", values_to = "zscore")

##Arrange data frames
d.p3 <- d.p2 %>% 
  mutate(order = factor(id.desc, levels = row.order)) %>%
  arrange(order) #%>%
  #mutate(id.desc = str_replace(id.desc, "-like", ""),
         #order = str_replace(order, "-like", "")
         
#)


mypalette<-rev(brewer.pal(11,"RdBu"))
lev<- c(paste0("s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact")

plot <- ggplot(d.p3, aes(x=factor(st, levels = lev), y=order)) +
  geom_tile(aes(fill= zscore)) +
  scale_fill_gradientn(colors=mypalette) +
  theme_cowplot()+
  theme(axis.text.x = element_text(angle= 90, hjust= 0.99, vjust = 0.95), panel.spacing.x = unit(0.5, "pt")) +
  ylab("Genes") +
  xlab("Stages") +
  labs(fill=nam)


ggsave(paste0("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/mixed-analysis/02_go/figures/go-genes-plot-mixed-reference-just-selected-genes-de-lfc0-", nam, ".pdf"), plot, height=20, width=10, limitsize=FALSE)

return(row.order)
})


g.to.wr<- l[[2]] 

write.table(g.to.wr, file=paste0("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/mixed-analysis/02_go/tables/go-genes-plot-mixed-reference-just-selected-genes-de-lfc0.txt"), row.names = FALSE, col.names= FALSE, 
            quote = FALSE, sep = '\t')
