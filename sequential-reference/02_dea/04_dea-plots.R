#number-de-genes-lfc", lfc,"_seq-ref.txt" DE gene plots
# Variety of plots from DEA output
# Brenda Pardo
# 2023-05-15


#Library
library(purrr)
library(tidyr)
library(dplyr)
library(stringr)
library(ggplot2)
library(cowplot)
library(openxlsx)
library(RColorBrewer)
library(here)
library(tidyverse)


##Plot number of DE genes across comparisons

lfc<- paste0("lfc", 0:2)
files <-here(paste0("sequential-reference/02_dea/tables/number-de-genes-", lfc, "_seq-ref.txt"))

de.ct <- map(files, function(x){
  read_tsv(x)
})
names(de.ct) <- lfc

order <- c("3dpa-1dpa", "6dpa-3dpa", "9dpa-6dpa", "12dpa-9dpa", "15dpa-12dpa", 
           "21dpa-15dpa", "28dpa-21dpa", "intact-28dpa")

p.count <- map(lfc, function(x){
  
  #Buil dataframe
  de.ct.p= de.ct[[x]] %>%
    mutate(tag=factor(tag, levels= de.ct[[x]]$tag)) %>%
    pivot_longer(cols=c(de.ct.up, de.ct.down), names_to = "direction", values_to = "number")
  
  #Arrange upregulated
  de.ct.p.up<- de.ct.p %>%
    filter(direction=="de.ct.up") %>%
    mutate(direction = str_replace(direction, "de.ct.up", "Upregulated genes")) %>%
    dplyr::rename("Direction" = "direction") %>%
    mutate(tag = str_replace_all(tag, "s_", "")) %>%
    mutate(tag = str_replace_all(tag, "-seq-ref", "")) %>%
    mutate(tag= factor(tag, levels=order)) 
  
  #Arrange downregulated
  de.ct.p.down<- de.ct.p %>%
    filter(direction=="de.ct.down") %>%
    mutate(direction = str_replace(direction, "de.ct.down", "Downregulated genes")) %>%
    dplyr::rename("Direction" = "direction") %>%
    mutate(tag = str_replace_all(tag, "s_", "")) %>%
    mutate(tag = str_replace_all(tag, "-seq-ref", "")) %>%
    mutate(tag= factor(tag, levels=order))
  
  #Plot upregulated
  p_up<- de.ct.p.up %>%
    ggplot(aes(fill=Direction, y=number, x=tag)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = "#ADD9F4") +
    ggtitle(paste0("Number of DE genes accross vs-intact comparisons, upregulated")) +
    xlab("Comparison") + 
    ylab("Number of DE genes") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=90),
          plot.title = element_text(size = 10)) +
    ylim(NA, max(c(de.ct.p.up$number, de.ct.p.down$number)))
  pdf(file=here(paste0("sequential-reference/02_dea/figures/number-de-genes-up-", x, "_seq-ref.pdf")))
  print(p_up)
  dev.off()
  
  #Plot downregulated
  p_down<- de.ct.p.down %>%
    ggplot(aes(fill=Direction, y=number, x=tag)) + 
    geom_bar(position="stack", stat="identity") +
    scale_fill_manual(values = "#984447") +
    ggtitle(paste0("Number of DE genes accross vs-intact comparisons, downregulated")) +
    xlab("Comparison") + 
    ylab("Number of DE genes") +
    theme_cowplot() +
    theme(axis.text.x = element_text(angle=90),
          plot.title = element_text(size = 10)) +
    ylim(NA, max(c(de.ct.p.up$number, de.ct.p.down$number)))
  
  
  pdf(file=here(paste0("sequential-reference/02_dea/figures/number-de-genes-down-", x, "_seq-ref.pdf")))
  print(p_down)
  dev.off()
  
})










##Pax Piwi Rho expression, lineplot

#Expression tables
tpm<- read_csv(here("02_processed-data/tpm-mean.csv.gz"))
tpm.l1p<- read_csv(here("02_processed-data/tpm-mean-log1p.csv.gz"))
tpm.z<- read_csv(here("02_processed-data/tpm-zscores.csv.gz"))
exp.tab<- list(tpm, tpm.l1p, tpm.z)
exp.nam<- c("tpm", "log1p-tpm", "zscore")

#Genes to plot
gn <- c("LOC112559942", "LOC112576458", "LOC112563208")
names(gn) <- c("pax6", "rho", "piwi")

#Stages
st<- c(paste0(c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "intact")

#Plot
#Non smothened line
ppr<- pmap(list(exp.tab, exp.nam), function(exp, nam){
  
  exp1 <- exp %>%
    filter(gene_id %in% gn) 
  
  exp.t <- exp1 %>% pivot_longer(starts_with("s_"), names_to = "comp", values_to = "val") %>%
    mutate(comp = str_replace(comp, "s_", "")) %>%
    mutate(gene = gene_id) %>%
    mutate(gene = str_replace(gene, "LOC112559942", "pax6")) %>%
    mutate(gene = str_replace(gene, "LOC112576458", "rho")) %>%
    mutate(gene = str_replace(gene, "LOC112563208", "piwi"))
  
  p <- exp.t %>%
    ggplot(aes(x=factor(comp, levels = st), y= val, group= gene)) +
    geom_line(aes(colour= gene)) +
    theme_classic() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab("") +
    ylab(nam) 
return(p)
})

#With log1p(tpm) and smothened line
exp1 <- tpm.l1p %>%
  filter(gene_id %in% gn) 

exp.t <- exp1 %>% pivot_longer(starts_with("s_"), names_to = "comp", values_to = "val") %>%
  mutate(comp = str_replace(comp, "s_", "")) %>%
  mutate(gene = gene_id) %>%
  mutate(gene = str_replace(gene, "LOC112559942", "pax6")) %>%
  mutate(gene = str_replace(gene, "LOC112576458", "rho")) %>%
  mutate(gene = str_replace(gene, "LOC112563208", "piwi"))

p <- exp.t %>%
  ggplot(aes(x=factor(comp, levels = st), y= val, group= gene)) +
  geom_smooth(aes(colour= gene), method = "lm", formula = y ~ poly(x, 3), se = FALSE) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=90)) +
  xlab("") +
  ylab("log1p-tpm (smoothened)") 

pdf(file=here(paste0("sequential-reference/02_dea/figures/pax-pw1-rho.pdf")), width = 5, height= 5)
ppr
p
dev.off()










##TF expression plots

#Import list of DE genes
####******************************** is this right!!!
de.list<- read.csv(here("1dpa-reference/02_dea/tables/de-up_lfc2.csv.gz"))
de.list<- de.list$gene_id

#Read TF list
tfs <- read.table(file=here("01_raw-data/putative_tfs-ejr.list"))
names(tfs)<- "gene_id"

#Keep tfs if they are DE
tfs<- tfs %>%
  filter(gene_id %in% de.list) 
tfs <- tfs$gene_id

chunk_length= 50
tfs_list= split(tfs, ceiling(seq_along(tfs) / chunk_length)) #split into smaller chunks

#Read logFC table
lfc.seq<-read.csv(here("sequential-reference/02_dea/tables/logFC-seq-ref.csv.gz"))
exp.tab<- list(lfc.seq, tpm.l1p, tpm.z)
exp.nam<- c("logFC", "log1p-tpm", "zscore")

#Reaf references
ref<- read_tsv(here("01_raw-data/gene-ref.txt.gz"))

ref$id.desc= paste0(ref$gene_id,"_", ref$description)
ref <- ref %>%
  select(gene_id, id.desc)

limscal<- list(c(4, -4), c(8.738117, 0), c(2, -2))

a<- pmap(list(exp.tab, exp.nam, limscal), function(exp.t, exp.n, lsc){
  
  b<- map(tfs_list, function(tfs.s){ 
    #Load exp
    file<- exp.t
    
    p <- file %>%
      inner_join(y= ref, by="gene_id") %>%
      filter(gene_id %in% tfs.s)
    
    #Matrix for clustering
    m= as.data.frame(p)
    rownames(m)<- m$id.desc
    m$gene_id<- m$id.desc<- NULL
    m <- as.matrix(m)
    
    #Cluster
    h <- hclust(dist(m), method = 'ward.D2')
    row.order <- rownames(m[h$order,])
    
    p2 <- p %>% pivot_longer(starts_with("s"), names_to = "comp", values_to = "val")
    #av <- max(p2$val) #of all tf list
    p3 <- p2 %>% mutate(order = factor(id.desc, levels = row.order)) %>%
      arrange(order)%>%
      mutate(val = replace(val, val >= lsc[1], lsc[1]))  %>%
      mutate(val = replace(val, val <= lsc[2], lsc[2]))
    
    mypalette<-rev(brewer.pal(11,"RdBu"))
    
    plot<-ggplot(p3, aes(y=order, x=factor(comp, levels = names(file[2:length(names(file))])))) +
      geom_tile(aes(fill= val)) +
      scale_fill_gradientn(colors=mypalette, name= exp.n) +
      theme(axis.text.x = element_text(angle=90))+
      xlab("Stages")+
      ylab("Genes")
    
    return(plot)
  })
  
  pdf(file=here(paste0("sequential-reference/02_dea/figures/tfs-", exp.n, ".pdf")), width = 20, height= 15)
  print(b)
  dev.off()
})











## Eye marker expression plot
eymk <- read.xlsx(here("01_raw-data/eye-markers-aai.xlsx"), 1, colNames = FALSE) %>%
  na.omit() 
names(eymk)<- "gene_id"

eymk.filt <- eymk %>%
  filter(gene_id %in% de.list) 
eymk.filt <- eymk.filt$gene_id
eymk <- eymk$gene_id

limscal<- list(c(4, -4), c(9.017505, 0), c(2, -2))

a<- pmap(list(exp.tab, exp.nam, limscal), function(exp.t, exp.n, lsc){
  #Load expression table
  file<- exp.t
  
  p <- file %>%
    left_join(y= ref, by="gene_id") %>%
    filter(gene_id %in% eymk.filt)
  
  #Matrix for clustering
  m= p
  rownames(m)<- m$id.desc
  m$gene_id<- m$id.desc<- NULL
  m <- as.matrix(m)
  
  #Cluster
  h <- hclust(dist(m), method = 'ward.D2')
  row.order <- rownames(m[h$order,])
  
  p2 <- p %>% pivot_longer(starts_with("s"), names_to = "comp", values_to = "val")
  
  #av <- max(p2$val)
  p3 <- p2 %>% mutate(order = factor(id.desc, levels = row.order)) %>%
    arrange(order) %>%
    mutate(val = replace(val, val >= lsc[1], lsc[1]))  %>%
    mutate(val = replace(val, val <= lsc[2], lsc[2]))
  
  mypalette<-rev(brewer.pal(11,"RdBu"))
  
  plot<-ggplot(p3, aes(y=order, x=factor(comp, levels = names(file[2:length(names(file))])))) +
    geom_tile(aes(fill= val)) +
    scale_fill_gradientn(colors=mypalette, name= exp.n) +
    theme(axis.text.x = element_text(angle=90)) +
    xlab("Stages") +
    ylab("Genes")
  
  pdf(file=here(paste0("sequential-reference/02_dea/figures/eye-mkrs-", exp.n, ".pdf")), width = 20, height= 15)
  print(plot)
  dev.off()

})









#### Eye marker boolean plot
comp <- paste0(c(paste0("s_", c(3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact"), "-s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa")
path<- here("sequential-reference/02_dea/tables/")
files<- paste0("de-up_", comp,"-seq-ref_")

boolean.marker<-map(lfc, function(lf){
  f <- paste0(path, files, lf, ".csv.gz")
  
  de <- map(f, function(f){
    read.csv(f)
  })
  #Comparison-names vector
  names(de) <- files
  
  #Evaluate if marker genes are in DE genes at any time point
  is.mkr <- map(de, function(x){
    list.de<- x$gene_id
    bol <- eymk %in% list.de %>%
      data.frame() 
    rownames(bol) <- eymk
    return(bol)
  })
  bol.mkr <- do.call("cbind", is.mkr)
  colnames(bol.mkr) <- names(de)
  bol.mkr$gene_id <- rownames(bol.mkr)
  
  #Up vsintact
  p <- bol.mkr %>%
    left_join(y= ref, by="gene_id") %>%
    pivot_longer(starts_with("de"), names_to = "comp", values_to = "bol") %>%
    filter(comp %in% files) %>%
    mutate(comp = str_replace(comp, "de-up_s_", ""),
           comp = str_replace(comp, "-seq-ref_", ""),
           comp = str_replace(comp, "s_", "")) %>%
    mutate(comp=factor(comp, levels=order))
  
  #Plot
  p<-ggplot(p,                                # Draw heatmap-like plot
            aes(y=id.desc, x=comp, fill = bol)) +
    geom_tile() +
    theme(axis.text.x = element_text(angle=90)) +
    xlab("Comparisons") +
    ylab("Genes")
  
  pdf(file=here(paste0("sequential-reference/02_dea/figures/eye-mkrs-boolean_", lf, ".pdf")), width = 20, height= 15)
  print(p)
  dev.off()
  
})










#### TF Wound healing, division and eye markers plot
#Load data
#List of genes
wound.div<- read.xlsx(here("01_raw-data/selection-woundhealing-celldivision-aai.xlsx"),  1, colNames = FALSE)
wound.div<- wound.div%>%
  rename(gene_id=X1)%>%
  filter(gene_id %in% de.list) 
wound.div<- wound.div$gene_id

tf<- read.xlsx(here("01_raw-data/selection-markers-tfs-aai.xlsx"), 1, colNames = FALSE)
tf<- tf%>%
  rename(gene_id=X1)%>%
  filter(gene_id %in% de.list) 
tf<- tf$gene_id

#Join zscore table with reference
zsc <- tpm.z %>%
  left_join(y=ref, by= "gene_id") %>%
  na.omit()

genes<- unique(c(wound.div, tf))

#Subset data for plotting 
d.p<- zsc %>%
  filter(gene_id %in% genes)

##Matrix for clustering
m<- d.p
rownames(m)<- m$id.desc
m$gene_id<- m$id.desc<- NULL
m <- as.matrix(m)
##Cluster
h <- hclust(dist(m), method = 'ward.D2')
row.order <- rownames(m[h$order,])

##Pivot dataframes
d.p2 <- d.p %>% pivot_longer(starts_with("s"), names_to = "st", values_to = "zscore")

##Arrange data frames
d.p3 <- d.p2 %>% 
  mutate(order = factor(id.desc, levels = row.order)) %>%
  arrange(order) 

##Plots 
mypalette<-rev(brewer.pal(11,"RdBu"))
###Plot for regeneration
lev<- c(paste0("s_", c(1, 3, 6, 9, 12, 15, 21, 28), "dpa"), "s_intact")

plot <- ggplot(d.p3, aes(y=order, x=factor(st, levels = lev))) +
  geom_tile(aes(fill= zscore)) +
  scale_fill_gradientn(colors=mypalette) +
  theme(axis.text.x = element_text(angle=90)) +
  ylab("Genes") +
  xlab("Stages")

ggsave(filename=here("sequential-reference/02_dea/figures/tf-woundhealing-div.pdf"), width=10, height= 20)

