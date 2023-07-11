# Run DE analysis with edgeR using sequential time-points as the reference.
# Brenda Pardo
# 2023-05-15

#Libraries
library(edgeR)
library(dplyr)
library(ggplot2)
library(tibble)
library(readr)
library(cowplot)
library(corrplot)
library(stringr)
library(ggrepel)
library(RColorBrewer)
library(purrr)
library(here)


##Generate DGE object as an input for edgeR

#Read raw counts
counts<- read.csv(here("01_raw-data/star_count.csv.gz"), row.names="gene_id")

#Reorder tibble
col_order<- c(paste0("s_1dpa_", 1:4), 
              paste0("s_3dpa_", 1:4), paste0("s_6dpa_", 1:4), 
              paste0("s_9dpa_", 1:4), paste0("s_12dpa_", 1:4), 
              paste0("s_15dpa_", 1:4), paste0("s_21dpa_", 1:4), 
              paste0("s_28dpa_", 1:4), paste0("s_intact_", 1:4))
counts<- counts[, col_order]

#Set groups
groups<- c("s_1dpa", "s_3dpa", "s_6dpa", "s_9dpa", "s_12dpa", "s_15dpa", "s_21dpa",
           "s_28dpa", "s_intact")
group<- factor(rep(groups, each=4), levels = groups)

#Create edgeR object
dge <- DGEList(counts=counts, group=group)





##Filtering by low expression

keep <- filterByExpr(dge) 
dge <- dge[keep, , keep.lib.sizes=FALSE]





##Boxplot counts distrib per sample, PRE-NORMALIZATION

#pdf("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/sequential-reference/02_dea/figures/boxplot-counts-prenorm.pdf")
pdf(here("sequential-reference/02_dea/figures/boxplot-counts-prenorm.pdf"))
boxplot(log2(dge$counts +1), las=2, main="Log2 read counts per sample, prenormalization")
dev.off()





##Normalization
#Calculate normalization factors
dge <- calcNormFactors(dge)

#Calculate gene dispersion
dge <- estimateDisp(dge)





##Boxplot counts distrib per sample, POST-NORMALIZATION
dge.norm<- dge$counts*dge$samples$norm.factors
#pdf("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/sequential-reference/02_dea/figures/boxplot-counts-postnorm.pdf")
pdf(here("sequential-reference/02_dea/figures/boxplot-counts-postnorm.pdf"))
boxplot(log2(dge.norm +1), las=2, main="Log2 read counts per sample, post-normalization")
dev.off()





##Read REFSEQ references table

#ref<- read.table("/home/bp2582/projects/eye-reg_rnaseq/github/SnailEyeReg_RNASeq-pipeline/01_raw-data/gene-ref.txt", sep="\t", quote = "", header= TRUE)
ref<- read.table(here("01_raw-data/gene-ref.txt.gz"), sep="\t", quote = "", header= TRUE)





##Differential Expression Analysis

#Comparison list
comp <- list(
  c("s_1dpa", "s_3dpa"),
  c("s_3dpa", "s_6dpa"),
  c("s_6dpa", "s_9dpa"),
  c("s_9dpa", "s_12dpa"),
  c("s_12dpa", "s_15dpa"),
  c("s_15dpa", "s_21dpa"),  
  c("s_21dpa", "s_28dpa"),
  c("s_28dpa", "s_intact")
)

#logFC list
lfc_list <- c(0:2)


#DEA
df.lfc <- map(lfc_list, function(lfc){
  df <- map(comp, function(x){
  
    #Exact test
    tag <- paste0(x[2],"-", x[1], "-seq-ref")
    et <- exactTest(dge, pair= x, dispersion="tagwise")
    et <- topTags(et, n = Inf) %>%
      data.frame() %>%
      rownames_to_column(var = "gene_id")
    
    #Write table with test values
    write_csv(et, 
              file=here("sequential-reference", "02_dea", "tables",paste0("eTest_",tag, ".csv.gz")))
    
    #UList of upregulated genes
    de.up <- et[which(et$FDR<= 1e-5 & et$logFC>= lfc), ]
    write_csv(de.up, 
              file=here("sequential-reference", "02_dea", "tables",paste0("de-up_",tag, "_lfc", lfc, ".csv.gz")))
    
    #List of downregulated genes
    de.down <- et[which(et$FDR<= 1e-5 & et$logFC<= -lfc), ]
    write_csv(de.down, 
              file=here("sequential-reference", "02_dea", "tables", paste0("de-down_",tag, "_lfc", lfc, ".csv.gz")))  
    
    #List of upregulated and downregulated genes
    de.ud <- et[which((et$FDR<= 1e-5 & et$FDR>= lfc) | 
                        (et$FDR<= 1e-5 & et$FDR<= -lfc)), ]    
    write_csv(de.ud, 
              file=here("sequential-reference", "02_dea", "tables", paste0("de-ud_",tag, "_lfc", lfc, ".csv.gz")) )
    
    
    
    
    
    #List of upregulated genes with rank and gene description
    de.up.desc= de.up %>%
      left_join(y= ref, by="gene_id") %>%
      mutate(rank = dense_rank(plyr::desc(logFC))) %>%
      arrange(plyr::desc(rank))

    #Write upregulated genes to a csv file
    write_csv(de.up.desc, 
              file=here("sequential-reference", "02_dea", "tables", paste0("de-up_",tag, "_lfc", lfc, "-desc.csv.gz")))
    
    #List of downregulated genes
    de.down.desc= de.down %>%
      left_join(y= ref, by="gene_id") %>%
      mutate(rank = dense_rank(logFC)) %>%
      arrange(plyr::desc(rank))
    #Write upregulated genes to a csv file
    write_csv(de.down.desc, 
              file=here("sequential-reference", "02_dea", "tables", paste0("de-down_",tag, "_lfc", lfc, "-desc.csv.gz")))

    
    
    
    
    ##DE genes plots
    
    #Volcano plot
    vm.p<- et
    vm.p$diffexpressed <- "NO"
    vm.p$diffexpressed[vm.p$logFC >= lfc & vm.p$FDR<= 1e-5] <- "UP"
    vm.p$diffexpressed[vm.p$logFC <= -lfc & vm.p$FDR<= 1e-5] <- "DOWN"
    
    volcano<- vm.p%>%
      ggplot(aes(x=logFC+ 1e-10000, y=-log10(PValue + 1e-10000), col=diffexpressed)) +
      geom_point(aes(alpha=0.05)) + 
      theme_cowplot() +
      scale_color_manual(values=c("blue", "black", "red")) +
      ylab("-log10(p-val)") +
      xlab("logFC") +
      ggtitle(paste0("Volcano plot. ", tag)) +
      geom_vline(xintercept=c(-lfc, lfc), col="grey")+
      theme(legend.position="none") 
    ggsave(filename= here(paste0("sequential-reference/02_dea/figures/volcano_", tag, "_lfc", lfc, ".pdf")), plot=volcano, width= 5, height=5)
  
    #MA plot
    ma <- vm.p %>%
      ggplot(aes(y=logFC, x=logCPM, col=diffexpressed)) +
      geom_point(aes(alpha=0.05)) +
      theme_cowplot() +
      ylab("logFC") +
      xlab("Average logCPM") +
      ggtitle(paste0("MA plot. ", tag)) +
      scale_color_manual(values=c("blue", "black", "red")) +
      geom_hline(yintercept=c(-lfc, lfc), col="grey")+
      theme(legend.position="none") 
    ggsave(filename=here( paste0("sequential-reference/02_dea/figures/maplot_", tag, "_lfc", lfc, ".pdf")), plot=ma, width= 5, height=5)
  
    
    
    
    
    #Count number of DE genes SAVE NUMBERS
    
    de.ct.up <- dim(de.up)[1]
    de.ct.down <- dim(de.down)[1]
    return(c(tag, de.ct.up, de.ct.down))
    
  })
  
  
  
  
  
  #Count number of DE genes BUILD TABLE
  
  de.ct <- do.call("rbind", df) %>%
    data.frame() 
  names(de.ct)=c("tag", "de.ct.up", "de.ct.down")
  write_csv(de.ct, 
            file=here(paste0("sequential-reference/02_dea/tables/number-de-genes-lfc", lfc,"_seq-ref.csv.gz")))
})





##Spearman correlation plot

#Select top 1000 genes
dge.top= as_tibble(dge$counts)
dge.top= dge.top %>%
  mutate(Mean= rowMeans(dge.top)) %>%
  arrange(dplyr::desc(Mean)) %>%
  top_n(n=1000)

#Calculate correlations
dge.top$Mean <- NULL
cor.sp= cor(dge.top, y=NULL, method= "spearman")

myPalette <- colorRampPalette(rev(c("#1f1248",  "#D25868", "#F2A979", "#FDFDC6")))
pdf(file=here("sequential-reference/02_dea/figures/correlation-spea-top1000.pdf"), width=20, height=20)
cor.spp<- corrplot(cor.sp, 
                   method= "color", 
                   main="Spearman correlation", 
                   col.lim = c(min(cor.sp), max(cor.sp)), 
                   is.corr = FALSE,
                   col=myPalette(100), 
                   tl.col="black",
                   tl.cex=2)
dev.off()





##Dimension reduction plot, MDS plot

#With all genes
mds <- plotMDS(dge,  plot=FALSE, top=1000)

#Order data
mds_df <- data.frame(replicate = rownames(mds$distance.matrix), x=mds$x, y=mds$y) %>%
  mutate(timepoint = str_replace(replicate, "s_", ""))%>%
  mutate(timepoint = str_replace(timepoint, "_.+$", "")) %>%
  mutate(rep1 = ifelse(str_detect(replicate, "_1$"), timepoint, NA)) 

#Plot
plot_order <- c("1dpa", "3dpa", "6dpa", "9dpa", "12dpa", "15dpa", "21dpa",
                "28dpa", "intact")
myPalette <- brewer.pal(9, "Paired")
names(myPalette) <- plot_order

#Plot with labels
mds_plot <- ggplot(mds_df, aes(x=x, y=y, color=factor(timepoint, level=plot_order), label=rep1)) +
  geom_point() +
  theme_cowplot() +
  geom_text_repel(aes(x = x, 
                      y = y, 
                      label = rep1), nudge_x=.5, segment.colour = NA) +
  theme(legend.position = "none")+
  scale_color_manual(values = myPalette)  +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  ggtitle(paste0("MDS PLOT: top 1000 genes")) +
  xlim(NA, 2.5)
ggsave(filename = here("sequential-reference/02_dea/figures/mds-eric-top1000.pdf"), width=6, height=5)

#Plot NO labels
mds_plot <- ggplot(mds_df, aes(x=x, y=y, color=factor(timepoint, level=plot_order))) +
  geom_point() +
  theme_cowplot() +
  theme(legend.position = "none")+
  scale_color_manual(values = myPalette)  +
  xlab("Dimension 1") +
  ylab("Dimension 2") +
  ggtitle(paste0("MDS PLOT: top 1000 genes")) +
  xlim(NA, 2.5)
ggsave(filename = here("sequential-reference/02_dea/figures/mds-eric-nlab-top1000.pdf"), width=6, height=5)





#List of all upregulated genes accross all comparisons
df.lfc <- map(lfc_list, function(lfc){
  df <- map(comp, function(x){
    tag <- paste0(x[2],"-", x[1], "-seq-ref")
    tab<- read.csv(here("sequential-reference/02_dea/tables/", paste0("de-up_",tag, "_lfc", lfc, ".csv.gz")))
    return(tab)
    
  })
  tab <- do.call("rbind", df) %>%
    data.frame() %>%
    select(gene_id) %>%
    unique()
  
  write_csv(tab, 
            file=here(paste0("sequential-reference/02_dea/tables/de-up_lfc", lfc, ".csv.gz")))
  
})
  

  
