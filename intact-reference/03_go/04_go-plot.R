# GO enrichment plot using intact time-point as the reference. Terms were manually
# selected based on GO bar plots (home/bp2582/projects/eye-reg_rnaseq
# /github/SnailEyeReg_RNASeq-pipeline/intact-reference/03_go/output). Made by 
# Alice Accorsi
# Code written by ejr
# Adapted by Brenda Pardo
# 2023-05-15

library(tidyverse)
library(cowplot)
library(RColorBrewer)
library(here)

# function to capitalize the first letter of string
CapStr <- function(y) {
  c <- strsplit(y, " ")[[1]]
  paste(toupper(substring(c, 1,1)), substring(c, 2),
        sep="", collapse=" ")
}

#logFC thresholds
thresholds <- c(0, 1, 2)

#Term list
both_list <- read_tsv(here("01_raw-data/selection-GO-73-aai.txt.gz"))

plot_order <- c(
  "1dpa-intact",
  "3dpa-intact",
  "6dpa-intact",
  "9dpa-intact",
  "12dpa-intact",
  "15dpa-intact",
  "21dpa-intact",
  "28dpa-intact"
)

map(thresholds, function(x) {
  # read in GO enrichment
  go <- read_tsv(here(paste0("intact-reference/03_go/output/go-breakdown_lfc", x, ".txt.gz"))) %>%
    mutate(file=cluster) %>%
    filter(str_detect(file, "de-up")) %>%
    dplyr::rename(go_id = ID) %>%
    mutate(file = str_replace(file, "de-", "")) %>%
    mutate(file = str_replace(file, "-go_s", "")) %>%
    mutate(file = str_replace(file, "-s_", "-")) %>%
    mutate(file = str_replace(file, ".txt", "")) %>%
    mutate(file = str_replace(file, "-intact-ref_lfc", "_lfc_")) %>%
    separate(file, c("direction", "timepoint","none", "lfc"), sep="_")
  
  # Create UP regulated figure
  label_levels <-   go %>%
    left_join(both_list) %>%
    filter(! is.na(order)) %>%
    mutate(label = paste0(go_id, ": ", Description)) %>%
    rowwise() %>%
    mutate(label2 = CapStr(Description))%>%
    filter(direction == "up") %>%
    select(label, label2, order) %>%
    distinct() %>%
    arrange(desc(order))
  
  go_fig <- go %>%
    left_join(both_list) %>%
    filter(! is.na(order)) %>%
    mutate(label = paste0(go_id, ": ", Description)) %>%
    mutate(label = factor(label, levels=label_levels$label)) %>%
    rowwise() %>%
    mutate(label2 = CapStr(Description)) %>%
    mutate(label2 = factor(label2, levels=label_levels$label2)) %>%
    filter(direction == "up") %>%
    mutate(`-log10(pval)` = log10(p.adjust) * -1) %>%
    rename(`# genes` = Count)
  
  mypalette<-rev(brewer.pal(11,"RdBu"))
  lim <- c(max(abs(go_fig$`-log10(pval)`)) * -1, max(abs(go_fig$`-log10(pval)`)))
  
  p2 <- ggplot(data = go_fig, aes(x = factor(timepoint, levels=plot_order), y = label, size = `# genes`, )) +
    theme_cowplot() +
    geom_hline(aes(yintercept = label), alpha = .1) +
    #    scale_colour_gradientn(colors=mypalette, limits=lim) +
    #scale_colour_distiller(palette = "YlGnBu") +
    #scale_colour_gradientn(colours = c("#3753a5", "#6e499e", "#bd1b8d", "#e01662", "#ec1b26"))+
    geom_point() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1)) +
    scale_x_discrete(limits = plot_order)
  
  ggsave(here(paste0("intact-reference/03_go/figures/select_GO_id_lfc", x, "-intact-ref.pdf")),p2, height=16, width=14)
  
  
  p3 <- ggplot(data = go_fig, aes(x = factor(timepoint, levels=plot_order), y = label2, size = `# genes`)) +
    theme_cowplot() +
    geom_hline(aes(yintercept = label), alpha = .1) +
    #    scale_colour_gradientn(colors=mypalette, limits=lim) +
    #scale_colour_distiller(palette = "YlGnBu") +
    #scale_colour_gradientn(colours = c("#3753a5", "#6e499e", "#bd1b8d", "#e01662", "#ec1b26"))+
    geom_point() +
    xlab("") +
    ylab("") +
    theme(axis.text.x = element_text(angle = 45, vjust = 1.1, hjust=1.1))+
    scale_x_discrete(limits = plot_order)
  
  ggsave(here(paste0("intact-reference/03_go/figures/select_GO_noid_lfc", x, "-intact-ref.pdf")),p3, height=16, width=14)
  
})

