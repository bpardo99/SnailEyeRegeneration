# Create tpm tables for development data set
# Brenda Pardo
# 2023-05-15
library(tidyverse)
library(here)
#Get tpms
d<-read_tsv(here("01_raw-data/development_cpms.tsv.gz"),  header=TRUE)

#Arrange table
st<- c(paste0("st", 2:10), "h")
cpm<- d %>%
  dplyr::rename(st2= X2d, st3=X3d, st4=X4d, st5=X5d, st6=X6d, st7=X7dE, 
                st8=X11d, st9=X13d, st10=X16d, h=X19d) %>%
  select("target_id", dplyr::starts_with(st)) %>%
  dplyr::rename(gene_id= target_id)

#Filter by low expression
cpm.max<- cpm %>%
  mutate(max.cpm = rowMaxs(as.matrix(cpm[,st])))

cpm.f<- cpm.max %>% 
  filter (max.cpm > 0.1) %>%
  select(!max.cpm)

rownames(cpm.f)= cpm.f$gene_id
cpm.f$gene_id <- NULL

#Calculate z-scores
z.cpm= t(cpm.f)
z.cpm=t(scale(z.cpm))
z.cpm=na.omit(z.cpm)

z.cpm= data.frame(z.cpm)
z.cpm$gene_id=rownames(z.cpm)

z.cpm<- z.cpm %>% 
  dplyr::select(gene_id, dplyr::starts_with(st))
#summary(rowSums(z.cpm[, st]))

#Write table
write_csv(z.cpm, here("02_processed-data/dev-cpm-zscores.csv.gz"))
