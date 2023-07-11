# Create tpm tables
# Brenda Pardo
# 2023-05-15

#Libraries
library(dplyr)
library(matrixStats)
library(readr)


#Load tpms
tpm=read_csv(here("01_raw-data/RSEM_TPM_table.csv.gz"))

#Arrange data frame 
tpm<- tpm %>%
  rename(gene_id=Gene_ID) %>%
  select(!Name)

#Time point vector
st<- c(paste0("s_", c(1,3,6,9,12,15,21,28), "dpa"), "s_intact")


##Mean tpm table
tpm.s= tpm %>%
  mutate(s_1dpa=rowMeans(tpm[, paste0("s_1dpa_", 1:4)]),
         s_3dpa=rowMeans(tpm[, paste0("s_3dpa_", 1:4)]),
         s_6dpa=rowMeans(tpm[, paste0("s_6dpa_", 1:4)]),
         s_9dpa=rowMeans(tpm[, paste0("s_9dpa_", 1:4)]),
         s_12dpa=rowMeans(tpm[, paste0("s_12dpa_", 1:4)]),
         s_15dpa=rowMeans(tpm[, paste0("s_15dpa_", 1:4)]),
         s_21dpa=rowMeans(tpm[, paste0("s_21dpa_", 1:4)]),
         s_28dpa=rowMeans(tpm[, paste0("s_28dpa_", 1:4)]),
         s_intact=rowMeans(tpm[, paste0("s_intact_", 1:4)])
  ) %>%
  dplyr::select(gene_id, ends_with(st))

#Filter by low expression
tpm.max<- tpm.s %>%
  mutate(max.tpm = rowMaxs(as.matrix(tpm.s[,st])))

tpm.f<- tpm.max %>% 
  filter (max.tpm > 0.1) %>%
  select(!max.tpm)

write_csv(tpm.f, here("02_processed-data/tpm-mean.csv.gz"))





##Z-score(tpm) table

z.tpm<- tpm.f
rownames(z.tpm)= z.tpm$gene_id
z.tpm$gene_id <- NULL

#Calculate z scores
z.tpm= t(z.tpm)
z.tpm=t(scale(z.tpm))
z.tpm=na.omit(z.tpm)

#Build table
z.tpm= data.frame(z.tpm)
z.tpm$gene_id=rownames(z.tpm)
z.tpm<- z.tpm %>% 
  dplyr::select(gene_id, all_of(st))

write_csv(z.tpm, here("02_processed-data/tpm-zscores.csv.gz"))





#Log1p (tpm) table
tpm.log<- tpm.f
rownames(tpm.log)<- tpm.log$gene_id
tpm.log$gene_id= NULL
tpm.log= log1p(tpm.log)
tpm.log$gene_id= rownames(tpm.log)

tpm.log <- tpm.log %>%
  dplyr::select(gene_id, s_1dpa, s_3dpa, s_6dpa, s_9dpa, s_12dpa, s_15dpa, 
                s_21dpa, s_28dpa, s_intact)

#Save table with collapsed replicates
write_csv(tpm.log, here("02_processed-data/tpm-mean-log1p.csv.gz"))
