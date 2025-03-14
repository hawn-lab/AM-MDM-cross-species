library(tidyverse)
library(edgeR)
library(limma)

### AM ####
#### Raw data ####
#Counts
count <- read_tsv("data_raw/human/AM-MDM/AM.MDM.featurecounts.paired.tsv",
                     skip=1) %>% 
  select(-c(Chr:Length)) %>% 
  rename_all(~gsub("project1/data/bam_filter_paired/|_filter_paired.bam","",.)) %>% 
  rename(ensembl_gene_id=Geneid)

#create metadata
meta <- data.frame(libID=colnames(count)[-1]) %>% 
  separate(libID, into=c("ptID","cell","mtb"), sep="_", remove = FALSE) %>% 
  mutate(mtb = recode(mtb, "none"="Media", "TB"="Mtb")) %>% 
  arrange(libID)

#### Protein-coding ####
ensembl_human <- biomaRt::useEnsembl(biomart="ensembl", 
                                     dataset="hsapiens_gene_ensembl",
                                     mirror = "useast")

## Format gene key
key_human <- biomaRt::getBM(attributes=c("hgnc_symbol", "ensembl_gene_id",
                                         "gene_biotype"), 
                            mart=ensembl_human) %>% 
  #Filter protein coding genes
  filter(gene_biotype == "protein_coding") %>% 
  rename(symbol=hgnc_symbol) %>% 
  arrange(symbol)

count_pc <- key_human %>% 
  inner_join(count, by="ensembl_gene_id")

dups <- count_pc %>% 
  filter(symbol != "") %>% 
  count(symbol) %>% 
  filter(n>1) %>% pull(symbol)

#also rename to symbol and deals with duplicate symbls
count_pc_rename <- count_pc %>% 
  mutate(symbol = case_when(symbol %in% dups ~ paste(symbol, ensembl_gene_id, sep="_"),
                            symbol=="" | is.na(symbol) ~ ensembl_gene_id,
                            TRUE~symbol)) %>% 
  select(-gene_biotype,-ensembl_gene_id)

key_human_rename <- key_human %>% 
  mutate(symbol = case_when(symbol %in% dups ~ paste(symbol, ensembl_gene_id, sep="_"),
                            symbol=="" | is.na(symbol) ~ ensembl_gene_id,
                            TRUE~symbol))

#### DGEList ####
count_pc_m <- count_pc_rename %>% 
  arrange(symbol) %>% 
  column_to_rownames("symbol") %>% 
  select(all_of(meta$libID)) %>% 
  as.matrix()
key <- key_human_rename %>% 
  filter(ensembl_gene_id %in% count_pc$ensembl_gene_id) %>% 
  filter(symbol %in% count_pc_rename$symbol) 

dat <- DGEList(counts=count_pc_m, samples=meta, genes=key)

#### Low abundance ####
#separate for am and mdm
dat_abund <- RNAetc::filter_rare(dat,
                                    min_CPM = 0.1, min.sample = 2,
                                    gene.var="symbol")
#3720 (18.67%) of 19921 genes removed.

#### Normalize ####
dat_abund_norm <- calcNormFactors(dat_abund, method = "TMM")

dat_abund_norm_voom <- voomWithQualityWeights(
  dat_abund_norm, plot = TRUE,
  design = model.matrix(~ cell*mtb, 
                        data = dat_abund_norm$samples))

voom_hum <- dat_abund_norm_voom

#### Save ####
save(voom_hum, 
     file = "data_clean/human_dat_clean.RData")
