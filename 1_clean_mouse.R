library(tidyverse)
library(readxl)
library(edgeR)
library(limma)
# library(BIGpicture)

### AM ####
#### Raw data ####
#Counts
count_am <- read_excel("data_raw/murine/AM_H37Rv_RNAseq.xlsx",
                       sheet="Raw counts")
count_bmdm <- readRDS("data_raw/murine/BMDM_allRawCounts.rds") %>% 
  select(gene, contains("gso30.Bl6.N.N.4"), contains("gso30.Bl6.N.Rv.4"))

count <- full_join(count_am, count_bmdm)
count[is.na(count)] <- 0

#create metadata
meta_am <- data.frame(libID=colnames(count_am)[-1]) %>% 
  separate(libID, into=c("treatment","mtb","ptID"), sep="[.]", remove = FALSE) %>% 
  mutate(mtb = recode(mtb, "none"="Media", "mtb"="Mtb")) %>% 
  mutate(treatment = recode(treatment, 
                            "cotb"="AM_coTB",
                            "unvax"="AM_control",
                            "scBCG"="AM_scBCG")) %>% 
  arrange(libID)

meta_bmdm <- data.frame(libID=colnames(count_bmdm)[-1]) %>% 
  separate(libID, into=c("study","breed","trash1",
                         "mtb","trash2","ptID"), 
           sep="[.]", remove=FALSE) %>% 
  select(-trash1, -trash2) %>% 
  mutate(mtb = fct_recode(mtb,"Media"="N","Mtb"="Rv"),
         treatment="BMDM")

meta <- full_join(meta_am, meta_bmdm)

#### Protein-coding ####
ensembl_mouse <- biomaRt::useEnsembl(biomart="ensembl", 
                                     dataset="mmusculus_gene_ensembl",
                                     mirror = "useast")

## Format gene key
key_mouse <- biomaRt::getBM(attributes=c("mgi_symbol", "gene_biotype"), 
                            mart=ensembl_mouse) %>% 
  #Filter protein coding genes
  filter(gene_biotype == "protein_coding") %>% 
  rename(symbol=mgi_symbol) %>% 
  arrange(symbol) %>% 
  filter(symbol %in% count$gene)

count_pc <- count %>% 
  filter(gene %in% key_mouse$symbol) %>% 
  column_to_rownames("gene") %>% 
  as.matrix()

#### DGEList ####
dat_am <- DGEList(counts=count_pc, samples=meta, genes=key_mouse)

#### Low abundance ####
# Remove only 0s
dat_abund <- RNAetc::filter_rare(dat_am,
                                 min_CPM = 0.1, min.sample = 2,
                                 gene.var="symbol")
# 37 (0.29%) of 12911 genes removed.

#### Normalize ####
dat_abund_norm <- calcNormFactors(dat_abund, method = "TMM")

dat_abund_norm_voom <- voomWithQualityWeights(
  dat_abund_norm, plot = TRUE,
  design = model.matrix(~ treatment*mtb, 
                        data = dat_abund_norm$samples))

voom_mur <- dat_abund_norm_voom

#### Save ####
save(voom_mur,
     file = "data_clean/mouse_dat_clean.RData")
