library(tidyverse)
library(SEARchways)

#### Data ####
load("results/lm_mouse.RData")
load("results/lme_human.RData")

#### FGSEA Mtb vs media - human ####
hum_key <- data.frame(
  contrast_ref = c("AM Media","MDM Media"),
  contrast_lvl = c("AM Mtb","MDM Mtb")
)

hum_fgsea <- model_hum$lme.contrast %>%  
  inner_join(hum_key) %>% 
  mutate(group = case_when(contrast_ref=="AM Media"~"human AM",
                           contrast_ref=="MDM Media"~"human MDM")) %>% 
  select(group, gene, estimate) %>% 
  BIGsea(gene_df = .,
         category = "H",
         species = "human", ID="SYMBOL",
         nperm = 1E5)

#### FGSEA Mtb vs media - Murine ####
mur_key <- data.frame(
  contrast_ref = c("AM_control Media","AM_coTB Media","AM_scBCG Media",
                   "BMDM Media"),
  contrast_lvl = c("AM_control Mtb","AM_coTB Mtb","AM_scBCG Mtb",
                   "BMDM Mtb")
)

mur_fgsea <- model_mur$lm.contrast %>% 
  inner_join(mur_key) %>% 
  mutate(group = case_when(contrast_ref=="AM_control Media"~"murine AM",
                           contrast_ref=="AM_coTB Media"~"murine AM coTB",
                           contrast_ref=="AM_scBCG Media"~"murine AM scBCG",
                           contrast_ref=="BMDM Media"~"murine BMDM")) %>% 
  select(group, gene, estimate) %>% 
  BIGsea(gene_df = .,
         category = "H",
         species = "mouse", ID="SYMBOL",
         nperm = 1E5)

#### Save ####
save(hum_fgsea, mur_fgsea, file = "results/fgsea.RData")

#### FGSEA AM vs MDM - human ####
hum_key <- data.frame(
  contrast_ref = c("AM Media","AM Mtb"),
  contrast_lvl = c("MDM Media","MDM Mtb")
)

hum_fgsea2 <- model_hum$lme.contrast %>%  
  inner_join(hum_key) %>% 
  mutate(group = case_when(contrast_ref=="AM Media"~"human media",
                           contrast_ref=="AM Mtb"~"human mtb")) %>% 
  select(group, gene, estimate) %>% 
  BIGsea(gene_df = .,
         category = "H",
         species = "human", ID="SYMBOL",
         nperm = 1E5)

#### FGSEA AM vs MDM - murine ####
mur_key <- data.frame(
  contrast_ref = c("AM_control Media","AM_coTB Media","AM_scBCG Media",
                   "AM_control Mtb","AM_coTB Mtb","AM_scBCG Mtb"),
  contrast_lvl = c("BMDM Media", "BMDM Media", "BMDM Media",
                   "BMDM Mtb", "BMDM Mtb", "BMDM Mtb")
)

mur_fgsea2 <- model_mur$lm.contrast %>% 
  inner_join(mur_key) %>% 
  mutate(group = case_when(contrast_ref=="AM_control Media"~"murine control media",
                           contrast_ref=="AM_coTB Media"~"murine coTB media",
                           contrast_ref=="AM_scBCG Media"~"murine scBCG media",
                           contrast_ref=="AM_control Mtb"~"murine control mtb",
                           contrast_ref=="AM_coTB Mtb"~"murine coTB mtb",
                           contrast_ref=="AM_scBCG Mtb"~"murine scBCG mtb"
                           )) %>% 
  select(group, gene, estimate) %>% 
  BIGsea(gene_df = .,
         category = "H",
         species = "mouse", ID="SYMBOL",
         nperm = 1E5)

#### Save ####
save(hum_fgsea2, mur_fgsea2, file = "results/fgsea_AMvMDM.RData")
