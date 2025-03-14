library(tidyverse)
library(openxlsx)

#### Table S1: GSEA ####
attach("results/fgsea.RData")
H_all_human <- hum_fgsea %>% 
  select(group, pathway, NES, pval, FDR, leadingEdge)
H_all_mouse <- mur_fgsea %>% 
  select(group, pathway, NES, pval, FDR, leadingEdge)
attach("results/fgsea_AMvMDM.RData")
H_all_human2 <- hum_fgsea2 %>% 
  select(group, pathway, NES, pval, FDR, leadingEdge)
H_all_mouse2 <- mur_fgsea2 %>% 
  select(group, pathway, NES, pval, FDR, leadingEdge)

gsea.ls <- list()

gsea.ls[["Human AM Mtb-Media"]] <- H_all_human %>% 
  filter(group == "human AM") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in human AM", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Human MDM Mtb-Media"]] <- H_all_human %>% 
  filter(group == "human MDM") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in human MDM", .before = 1) %>% 
  arrange(FDR)
###
gsea.ls[["Murine AM coMtb Mtb-Media"]] <- H_all_mouse %>% 
  filter(group == "murine AM coTB") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in murine AM coMtb", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Murine AM scBCG Mtb-Media"]] <- H_all_mouse %>% 
  filter(group == "murine AM scBCG") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in murine AM scBCG", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Murine AM control Mtb-Media"]] <- H_all_mouse %>% 
  filter(group == "murine AM") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in murine AM control", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Murine BMDM Mtb-Media"]] <- H_all_mouse %>% 
  filter(group == "murine BMDM") %>% 
  select(-group) %>% 
  mutate(foldChange = "Mtb vs media in murine BMDM", .before = 1) %>% 
  arrange(FDR)
###
gsea.ls[["Human Media MDM-AM"]] <- H_all_human2 %>% 
  filter(group == "human media") %>% 
  select(-group) %>% 
  mutate(foldChange = "MDM vs AM in human media", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Human Mtb MDM-AM"]] <- H_all_human2 %>% 
  filter(group == "human mtb") %>% 
  select(-group) %>% 
  mutate(foldChange = "MDM vs AM in human Mtb", .before = 1) %>% 
  arrange(FDR)
###
gsea.ls[["Murine Media BMDM-AMcontrol"]] <- H_all_mouse2 %>% 
  filter(group == "murine control media") %>% 
  select(-group) %>% 
  mutate(foldChange = "BMDM vs AM control in murine media", .before = 1) %>% 
  arrange(FDR)

gsea.ls[["Murine Mtb BMDM-AMcontrol"]] <- H_all_mouse2 %>% 
  filter(group == "murine control mtb") %>% 
  select(-group) %>% 
  mutate(foldChange = "BMDM vs AM control in murine Mtb", .before = 1) %>% 
  arrange(FDR)

write.xlsx(gsea.ls, file="publication/TableS2_gsea.xlsx")

#### Table S2: Linear models ####
load("results/lme_human.RData")
load("results/lm_mouse.RData")

lm.ls <- list()

lm.ls[["Human AM Mtb-Media"]] <- model_hum$lme.contrast %>% 
  filter(contrast_ref== "AM Media" & contrast_lvl=="AM Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in human AM", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Human MDM Mtb-Media"]] <- model_hum$lme.contrast %>% 
  filter(contrast_ref== "MDM Media" & contrast_lvl=="MDM Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in human MDM", .before = 1) %>% 
  arrange(FDR)
###
lm.ls[["Murine AM coMtb Mtb-Media"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "AM_coTB Media" & contrast_lvl=="AM_coTB Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in murine AM coMtb", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Murine AM scBCG Mtb-Media"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "AM_scBCG Media" & contrast_lvl=="AM_scBCG Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in murine AM scBCG", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Murine AM control Mtb-Media"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "AM_control Media" & contrast_lvl=="AM_control Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in murine AM control", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Murine BMDM Mtb-Media"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "BMDM Media" & contrast_lvl=="BMDM Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "Mtb vs media in murine BMDM", .before = 1) %>% 
  arrange(FDR)
###
lm.ls[["Human Media MDM-AM"]] <- model_hum$lme.contrast %>% 
  filter(contrast_ref== "AM Media" & contrast_lvl=="MDM Media") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "MDM vs AM in human media", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Human Mtb MDM-AM"]] <- model_hum$lme.contrast %>% 
  filter(contrast_ref== "AM Mtb" & contrast_lvl=="MDM Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "MDM vs AM in human Mtb", .before = 1) %>% 
  arrange(FDR)
###
lm.ls[["Murine Media BMDM-AMcontrol"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "AM_control Media" & contrast_lvl=="BMDM Media") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "BMDM vs AM control in murine media", .before = 1) %>% 
  arrange(FDR)

lm.ls[["Murine Mtb BMDM-AMcontrol"]] <- model_mur$lm.contrast %>% 
  filter(contrast_ref== "AM_control Mtb" & contrast_lvl=="BMDM Mtb") %>% 
  select(gene, variable, estimate, pval, FDR) %>% 
  mutate(foldChange = "BMDM vs AM control in murine Mtb", .before = 1) %>% 
  arrange(FDR)

write.xlsx(lm.ls, file="publication/TableS1_linear_models.xlsx")
