library(tidyverse)
library(viridis)
library(limma)
library(ComplexHeatmap)
library(patchwork)
signif.pw <- 0.01
signif.gene <- 0.1
source("publication/plot_fxns.R")

#### Data ####
load("results/fgsea.RData")

H_all <- bind_rows(hum_fgsea, mur_fgsea) %>% 
  mutate(group2 = recode(group,
                         "human AM"="Human\nAM",
                         "human MDM"="Human\nMDM",
                         "murine AM scBCG"="scBCG\nMurine AM",
                         "murine AM coTB"="coMtb\nMurine AM",
                         "murine BMDM"="Murine\nBMDM",
                         "murine AM"="control\nMurine AM"))%>% 
  mutate(fdr.group = case_when(FDR < 0.01 ~ "FDR < 0.01",
                               FDR < 0.1 ~ "FDR < 0.1",
                               TRUE~"NS"),
         fdr.group = factor(fdr.group, levels=rev(c("FDR < 0.01",
                                                    "FDR < 0.1",
                                                    "NS"))))

#### Define gene lists ####
# Signif in both
MDM_signif <- H_all %>% 
  filter(grepl("MDM", group2)) %>% 
  mutate(species = case_when(grepl("Human",group2)~"Human",
                             grepl("Murine",group2)~"Murine")) %>% 
  select(species, FDR, pathway) %>% 
  pivot_wider(names_from = species, values_from = FDR) %>% 
  filter(Human < signif.pw & Murine < signif.pw) %>% 
  pull(pathway) %>% unique()

#Relax fdr, signif in at least 1
MDM_signif2 <- H_all %>% 
  filter(grepl("MDM", group2)) %>% 
  mutate(species = case_when(grepl("Human",group2)~"Human",
                             grepl("Murine",group2)~"Murine")) %>% 
  select(species, FDR, pathway) %>% 
  pivot_wider(names_from = species, values_from = FDR) %>% 
  filter(Human < 0.1 | Murine < 0.1) %>% 
  pull(pathway) %>% unique()

#Signif and different directions in human vs murine
unique_MDM <- H_all %>% 
  filter(grepl("MDM", group2) & pathway %in% MDM_signif2) %>% 
  select(group2, pathway, NES) %>% 
  pivot_wider(names_from = group2, values_from = NES) %>% 
  filter((`Human\nMDM` < 0 & `Murine\nBMDM`> 0) | (`Human\nMDM` > 0 & `Murine\nBMDM` < 0))

#Same direction in human and murine 
shared_MDM <- H_all %>% 
  filter(grepl("MDM", group2) & pathway %in% MDM_signif) %>% 
  select(group2, pathway, NES) %>% 
  pivot_wider(names_from = group2, values_from = NES) %>% 
  filter((`Human\nMDM` < 0 & `Murine\nBMDM` < 0) | (`Human\nMDM` > 0 & `Murine\nBMDM` > 0)) 

#### Filter data ####
dat_shar <- H_all %>% 
  mutate(fdr.group = factor(fdr.group)) %>% 
  filter(pathway %in% shared_MDM$pathway) %>% 
  filter(grepl("MDM", group2)) %>%   
  #Clean pathway names
  mutate(pathway.orig=pathway,
         pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway),
         pathway = gsub("interferon gamma", "IFNG", pathway),
         pathway = gsub("interferon alpha", "IFNA", pathway),
         pathway = gsub("tnfa", "TNFA", pathway),
         pathway = gsub("nfkb", "NF-kB", pathway),
         pathway = gsub("il2 stat5", "IL2 STAT5", pathway),
         pathway = gsub("il6 jak stat", "IL6 JAK STAT", pathway),
         pathway = gsub("uv", "UV", pathway),
         pathway = gsub("dna", "DNA", pathway),
         pathway = gsub("wnt beta catenin", "Wnt/B-catenin", pathway),
         pathway = gsub("notch", "NOTCH", pathway),
         pathway = gsub("pi3k akt mtor", "PI3K/AKT/mTOR", pathway),
         pathway = gsub("kras", "KRAS", pathway),
         pathway = gsub("mtorc1", "MTORC1", pathway),
         pathway = gsub("tgf beta", "TGF-B", pathway)
  ) %>% 
  mutate(pathway = str_wrap(pathway, width = 23))

dat_unq <- H_all %>% 
  filter(pathway %in% unique_MDM$pathway) %>% 
  # filter(group2 %in% c("Human AM","Murine AM scBCG")) %>% 
  filter(grepl("MDM", group2)) %>% 
  #Clean pathway names
  mutate(pathway.orig=pathway,
         pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway),
         pathway = gsub("myc", "MYC", pathway),
         pathway = gsub("g2m", "G2M", pathway)
  )%>% 
  mutate(pathway = str_wrap(pathway, width = 8))

#### GSEA plot: shared ####
#pathway order
path.ord <- dat_shar %>% 
  filter(group2=="Human\nMDM") %>% 
  arrange(NES)

p1 <- dat_shar %>% 
  mutate(group2 = factor(group2, levels=c("Human\nMDM",
                                          "Murine\nBMDM"))) %>% 
  mutate(pathway = factor(pathway, levels=c(path.ord$pathway))) %>% 
  ggplot(aes(x=pathway, y=group2, 
             color=NES, size=fdr.group)) +
  geom_point(show.legend = TRUE) +
  theme_classic(base_size = 10) +
  coord_flip() + 
  labs(y="", x="Hallmark pathway",
       fill="", size="", 
       color = "+Mtb\n\n\nNES\n\n\n-Mtb",
       title="Shared responses") +
  scale_color_gradient2(low="darkblue", mid = "white", high = "red",
                        midpoint = 0, 
                        limits = c(min(c(dat_shar$NES,dat_unq$NES)), 
                                   max(c(dat_shar$NES,dat_unq$NES)))) +
  guides(color = guide_colorbar(title.position = "right",
                                order = 1),
         size = guide_legend(order = 2)) +
  scale_size_manual(values=c(1,3,4), drop=FALSE) +
  theme(legend.title = element_text(size=9),
        legend.margin = margin(-0.2,0,0,0, unit="cm"),
        legend.key.height = unit(0.8, "lines")
  )
# p1

#### Gene heatmap: shared ####
#Hallmark genes
db <- combine_msigdb(category="H")

#List pathways
pw_to_plot <- path.ord %>% 
  slice_max(NES, n=3) %>% 
  pull(pathway.orig)

genes.temp <- select_genes(cell="MDM", db=db, 
                           pw=pw_to_plot,
                           FDR = 0.1, topLE=20)
hm.ls <- gene_heatmaps(cell="MDM", genes=genes.temp, 
                       pw=pw_to_plot,
                       category="H", 
                       font_size = 8, row_splits = 4,
                       cluster_by = "FC",
                       pw_anno = "all")


#### GSEA plot: unique ####
#pathway order
path.ord2 <- dat_unq %>% 
  filter(group2=="Human\nMDM") %>% 
  arrange(NES)

p3 <- dat_unq %>% 
  mutate(group2 = factor(group2, levels=c("Human\nMDM",
                                          "Murine\nBMDM"))) %>% 
  mutate(pathway = factor(pathway, levels=c(path.ord2$pathway))) %>% 
  ggplot(aes(x=pathway, y=group2, 
             color=NES, size=fdr.group)) +
  geom_point(show.legend = TRUE) +
  theme_classic(base_size = 10) +
  coord_flip() + 
  labs(y="", x="Hallmark pathway",
       fill="", size="", 
       color = "+Mtb\n\n\nNES\n\n\n-Mtb",
       title="Differential responses") +
  scale_color_gradient2(low="darkblue", mid = "white", high = "red",
                        midpoint = 0, 
                        limits = c(min(c(dat_shar$NES,dat_unq$NES)), 
                                   max(c(dat_shar$NES,dat_unq$NES)))) +
  guides(color = guide_colorbar(title.position = "right",
                                order = 1),
         size = guide_legend(order = 2)) +
  scale_size_manual(values=c(1,3,4), drop=FALSE) +
  theme(legend.position="none"
  )
# p3

#### Gene heatmap: unique ####
pw_to_plot2 <- path.ord2$pathway.orig[c(1,3)]

genes.temp2 <- select_genes(cell="MDM", db=db, 
                           pw=pw_to_plot2,
                           FDR = 0.1, topLE=20)
hm.ls2 <- gene_heatmaps(cell="MDM", genes=genes.temp2, 
                        pw=pw_to_plot2,
                        category="H", 
                        font_size = 8, row_splits = 3,
                        cluster_by = "FC",
                        pw_anno = "all")

#### Save Figure 3 ####
lo <- "
AC
AC
BD
BD
BD
BD
"
p_all <- p1+hm.ls[2]+p3+hm.ls2[2] +
  plot_layout(design = lo) +
  # plot_layout(widths=c(1,1.2,1.2)) +
  plot_annotation(tag_levels = "A")
# p_all

ggsave(p_all, filename="publication/temp/Fig3_MDM_fgsea.pdf",
       width=6, height=8)

#### Save Figure S3 ####
p_all2 <- plot_spacer() +hm.ls[1]+hm.ls2[1] +
  plot_layout(widths=c(0,1,1)) +
  plot_annotation(tag_levels = list(c("A  Shared responses","B  Differential responses")))
# p_all2
ggsave(p_all2, filename="publication/temp/FigS3_MDM_fgsea_express.pdf",
       width=5.5, height=6)

