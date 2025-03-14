library(tidyverse)
library(limma)
library(viridis)
library(ComplexHeatmap)
library(patchwork)
signif.pw <- 0.01
signif.gene <- 0.1
source("publication/plot_fxns.R")
#
#### Data ####
attach("results/fgsea.RData")

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
#Same direction in human and mouse BCG and coTB AM
#Human + at least 1 murine signif
AM_signif <- H_all %>% 
  filter(grepl("AM", group2)) %>% 
  mutate(species = case_when(grepl("Human",group2)~"Human",
                             grepl("Murine",group2)~"Murine")) %>% 
  group_by(species, pathway) %>% 
  summarise(minFDR = min(FDR, na.rm=TRUE)) %>% 
  pivot_wider(names_from = species, values_from = minFDR) %>% 
  filter(Human < signif.pw & Murine < signif.pw) %>% 
  pull(pathway) %>% unique()

#Relax fdr
AM_signif2 <- H_all %>% 
  filter(grepl("AM", group2)) %>% 
  mutate(species = case_when(grepl("Human",group2)~"Human",
                             grepl("Murine",group2)~"Murine")) %>% 
  group_by(species, pathway) %>% 
  summarise(minFDR = min(FDR, na.rm=TRUE)) %>% 
  pivot_wider(names_from = species, values_from = minFDR) %>% 
  filter(Human < 0.1 & Murine < 0.1) %>% 
  pull(pathway) %>% unique()

#Signif and different directions in human vs at least 1 mouse AM
unique_AM <- H_all %>% 
  filter(grepl("AM", group2) & pathway %in% AM_signif2) %>% 
  select(group2, pathway, NES) %>% 
  pivot_wider(names_from = group2, values_from = NES) %>% 
  rowwise() %>% 
  mutate(mouse_min = min(`scBCG\nMurine AM`,`coMtb\nMurine AM`,
                         na.rm = TRUE),
         mouse_max = max(`scBCG\nMurine AM`,`coMtb\nMurine AM`,
                         na.rm = TRUE)) %>% 
  filter((`Human\nAM` < 0 & mouse_max > 0) | (`Human\nAM` > 0 & mouse_min < 0))

#Same direction in human and murine scBCG or coTB
shared_AM <- H_all %>% 
  filter(grepl("AM", group2) & pathway %in% AM_signif) %>% 
  select(group2, pathway, NES) %>% 
  pivot_wider(names_from = group2, values_from = NES) %>% 
  rowwise() %>% 
  mutate(mouse_min = min(`scBCG\nMurine AM`,`coMtb\nMurine AM`),
         mouse_max = max(`scBCG\nMurine AM`,`coMtb\nMurine AM`)) %>% 
  filter((`Human\nAM` < 0 & mouse_max < 0) | (`Human\nAM` > 0 & mouse_min > 0)) 

#### Filter data ####
dat_shar <- H_all %>% 
  filter(pathway %in% shared_AM$pathway) %>% 
  filter(grepl("AM", group2))  %>%   
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
         pathway = gsub("kras", "KRAS", pathway)
  ) %>% 
  mutate(pathway = str_wrap(pathway, width = 23))

dat_unq <- H_all %>% 
  filter(pathway %in% unique_AM$pathway) %>% 
  # filter(group2 %in% c("Human AM","Murine AM scBCG")) %>% 
  filter(grepl("AM", group2)) %>% 
  #Clean pathway names
  mutate(pathway.orig=pathway,
         pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway)
  ) %>% 
  mutate(pathway = str_wrap(pathway, width = 10))

#### GSEA plot: shared ####
#pathway order
path.ord1 <- dat_shar %>% 
  filter(group2=="Human\nAM") %>% 
  arrange(NES)

p1 <- dat_shar %>% 
  mutate(group2 = factor(group2, levels=c("Human\nAM",
                                          "coMtb\nMurine AM",
                                          "scBCG\nMurine AM",
                                          "control\nMurine AM"))) %>% 
  mutate(pathway = factor(pathway, levels=c(path.ord1$pathway))) %>% 
  ggplot(aes(x=pathway, y=group2, 
             color=NES, size=fdr.group)) +
  geom_point() +
  theme_classic(base_size = 10) +
  coord_flip() + 
  labs(y="", x="Hallmark pathway",
       fill="", size="", 
       color = "+Mtb\n\n\nNES\n\n\n-Mtb",
       title="Shared responses") +
  scale_color_gradient2(low="darkblue", mid = "white", high = "red",
                        midpoint = 0, 
                        limits = c(min(c(dat_shar$NES, dat_unq$NES)), 
                                   max(c(dat_shar$NES, dat_unq$NES)))) +
  
  guides(color = guide_colorbar(title.position = "right",
                                order = 1),
         size = guide_legend(order = 2)) +
  scale_size_manual(values=c(1,3,4)) +
  theme(legend.title = element_text(size=9),
        legend.margin = margin(-0.2,0,0,0, unit="cm"),
        legend.key.height = unit(0.8, "lines"),
        legend.key.width= unit(0.3, 'cm')
  ) 
# p1

#### Gene heatmap: shared ####
#Hallmark genes
db <- combine_msigdb(category="H")

pw_to_plot <- path.ord1 %>% 
  slice_max(NES, n=4) %>% 
  pull(pathway.orig)
pw_to_plot <- pw_to_plot[c(1,3,2,4)]

genes.temp <- select_genes(cell="AM", db=db, 
                           pw=pw_to_plot,
                           FDR = 0.1, topLE=20)
hm.ls <- gene_heatmaps(cell="AM", genes=genes.temp, 
                       pw=pw_to_plot,
                       category="H", 
                       font_size = 8, row_splits = 3,
                       cluster_by = "FC",
                       pw_anno = "all")

#### GSEA plot: unique ####
#pathway order
path.ord2 <- dat_unq %>% 
  filter(group2=="Human\nAM") %>% 
  arrange(NES)

p3 <- dat_unq %>% 
  mutate(group2 = factor(group2, levels=c("Human\nAM",
                                          "coMtb\nMurine AM",
                                          "scBCG\nMurine AM",
                                          "control\nMurine AM"))) %>% 
  mutate(pathway = factor(pathway, levels=c(path.ord2$pathway))) %>% 
  ggplot(aes(x=pathway, y=group2, 
             color=NES, size=fdr.group)) +
  geom_point() +
  theme_classic(base_size = 10) +
  coord_flip() + 
  labs(y="", x="Hallmark pathway",
       fill="", size="", 
       color = "Up +Mtb     NES     Down +Mtb",
       title="Differential responses") +
  scale_color_gradient2(low="darkblue", mid = "white", high = "red",
                        midpoint = 0, 
                        limits = c(min(c(dat_shar$NES, dat_unq$NES)), 
                                   max(c(dat_shar$NES, dat_unq$NES)))) +
  scale_size_manual(values=c(1,3,4)) +
  theme(legend.position="none"
  ) 
# p3

#### Gene heatmap: unique ####
pw_to_plot2 <- path.ord2 %>% 
  pull(pathway.orig)
genes.temp2 <- select_genes(cell="AM", db=db, 
                            pw=pw_to_plot2[2],
                            FDR = 0.1, topLE=20)
hm.ls2 <- gene_heatmaps(cell="AM", genes=genes.temp2, 
                        pw=pw_to_plot2[2],
                        category="H", 
                        font_size = 8, row_splits = 3,
                        cluster_by = "FC",
                        pw_anno = "none")

#Separate suppl for androgen response
genes.temp3 <- select_genes(cell="AM", db=db, 
                            pw=pw_to_plot2[1],
                            FDR = 0.1, topLE=20)
hm.ls3 <- gene_heatmaps(cell="AM", genes=genes.temp2, 
                        pw=pw_to_plot2[1],
                        category="H", 
                        font_size = 8, row_splits = 3,
                        cluster_by = "FC",
                        pw_anno = "none")

#### Save Fig 2 ####
lo <- "
AC
AC
BD
BD
BD
"
p_all <- p1+hm.ls[2]+p3+hm.ls2[2] +
  plot_layout(design = lo) +
  # plot_layout(widths=c(1,1.2,1.2)) +
  plot_annotation(tag_levels = "A")
p_all

# ggsave(p_all, filename="publication/Fig2_AM_fgsea.png",
#        width=7.08, height=6)

ggsave(p_all, filename="publication/temp/Fig2_AM_fgsea.pdf",
       width=7.08, height=6)

#### Save Figure S1 ####
p_all2 <- plot_spacer() +hm.ls[1]+hm.ls2[1] +
  plot_layout(widths=c(0,1,1)) +
  plot_annotation(tag_levels = list(c("A  Shared responses","B  Differential responses")))
# p_all2
ggsave(p_all2, filename="publication/temp/FigS1_AM_fgsea_express.pdf",
       width=7.08, height=3.6)
