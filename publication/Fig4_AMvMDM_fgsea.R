library(tidyverse)
# library(limma)
# library(viridis)
# library(ComplexHeatmap)
library(patchwork)
signif.pw <- 0.01
signif.gene <- 0.1

#### Data ####
attach("results/fgsea_AMvMDM.RData")

#human MDM-AM
#murine BMDM-AM

H_all <- bind_rows(hum_fgsea2, mur_fgsea2) %>% 
  separate(group, into=c("species","tb"), sep=" ", extra = "merge") %>% 
  filter(tb %in% c("media","mtb","control media","control mtb")) %>% 
  mutate(species = recode(species, "human"="Human", "murine"="Murine"),
         tb = ifelse(grepl("media",tb), "Media", "+Mtb"),
         tb = factor(tb, levels=c("Media","+Mtb"))) %>% 
  mutate(fdr.group = case_when(FDR < 0.01 ~ "FDR < 0.01",
                               FDR < 0.1 ~ "FDR < 0.1",
                               TRUE~"NS"),
         fdr.group = factor(fdr.group, levels=rev(c("FDR < 0.01",
                                                    "FDR < 0.1",
                                                    "NS"))))
#### Determine overlap ####
## MEDIA
#Same direction and signif in human and mouse control
media_signif <- H_all %>% 
  filter(tb=="Media") %>% 
  select(species, pathway, FDR) %>% 
  pivot_wider(names_from = species, values_from = FDR) %>% 
  filter(Human < signif.pw & Murine < signif.pw) %>% 
  pull(pathway) %>% unique()

#Same direction in human and murine control
shared_media <- H_all %>% 
  filter(tb=="Media" & pathway %in% media_signif) %>% 
  select(species, pathway, NES) %>% 
  pivot_wider(names_from = species, values_from = NES) %>% 
  filter((Human < 0 & Murine < 0) |
           (Human > 0 & Murine > 0))

## MTB
#Same direction and signif in human and mouse control
mtb_signif <- H_all %>% 
  filter(tb=="+Mtb") %>% 
  select(species, pathway, FDR) %>% 
  pivot_wider(names_from = species, values_from = FDR) %>% 
  filter(Human < signif.pw & Murine < signif.pw) %>% 
  pull(pathway) %>% unique()

#Same direction in human and murine control
shared_mtb <- H_all %>% 
  filter(tb=="+Mtb" & pathway %in% mtb_signif) %>% 
  select(species, pathway, NES) %>% 
  pivot_wider(names_from = species, values_from = NES) %>% 
  filter((Human < 0 & Murine < 0) |
           (Human > 0 & Murine > 0))

dat <- H_all %>% 
  filter(pathway %in% c(shared_media$pathway, shared_mtb$pathway)) %>% 
  #Clean pathway names
  mutate(pathway.orig=pathway,
         pathway = gsub("HALLMARK_", "", pathway),
         pathway = gsub("_", " ", pathway),
         pathway = tolower(pathway),
         pathway = gsub("tnfa", "TNFA", pathway),
         pathway = gsub("nfkb", "NF-kB", pathway),
         pathway = gsub("e2f", "E2F", pathway),
         pathway = gsub("g2m", "G2M", pathway)
  ) 

#### GSEA plot ####
#pathway order
# path.ord <- dat %>%
#   filter(group2=="Human\nMedia") %>%
  # arrange(NES)

path.ord <- rev(c("E2F targets","G2M checkpoint","allograft rejection",
              "inflammatory response","TNFA signaling via NF-kB",
              "myogenesis"))
p1 <- dat %>% 
  mutate(pathway = factor(pathway, levels=path.ord)) %>% 
  ggplot(aes(x=pathway, y=species, 
             color=NES, size=fdr.group)) +
  geom_point() +
  theme_classic(base_size = 10) +
  coord_flip() + 
  labs(y="", x="Hallmark pathway",
       fill="", size="", 
       color = "Up in MDM\n\n\nNES\n\n\nUp in AM") +
  scale_color_gradient2(low="#762a83", mid = "white", high = "#1b7837",
                        midpoint = 0, limits = c(min(dat$NES), max(dat$NES)),) +
  guides(color = guide_colorbar(title.position = "right",
                                order = 1),
         size = guide_legend(order = 2)) +
  scale_size_manual(values=c(1,3,4)) +
  theme(legend.title = element_text(size=9),
        legend.margin = margin(-0.2,0,0,0, unit="cm"),
        legend.key.height = unit(0.8, "lines"),
        panel.border = element_rect(fill=NA)
  ) +
  facet_wrap(~tb, scales="free_x")
p1

#### Save ####
ggsave(p1, filename="publication/Fig4_AMvMDM_fgsea.png",
       width=5.5, height=2)

ggsave(p1, filename="publication/Fig4_AMvMDM_fgsea.pdf",
       width=5.5, height=2)
