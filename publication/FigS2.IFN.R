library(tidyverse)
library(patchwork)

#### Fold change ####
attach("data_clean/human_dat_clean.RData")
attach("results/lme_human.RData")
FC_human <- model_hum$lme.contrast %>% 
  filter(grepl("Media",contrast_ref) & grepl("Mtb", contrast_lvl)) %>% 
  separate(contrast_ref, into=c("cell"), sep=" ", extra="drop") %>% 
  separate(contrast_lvl, into=c("cell2","variable"), sep=" ") %>% 
  filter(cell==cell2) %>% 
  mutate(species="Human") %>% 
  select(model,species,cell,gene,variable,estimate,FDR)

attach("data_clean/mouse_dat_clean.RData")
attach("results/lm_mouse.RData")

FC_mouse <- model_mur$lm.contrast %>% 
  filter(grepl("Media",contrast_ref) & grepl("Mtb", contrast_lvl)) %>% 
  separate(contrast_ref, into=c("cell"), sep=" ", extra="drop") %>% 
  separate(contrast_lvl, into=c("cell2","variable"), sep=" ") %>% 
  filter(cell==cell2) %>% 
  mutate(species="Murine") %>% 
  select(model,species,cell,gene,variable,estimate,FDR)

FC_all <- bind_rows(FC_human, FC_mouse) %>% 
  mutate(gene = toupper(gene))

#### IFN ####
ifn_genes <- unique(FC_all$gene)
ifn_genes <- ifn_genes[grepl("^IFN",ifn_genes)] #IFN

FC_ifn <- FC_all %>% 
  filter(gene %in% ifn_genes) 

#### Determine NS vs missing ####
# all IFN genes in human
ensembl_human <- biomaRt::useEnsembl(biomart="ensembl", 
                                     dataset="hsapiens_gene_ensembl")
ifn_human <- biomaRt::getBM(attributes=c("hgnc_symbol", "gene_biotype"), 
                            mart=ensembl_human) %>% 
  mutate(symbol=hgnc_symbol) %>% 
  filter(grepl("^IFN", hgnc_symbol) & gene_biotype=="protein_coding") %>% 
  pull(symbol) %>% unique()

# all IFN genes in mouse
ensembl_mouse <- biomaRt::useEnsembl(biomart="ensembl", 
                                     dataset="mmusculus_gene_ensembl")
ifn_mouse <- biomaRt::getBM(attributes=c("mgi_symbol", "hgnc_symbol","gene_biotype"), 
                            mart=ensembl_mouse) %>% 
  mutate(symbol=toupper(mgi_symbol)) %>% 
  filter(gene_biotype=="protein_coding") %>% 
  filter(grepl("^IFN", symbol) |  grepl("^IFN", hgnc_symbol) ) %>% 
  pull(symbol, hgnc_symbol) %>% unique()

#shared
ifn_shared <- intersect(ifn_human, ifn_mouse)
#unique
ifn_hum_only <- ifn_human[!ifn_human %in% ifn_shared]
ifn_mur_only <- ifn_mouse[!ifn_mouse %in% ifn_shared]
#in expression data
ifn_shared <- ifn_shared[ifn_shared %in% unique(FC_ifn$gene)]
ifn_hum_only <- ifn_hum_only[ifn_hum_only %in% unique(FC_ifn$gene)]
ifn_mur_only <- ifn_mur_only[ifn_mur_only %in% unique(FC_ifn$gene)]

#Add groups for NS, missing, or not measured
FC_miss <- FC_ifn %>% 
  select(species,cell,variable,estimate, gene) %>% 
  mutate(estimate="Y") %>% 
  pivot_wider(names_from = gene, values_from = estimate) %>% 
  mutate(across(any_of(c(ifn_shared)), ~replace(.,is.na(.),"N/D"))) %>% 
  mutate(across(any_of(c(ifn_hum_only, ifn_mur_only)), 
                ~replace(.,is.na(.),"MISS"))) %>% 
  pivot_longer(-c(species,cell,variable), names_to = "gene", 
               values_to = "miss") %>% 
  filter(miss%in% c("N/D","MISS")) %>% 
  mutate(fdr.group = "NS")

FC_ifn_miss <- FC_ifn %>% 
  mutate(miss="D") %>% 
  mutate(miss=" ") %>% 
  mutate(fdr.group = case_when(FDR < 0.01 ~ "FDR < 0.01",
                               FDR < 0.1 ~ "FDR < 0.1",
                               TRUE~"NS"))%>%
  bind_rows(FC_miss)

#### Split Type I / II / III ####
FC_ifn_miss_type <- FC_ifn_miss %>% 
  mutate(type = case_when(grepl("^IFNA", gene)~"Type I",
                          grepl("^IFNB", gene)~"Type I",
                          grepl("^IFNE", gene)~"Type I",
                          grepl("^IFNK", gene)~"Type I",
                          
                          grepl("^IFNG", gene)~"Type II",
                          
                          grepl("^IFNL", gene)~"Type III",
  )) %>% 
  #order IFN
  mutate(gene = factor(gene, levels=sort(ifn_genes)[c(1,3:7,2,8:length(ifn_genes))])) %>% 
  mutate(fdr.group = factor(fdr.group, levels=rev(c("FDR < 0.01",
                                                    "FDR < 0.1",
                                                    "NS")))) %>% 
  #x label
  separate(cell, into=c("cell","treatment"), sep="_") %>% 
  mutate(x=case_when(species == "Human" ~ species,
                     species == "Murine" & cell=="BMDM" ~ species,
                     species == "Murine" & cell!="BMDM" ~ paste(species,treatment, sep="\n"))) %>% 
  mutate(cell = recode(cell, "BMDM"="MDM")) %>% 
  mutate(x = recode(x, "Murine"="Murine\ncontrol", "Murine\ncoTB"="Murine\ncoMtb")) %>% 
  mutate(x=factor(x, levels=c("Human",
                              "Murine\ncoMtb",
                              "Murine\nscBCG",
                              "Murine\ncontrol")))

#### Plot shared IFN ####
p1a <- FC_ifn_miss_type %>%
  filter(gene %in% ifn_shared) %>% 
  filter(type=="Type I") %>% 
  
  ggplot(aes(x=x, y=cell)) +
  geom_point(aes(size=fdr.group, fill=estimate, shape=miss)) +
  theme_classic(base_size = 10) +
  scale_fill_gradient2(low="darkblue", mid = "white", high = "red",
                       midpoint = 0, limits = c(min(FC_ifn_miss_type$estimate, na.rm = TRUE),
                                                max(FC_ifn_miss_type$estimate, na.rm = TRUE))) +
  scale_shape_manual(values=c(21,4)) +
  scale_size_manual(drop = FALSE, values = c(1,3,4)) +
  labs(x="",y="",fill="+Mtb - Media\nlog2 fold change",
       size="", shape="",title="Type I IFN") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(alpha = c(0, 1))), 
         size = guide_legend(order = 2)) +
  facet_wrap(~gene, scales="free_x", ncol=6) 
p1a

p1b <- FC_ifn_miss_type %>%
  filter(gene %in% ifn_shared) %>% 
  filter(type!="Type I") %>% 
  
  ggplot(aes(x=x, y=cell)) +
  geom_point(aes(size=fdr.group, fill=estimate, shape=miss)) +
  theme_classic(base_size = 10) +
  scale_fill_gradient2(low="darkblue", mid = "white", high = "red",
                       midpoint = 0, limits = c(min(FC_ifn_miss_type$estimate, 
                                                    na.rm = TRUE),
                                                max(FC_ifn_miss_type$estimate, 
                                                    na.rm = TRUE))) +
  scale_shape_manual(values=c(21,4), drop=FALSE) +
  scale_size_manual(drop = FALSE, values = c(1,3,4)) +
  labs(x="",y="",fill="+Mtb - Media\nlog2 fold change",
       size="", shape="",title="Type II and III IFN") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(alpha = c(0))),
         size = guide_legend(order = 2)) +
  facet_wrap(type~gene, scales="free_x", ncol=4) 
# p1b

#### Plot human only IFN ####
p1c <- FC_ifn_miss_type %>%
  filter(gene %in% ifn_hum_only) %>% 
  filter(grepl("Human",x)) %>% 
  
  ggplot(aes(x=x, y=cell)) +
  geom_point(aes(size=fdr.group, fill=estimate, shape=miss)) +
  theme_classic(base_size = 10) +
  scale_fill_gradient2(low="darkblue", mid = "white", high = "red",
                       midpoint = 0, limits = c(min(FC_ifn_miss_type$estimate, na.rm = TRUE),
                                                max(FC_ifn_miss_type$estimate, na.rm = TRUE))) +
  scale_shape_manual(values=c(21,4)) +
  scale_size_manual(values = c(1,4)) +
  labs(x="",y="",fill="+Mtb - Media\nlog2 fold change",
       size="", title="Human only IFN") +
  guides(shape = guide_legend(order = 1,
                              override.aes = list(alpha = c(0))),
         size = guide_legend(order = 2)) +
  facet_wrap(type~gene, scales="free_x", ncol=4) +
  theme(legend.position ="none")
  
p1c

#### Plot murine only IFN ####
#None
ifn_mur_only

#### Save ####
lo <- "
AAA
AAA
AAA
BBC
"

p_all <- p1a+p1b+p1c + plot_annotation(tag_levels = "A") + plot_layout(design=lo, guides = "collect")
p_all

ggsave("publication/FigS2.IFN.png", p_all, height=5,width=13)
ggsave("publication/FigS2.IFN.pdf", p_all, height=5,width=13)
