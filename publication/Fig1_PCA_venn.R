library(tidyverse)
library(plotly)
library(ggvenn)
library(patchwork)

#### Data ####
load("data_clean/human_dat_clean.RData")
load("data_clean/mouse_dat_clean.RData")

#Recode murine genes to human names
source("publication/plot_fxns.R")
#Hallmark genes
dbH <- combine_msigdb(category="H")
db2 <- combine_msigdb(category="C2")
db5 <- combine_msigdb(category="C5")
db7 <- combine_msigdb(category="C7")

db <- bind_rows(dbH, db2, db5, db7) %>% 
  distinct()
rm(dbH, db2, db5, db7)

voom_mur_rename <- voom_mur
voom_mur_rename$genes <- voom_mur_rename$genes %>% 
  rename(mouse_gene_symbol=symbol) %>% 
  left_join(db) %>% 
  rename(symbol=human_gene_symbol)
voom_mur_rename$E <- as.data.frame(voom_mur_rename$E) %>% 
  rownames_to_column("mouse_gene_symbol") %>% 
  left_join(db) %>% 
  select(-mouse_gene_symbol) %>% 
  drop_na(human_gene_symbol) %>% 
  #Mean expression of human genes with multiple mouse names
  group_by(human_gene_symbol) %>% 
  summarise(across(where(is.numeric), ~mean(., na.rm=TRUE))) %>% 
  ungroup() %>% 
  column_to_rownames("human_gene_symbol")

#Fix metadata to match human
voom_mur_rename$targets <- voom_mur_rename$targets %>% 
  separate(treatment, into = c("cell","treatment"), sep="_")

#Combine mouse and human
voom_all <- voom_hum
voom_all$genes <- full_join(voom_hum$genes, voom_mur_rename$genes)
voom_all$targets <- bind_rows(voom_hum$targets %>% mutate(species="Human"),
                              voom_mur_rename$targets %>% mutate(species="Murine"))
voom_all$E <- full_join(as.data.frame(voom_hum$E) %>% 
                          rownames_to_column(),
                        as.data.frame(voom_mur_rename$E)%>% 
                          rownames_to_column()) %>% 
  column_to_rownames() %>% 
  #remove NA
  drop_na()

dim(voom_all$E)

#### Fold change ####
voom_FC <- voom_all$E %>% 
  rownames_to_column() %>% 
  pivot_longer(-rowname, names_to = "libID") %>% 
  full_join(voom_all$targets %>% select(libID, ptID, cell,
                                        mtb, species, treatment)) %>% 
  select(-libID) %>% 
  pivot_wider(names_from = mtb) %>% 
  mutate(value=Mtb-Media) %>% 
  mutate(name = paste(ptID, cell, species, treatment)) %>% 
  select(rowname,name,value) %>% 
  pivot_wider() %>% 
  column_to_rownames()

#### PCA data ####
set.seed(8456)
PCA_FC <- prcomp(t(voom_FC), scale. = TRUE, center = TRUE)
PC1.FC <- paste("PC", 1, " (", 
                round(summary(PCA_FC)$importance[2,1], 2) * 100, 
                "%)", sep = "")
PC2.FC <- paste("PC", 2, " (", 
                round(summary(PCA_FC)$importance[2, 2], 2) * 100, 
                "%)", sep = "")
PC3.FC <- paste("PC", 3, " (", 
                round(summary(PCA_FC)$importance[2, 3], 2) * 100, 
                "%)", sep = "")

pca.dat.fc <- as.data.frame(PCA_FC$x) %>% 
  rownames_to_column("libID") %>%
  separate(libID, into=c("ptID", "Cell", "Species", "Treatment"), sep=" ") %>% 
  mutate(group = paste(Species, Cell)) %>% 
  mutate(Treatment = ifelse(Treatment%in%c("","NA")," ", Treatment),
         Treatment = gsub("coTB","coMtb",Treatment),
         Treatment = factor(Treatment, 
                            levels = c(" ","coMtb","scBCG","control")))

#### PCA Plot ####
# p1 <- pca.dat.fc %>%
#   ggplot(aes(PC1, PC2, color = group, shape = Treatment)) +
#   geom_point(size = 3) +
#   theme_classic() +
#   labs(x = PC1.FC, y = PC2.FC, color = "Cell") +
#   coord_fixed(ratio = 1) +
#   scale_color_manual(values=c("#332288","#DDCC77","#88CCEE","#CC6677"))
# p1

#3D plot
scene = list(camera = list(eye = list(x = -0.75, 
                                      y = 1.25,
                                      z = 1)),
             xaxis = list(title = PC1.FC, 
                          showticklabels = FALSE, ticks = ""),
             yaxis = list(title = PC2.FC, 
                          showticklabels = FALSE, ticks = ""),
             zaxis = list(title = PC3.FC, 
                          showticklabels = FALSE, ticks = ""))

fig <- plot_ly(pca.dat.fc, 
               x = ~PC1,
               y = ~PC2, 
               z = ~PC3, 
               color = ~group, 
               colors = c("#332288","#DDCC77","#88CCEE","#CC6677"),
               symbol = ~Treatment,
               symbols = c("circle","square","diamond"),
               marker = list(opacity = 0.6)) %>%
  add_markers() %>% 
  layout(scene = scene)

# fig

# reticulate::py_run_string("import sys")
save_image(fig, "publication/Fig1_PCA.pdf",
           width=600, height=400, scale=3)
# save_image(fig, "publication/Fig1_PCA.png",
#            width=1000, height=300, scale=2)

#### Venn data ####
# human
load("results/lme_human.RData")
key_hum <- data.frame(
  contrast_ref = c("AM Media","MDM Media"),
  contrast_lvl = c("AM Mtb",  "MDM Mtb")
)
model_hum <- model_hum$lme.contrast %>% 
  inner_join(key_hum) %>% 
  rename(human_gene_symbol=gene) %>% 
  left_join(db) %>%
  filter(FDR<0.1)

#mouse
load("results/lm_mouse.RData")
key_mur <- data.frame(
  contrast_ref = c("AM_control Media","AM_coTB Media", "AM_scBCG Media",
                   "BMDM Media"),
  contrast_lvl = c("AM_control Mtb",  "AM_coTB Mtb", "AM_scBCG Mtb",
                   "BMDM Mtb")
)
model_mur <- model_mur$lm.contrast %>% 
  inner_join(key_mur) %>% 
  rename(mouse_gene_symbol=gene) %>% 
  left_join(db) %>% 
  mutate(human_orig = human_gene_symbol) %>% 
  mutate(human_gene_symbol = ifelse(is.na(human_gene_symbol),
                                    mouse_gene_symbol,human_gene_symbol)) %>% 
  filter(FDR<0.1)

#### MDM ####
venn.ls2 <- list()
venn.ls2[["Human\nMDM"]] <- model_hum %>%
  filter(contrast_ref == "MDM Media") %>% 
  pull(human_gene_symbol) %>% unique()
venn.ls2[["Murine\nBMDM"]] <- model_mur %>% 
  filter(contrast_ref == "BMDM Media") %>% 
  pull(human_gene_symbol) %>% unique()
## Human signif w/o mouse counterpart
model_hum %>%
  filter(contrast_ref == "MDM Media") %>% 
  filter(is.na(mouse_gene_symbol)) %>% 
  pull(human_gene_symbol) %>% unique() %>% length()
## mouse signif w/o human counterpart
model_mur %>%
  filter(contrast_ref == "BMDM Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(mouse_gene_symbol) %>% unique() %>% length()

#### AM ####
#lists
venn.ls <- list()
venn.ls[["Human AM"]] <- model_hum %>% 
  filter(contrast_ref == "AM Media") %>% 
  pull(human_gene_symbol) %>% unique()
venn.ls[["Murine AM\ncoMtb"]] <- model_mur %>% 
  filter(contrast_ref == "AM_coTB Media") %>% 
  pull(human_gene_symbol) %>% unique()
venn.ls[["Murine AM\ncontrol"]] <- model_mur %>% 
  filter(contrast_ref == "AM_control Media") %>% 
  pull(human_gene_symbol) %>% unique()
venn.ls[["Murine AM\nscBCG"]] <- model_mur %>% 
  filter(contrast_ref == "AM_scBCG Media") %>% 
  pull(human_gene_symbol) %>% unique()

## Human signif w/o mouse counterpart
model_hum %>%
  filter(contrast_ref == "AM Media") %>% 
  filter(is.na(mouse_gene_symbol)) %>% 
  pull(human_gene_symbol) %>% unique() %>% length()
## mouse signif w/o human counterpart
model_mur %>%
  filter(contrast_ref == "AM_coTB Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(mouse_gene_symbol) %>% unique() %>% length()
model_mur %>%
  filter(contrast_ref == "AM_control Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(mouse_gene_symbol) %>% unique() %>% length()
model_mur %>%
  filter(contrast_ref == "AM_scBCG Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(mouse_gene_symbol) %>% unique() %>% length()

#no human counterpart venn
venn.ls3 <- list()
venn.ls3[["Human AM"]] <- c('test')
venn.ls3[["Murine AM\ncoMtb"]] <- model_mur %>% 
  filter(contrast_ref == "AM_coTB Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(human_gene_symbol) %>% unique()
venn.ls3[["Murine AM\ncontrol"]] <- model_mur %>% 
  filter(contrast_ref == "AM_control Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(human_gene_symbol) %>% unique()
venn.ls3[["Murine AM\nscBCG"]] <- model_mur %>% 
  filter(contrast_ref == "AM_scBCG Media") %>% 
  filter(is.na(human_orig)) %>%
  pull(human_gene_symbol) %>% unique()

### Venn ####
#AM
#Total unique
v1 <- ggvenn(venn.ls, show_percentage = FALSE,
             fill_color = c("#332288","#88CCEE","#88CCEE","#88CCEE"),
             set_name_size = 2.5,text_size = 2.5, stroke_size = 0.25) +
  # theme_classic() +
  lims(x=c(-2.5,2.5), y=c(-1.7,1.7))
# v1

#AM
#Missing ortholog
v3 <- ggvenn(venn.ls3, show_percentage = FALSE, 
             fill_color = c("#332288","#88CCEE","#88CCEE","#88CCEE"),
             set_name_size = 2.5,text_size = 2.5, stroke_size = 0.25)
# v3

#MDM
v2 <- ggvenn(venn.ls2, show_percentage = FALSE, 
             fill_color = c("#DDCC77","#CC6677"),
             set_name_size = 2.5,text_size = 2.5, stroke_size = 0.25)+
  # theme_classic() +
  lims(x=c(-2,2), y=c(-1.3,1.7))
# v2



v_all <- v1+v2 + #plot_annotation(tag_levels = list(c("C"))) +
  plot_layout(widths = c(1.5,1))
# v_all

ggsave(v_all, filename = "publication/temp/Fig1_venn.pdf", width=4, height=2)

#### DEG numbers for text ####
#### total DEGs ####
##human
model_hum %>%
  filter(contrast_ref == "MDM Media") %>% 
  pull(human_gene_symbol) %>% unique() %>% length()
model_hum %>%
  filter(contrast_ref == "AM Media") %>% 
  pull(human_gene_symbol) %>% unique() %>% length()

##Murine
model_mur %>% 
  filter(contrast_ref == "BMDM Media") %>% 
  pull(mouse_gene_symbol) %>% unique() %>% length()
model_mur %>% 
  filter(contrast_ref == "AM_coTB Media") %>% 
  pull(mouse_gene_symbol) %>% unique() %>% length()
model_mur %>% 
  filter(contrast_ref == "AM_control Media") %>% 
  pull(mouse_gene_symbol) %>% unique() %>% length()
model_mur %>% 
  filter(contrast_ref == "AM_scBCG Media") %>% 
  pull(mouse_gene_symbol) %>% unique() %>% length()

#### Orthologs ####
mdm_h <- model_hum %>% 
  filter(contrast_ref == "MDM Media") %>% 
  drop_na(human_gene_symbol) %>% 
  pull(human_gene_symbol) %>% unique()
model_mur %>% 
  filter(contrast_ref == "BMDM Media") %>% 
  drop_na(human_gene_symbol) %>% 
  filter(human_gene_symbol %in% mdm_h) %>% 
  pull(human_gene_symbol) %>% unique() %>% length()


am_h <- model_hum %>% 
  filter(contrast_ref %in% ("AM Media")) %>% 
  drop_na(human_gene_symbol) %>% 
  pull(human_gene_symbol) %>% unique()

#Murine no orth
am_m_1 <- model_mur %>% 
  filter(contrast_ref %in% c("AM_coTB Media","AM_scBCG Media",
                             "AM_control Media")) %>% 
  filter(is.na(human_orig)) %>% 
  pull(mouse_gene_symbol) %>% unique()
am_m_2 <- model_mur %>% 
  filter(contrast_ref %in% c("AM_coTB Media","AM_scBCG Media",
                             "AM_control Media")) %>% 
  filter(!is.na(human_orig)) %>% 
  pull(human_gene_symbol) %>% unique() 
length(am_m_1)
length(am_m_2)

am_m <- model_mur %>% 
  filter(contrast_ref %in% c("AM_coTB Media","AM_scBCG Media",
                            "AM_control Media")) %>% 
  drop_na(human_gene_symbol) %>% 
  filter(human_gene_symbol %in% am_h) %>% 
  pull(human_gene_symbol) %>% unique() %>% length()
