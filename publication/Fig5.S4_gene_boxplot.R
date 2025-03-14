library(tidyverse)
library(patchwork)
library(ggpubr)

#### Define genes ####
geneOI1 <- c("MAF","IL10")
geneOI2 <- list()
geneOI2[["Antimicrobial"]] <- c("NOS1","NOS2","TNF","IL1B","CAMP")
geneOI2[["Reactive oxygen species"]] <- c("CYBB")
# geneOI2[["Reactive nitrogen intermediates"]] <- c("ARG1","SLC7A2")
geneOI2[["Antimicrobial peptides"]] <- c("DEFB4A","LYZ","LYZ2","PLA2G2A")
# geneOI2[["Autophagy"]] <- c("MAP1LC3B","SQSTM1","ATG5","BECN1","TBK1")
geneOI2[["Other"]] <- c("NRF2","NRAMP1","HIF1A","NLRP3","CD14")
geneOI2[["Mendelian Susceptibility to Mycobacterial Disease"]] <- c("IL12p40","IL12R","STAT1","GATA2","IRF8","ISG15")

#### Gene expression ####
attach("data_clean/human_dat_clean.RData")
attach("data_clean/mouse_dat_clean.RData")

# combine
voom_hum <- RNAetc::collapse_voom(voom_hum, geneID = "symbol") %>% 
  mutate(species="Human")
voom_mur <- RNAetc::collapse_voom(voom_mur, geneID = "symbol") %>% 
  mutate(species="Murine")

#### Linear models ####
load("results/lme_human.RData")
load("results/lm_mouse.RData")

key_hum <- data.frame(
  contrast_ref = c("AM Media","MDM Media", "AM Media"),
  contrast_lvl = c("AM Mtb",  "MDM Mtb", "MDM Media")
)

model_hum <- model_hum$lme.contrast %>% 
  inner_join(key_hum) %>% 
  mutate(species="Human") %>% 
  filter(FDR < 0.1)

key_mur <- data.frame(
  contrast_ref = c("AM_control Media","AM_coTB Media", "AM_scBCG Media",
                   "BMDM Media",
                   "AM_control Media","AM_coTB Media", "AM_scBCG Media"),
  contrast_lvl = c("AM_control Mtb",  "AM_coTB Mtb", "AM_scBCG Mtb",
                   "BMDM Mtb",
                   "BMDM Media","BMDM Media","BMDM Media")
)

model_mur <- model_mur$lm.contrast %>% 
  inner_join(key_mur) %>% 
  mutate(species="Murine") %>% 
  filter(FDR < 0.1)

#### Combine human and mouse ####
#Expression
dat <- voom_mur %>% 
  #Merge murine and human
  mutate(symbol=toupper(symbol)) %>% 
  filter(symbol %in% c(geneOI1, unlist(geneOI2))) %>% 
  separate(treatment, into=c("cell","treatment"), sep="_") %>% 
  bind_rows(voom_hum) %>% 
  filter(symbol %in% c(geneOI1, unlist(geneOI2))) %>% 
  #Beautify labels
  mutate(mtb = recode(mtb, "Mtb"="+Mtb"),
         mtb = factor(mtb, levels=c("Media","+Mtb"))) %>% 
  mutate(group = paste(ptID,species, cell, treatment)) %>% 
  mutate(x=paste(treatment, cell, sep=" "),
         x=paste(mtb, x, sep="\n"),
         x=gsub("\nNA ","\n",x),
         x=gsub("coTB","coMtb",x)) %>% 
  #Blank panel of NOS2 human
  add_row(group = "AM10 Human AM NA", species="Human", symbol="NOS2",
          x="Media\nAM", value=NULL, mtb="Media") %>% 
  add_row(group = "AM10 Human AM NA", species="Human", symbol="NOS2",
          x="+Mtb\nAM", value=NULL, mtb="+Mtb") %>% 
  add_row(group = "AM10 Human AM NA", species="Human", symbol="NOS2",
          x="Media\nMDM", value=NULL, mtb="Media") %>% 
  add_row(group = "AM10 Human AM NA", species="Human", symbol="NOS2",
          x="+Mtb\nMDM", value=NULL, mtb="+Mtb") %>% 
  #Order x variable
  mutate(x=factor(x, levels=c("Media\nAM","+Mtb\nAM",
                              "Media\ncoMtb AM","+Mtb\ncoMtb AM",
                              "Media\nscBCG AM","+Mtb\nscBCG AM",
                              "Media\ncontrol AM","+Mtb\ncontrol AM",
                              "Media\nMDM","+Mtb\nMDM",
                              "Media\nBMDM","+Mtb\nBMDM")))

fc_color <- dat %>% 
  filter(symbol!="NOS2" | species != "Human") %>%  #remove blank place holder
  select(ptID, symbol, cell, species, treatment, mtb, value) %>% 
  pivot_wider(names_from = mtb) %>% 
  mutate(fc=`+Mtb`-Media,
         fc = ifelse(fc<0,"down","up")) %>% 
  select(ptID, symbol, cell,species, treatment, fc) 

dat <- left_join(dat, fc_color) %>% 
  select(group, species, symbol, libID, x, value, fc) %>% 
  #get rid of unpaired sample lines
  mutate(group = ifelse(grepl("Murine",group),NA, group))

# Linear models
fdr <- model_mur %>% 
  mutate(gene=toupper(gene)) %>% 
  bind_rows(model_hum) %>% 
  filter(gene %in% c(geneOI1, unlist(geneOI2))) %>% 
  rename(symbol=gene) %>% 
  separate(contrast_ref, into=c("cell1","mtb1"), sep=" ") %>% 
  separate(cell1, into=c("cell1","treatment1"), sep="_") %>% 
  separate(contrast_lvl, into=c("cell2","mtb2"), sep=" ") %>% 
  separate(cell2, into=c("cell2","treatment2"), sep="_") %>% 
  mutate(across(c(mtb1,mtb2), ~ recode(., "Mtb"="+Mtb"))) %>% 
  # start and end positions
  mutate(group1 = paste(treatment1, cell1, sep=" "),
         group1 = paste(mtb1, group1, sep="\n"),
         group1 = gsub("\nNA ","\n",group1),
         group1 = gsub(" NA","",group1),
         group1 = gsub("coTB","coMtb",group1),
         group2 = paste(treatment2, cell2, sep=" "),
         group2 = paste(mtb2, group2, sep="\n"),
         group2 = gsub("\nNA ","\n",group2),
         group2 = gsub(" NA","",group2),
         group2 = gsub("coTB","coMtb",group2)) %>%
  select(species, symbol, group1, group2, FDR) %>% 
  #Y position
  arrange(symbol, species, group1, group2) %>% 
  distinct(group1, group2, FDR, symbol, species) %>% 
  mutate(FDR = signif(FDR, digits=2)) %>% 
  mutate(signif = case_when(FDR < 0.01~"**",
                            FDR < 0.1~"*",
                            TRUE~"ns"))

#Mean murine lines
dat_mean <- dat %>% 
  filter(species=="Murine") %>% 
  group_by(species, symbol, x) %>% 
  summarise(value=mean(value)) %>% 
  ungroup() %>% 
  separate(x, into=c("tb","x"), sep="\n", extra = "merge") %>% 
  pivot_wider(names_from = tb) %>% 
  mutate(fc=`+Mtb`-Media,
         fc = ifelse(fc<0,"down","up"),
         fc=factor(fc)) %>% 
  pivot_longer(c(Media,`+Mtb`)) %>% 
  mutate(group=x) %>% 
  mutate(x=paste(name, x, sep="\n"))

#### Plot: MAF and IL10 ####
fdr_temp <- fdr %>% 
  filter(symbol%in%geneOI1) %>% 
  mutate(y.position=c(9.5,9,1.5,7.5,9.5,1.5,8.5,7,
                      7.5,10.5,10,5,10.5,12.5,11.5,9.8))

plot.ls <- list()
for(g in geneOI1){
  plot.ls[[g]] <-
    dat %>% filter(symbol==g) %>% 
    
    ggplot(aes(x=x,y=value)) +
    geom_point(color="black") +
    geom_line(data=dat %>% filter(species=="Human" & symbol==g),
              aes(group = group, color=fc)) +
    geom_line(data=dat_mean %>% filter(symbol==g),
              aes(group = group, color=fc)) +
    stat_pvalue_manual(data=fdr_temp %>% filter(symbol==g),
                       label = "signif", tip.length = 0.01, size=3) +
    scale_color_manual(values=c("down"="blue","up"="red"), drop=FALSE) +
    theme_bw(base_size = 10) +
    facet_grid(~symbol+species, scales="free_x", space="free") +
    labs(color="Fold change", y="Log2 CPM", x="") +
    theme(strip.background=element_rect(fill="white"),
          legend.position = "bottom",
          legend.direction = "horizontal")
}

p_all <- wrap_plots(plot.ls) + 
  plot_annotation(tag_levels = "A") + 
  plot_layout(guides = "collect", ncol=1)  &
  theme(legend.position='bottom')
p_all

ggsave("publication/temp/Fig5_gene_boxplot.pdf",
       p_all,width=5, height=5.5)

#### Other genes ####
#Y position of FDR
dat_lims <- dat %>% 
  filter(symbol %in% unlist(geneOI2)) %>% 
  group_by(species, symbol) %>% 
  summarise(minE=min(value),
            maxE=max(value),
            range=maxE-minE) %>% 
  ungroup()

fdr_temp <- fdr %>% 
  filter(symbol %in% unlist(geneOI2)) %>% 
  left_join(dat_lims) %>% 
  #rank FDR lines
  group_by(species, symbol) %>% 
  mutate(rk = as.numeric(as.factor(paste(group1,group2)))) %>% 
  ungroup() %>% 
  mutate(y.position = maxE+(rk*range/10)) %>% 
  mutate(symbol2 = recode(symbol,"LYZ"="LYZ / LYZ2",
                         "LYZ2"="LYZ / LYZ2"))

plot.ls.combo <- list()
for(group in names(geneOI2)){
  # print(group)
  gene.temp <- sort(geneOI2[[group]])
  gene.temp <- gene.temp[gene.temp %in% dat$symbol]
  plot.ls2 <- list()
  for(g in gene.temp){
    if(g=="LYZ2") { next }
    if(g=="LYZ"){ g <- c("LYZ","LYZ2") }
    p_temp <- dat %>% filter(symbol%in%g) %>% 
      mutate(symbol = recode(symbol,"LYZ"="LYZ / LYZ2",
                             "LYZ2"="LYZ / LYZ2")) %>% 
      
      ggplot(aes(x=x,y=value)) +
      geom_point(color="black") +
      geom_line(data=dat %>% filter(species=="Human" & symbol%in%g),
                aes(group = group, color=fc)) +
      geom_line(data=dat_mean %>% filter(symbol%in%g),
                aes(group = group, color=fc)) +
      stat_pvalue_manual(data=fdr_temp %>% filter(symbol%in%g),
                         label = "signif", tip.length = 0.01, size=3) +
      scale_color_manual(values=c("down"="blue","up"="red"), drop=FALSE) +
      theme_bw(base_size = 8) +
      facet_grid(~symbol2+species, scales="free_x", space = "free") +
      theme(strip.background=element_rect(fill="white"),
            legend.position = "bottom",
            legend.direction = "horizontal")
    
    # if(g[1]==gene.temp[1]){
      plot.ls2[[paste(group,g[1])]] <- p_temp +
        labs(color="Fold change", y="Log2 CPM", x="",
             title = group)
    # } else{
      # plot.ls2[[paste(group,g[1])]] <- p_temp +
      #   labs(color="Fold change", y="", x="")
    # }
  }
  # plot.ls.combo[[group]] <- wrap_plots(plot.ls2,nrow = 1)
  plot.ls.combo[[group]] <- plot.ls2
  
}

p_all2 <- plot.ls.combo[[1]][[1]] + plot.ls.combo[[1]][[2]] + plot.ls.combo[[1]][[3]] +
  plot.ls.combo[[1]][[4]] + plot.ls.combo[[1]][[5]] + 
  plot.ls.combo[[3]][[1]] + 
  plot.ls.combo[[5]][[1]] + plot.ls.combo[[5]][[2]] + plot.ls.combo[[5]][[3]] + 
  plot.ls.combo[[5]][[4]] + plot.ls.combo[[2]][[1]] + guide_area() +
  plot.ls.combo[[4]][[1]] + plot.ls.combo[[4]][[2]] + plot.ls.combo[[4]][[3]] + 
  # plot_annotation(tag_levels = list(c("A","","","","",
  #                                     "C","D","","B","E",
  #                                     "","","","","G",
  #                                     "","","","F",""))) +
  plot_layout(guides = "collect", ncol=3)  &
  theme(legend.position='bottom')
# p_all2

ggsave("publication/temp/FigS4_other_boxplots.pdf",
       p_all2,width=7.5, height=10)
