#Combine human and mouse msigDB
combine_msigdb <- function(category=NULL, subcategory=NULL){
  require(tidyverse)
  require(msigdbr)
  
  #Recode mouse to hgnc when available
  db.mouse <- msigdbr(species="mouse", 
                      category=category, subcategory=subcategory) %>% 
    distinct(gene_symbol, human_gene_symbol) %>% 
    rename(mouse_gene_symbol=gene_symbol)
  
  db.human <- msigdbr(species="human",
                      category=category, subcategory=subcategory) %>% 
    distinct(gene_symbol) %>% 
    rename(human_gene_symbol=gene_symbol)
  
  dat_all <- full_join(db.mouse, db.human, 
                       by = join_by(human_gene_symbol)) 
  
  return(dat_all)
}

#Select genes that are DEGs and GSEA leading edge
select_genes <- function(cell=NULL, db=NULL, pw=NULL,
                         FDR = 0.1, topLE=25){
  require(tidyverse)
  FDR.cutoff<-FDR
  
  #### GSEA leading edge ####
  cell_lower <- tolower(cell)
  
  allLE <- H_all %>% 
    filter(pathway%in%pw) %>% 
    filter(group2 != "control\nMurine AM") %>% 
    filter(grepl(cell, group2)) %>% 
    unnest(leadingEdge) %>% 
    pull(leadingEdge) %>% unique()
  
  #### Linear models ####
  load("results/lme_human.RData")
  load("results/lm_mouse.RData")
  
  allLE.DEG.H <- model_hum$lme.contrast %>%
    filter(grepl(cell,contrast_ref) & grepl(cell, contrast_lvl)) %>% 
    filter(FDR < FDR.cutoff) %>% 
    filter(gene %in% allLE)
  
  if(topLE<nrow(allLE.DEG.H)){
    allLE.DEG.H <- allLE.DEG.H  %>% 
      slice_min(order_by = FDR, n=topLE) %>% 
      pull(gene)
  } else{
    allLE.DEG.H <- allLE.DEG.H %>% 
      pull(gene)
  }
  
  if(cell=="AM"){
    allLE.DEG.M <- model_mur$lm.contrast %>% 
      filter((contrast_ref=="AM_coTB Media" & contrast_lvl=="AM_coTB Mtb") |
               (contrast_ref=="AM_scBCG Media" & contrast_lvl=="AM_scBCG Mtb")) %>% 
      filter(FDR < FDR.cutoff) %>% 
      filter(gene %in% allLE)
    
    if(topLE<nrow(allLE.DEG.M)){
      allLE.DEG.M <- allLE.DEG.M %>% 
        slice_min(order_by = FDR, n=topLE) 
    }else{
      allLE.DEG.M <- allLE.DEG.M %>% 
        pull(gene)
    }
    
  } else if(cell=="MDM"){
    allLE.DEG.M <- model_mur$lm.contrast %>% 
      filter(grepl(cell,contrast_ref) & grepl(cell, contrast_lvl)) %>% 
      filter(FDR < FDR.cutoff) %>% 
      filter(gene %in% allLE)
    
    if(topLE<nrow(allLE.DEG.M)){
      allLE.DEG.M <- allLE.DEG.M  %>% 
        slice_min(order_by = FDR, n=topLE) %>% 
        pull(gene)
    } else{
      allLE.DEG.M <- allLE.DEG.M %>% 
        pull(gene)
    }
  }
  
  
  #all selected genes
  allLE.DEG <- db %>% 
    filter(mouse_gene_symbol %in% allLE.DEG.M |
             human_gene_symbol %in% allLE.DEG.H)
  return(allLE.DEG)
}

#Plot heatmaps of select genes
gene_heatmaps <- function(cell=NULL, genes=NULL, pw=NULL,
                          category=NULL, subcategory=NULL,
                          font_size = 8, row_splits = 4,
                          cluster_by = "E",
                          pw_anno = "none"){
  require(tidyverse)
  require(limma)
  require(ComplexHeatmap)
  require(viridis)
  cell.upper <- toupper(cell)
  hm.ls <- list()
  cellOI <- cell
  
  #### Gene expression ####
  load("data_clean/human_dat_clean.RData")
  load("data_clean/mouse_dat_clean.RData")
  
  #Human expression
  dat.human <- as.data.frame(voom_hum$E) %>%
    rownames_to_column("human_gene_symbol") %>%
    left_join(db, by="human_gene_symbol") %>% 
    filter(human_gene_symbol %in% genes$human_gene_symbol | 
             mouse_gene_symbol %in% genes$mouse_gene_symbol) %>% 
    pivot_longer(-c(human_gene_symbol,mouse_gene_symbol),
                 names_to = "libID") %>%
    left_join(voom_hum$targets, by="libID") %>%
    filter(cell==cellOI) %>% 
    select(human_gene_symbol, mouse_gene_symbol,
           libID, cell, mtb, value) %>%
    mutate(species="Human") %>%
    mutate(mtb = recode(mtb, "Mtb"="+Mtb"),
           mtb = factor(mtb, levels=c("Media","+Mtb"))) 
  
  #Mouse expression
  dat.mur <- as.data.frame(voom_mur$E) %>%
    rownames_to_column("mouse_gene_symbol") %>%
    left_join(db, by="mouse_gene_symbol") %>%
    filter(human_gene_symbol %in% genes$human_gene_symbol |
             mouse_gene_symbol %in% genes$mouse_gene_symbol) %>%
    pivot_longer(-c(human_gene_symbol,mouse_gene_symbol),
                 names_to = "libID") %>%
    left_join(voom_mur$targets, by="libID") %>%
    filter(grepl(cellOI, treatment)) %>% 
    mutate(cell=cell.upper)%>%  
    select(human_gene_symbol, mouse_gene_symbol,
           libID, cell, treatment, mtb, value) %>%
    mutate(species="Murine") %>%
    mutate(mtb = recode(mtb, "Mtb"="+Mtb"),
           mtb = factor(mtb, levels=c("Media","+Mtb"))) %>% 
    mutate(treatment = gsub("AM_","",treatment))
  
  #Combine human and mouse
  dat_all <- bind_rows(dat.human, dat.mur) %>% 
    mutate(group2 = case_when(
      is.na(treatment) ~ paste(species, cell, sep="\n"),
      treatment=="BMDM" ~ paste(species, treatment, sep="\n"),
      !is.na(treatment) ~ paste0(treatment, "\n", species, " ", cell)))
  
  #### Format data ####
  dat.E.temp <- dat_all %>%
    drop_na(human_gene_symbol, mouse_gene_symbol) %>% 
    mutate(all_symbol = ifelse(
      human_gene_symbol==toupper(mouse_gene_symbol), 
      human_gene_symbol, 
      paste(human_gene_symbol,mouse_gene_symbol,sep="/"))) %>% 
    group_by(group2, all_symbol, species, cell, treatment, mtb) %>% 
    summarise(meanE=mean(value, na.rm = TRUE),
              .groups = "drop") %>% 
    ungroup() %>% 
    arrange(all_symbol) %>% 
    mutate(name = paste(group2,mtb, sep="\n"))
  
  #expression
  dat.E.mat <- dat.E.temp %>% 
    select(all_symbol, name, meanE) %>%
    # mutate(name = gsub(" ","\n",name)) %>% 
    mutate(name = factor(name, levels=c("Human\nAM\nMedia", 
                                        "Human\nAM\n+Mtb",
                                        "coTB\nMurine AM\nMedia",
                                        "coTB\nMurine AM\n+Mtb",
                                        "scBCG\nMurine AM\nMedia", 
                                        "scBCG\nMurine AM\n+Mtb",
                                        "control\nMurine AM\nMedia",
                                        "control\nMurine AM\n+Mtb",
                                        "Human\nMDM\nMedia", 
                                        "Human\nMDM\n+Mtb",
                                        "Murine\nBMDM\nMedia", 
                                        "Murine\nBMDM\n+Mtb"
    ))) %>% 
    arrange(name) %>% 
    pivot_wider(names_from = name, values_from = meanE) %>% 
    #remove non-shared genes
    drop_na() %>% 
    arrange(all_symbol) %>% 
    column_to_rownames("all_symbol") %>% 
    as.matrix()
  
  #Calculate fold change
  dat.FC.mat <- dat.E.temp %>% 
    select(all_symbol, group2, mtb, meanE) %>%
    # arrange(s) %>% 
    #Calculate fold change
    pivot_wider(names_from=mtb, values_from = meanE) %>% 
    mutate(FC=`+Mtb`-Media) %>% 
    select(-`+Mtb`,-Media) %>% 
    #Pivot wider
    mutate(group2 = factor(group2, levels=c("Human\nAM",
                                            "Human\nMDM",
                                            "coTB\nMurine AM",
                                            "scBCG\nMurine AM",
                                            "control\nMurine AM",
                                            "Murine\nBMDM"))) %>%
    arrange(group2) %>% 
    pivot_wider(names_from = group2, values_from = FC) %>% 
    #remove non-shared genes
    drop_na() %>% 
    arrange(all_symbol) %>% 
    column_to_rownames("all_symbol") %>% 
    as.matrix()
  
  #### Heatmap params ####
  anno_ls <- strsplit(colnames(dat.E.mat), split = "\n")
  anno_ls2 <- strsplit(colnames(dat.FC.mat), split = "\n")
  if(cell.upper=="AM"){
    # anno_ls[[1]][4] <- anno_ls[[1]][3]
    # anno_ls[[1]][3] <- " "
    # anno_ls[[2]][4] <- anno_ls[[2]][3]
    # anno_ls[[2]][3] <- " "
    # anno_ls2[[1]][4] <- anno_ls2[[1]][3]
    # anno_ls2[[1]][3] <- " "
    
    column_names <- paste(unlist(lapply(anno_ls, `[[`, 1)),
                          unlist(lapply(anno_ls, `[[`, 2)),
                          unlist(lapply(anno_ls, `[[`, 3)),sep="\n")
    column_names2 <- paste(unlist(lapply(anno_ls2, `[[`, 1)),
                           unlist(lapply(anno_ls2, `[[`, 2)),sep="\n")
    #Order column names
    column_names2 <- factor(column_names2, levels=c(
      "Human\nAM", "Human\nMDM",
      "coTB\nMurine AM", "scBCG\nMurine AM",
      "control\nMurine AM","Murine\nMDM"))
    column_names <- gsub("\nMedia|\n[+]Mtb", "", column_names)
    column_names <- factor(column_names, levels=c(
      "Human\nAM", "Human\nMDM",
      "coTB\nMurine AM", "scBCG\nMurine AM",
      "control\nMurine AM","Murine\nMDM"))
    
    treat_label <- 3
    cell_width <- 5
  } else if(cell.upper=="MDM"){
    column_names <- paste(unlist(lapply(anno_ls, `[[`, 1)),
                          unlist(lapply(anno_ls, `[[`, 2)),sep="\n")
    column_names2 <- paste(unlist(lapply(anno_ls2, `[[`, 1)),
                           unlist(lapply(anno_ls2, `[[`, 2)),sep="\n")
    treat_label <- 3
    cell_width <- 2.5
  }
  
  col_ha <- HeatmapAnnotation(
    Infection = factor(unlist(lapply(anno_ls, `[[`, treat_label)), 
                       level=c("Media","+Mtb")),
    col = list(Infection = c("Media"="grey70",
                             "+Mtb"="grey40")),
    # show_legend = legs,
    annotation_name_gp = gpar(fontsize = font_size), 
    annotation_legend_param = list(
      Infection = list(nrow = 2, title_position = "lefttop",
                       title_gp = grid::gpar(fontsize = font_size),
                       legend_label_gp = grid::gpar(fontsize = font_size),
                       labels_gp = gpar(fontsize = font_size))),
    simple_anno_size = unit(0.2, "cm"),
    annotation_height = unit(0.1, "cm"),
    annotation_width = unit(0.1, "cm"))
  
  plot_title <- "log2 CPM"
  plot_title2 <- "log2 fold change\nMtb - Media"
  
  #### row annotation ####
  if(pw_anno!="none"){
    require(msigdbr)
    db.mouse <- msigdbr(species="mouse", 
                        category=category, subcategory=subcategory) %>% 
      filter(gs_name %in% pw) %>% 
      distinct(gs_name, gene_symbol, human_gene_symbol) %>% 
      rename(mouse_gene_symbol=gene_symbol)
    
    db.human <- msigdbr(species="human",
                        category=category, subcategory=subcategory) %>% 
      filter(gs_name %in% pw) %>% 
      distinct(gs_name, gene_symbol) %>% 
      rename(human_gene_symbol=gene_symbol)
    
    db.all <- full_join(db.mouse, db.human, 
                        by = join_by(gs_name, human_gene_symbol)) %>% 
      #clean db names
      mutate(gs_name_orig = factor(gs_name, levels=pw),
             gs_name = gsub("HALLMARK_|GOBP_","",gs_name),
             gs_name = gsub("_"," ", gs_name),
             gs_name = tolower(gs_name),
             gs_name = recode(gs_name,
                              "tnfa signaling via nfkb"="TNFA",
                              "interferon gamma response"="IFNG",
                              "interferon alpha response"="IFNA",
                              "inflammatory response"="inflammatory",
                              "il6 jak stat3 signaling"="IL6",
                              "g2m checkpoint"="G2M",
                              "myc targets v1"="MYC"
                              ),
             gs_name = fct_reorder(gs_name, gs_name_orig, 
                                   .fun = unique),
             gs_name = factor(gs_name, levels=c(levels(gs_name), "z"))) %>% 
      mutate(value=gs_name) %>% 
      select(-gs_name_orig) %>% 
      arrange(gs_name) %>% 
      # mutate(gs_name = as.numeric(as.factor(gs_name))) %>% 
      pivot_wider(names_from = gs_name) %>% 
      mutate(all_symbol = ifelse(
        human_gene_symbol==toupper(mouse_gene_symbol),
        human_gene_symbol,
        paste(human_gene_symbol,mouse_gene_symbol,sep="/"))) %>% 
      select(-mouse_gene_symbol, -human_gene_symbol) %>% 
      select(all_symbol, everything())
    
    col_list <- list()
    col_vec <- c("#332288","#88CCEE","#DDCC77","#CC6677","#AA4499")
    # col_vec <- rep("grey10",10)
    for (i in 1:(ncol(db.all)-1)){
      vec.temp <- c(col_vec[i], "white")
      pw_name <- unique(db.all[,i+1])
      pw_name <- pw_name[!is.na(pw_name)]
      # print(pw_name)
      names(vec.temp) <- c(as.character(pw_name), "z")
      col_list[[as.character(pw_name)]] <- vec.temp
    }
    
  } else{
    db.all <- NULL
  }
  
  if(cluster_by=="FC"){
    #### Heatmaps order by FC ####
    ##### Fold change heatmap #####
    limit <- max(abs(min(dat.FC.mat)),abs(max(max(dat.FC.mat))))
    
    if(pw_anno %in% c("all","FC")){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.FC.mat)) %>% 
        arrange(match(all_symbol, rownames(dat.FC.mat))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
    } else{ row_anno <- NULL }
    
    if(pw_anno == "none") {
      pretty_pw <- gsub("HALLMARK_|GOBP_","",pw)
      pretty_pw <- gsub("_"," ", pretty_pw)
      pretty_pw <- str_to_title(pretty_pw)
      rtitle <- pretty_pw
    } else {rtitle <- NULL}
    
    hm2 <- Heatmap(dat.FC.mat, name = plot_title2,
                   col = circlize::colorRamp2(
                     c(-limit, 0, limit),
                     c("darkblue", "white", "red")),
                   heatmap_legend_param = list(
                     direction = "horizontal",
                     title_position = "lefttop",
                     title_gp = grid::gpar(fontsize = font_size),
                     legend_label_gp = grid::gpar(fontsize = font_size),
                     labels_gp = gpar(fontsize = font_size)
                   ),
                   #rows
                   row_title=rtitle, 
                   row_title_gp = grid::gpar(fontsize = font_size+2),
                   cluster_rows = TRUE,
                   show_row_dend = FALSE,
                   row_split=row_splits,
                   row_names_gp = grid::gpar(fontsize = font_size),
                   left_annotation = row_anno,
                   width = unit(cell_width, "cm"),
                   #columns
                   column_title_side = "bottom",
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   # bottom_annotation = col_ha,
                   column_split = column_names2,
                   column_title_gp = grid::gpar(fontsize = font_size))
    
    hm2_temp <- draw(hm2)
    
    #### Expression heatmap ####
    # Reorder by fold change heatmap
    row_hm2 <- row_order(hm2_temp)
    dat.E.mat.ord <- dat.E.mat[unlist(row_hm2),]
    row_splits <- c()
    for(i in length(row_hm2):1){
      row_splits <- c(rep(i,length(row_hm2[[i]])), row_splits)
    }
    
    if(pw_anno %in% c("all","E")){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.E.mat.ord)) %>% 
        arrange(match(all_symbol, rownames(dat.E.mat.ord))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
    } else{ row_anno <- NULL }
    
    hm <- Heatmap(dat.E.mat.ord, name = plot_title, 
                  col = circlize::colorRamp2(
                    c(min(dat.E.mat), mean(dat.E.mat), max(dat.E.mat)),
                    c("#440154FF", "#21908CFF", "#FDE725FF")),
                  heatmap_legend_param = list(
                    direction = "horizontal",
                    title_position = "lefttop",
                    title_gp = grid::gpar(fontsize = font_size),
                    legend_label_gp = grid::gpar(fontsize = font_size),
                    labels_gp = gpar(fontsize = font_size)
                  ),
                  #rows
                  row_title=rtitle,
                  row_title_gp = grid::gpar(fontsize = font_size+2),
                  left_annotation = row_anno,
                  cluster_rows = FALSE,
                  show_row_dend = FALSE,
                  row_split=row_splits,
                  row_names_gp = grid::gpar(fontsize = font_size),
                  width = unit(cell_width, "cm"),
                  #columns
                  column_title_side = "bottom",
                  show_column_names = FALSE,
                  cluster_columns = FALSE,
                  bottom_annotation = col_ha,
                  column_split = column_names,
                  column_title_gp = grid::gpar(fontsize = font_size)) 
    
    if(pw_anno=="separate"){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.E.mat.ord)) %>% 
        arrange(match(all_symbol, rownames(dat.E.mat.ord))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        labels = anno_text(rownames(dat.E.mat.ord), which="row", 
                           show_name=TRUE,
                           gp=gpar(fontsize = font_size)),
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
      # hm_anno <- draw(row_anno)
    } 
  } else if(cluster_by=="E"){
    #### Heatmaps order by E ####
    ##### Expression heatmap #####
    if(pw_anno %in% c("all","E")){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.E.mat)) %>% 
        arrange(match(all_symbol, rownames(dat.E.mat))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
    } else{ row_anno <- NULL }
    
    if(pw_anno == "none") {
      pretty_pw <- gsub("HALLMARK_|GOBP_","",pw)
      pretty_pw <- gsub("_"," ", pretty_pw)
      pretty_pw <- str_to_title(pretty_pw)
      rtitle <- pretty_pw
    } else {rtitle <- NULL}
    
    hm <- Heatmap(dat.E.mat, name = plot_title, 
                  col = circlize::colorRamp2(
                    c(min(dat.E.mat), mean(dat.E.mat), max(dat.E.mat)),
                    c("#440154FF", "#21908CFF", "#FDE725FF")),
                  heatmap_legend_param = list(
                    direction = "horizontal",
                    title_position = "lefttop",
                    title_gp = grid::gpar(fontsize = font_size),
                    legend_label_gp = grid::gpar(fontsize = font_size),
                    labels_gp = gpar(fontsize = font_size)
                  ),
                  #rows
                  row_title=rtitle,
                  row_title_gp = grid::gpar(fontsize = font_size+2),
                  left_annotation = row_anno,
                  cluster_rows = TRUE,
                  show_row_dend = FALSE,
                  row_split=row_splits,
                  row_names_gp = grid::gpar(fontsize = font_size),
                  width = unit(cell_width, "cm"),
                  #columns
                  column_title_side = "bottom",
                  show_column_names = FALSE,
                  cluster_columns = FALSE,
                  bottom_annotation = col_ha,
                  column_split = column_names,
                  column_title_gp = grid::gpar(fontsize = font_size)) 
    
    hm_temp <- draw(hm)
    
    ##### Fold change heatmap #####
    # Reorder by expression heatmap
    row_hm <- row_order(hm_temp)
    dat.FC.mat.ord <- dat.FC.mat[unlist(row_hm),]
    row_splits <- c()
    for(i in length(row_hm):1){
      row_splits <- c(rep(i,length(row_hm[[i]])), row_splits)
    }
    
    if(pw_anno %in% c("all","FC")){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.FC.mat.ord)) %>% 
        arrange(match(all_symbol, rownames(dat.FC.mat.ord))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
    } else{ row_anno <- NULL }
    
    limit <- max(abs(min(dat.FC.mat)),abs(max(max(dat.FC.mat))))
    hm2 <- Heatmap(dat.FC.mat.ord, name = plot_title2, 
                   col = circlize::colorRamp2(
                     c(-limit, 0, limit),
                     c("darkblue", "white", "red")),
                   heatmap_legend_param = list(
                     direction = "horizontal",
                     title_position = "lefttop",
                     title_gp = grid::gpar(fontsize = font_size),
                     legend_label_gp = grid::gpar(fontsize = font_size),
                     labels_gp = gpar(fontsize = font_size)
                   ),
                   #rows
                   row_title=rtitle,
                   row_title_gp = grid::gpar(fontsize = font_size+2),
                   left_annotation = row_anno,
                   cluster_rows = FALSE,
                   show_row_dend = FALSE,
                   row_split=row_splits,
                   row_names_gp = grid::gpar(fontsize = font_size),
                   width = unit(cell_width, "cm"),
                   #columns
                   column_title_side = "bottom",
                   show_column_names = FALSE,
                   cluster_columns = FALSE,
                   # bottom_annotation = col_ha,
                   column_split = column_names2,
                   column_title_gp = grid::gpar(fontsize = font_size))
    
    if(pw_anno=="separate"){
      db.ord <- db.all %>% 
        filter(all_symbol %in% rownames(dat.FC.mat.ord)) %>% 
        arrange(match(all_symbol, rownames(dat.FC.mat.ord))) %>% 
        column_to_rownames("all_symbol")
      db.ord[is.na(db.ord)] <- "z"
      
      row_anno <- rowAnnotation(
        df = db.ord,
        col = col_list,
        labels = anno_text(rownames(dat.FC.mat.ord), which="row", 
                           show_name=TRUE,
                           gp=gpar(fontsize = font_size)),
        show_legend = FALSE,
        show_annotation_name = TRUE, 
        annotation_name_gp = gpar(fontsize = font_size),
        annotation_legend_param = list(
          title_gp = gpar(fontsize = font_size),  
          labels_gp = gpar(fontsize = font_size)  
        ))
      # hm_anno <- draw(row_anno)
    } 
  } else{
    stop("cluster_by must be one of 'E' for expression or 'FC' for fold change")
  }
  
  #### Save ####
  hm.ls[[paste(pw, collapse=" ")]] <- grid.grabExpr(draw(
    hm, 
    # column_title="Expression",
    column_title_gp = grid::gpar(fontsize = font_size+2),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends=TRUE))
  
  hm.ls[[paste0(paste(pw, collapse=" "),"_FC")]] <- grid.grabExpr(draw(
    hm2, 
    # column_title="Fold Change\n",
    column_title_gp = grid::gpar(fontsize = font_size+2),
    heatmap_legend_side = "bottom",
    annotation_legend_side = "bottom",
    merge_legends=TRUE))
  
  if(pw_anno=="separate"){
    hm.ls[[paste0(pw,"_anno")]] <- grid.grabExpr(draw(
      row_anno))
  }
  return(hm.ls)
}
