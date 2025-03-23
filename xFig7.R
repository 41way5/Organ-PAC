#xFig7------------------------------------------------------------------------------------------------
# Corr between PAC and brain structure

{
  ukbb_index <- read_excel("Project/Proteomic aging/Manuscript/Supplementary File.xlsx",sheet="ST5") %>% 
    mutate(Category = if_else(Category == "", NA_character_, Category)) %>%
    fill(Category)
  
  ukbb_trait_list_imaging <- ukbb_index %>% filter(Modality=="Imaging visit") %>%
    filter(Category %in% c("Brain MRI")) %>% 
    mutate(fid=paste0(`UKB Field ID`,"-2.0")) %>% filter(!is.na(Abbr))
  
  ukbb_trait_list_imaging_wrt <- ukbb_index %>% filter(Modality=="Imaging visit") %>%
    filter(Category %in% c("Brain MRI")) %>% 
    mutate(fid=paste0(`UKB Field ID`,"-2.0")) %>% filter(is.na(Abbr)) %>% filter(!is.na(Note)) %>% 
    filter(!grepl("L1|L2|L3", Phenotype)) %>% 
    mutate(wrt=str_extract(Phenotype, "(?<=Weighted-mean\\s)\\w+"),
           lr=str_extract(Phenotype, "(?<=\\()\\w+"),
           lr=ifelse(lr=="left","L","R"),
           region=str_extract(Phenotype, "(?<=tract\\s)([^()]+|$)"),
           rg=toupper(gsub("\\b(\\w)\\w*\\s*", "\\1", region)),
           rg=case_when(rg=="CGPOC"~"CGC",
                        rg=="IF-OF"~"IFO",
                        T~rg
           ),
           rg=ifelse(region=="forceps major","FMA",rg),rg=ifelse(region=="forceps minor","FMI",rg),
           tract=paste0(rg,"_",lr),
           tract=ifelse(is.na(lr),rg,tract),
           full=paste0(wrt,"_",tract)
    ) %>% filter(rg!="PPOC") %>% group_by(wrt,tract)
  
  
  
  #Baseline trait name list (imaging)
  trait_base_imaging <- ukbb_full %>% dplyr::select(eid,ukbb_trait_list_imaging$fid) %>%
    rename_with(~ as.character(ukbb_trait_list_imaging$Abbr), .cols = -eid) %>% dplyr::select(-eid) %>% colnames()
  
  trait_base_imaging2 <- ukbb_full %>% dplyr::select(eid,ukbb_trait_list_imaging$fid) %>%
    rename_with(~ as.character(ukbb_trait_list_imaging$Phenotype), .cols = -eid) %>% dplyr::select(-eid) %>% colnames()
  
  
  ukbb_trait_base_imaging <- ukbb_full %>% dplyr::select(eid,ukbb_trait_list_imaging$fid) %>% 
    mutate_at(vars(-eid), scale) %>% 
    mutate_at(vars(-eid), as.numeric)
  # mutate(`20023-0.0`=as.numeric(scale(`20023-0.0`)),
  #        `4282-0.0`=as.numeric(scale(`4282-0.0`)),
  #        `20016-0.0`=as.numeric(scale(`20016-0.0`)),
  #        `399-0.1`=as.numeric(scale(`399-0.1`))
  # )
  
  ukbb_trait_imaging <-  ukb_olink %>% dplyr::select(eid,age_at_recruitment,c("age","sex","ethnic","edu","tdi","center","pa"),organ_pac,"phenoage_res","kdm_res") %>% 
    left_join(ukbb_trait_base_imaging,by="eid")
  
  
  
  #Baseline trait name list (imaging) for White tract region
  trait_base_imaging_wrt <- ukbb_full %>% dplyr::select(eid,ukbb_trait_list_imaging_wrt$fid) %>%
    rename_with(~ as.character(ukbb_trait_list_imaging_wrt$full), .cols = -eid) %>% dplyr::select(-eid) %>% colnames()
  
  ukbb_trait_base_imaging_wrt <- ukbb_full %>% dplyr::select(eid,ukbb_trait_list_imaging_wrt$fid) %>% 
    mutate_at(vars(-eid), scale) %>% 
    mutate_at(vars(-eid), as.numeric)
  # mutate(`20023-0.0`=as.numeric(scale(`20023-0.0`)),
  #        `4282-0.0`=as.numeric(scale(`4282-0.0`)),
  #        `20016-0.0`=as.numeric(scale(`20016-0.0`)),
  #        `399-0.1`=as.numeric(scale(`399-0.1`))
  # )
  
  ukbb_trait_imaging_wrt <-  ukb_olink %>% dplyr::select(eid,age_at_recruitment,c("age","sex","ethnic","edu","tdi","center","pa"),organ_pac,"phenoage_res","kdm_res") %>% 
    left_join(ukbb_trait_base_imaging_wrt,by="eid")
  
  
  # qc_id$FID
  # 
  # qc_id %>%
  #   group_by(FID) %>%
  #   filter(n() > 1) %>%
  #   distinct(FID)
  # 
  # 
  # qc_id %>%
  #   distinct(FID, .keep_all = TRUE)
  # 
  # write.table(qc_id %>%
  #               distinct(FID, .keep_all = TRUE), file = "./Data/pac.phe2", row.names = FALSE, quote = FALSE, sep = "\t")
  # 
  # colSums(is.na(qc_id))
  
  
  describe(ukbb_trait_imaging$`25851-2.0`)
  describe(ukbb_trait_imaging_wrt$`25703-2.0`)
  
  
  #Grey matter
  beta_matrix_imaging <- calculate_beta_matrix(ukbb_trait_imaging, 
                                               var_list1=c(organ_pac,"phenoage_res","kdm_res"), 
                                               var_list2=ukbb_trait_list_imaging$fid, 
                                               covariates = c("age","sex","ethnic","edu","tdi","center","pa")
  )
  
  
  colnames(beta_matrix_imaging$beta_matrix) <- trait_base_imaging
  rownames(beta_matrix_imaging$beta_matrix) <- c(organ_list,"Phenoage gap","KDMage gap")
  colnames(beta_matrix_imaging$fdr_p_values) <- trait_base_imaging
  rownames(beta_matrix_imaging$fdr_p_values) <- c(organ_list,"Phenoage gap","KDMage gap")
  
  
  
  #p1
  mat_beta1 <- (beta_matrix_imaging$beta_matrix)
      colnames(mat_beta1) <- trait_base_imaging2
      rownames(mat_beta1) <- c(organ_list,"Phenoage gap","KDMage gap")
  
  mat_beta1<- mat_beta1[organ_list, 1:96]
  
  
  mat_p1 <- (beta_matrix_imaging$fdr_p_values)
    colnames(mat_p1) <- trait_base_imaging2
    rownames(mat_p1) <- c(organ_list,"Phenoage gap","KDMage gap")
  mat_p1<- mat_p1[organ_list, 1:96]
  
  
  #p2
  mat_beta <- (beta_matrix_imaging$beta_matrix)
  #Subset the summary results
  mat_beta <- mat_beta[organ_list, 97:144]
  
  mat_p <- (beta_matrix_imaging$fdr_p_values)
  mat_p <- mat_p[organ_list, 97:144]
  
  
  
  
  
  
  #White matter
  beta_matrix_imaging_wrt <- calculate_beta_matrix(ukbb_trait_imaging_wrt, 
                                                   var_list1=c(organ_pac,"phenoage_res","kdm_res"), 
                                                   var_list2=ukbb_trait_list_imaging_wrt$fid, 
                                                   covariates = c("age","sex","ethnic","edu","tdi","center","pa")
  )
  
  
  colnames(beta_matrix_imaging_wrt$beta_matrix) <- trait_base_imaging_wrt
  rownames(beta_matrix_imaging_wrt$beta_matrix) <- c(organ_list,"Phenoage gap","KDMage gap")
  colnames(beta_matrix_imaging_wrt$fdr_p_values) <- trait_base_imaging_wrt
  rownames(beta_matrix_imaging_wrt$fdr_p_values) <- c(organ_list,"Phenoage gap","KDMage gap")
  
  #Brain
  mat_beta_wrt <- (beta_matrix_imaging_wrt$beta_matrix)
  mat_beta_wrt <- mat_beta_wrt["Brain",] %>% t()
  mat_beta_wrt <- as.data.frame(mat_beta_wrt) %>%  
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>% 
    separate(variable, into = c("v1", "v2"), sep = "_", extra = "merge")
  
  mat_beta_wrt <- xtabs(value ~ v1 + v2, mat_beta_wrt) %>% as.matrix()
  mat_beta_wrt_brain <- mat_beta_wrt[-5,]
  
  
  mat_p_wrt <- (beta_matrix_imaging_wrt$fdr_p_values)
  mat_p_wrt <- mat_p_wrt["Brain",] %>% t()
  mat_p_wrt <- as.data.frame(mat_p_wrt) %>%  
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>% 
    separate(variable, into = c("v1", "v2"), sep = "_", extra = "merge")
  
  mat_p_wrt <- xtabs(value ~ v1 + v2, mat_p_wrt) %>% as.matrix()
  mat_p_wrt_brain <- mat_p_wrt[-5,]
  
  
  #Organismal
  mat_beta_wrt <- (beta_matrix_imaging_wrt$beta_matrix)
  mat_beta_wrt <- mat_beta_wrt["Organismal",] %>% t()
  mat_beta_wrt <- as.data.frame(mat_beta_wrt) %>%  
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>% 
    separate(variable, into = c("v1", "v2"), sep = "_", extra = "merge")
  
  mat_beta_wrt <- xtabs(value ~ v1 + v2, mat_beta_wrt) %>% as.matrix()
  mat_beta_wrt_organ <- mat_beta_wrt[-5,]
  
  
  mat_p_wrt <- (beta_matrix_imaging_wrt$fdr_p_values)
  mat_p_wrt <- mat_p_wrt["Organismal",] %>% t()
  mat_p_wrt <- as.data.frame(mat_p_wrt) %>%  
    pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>% 
    separate(variable, into = c("v1", "v2"), sep = "_", extra = "merge")
  
  mat_p_wrt <- xtabs(value ~ v1 + v2, mat_p_wrt) %>% as.matrix()
  mat_p_wrt_organ <- mat_p_wrt[-5,]
  
  

  
  p_imaging_cortical <-grid.grabExpr(draw( 
    Heatmap(mat_beta1,
            name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
            col = colorRamp2(c(-0.2,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
            show_row_names = TRUE, 
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 3.5),
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot=45,
            #column_split=ukbb_trait_list$Categroy_short,
            column_title_gp = gpar(fontsize = 5),
            #row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
            # column_split = factor(c(rep("Subcortical regions GMV", each = 15),
            #                         rep("Cerebellar regions GMV", each = 28),
            #                         rep("Total brain", each = 5)
            # ),
            # levels = c("Subcortical regions GMV","Cerebellar regions GMV","Total brain")),
            column_title_side = "top",
            row_names_side = "right",
            row_title = NULL,
            #column_title = NULL,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if((mat_p1 )[i, j] < 0.001) {
                grid.text("***", x, y,gp = gpar(fontsize = 3))
              } else if((mat_p1 )[i, j] < 0.01) {
                grid.text("**", x, y,gp = gpar(fontsize = 3))
              } else if((mat_p1 )[i, j] < 0.05) {
                grid.text("*", x, y,gp = gpar(fontsize = 3))
              }
            },
            heatmap_legend_param = list(
              #at = c(-2, 0, 2),
              #labels = c("0", "10%", "20%", "30%", "40%"),
              title = "",
              legend_height = unit(1, "cm"),
              grid_width = unit(0.5, "mm"),
              #legend_width =unit(1, "mm"),
              #title_position = "lefttop-rot",
              legend_title_gp = gpar(fontsize = 5), 
              labels_gp = gpar(col = "black", fontsize = 5)
              
            ),
            width = ncol(mat_beta)*unit(2.8, "mm"), 
            height =nrow(mat_beta)*unit(2.2, "mm")
            
    ),padding = unit(c(0, 0, 0, 0), "mm")
  )
  )
  
  
  plot_grid(p_imaging_cortical)
  
  
  p_imaging <-grid.grabExpr(draw( 
    Heatmap(mat_beta,
            name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
            col = colorRamp2(c(-0.2,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
            show_row_names = TRUE, 
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 5),
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot=45,
            #column_split=ukbb_trait_list$Categroy_short,
            column_title_gp = gpar(fontsize = 5),
            #row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
            column_split = factor(c(rep("Subcortical regions GMV", each = 15),
                                    rep("Cerebellar regions GMV", each = 28),
                                    rep("Total brain", each = 5)
            ),
            levels = c("Subcortical regions GMV","Cerebellar regions GMV","Total brain")),
            column_title_side = "top",
            row_names_side = "right",
            row_title = NULL,
            #column_title = NULL,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if((mat_p)[i, j] < 0.001) {
                grid.text("***", x, y,gp = gpar(fontsize = 4))
              } else if((mat_p)[i, j] < 0.01) {
                grid.text("**", x, y,gp = gpar(fontsize = 4))
              } else if((mat_p)[i, j] < 0.05) {
                grid.text("*", x, y,gp = gpar(fontsize = 4))
              }
            },
            heatmap_legend_param = list(
              #at = c(-2, 0, 2),
              #labels = c("0", "10%", "20%", "30%", "40%"),
              title = "",
              legend_height = unit(1, "cm"),
              grid_width = unit(0.5, "mm"),
              #legend_width =unit(1, "mm"),
              #title_position = "lefttop-rot",
              legend_title_gp = gpar(fontsize = 5), 
              labels_gp = gpar(col = "black", fontsize = 5)
              
            ),
            width = ncol(mat_beta)*unit(2.8, "mm"), 
            height =nrow(mat_beta)*unit(2.2, "mm")
            
    ),padding = unit(c(0, 0, 0, 0), "mm")
  )
  )
  
  plot_grid(p_imaging)
  
  
  
  
  #Panel for white matter tracts
  p_imaging_wrt <-
    Heatmap(mat_beta_wrt_brain,
            name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
            col = colorRamp2(c(-0.1,0,0.1), c("#4DBBD5B2","white","#DC0000B2")),
            show_row_names = TRUE, 
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 5),
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot=45,
            #column_split=ukbb_trait_list$Categroy_short,
            column_title_gp = gpar(fontsize = 5),
            #row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
            # column_split = factor(c(rep("Subcortical regions GMV", each = 15),
            #                         rep("Cerebellar regions GMV", each = 28),
            #                         rep("Total brain", each = 5)),
            # levels = c("Subcortical regions GMV","Cerebellar regions GMV","Total brain")),
            column_title_side = "top",
            row_names_side = "right",
            row_title = NULL,
            column_title = "Brain age gap",
            cell_fun = function(j, i, x, y, w, h, fill) {
              if((mat_p_wrt_brain)[i, j] < 0.001) {
                grid.text("***", x, y,gp = gpar(fontsize = 5))
              } else if((mat_p_wrt_brain)[i, j] < 0.01) {
                grid.text("**", x, y,gp = gpar(fontsize = 5))
              } else if((mat_p_wrt_brain)[i, j] < 0.05) {
                grid.text("*", x, y,gp = gpar(fontsize = 5))
              }
            },
            heatmap_legend_param = list(
              #at = c(-2, 0, 2),
              #labels = c("0", "10%", "20%", "30%", "40%"),
              title = NULL,
              legend_height = unit(1, "cm"),
              grid_width = unit(0.5, "mm"),
              #legend_width =unit(1, "mm"),
              title_position = "leftcenter-rot",
              title_gp=gpar(col = "black", fontsize = 5),
              labels_gp = gpar(col = "black", fontsize = 5),
              show_legend=F
              
            ),
            width = ncol(mat_beta)*unit(1.5, "mm"), 
            height =nrow(mat_beta)*unit(1.5, "mm")
            
    )
  
  p_imaging_wrt2 <-
    Heatmap(mat_beta_wrt_organ,
            name = "beta2",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
            col = colorRamp2(c(-0.1,0,0.1), c("#4DBBD5B2","white","#DC0000B2")),
            show_row_names = TRUE, 
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 5),
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot=45,
            #column_split=ukbb_trait_list$Categroy_short,
            column_title_gp = gpar(fontsize = 5),
            #row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
            # column_split = factor(c(rep("Subcortical regions GMV", each = 15),
            #                         rep("Cerebellar regions GMV", each = 28),
            #                         rep("Total brain", each = 5)),
            # levels = c("Subcortical regions GMV","Cerebellar regions GMV","Total brain")),
            column_title_side = "top",
            row_names_side = "right",
            row_title = NULL,
            column_title = "Organismal age gap",
            cell_fun = function(j, i, x, y, w, h, fill) {
              if((mat_p_wrt_organ)[i, j] < 0.001) {
                grid.text("***", x, y,gp = gpar(fontsize = 5))
              } else if((mat_p_wrt_organ)[i, j] < 0.01) {
                grid.text("**", x, y,gp = gpar(fontsize = 5))
              } else if((mat_p_wrt_organ)[i, j] < 0.05) {
                grid.text("*", x, y,gp = gpar(fontsize = 5))
              }
            },
            heatmap_legend_param = list(
              #at = c(-2, 0, 2),
              #labels = c("0", "10%", "20%", "30%", "40%"),
              title = NULL,
              legend_height = unit(1, "cm"),
              grid_width = unit(0.5, "mm"),
              #legend_width =unit(1, "mm"),
              title_position = "leftcenter-rot",
              title_gp=gpar(col = "black", fontsize = 5),
              labels_gp = gpar(col = "black", fontsize = 5)
              
              
            ),
            width = ncol(mat_beta)*unit(1.5, "mm"), 
            height =nrow(mat_beta)*unit(1.5, "mm")
            
    )
  
  
  
  library(ggseg)
  library(ggseg3d)
  library(ggsegHO)
  # ggsegHO::hoCort
  
  # plot(hoCort) +
  #   theme(plot.title = element_blank(),
  #         legend.position = "bottom",
  #         legend.text = element_text(size = 6),
  #         plot.margin = margin(0, 0, 0, 0, "cm")) +
  #   guides(fill = guide_legend(ncol = 5))+theme_brain2()
  
  ho <- as.data.frame(ggsegHO::hoCort) %>% arrange(hemi,region)
  
  ho_data <- as.data.frame(t(beta_matrix_imaging$fdr_p_values)) %>% 
    mutate( n=rowSums(across(everything(), ~ . < 0.05)),
            label=row.names(.)) %>% 
    filter(!grepl("Volume", label)) %>% 
    dplyr::select(label,n)
  
  ho_full <- ho %>% left_join(ho_data,by="label")
  
  brain_region_num <- ho_data %>%
    ggplot() + geom_brain(atlas = hoCort, 
                          #position = position_brain(hemi ~ side),
                          aes(fill = n),inherit.aes = T
    )+theme_minimal(base_line_size = 0.2,base_size = 5,base_family = "Arial")+
    theme(axis.title.x  = NULL,
          axis.title.y  = NULL,
          axis.text.x   = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          legend.key.height = unit(0.1, "cm"),
          legend.position = "bottom",
          legend.direction = "horizontal",
          legend.title = element_text(size = 5, family = "ARIAL"),
          legend.text = element_text(size = 5, family = "ARIAL"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          plot.title = element_text(size = 6, family = "ARIAL",hjust = 0.5),
    )+scale_fill_gradient(high = "#2166ac", low = "white",na.value = "lightgrey",
                          breaks=seq(0, 9, by = 3),
                          name = "No. of significant\nassociations"
    )+
    guides(fill = guide_colorbar(title.position = "bottom",title.hjust = 0.5))+
    labs(title="Grey matter volumes",hjust = 0.5)
  
  # summary(ukbb_trait_imaging$`25005-2.0`)
  # 
  # ggplot(ukbb_trait_imaging %>% filter(!is.na(`25005-2.0`)), aes(x = Organismal_res, y = `25005-2.0`)) +
  #   geom_point(size = 2, alpha = 0.2,color="#4dbbd5ff") +  # Add points
  #   labs(title = "",x = "Chronological age",y = "PhenoAge"
  #   )+
  #   theme(plot.caption = element_text(hjust = 1, vjust = 1, color = "black"))+
  #   sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2,size=0.5)+theme_classic(base_size = 6,base_line_size = 0.2)
  # 
  
  
  library(ggseg)
  library(ggplot2)
  library(ggseg3d)
  library(ggsegTracula)
  library(dplyr)
  
  plot(tracula) +
    theme(legend.position = "bottom",
          legend.text = element_text(size = 9)) +
    guides(fill = guide_legend(ncol = 3))
  
  wm_tract <- 
    ggseg(atlas = tracula, mapping = aes(fill = region))+theme_minimal(base_line_size = 0.2,base_size = 5,base_family = "Arial")+
    theme(axis.title.x  = NULL,
          axis.title.y  = NULL,
          axis.text.x   = element_blank(),
          axis.text.y  = element_blank(),
          panel.grid = element_blank(),
          legend.key.size = unit(0.2, "cm"),
          legend.key.height = unit(0.1, "cm"),
          legend.position = "right",
          legend.direction = "horizontal",
          legend.title = element_text(size = 5, family = "ARIAL"),
          legend.text = element_text(size = 5, family = "ARIAL"),
          plot.margin = margin(0, 0, 0, 0, "cm"),
          plot.title = element_text(size = 6, family = "ARIAL",hjust = 0.5),
    ) +xlab(NULL) +ylab(NULL)+scale_fill_brewer(palette = "Set1", na.value = "lightgrey")+
    guides(fill = guide_legend(ncol = 1,title=NULL))
  
  
  xf7 <- plot_grid(brain_region_num, p_imaging_cortical, p_imaging, wm_tract,
                   plot_grid(grid.grabExpr(draw(p_imaging_wrt + p_imaging_wrt2),show_legend=c(T,F)), labels="e",nrow = 1,label_size = 8), 
                   nrow = 5,labels=c("a","b","c","d",""),label_size = 8,rel_heights = c(1.2,0.9,1,0.8,1)
                   
  )      

  ggsave(filename = "Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Figure/Final/xFig.7.pdf",
         plot = xf7,
         device = cairo_pdf,width =7,height=8.2,units="in")
  
}
