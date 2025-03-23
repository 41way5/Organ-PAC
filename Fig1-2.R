# Lib===================================================
{
  library(lubridate)
  library(survival)
  library(patchwork)
  library(gridExtra)
  library(cowplot)
  library(readxl)
  library(ComplexHeatmap)
  library(circlize)
  library(ggsci)
  library(scales)
  library(magick)
  library(tidyverse)
  library(data.table)
  library(Hmisc)
  library(dplyr)
  library(tidyr)
  library(DESeq2)
  library(grid)
  library(ggseg)
  library(ggplot2)
  library(ggseg3d)
  library(ggsegTracula)
  library(ggrepel)
  library(smplot2)
  setwd("E:\\Tutor\\UK Biobank\\")
  select <- dplyr::select
  
}
load(file = "./Project/Proteomic aging/Analyses/ukb_pool.Rdata")

#-----*****Figure 1*****---------
  #Panel A-C----
    #Organ age gap distribution
  ukb_pred <- ukb_olink %>% select(ends_with("_pAD")) %>% mutate(source="UKB")
  ckb_pred <-read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/CKB/validation/ckb_preds_df.csv")
  nhs1_pred <- read.csv("E:\\Tutor\\UK Biobank\\Project\\Proteomic aging\\Submission\\Submission\\Nat Aging\\Revision\\2-16 R1\\Validation\\NHS\\nhs1\\nhs1\\ckb_preds_df.csv")
  
  pred <- do.call(rbind, list(ukb_pred,
                              ckb_pred %>% select(ends_with("_pAD")) %>% mutate(source="CKB"),
                              nhs1_pred %>% select(ends_with("_pAD")) %>% mutate(source="NHS")
                              # nhs2_pred %>% select(ends_with("_pAD")) %>% mutate(source="NHS II")
  )
  )
  colnames(pred) <- gsub("_pAD","",colnames(pred))
  
  
  pred_long <- pred %>%
    pivot_longer(cols = -source, names_to = "Organ", values_to = "OrganAgeGap") %>% 
    mutate(Organ = factor(Organ,levels= c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas")),
           source=factor(source, levels = c("UKB", "CKB", "NHS"))
    )
  
  #A----
  age_distri <- ggplot(pred_long, aes(x = Organ, y = OrganAgeGap,fill = source)) + 
    geom_violin(trim = FALSE, size = 0.2, alpha = 0.3, position = position_dodge(width = 0.6)) +  # Violin plot with alignment
    geom_boxplot(aes(group = interaction(Organ,source)),width = 0.1,size = 0.1, color="black", alpha = 0.5, position = position_dodge(width = 0.6),show.legend = F,
                 outlier.shape = NA) +  # Boxplot aligned with violin
    # stat_summary(fun = "median", geom = "crossbar", width = 0.09, color = "white", size = 0.1, 
    #              position = position_dodge(width = 0.6),show.legend = F) +  # White median line
    scale_fill_manual(values = c("UKB" = "#2171b5", "CKB" = "#6baed6", "NHS" = "#bdd7e7"))+  # Customize the colors for different cohorts
    labs(x = NULL, y = "ProtAgeGap (years)", fill = "Cohort") +  # Labels
    theme_bw(base_line_size = 0.2,base_size = 7) +
    theme(
      panel.grid = element_blank(),
      legend.direction = "horizontal",
      legend.spacing.x = unit(0.5, "cm"),  # Adjust space between legend items
      legend.key.size = unit(0.2, "cm"),  # Adjust the size of the legend keys (the boxes)
      legend.position = c(0.02, 0.98),  # Place legend inside the plot at the top-left corner
      legend.justification = c("left", "top"),  # Align legend to the top-left corner
      legend.box.just = "left",  # Ensure the legend is aligned to the left
      legend.title = element_text(size = 6),
      text = element_text(size = 7),
      axis.title = element_text(size = 7),
      axis.text = element_text(size = 7),
      axis.text.x = element_text(color = "black") , # Rotate x-axis labels for better readability
      axis.text.y = element_text(color = "black")  
      
    ) 
  #B&C----
  #Correlation between protage and chrono age
  mat_corr_age <- matrix(NA, nrow = length(organ_list), ncol = 3)
  for (i in 1:length(organ_list)) {
    mat_corr_age[i, 1] <- cor(olink_pred[[organ_list[i]]], olink_pred[["age_at_recruitment"]])
    mat_corr_age[i, 2] <- cor(ckb_pred[[organ_list[i]]],  ckb_pred[["age_at_recruitment"]])
    mat_corr_age[i, 3] <- cor(nhs1_pred[[organ_list[i]]], nhs1_pred[["age_at_recruitment"]])
    # mat_corr_age[i, 4] <- cor(nhs2_pred[[organ_list[i]]], nhs2_pred[["age_at_recruitment"]])
    
  }
  
  rownames(mat_corr_age) <- c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas")
  colnames(mat_corr_age) <- c("UKB","CKB","NHS")
  
  data_long <- as.data.frame(mat_corr_age) %>%
    rownames_to_column("Organ") %>%
    gather(key = "Dataset", value = "Correlation", -Organ)
  
  data_long$Shape <- rep(c(17, 3, 5), length.out = nrow(data_long))
  data_long$Organ <- factor(data_long$Organ,levels=rev(organ_list))
  data_long$Dataset <- factor(data_long$Dataset,levels=c("UKB","CKB","NHS"))
  
  
  r <- ggplot(data_long, aes(x = Correlation, y = Organ, color  = Dataset)) +
    geom_point(size = 1.5,aes(shape=Dataset), stroke = 0.1) +  # Use different shapes for each dataset
    scale_shape_manual(values = c(17, 3, 5)) +  # Assign specific shapes (17, 3, 5 for each dataset)
    scale_color_manual(
      values = rep(c("black"),3) #c("UKB" = "#2171b5", "CKB" = "#6baed6", "NHS I" = "#bdd7e7")
      
    )+  # Customize the colors for different cohorts
    labs(x = "ProtAge vs chronological age\n correlation across cohorts", y = NULL) +
    theme_bw(base_line_size = 0.2,base_size = 6) +
    theme(
      panel.grid = element_blank(),
      legend.direction = "vertical",
      legend.spacing.x = unit(0.25, "cm"),  # Adjust space between legend items
      legend.key.size = unit(0.1, "cm"),  # Adjust the size of the legend keys (the boxes)
      legend.position = c(0.98, 0.03),  # Place legend inside the plot at the top-left corner
      legend.justification = c("right", "bottom"),  # Align legend to the top-left corner
      legend.box.just = "left",  # Ensure the legend is aligned to the left
      legend.title = element_blank(),
      text = element_text(size = 6),
      axis.title = element_text(size = 6),
      axis.text = element_text(size = 6),
      axis.text.x = element_text(color = "black") , # Rotate x-axis labels for better readability
      axis.text.y = element_text(color = "black")  
      
    ) 
  
  r12 <- ggplot(as.data.frame(mat_corr_age), aes(x = UKB, y = CKB)) +
    geom_point(aes(color=row.names(mat_corr_age))) +
    labs(x = "ProtAge vs chronological age\n correlation (UKB discovery)", 
         y = "ProtAge vs chronological age\n correlation (CKB validation)"
         
    ) +
    geom_text_repel(aes(label = rownames(mat_corr_age)), size = 2, max.overlaps = 20, direction= "both",segment.size = 0.2)+
    sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2.1, size=0.5)+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    scale_color_manual(values=color_organ)+
    # scale_y_continuous(limits = c(0.2, 1.0), breaks = seq(0.2, 1.0, by = 0.2))+
    theme(legend.position = "none")
  
  r13 <- ggplot(as.data.frame(mat_corr_age), aes(x = UKB, y = NHS)) +
    geom_point(aes(color=row.names(mat_corr_age))) +
    labs(x = "ProtAge vs chronological age\n correlation (UKB discovery)", 
         y = "ProtAge vs chronological age\n correlation (NHS validation)"
         
    ) +
    geom_text_repel(aes(label = rownames(mat_corr_age)), size = 2, max.overlaps = 20, direction= "both",segment.size = 0.2)+
    sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2.1, size=0.5)+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    scale_color_manual(values=color_organ)+
    # scale_y_continuous(limits = c(0.2, 1.0), breaks = seq(0.2, 1.0, by = 0.2))+
    theme(legend.position = "none")
  
  r23 <-ggplot(as.data.frame(mat_corr_age), aes(x = CKB, y = NHS)) +
    geom_point(aes(color=row.names(mat_corr_age))) +
    labs(x = "ProtAge vs chronological age\n correlation (CKB validation)", 
         y = "ProtAge vs chronological age\n correlation (NHS validation)"
         
    ) +
    geom_text_repel(aes(label = rownames(mat_corr_age)), size = 2, max.overlaps = 20, direction= "both",segment.size = 0.2)+
    sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2.1, size=0.5)+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    scale_color_manual(values=color_organ)+
    # scale_y_continuous(limits = c(0.2, 1.0), breaks = seq(0.2, 1.0, by = 0.2))+
    theme(legend.position = "none")
  
  
  r123 <- plot_grid(r,r12,r13,r23,nrow=1)
  
  ###A-C
  fig1_ac <- plot_grid(age_distri,
            plot_grid(r,r12,r13,r23,nrow=1,labels = c("b","c",""),label_size = 8,vjust = 0.2, hjust = 0),
            # plot_grid(pac_corr123,v1,rel_widths = c(3,1)),
            # rel_heights = c(1,1.5,1),
            labels = c("a",""),label_size = 8,
            ncol=1)    
  ###D----
  pac_corr_ht <- function(data,legend=0,column_title=NULL){
    
    mat_corr <- data %>% dplyr::select(paste0(organ_list, "_pAD")) %>% cor()
    rownames(mat_corr) <- c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                            "Lung","Muscle","Pancreas")
    colnames(mat_corr) <- c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                            "Lung","Muscle","Pancreas")
    mat_corr[upper.tri(mat_corr, diag = FALSE)] <- NA
    print(mean(mat_corr[-1, -1][lower.tri(mat_corr[-1, -1])]))
    # grid.grabExpr(draw((mat_corr %>% cor_ht(.,cut = 0.5))))
    mat_corr %>% cor_ht(.,cut = 0.5,legend=legend,column_title = column_title)
  }
  
  pac_corr_ht(ukb_pred)
  pac_corr_ht(ckb_pred)
  pac_corr_ht(nhs1_pred)
  
  pac_corr123 <- grid.grabExpr(draw(
    pac_corr_ht(ckb_pred,column_title = expression(paste("UKB (mean ", italic("r"), " = 0.16)")))+
      pac_corr_ht(ukb_pred,column_title = expression(paste("CKB (mean ", italic("r"), " = 0.19)")))+
      pac_corr_ht(nhs1_pred,column_title = expression(paste("NHS (mean ", italic("r"), " = 0.10)")),legend = 1),auto_adjust = T)) %>% plot_grid()
  
  
  
  ###E----
  #Brouta selected protein
  prot_list <- read.csv(file = "Data/Proteomic/Data/boruta_protein_df.csv") %>% filter_all(any_vars(. != ""))
  x <- lapply(prot_list[, c("Organismal", "Brain", "Artery", "Heart")], function(x) x[x != ""])
  
  #Venn diagram
  display_venn <- function(x, ...){
    library(VennDiagram)
    grid.newpage()
    venn_object <- venn.diagram(x, filename = NULL, ...)
  }
  
  v1 <- display_venn(x[1:4],
                     category.names = c("Organismal","Brain","Artery","Heart"),
                     fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"),
                     lwd = 0.2,
                     cat.fontfamily="arial",
                     fontfamily="arial",cex=0.5,cat.cex = 0.55,margin = 0.1)
  
  plot_grid(v1)
  
  fig1_de <-plot_grid(pac_corr123,v1,rel_widths = c(3,1),labels = c("d","e"),label_size = 8)

  
  ###F-----
  #Association of PAC with INCIDENT disease and death
  length(organ_pac)
  #Run the 11 organ PAC with 10 outcomes
  result_list <- cox(data=ukb_olink,exposures=c(organ_pac),
                     outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death","copd","mi","stroke","hyp","ckd","ld","dm"), 
                     covar=c("age","sex","ethnic","edu","tdi","center","smok","pa"))
  
  ukb_olink_inci <- result_list$data
  
  result_disease_fomrat <-function(result_list) {
    
    result_disease <- result_list$result_list %>% 
      
      mutate(Outcome_label = case_when(
        Outcome == "acd" ~ "All-cause dementia",
        Outcome == "ad" ~ "Alzheimer's disease",
        Outcome == "vd" ~ "Vascular dementia",
        Outcome == "pd" ~ "Parkinson's disease",
        Outcome == "ms" ~ "Multiple sclerosis",
        Outcome == "scz" ~ "Schizophrenia",
        Outcome == "dep" ~ "Depression",
        Outcome == "anx" ~ "Anxiety",
        Outcome == "sd" ~ "Sleep disorder",
        Outcome == "sud" ~ "Substance use disorder",
        Outcome == "death" ~ "All-cause mortality",
        Outcome == "copd" ~ "COPD",
        Outcome == "mi" ~ "Myocardial infarction",
        Outcome == "stroke" ~ "Stroke",
        Outcome == "ckd" ~ "Chronic kidney disease",
        Outcome == "ld" ~ "Chronic liver disease",
        Outcome == "dm" ~ "Type 2 diabetes",
        Outcome == "hyp" ~ "Hypertension",
        TRUE ~ as.character(Outcome)
      ),
      Exposure_label = sub("_res", "", Exposure),
      Exposure_label = sub("_super$", "", Exposure_label),
      Exposure_label = sub("_extra$", "", Exposure_label),
      Exposure_label = case_when(
        Exposure == "phenoage_res" ~ "PhenoAge",
        Exposure == "kdm_res" ~ "KDMage",
        Exposure == "Telomere_length" ~ "Telomere length",
        TRUE ~ as.character(Exposure_label)
      )
      
      ) %>% 
      mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                           "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","All-cause mortality",
                                                           "Hypertension",
                                                           "Myocardial infarction","Stroke","COPD","Chronic kidney disease","Chronic liver disease","Type 2 diabetes"
      )),
      Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas"))
      )
    
    return(result_disease)
  }
  
  result_disease <- result_disease_fomrat(result_list)
  
  #Organ contribution
  result_list_organ<-cox_prot(data=ukb_olink,
                              exposures=organ_pac,
                              outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death",
                                         "copd","mi", "stroke","hyp","ckd","ld","dm"),
                              covar=c("age","sex","ethnic","tdi","center","smok","pa",organ_pac))
  
  result_disease_organ_pac <- result_list_organ$result_list %>%   
    mutate(Outcome_label = case_when(
      Outcome == "acd" ~ "All-cause dementia",
      Outcome == "ad" ~ "Alzheimer's disease",
      Outcome == "vd" ~ "Vascular dementia",
      Outcome == "pd" ~ "Parkinson's disease",
      Outcome == "ms" ~ "Multiple sclerosis",
      Outcome == "scz" ~ "Schizophrenia",
      Outcome == "dep" ~ "Depression",
      Outcome == "anx" ~ "Anxiety",
      Outcome == "sd" ~ "Sleep disorder",
      Outcome == "sud" ~ "Substance use disorder",
      Outcome == "death" ~ "All-cause mortality",
      Outcome == "copd" ~ "COPD",
      Outcome == "mi" ~ "Myocardial infarction",
      Outcome == "stroke" ~ "Stroke",
      Outcome == "ckd" ~ "Chronic kidney disease",
      Outcome == "ld" ~ "Chronic liver disease",
      Outcome == "dm" ~ "Type 2 diabetes",
      Outcome == "hyp" ~ "Hypertension",
      TRUE ~ as.character(Outcome)
    ),
    Exposure = gsub("_res","",Exposure)
    ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                               "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","All-cause mortality",
                                                               "Hypertension",
                                                               "Myocardial infarction","Stroke","COPD","Chronic kidney disease","Chronic liver disease","Type 2 diabetes")
    ),
    
    Exposure=factor(Exposure, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                                       "Lung","Muscle","Pancreas"))
    ) %>% 
    mutate(z_score_trans=ifelse(P_Value>0.05,0,z_score))
  
  
  #Disease association with individual organ PAC

  mat <- xtabs(HR ~ Outcome_label + Exposure_label, result_disease)
  mat_p <- xtabs(FDR_P_Value ~ Outcome_label + Exposure_label, result_disease)
  
  mat_case <- result_disease %>%
    select(Outcome_label, N_Cases_base, N_Cases_inci) %>%
    distinct(Outcome_label, .keep_all = TRUE) %>%
    arrange(factor(Outcome_label, levels = levels(result_disease$Outcome_label))) %>%
    mutate(
      N_Cases_base = ifelse(Outcome_label == "Hypertension", N_Cases_base / 4, N_Cases_base),
      N_Cases_inci = ifelse(Outcome_label == "Hypertension", N_Cases_inci / 4, N_Cases_inci)
    ) %>% 
    # mutate(N_Cases_inci=N_Cases_inci-N_Cases_base ) %>% 
    select(-Outcome_label)
  
  mat_disease_gap <- xtabs(Mean_diff ~ Outcome_label + Exposure_label, result_disease)
  
  mat_prot_trans <- xtabs(z_score_trans ~ Outcome_label + Exposure, result_disease_organ_pac) %>% as.data.frame() %>% 
    group_by(Outcome_label) %>% 
    mutate(Freq=ifelse(Freq<0,-Freq,Freq)) %>% 
    mutate(Freq_pct=Freq/sum(Freq))
  
  mat_prot_trans <- xtabs(Freq_pct ~ Outcome_label + Exposure, mat_prot_trans)
  mat_prot_trans_bar <- mat_prot_trans[]
  # mat_prot <- mat_prot[order(rowMeans(mat_prot)), ]
  # mat_prot_p <- mat_prot_p[order(rowMeans(mat_prot)), ]
  
  col_fun = colorRamp2(c(0,1,1.5), c("#4DBBD5B2","white","#DC0000B2"))
  # col_fun = colorRamp2(c(-10,0,20), c("#3288bd","white","#d53e4f"))
  col_panel_prot=c("black",c( "#A52A2A", "#4B0082","#FF7F50", "#FFFF00", "#808080", "#808000",
                                       "#D2B48C", "#00FFFF", "#0000FF","#008000"))
  
  ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "mm")
  
      fig1_f <-grid.grabExpr(draw( Heatmap(mat, col = col_fun, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                cluster_rows = F,cluster_columns = F,
                row_gap = unit(0.75, "mm"),
                width = unit(4, "cm"), height = unit(3.5, "cm"),
                row_names_side = "left",row_names_gp = gpar(fontsize = 6),column_names_rot=45,
                column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 6),
                row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",5),rep("All-cause mortality",1),rep("Others",7)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality","Others")),
                # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                column_title_gp = gpar(fontsize = 6),
                column_title_side = "bottom",
                cell_fun = function(j, i, x, y, width, height, fill) {if(mat_p[i, j] < 0.05)
                  grid.text(sprintf("%.1f", mat[i, j]), x, y, gp = gpar(fontsize = 5))
                },
                heatmap_legend_param = list(
                  #at = c(-2, 0, 2),
                  #labels = c("0", "10%", "20%", "30%", "40%"),
                  title = NULL,
                  legend_height = unit(1, "cm"),
                  grid_width = unit(0.5, "mm"),
                  #legend_width =unit(1, "mm"),
                  #title_position = "lefttop-rot",
                  labels_gp = gpar(col = "black", fontsize = 5)
                  
                ),
                bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                                                    col= list(" "=setNames(col_panel_prot, levels(result_disease_organ_pac$Exposure))),
                                                    simple_anno_size = unit(1.5, "mm"),
                                                    show_legend = FALSE), 
                
                right_annotation=rowAnnotation(
                  "Disease cases"= anno_barplot(mat_case, border=F, ylim=c(0,6000), gp = gpar(fill =c("#3c5488ff","#00a087ff"), col = c("#3c5488ff","#00a087ff")),
                                                axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,3000,6000),labels = c("0","3000", "6000"),labels_rot=0),
                  ),
                  
                  "Organ Age gap for each outcome" = anno_points(mat_disease_gap, ylim = c(-0.3,0.6),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                                                           pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                  
                  "Independent contribution of organ aging"= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                                    gp = gpar(fill =col_panel_prot, 
                                              col = col_panel_prot,lwd=0.2),bar_width = 0.8,
                                    axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                  ),
                  
                  width = unit(9, "cm"), annotation_name_gp= gpar(fontsize = 6),
                  annotation_width=unit(c(0.5, 1, 1), c("null", "null","null")),
                  gap = unit(1, "mm")
                  
                ),
                
                row_title =NULL,column_title = NULL
                
  )
      ), padding = unit(c(0, 0, 0, 0), "mm")
  
  )
       
       plot_grid(fig1_f)
  ###G-H-----
       #panel for Prevalent disease and age gap----
       
       add_data <- ukbb_full %>% 
         dplyr::select(c(eid,`26214-0.0`,`26275-0.0`,
                         `2050-0.0`,`2060-0.0`,`2070-0.0`,`2080-0.0`,#PHQ4 at baseline
                         `20518-0.0`,`20510-0.0`,`20507-0.0`,`20519-0.0`,`20514-0.0`,`20511-0.0`,`20513-0.0`,`20508-0.0`,`20517-0.0`,#PHQ9_fu
                         `20505-0.0`,`20512-0.0`,`20506-0.0`,`20509-0.0`,`20516-0.0`,`20515-0.0`,`20520-0.0`,#GAD7_fu
                         `20500-0.0`,`20480-0.0`,`20483-0.0`
         )) %>% 
         dplyr::rename(bd_prs=`26214-0.0`,scz_prs=`26275-0.0`) %>% 
         mutate(
           #PHQ4 and defined dep and anx AT BASELINE
           across(all_of(c("2050-0.0","2060-0.0","2070-0.0","2080-0.0")), ~ if_else(. < 0, NA, .)),
           across(all_of(c("2050-0.0","2060-0.0","2070-0.0","2080-0.0")), ~ . - 1),
           phq4=rowSums(across(all_of(c("2050-0.0","2060-0.0","2070-0.0","2080-0.0"))), na.rm = F),
           phq4_1=rowSums(across(all_of(c("2050-0.0","2060-0.0"))), na.rm = F),
           phq4_2=rowSums(across(all_of(c("2070-0.0","2080-0.0"))), na.rm = F),
           dep_bas=ifelse(phq4_1>=3,1,0),
           anx_bas=ifelse(phq4_2>=3,1,0),#Cutoff >=3
           
           #PHQ9/GAD7 and defined dep and anx AT FOLLOW-UP
           #PHQ9
           across(all_of(c("20518-0.0", "20510-0.0", "20507-0.0", "20519-0.0", 
                           "20514-0.0", "20511-0.0", "20513-0.0", "20508-0.0", "20517-0.0")), ~ if_else(. < 0, NA, .)),
           across(all_of(c("20518-0.0", "20510-0.0", "20507-0.0", "20519-0.0", 
                           "20514-0.0", "20511-0.0", "20513-0.0", "20508-0.0", "20517-0.0")), ~ . - 1),
           
           phq9_fu=rowSums(across(all_of(c("20518-0.0", "20510-0.0", "20507-0.0", "20519-0.0", 
                                           "20514-0.0", "20511-0.0", "20513-0.0", "20508-0.0", "20517-0.0"))), na.rm = F),
           dep_fu=ifelse(phq9_fu>=10,1,0),
           #GAD7
           across(all_of(c("20505-0.0", "20512-0.0", "20506-0.0", "20509-0.0", "20516-0.0", "20515-0.0", "20520-0.0")), ~ if_else(. < 0, NA, .)),
           across(all_of(c("20505-0.0", "20512-0.0", "20506-0.0", "20509-0.0", "20516-0.0", "20515-0.0", "20520-0.0")), ~ . - 1),
           
           gad7_fu=rowSums(across(all_of(c("20505-0.0", "20512-0.0", "20506-0.0", "20509-0.0", "20516-0.0", "20515-0.0", "20520-0.0"))), na.rm = F),
           anx_fu=ifelse(gad7_fu>=10,1,0),
           distress=ifelse(`20500-0.0`==1,1,0),
           self_harm=ifelse(`20480-0.0`==1,1,0),
           suicide=ifelse(`20483-0.0`==1,1,0),
         )
       
       #Trait list including FU data
       ukbb_trait_list_fu <- ukbb_index %>% 
         filter(Category %in% c("Physical measures and age-related phenotype",
                                "Cognitive function",
                                "Mental well-being",
                                "Blood biochemistry"
                                #,"Neurobiomarkers"
         )) %>% 
         mutate(fid=paste0(`UKB Field ID`,"-0.0")) %>% 
         mutate(fid=ifelse(`UKB Field ID` %in% c("31040","31041","31042","31043","31044"),paste0(`UKB Field ID`,"-2.0"),fid)) %>%
         mutate(fid=ifelse(`UKB Field ID` %in% c("399"),paste0(`UKB Field ID`,"-0.1"),fid)) %>% 
         mutate(Categroy_short=case_when(Category=="Physical measures and age-related phenotype"~"Age-related phenotype",
                                         Category=="Cognitive function"~"Cognition",
                                         Category=="Mental well-being"~"Mental",
                                         TRUE ~ Category)) %>% filter(!is.na(Abbr))
       
       
       cognition <- ukb_olink %>% left_join(ukbb_full %>% dplyr::select(eid,(ukbb_trait_list_fu %>%filter(Category == "Cognitive function"))$fid),by="eid") %>% 
         mutate(cog_flag=ifelse(rowSums(!is.na(.[c("4282-0.0", "20016-0.0","399-0.1")])) == 3, 1, 0)) %>% 
         mutate(rt=as.numeric(scale(`20023-0.0`)),
                nm=as.numeric(scale(`4282-0.0`)),nm_fu=as.numeric(scale(`20240-0.0`)),
                fi=as.numeric(scale(`20016-0.0`)),fi_fu=as.numeric(scale(`20191-0.0`)),
                vm=-as.numeric(scale(`399-0.1`)),vm_fu=-as.numeric(scale(`20132-0.0`))
         ) %>% mutate_at(vars(ends_with("_date")), ymd) %>% 
         group_by(edu) %>%
         mutate(vm_cut = mean(`399-0.1`, na.rm = TRUE)+1.5 * sd(`399-0.1`, na.rm = TRUE)) %>% ungroup() %>% 
         mutate(vm_change=ifelse(!is.na(vm)&!is.na(vm_fu)&`20132-0.0`>=vm_cut,1,0)) %>% 
         mutate(vm_change=factor(vm_change,levels=c(1,0))) %>% left_join(add_data,by="eid")
       
       
       #Trajectoris of age gap with time to diagnosis
       disease_prefixes <- c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","hyp",
                             "copd","mi","stroke","ckd","ld","dm")
       
       for (prefix in disease_prefixes) {
         interval_col <- paste0(prefix, "_interval")
         disease_col <- paste0(prefix, "_con")
         base_col <- paste0(prefix, "_base")
         inci_col <- paste0(prefix, "_inci")
         
         cognition <- cognition %>% 
           mutate(!!interval_col := if_else(!is.na(!!sym(paste0(prefix, "_date"))),
                                            as.numeric(difftime(.$recruit_date, !!sym(paste0(prefix, "_date")), units = "days")),0),
                  
                  !!disease_col :=if_else(!is.na(!!sym(paste0(prefix, "_date"))),1,0) %>% as.factor(),
                  
                  !!base_col :=if_else(!!sym(interval_col)>0,1,0) %>% as.factor(),
                  !!inci_col :=if_else(!!sym(interval_col)<0,1,0) %>% as.factor(),
                  
           )
       }
       
       
       progression <- cognition %>%  mutate(across(contains("_interval"), ~ ifelse(. / 365.25 >= 20, 20, . / 365.25)))
       
       # Match people for those with incident disease
       case_control_match <- function(data, disease_status, disease_interval) {
         # Subset cases (individuals where disease_status_var == 1)
         cases <- data[data[[disease_status]] == 1, ]
         
         # Initialize an empty list to store matched controls
         controls <- list()
         
         # Loop through each case
         for (i in 1:nrow(cases)) {
           # Extract the case information
           case <- cases[i, ]
           age_range <- case$age + c(-2, 2)
           sex <- case$sex
           
           # Subset potential controls based on age and sex
           potential_controls <- data[data[[disease_status]] == 0 & 
                                        data$age >= age_range[1] & 
                                        data$age <= age_range[2] & 
                                        data$sex == sex, ]
           
           # Randomly select 5 controls from potential_controls
           matched_controls <- potential_controls[sample(nrow(potential_controls), 5), ]
           
           # Assign the time_to_disease of cases to matched controls
           
           matched_controls[[disease_interval]] <- case[[disease_interval]] 
           # Add matched controls to the list
           controls[[i]] <- matched_controls
         }
         
         # Combine all matched controls into a single data frame
         matched_controls_df <- do.call(rbind, controls)
         
         # Return the matched controls data frame
         return(matched_controls_df)
       }
       
       
       control_acd <- case_control_match(progression, "acd_inci","acd_interval")
       control_ad <- case_control_match(progression, "ad_inci","ad_interval")
       control_ckd <- case_control_match(progression, "ckd_inci","ckd_interval")
       control_dep <- case_control_match(progression, "dep_inci","dep_interval")
       
       
       matched_ad <- rbind(progression %>% filter(ad_inci==1|ad_base==1),control_ad) %>% 
         mutate(ad_con=factor(ad_con,levels = c(0, 1),labels = c("Control", "Dementia")))
       
       matched_acd <- rbind(progression %>% filter(acd_inci==1|acd_base==1),control_acd) %>% 
         mutate(acd_con=factor(acd_con,levels = c(0, 1),labels = c("Control", "Dementia")))
       
       matched_ckd <- rbind(progression %>% filter(ckd_inci==1|ckd_base==1),control_ckd) %>% 
         mutate(ckd_con=factor(ckd_con,levels = c(0, 1),labels = c("Control", "CKD")))
       
       matched_dep <- rbind(progression %>% filter(dep_inci==1|dep_base==1),control_dep) %>% 
         mutate(dep_con=factor(dep_con,levels = c(0, 1),labels = c("Control", "Depression")))
       
       
       
       
       acd_curve <- ggplot(matched_acd, aes(x = acd_interval, y = Brain_res, color=acd_con, fill=acd_con, group=acd_con)) +
         geom_smooth(method = "loess", se = TRUE, size = 0.2) +  # Linear regression lines with confidence intervals
         labs(
           x = "Time since diagnosis",
           y = "Brain age gap",
           color = NULL,  # Remove legend title
           fill = NULL  # Remove fill legend title
         ) +
         scale_color_manual(values = c("black", "black")) +
         scale_fill_manual(values = alpha(c("#4292c6", "#fc9272"), 0.3)) +
         theme_classic(base_size = 6,base_line_size = 0.1) +
         theme(
           legend.position = c(0.3, 0.92),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, "cm")
           
         )+
         guides(fill = guide_legend(nrow = 2))
       
       
       ckd_curve <- ggplot(matched_ckd, aes(x = ckd_interval, y = Brain_res, color=ckd_con, fill=ckd_con, group=ckd_con)) +
         geom_smooth(method = "loess", se = TRUE, size = 0.2) +  # Linear regression lines with confidence intervals
         labs(
           x = "Time since diagnosis",
           y = "Brain age gap",
           color = NULL,  # Remove legend title
           fill = NULL  # Remove fill legend title
         ) +
         scale_color_manual(values = c("black", "black")) +
         scale_fill_manual(values = alpha(c("#4292c6", "#fc9272"), 0.3)) +
         theme_classic(base_size = 6,base_line_size = 0.1) +
         theme(
           legend.position = c(0.3, 0.92),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, "cm")
         )+
         guides(fill = guide_legend(nrow = 2))
       
       
       dep_curve <- ggplot(matched_dep, aes(x = dep_interval, y = Brain_res, color=dep_con, fill=dep_con, group=dep_con)) +
         geom_smooth(method = "loess", se = TRUE, size = 0.2) +  # Linear regression lines with confidence intervals
         labs(
           x = "Time since diagnosis",
           y = "Brain age gap",
           color = NULL,  # Remove legend title
           fill = NULL  # Remove fill legend title
         ) +
         scale_color_manual(values = c("black", "black")) +
         scale_fill_manual(values = alpha(c("#4292c6", "#fc9272"), 0.3)) +
         theme_classic(base_size = 6,base_line_size = 0.1) +
         theme(
           legend.position = c(0.3, 0.92),
           legend.title = element_blank(),
           legend.key.size = unit(0.2, "cm")
         )+
         guides(fill = guide_legend(nrow = 2))
       
       
       
       #Exclude those with incident disease (for analysis of prevalent cases)
       filtered_data <- progression %>%  mutate(across(ends_with("_interval"), ~ ifelse(. < 0, 0, .)))
       

       beta_matrix_prevalent <- calculate_beta_matrix(filtered_data,
                                                      var_list1=c(organ_pac), 
                                                      var_list2=paste0(c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","hyp",
                                                                         "mi","stroke","copd","ckd","ld","dm"),"_interval"),
                                                      covariates = c("age","sex","ethnic","edu","tdi","center","smok","pa"))
       
       rownames(beta_matrix_prevalent$beta_matrix) <- c(organ_list)
       colnames(beta_matrix_prevalent$beta_matrix) <- c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                        "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","Hypertension",
                                                        "Myocardial infarction","Stroke","COPD","Chronic kidney disease","Chronic liver disease","Type 2 diabetes")
       
       rownames(beta_matrix_prevalent$fdr_p_values) <- c(organ_list)
       colnames(beta_matrix_prevalent$fdr_p_values) <- c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                         "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","Hypertension",
                                                         "Myocardial infarction","Stroke","COPD","Chronic kidney disease","Chronic liver disease","Type 2 diabetes")
       
       
       
       mat_beta <- (beta_matrix_prevalent$beta_matrix) %>% t() 
       mat_p <- (beta_matrix_prevalent$fdr_p_values)%>% t() 
       
       fig1_gh <- grid.grabExpr(draw( 
         Heatmap(mat_beta,
                 name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                 col = colorRamp2(c(-0.4,0,0.4), c("#4DBBD5B2","white","#DC0000B2")),
                 show_row_names = TRUE, 
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 6),
                 column_names_gp = gpar(fontsize = 6),
                 
                 cluster_rows = T,cluster_columns = T,
                 row_dend_width = unit(0.3, "cm"),
                 column_dend_gp = gpar(lwd = 0.2),
                 row_dend_gp = gpar(lwd = 0.2),
                 show_column_dend = T, show_row_dend = F,
                 column_dend_height  = unit(0.1, "cm"),
                 column_names_rot=45,
                 # column_split=ukbb_trait_list$Categroy_short,
                 column_title_gp = gpar(fontsize = 6),
                 # row_split = factor(c(rep("Organ-specific aging", each = 11)),levels = c("Organ-specific aging")),
                 column_title_side = "top",
                 row_names_side = "left",
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
                   at = c(-0.4, 0, 0.4),
                   # labels = c(-0.1, 0, 0.5),
                   title = "r",
                   title_gp=gpar(col = "black", fontsize = 5),
                   legend_height = unit(1, "cm"),
                   grid_width = unit(0.5, "mm"),
                   legend_width =unit(1, "mm"),
                   # title_position = "lefttop-rot",
                   labels_gp = gpar(col = "black", fontsize = 5)
                   
                 ),
                 width = unit(4, "cm"), 
                 height =unit(3, "cm")
                 
         ),padding = unit(c(2, 1, 1, 0), "mm")
       )
       ) %>% plot_grid(.,plot_grid(ckd_curve,acd_curve,dep_curve,nrow=1,scale = 1),nrow=1,labels = c("g","h","",""),label_size = 8,rel_widths = c(1,1.5))
       
      
 #Fig. 1 OUTPUT--------------
       #7x7.8 inch
  plot_grid(fig1_ac+theme(plot.margin = margin(0, 0, 0, 0, "cm")),
            fig1_de+theme(plot.margin = margin(0, 0, 0, 0, "cm")),
            fig1_f,
            fig1_gh,
            ncol = 1,
            rel_heights = c(1,0.6,0.75,0.6))
                                             
                                            
 #----------*****Figure 2*****---------------      
       #Fig 2: Ageotypes and disease outcomes
       
       result_list_extra <- cox(data=ukb_olink,exposures=c(organ_extra),
                                outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death","copd","mi","stroke","ckd","ld","dm"), 
                                covar=c("age","sex","ethnic","edu","tdi","center","pa","smok"))
       
       result_list_super <- cox(data=ukb_olink,exposures=c(organ_super),
                                outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death","copd","mi","stroke","ckd","ld","dm"), 
                                covar=c("age","sex","ethnic","edu","tdi","center","pa","smok"))
       
       
       #Format Cox output
       result_disease_super <- result_disease_fomrat(result_list_super)
       result_disease_extra <- result_disease_fomrat(result_list_extra)    
       
       ##Forest plot For a specific outcome
       multi_pac_disease <- function(data,exposure_list=organ_list,outcome_list,Pval=T){
         
         organ_list <- c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver","Lung","Muscle","Pancreas")
         organ_pac <- paste0(organ_list, "_res")
         organ_extra <- paste0(organ_pac, "_extra")
         organ_super<- paste0(organ_pac, "_super")
         
         data <- data %>% filter(Exposure_label %in% exposure_list)
         
         plot_list <- list()
         
         for (outcome in outcome_list) {
           forest <- ggplot(data %>% filter(Outcome==outcome), aes(x = HR, y = reorder(Exposure_label, HR))) +
             geom_point(size = 1,shape=15)+
             geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.1) +
             geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
             # scale_x_continuous(limits = c(0.5, 2.6), 
             #                    breaks = c(1,2,2.5)
             #                    ) +
             scale_y_discrete()+
             labs(x = "HR (95% CI)", y = ifelse(grepl("_super", data$Exposure), "ProtAgeSuper",
                                                ifelse(grepl("_extra", data$Exposure), "ProtAgeExtra","ProtAgeAccel")), title =(data %>% filter(Outcome==outcome))$Outcome_label) +
             theme_classic(base_size = 6,base_line_size = 0.1)+
             theme(axis.text.x = element_text(size = 6,colour = "black"),
                   axis.text.y = element_text(size = 6,colour = "black"),
                   axis.title = element_text(size = 6,colour = "black"),
                   plot.margin = margin(0.1, 0, 0.1, 0.1, "cm"),
                   plot.title = element_text(hjust = 0.5, face = "bold",size=6) 
             )
           if(Pval==T){
             tab <- ggplot(data %>% filter(Outcome==outcome), aes(y=reorder(Exposure_label, HR)))+
               geom_text(aes(x = 0, label = ifelse(FDR_P_Value < 0.01, format(FDR_P_Value, scientific = TRUE, digits=3), sprintf("%0.2f", FDR_P_Value))), size = 2)+
               scale_x_continuous(
                 breaks = 0, labels = c('FDR p-value'), 
                 position = 'top', expand = expansion(add = 0.5)) +
               theme_void()+ theme(axis.text.x = element_text(size=6),
                                   plot.margin = margin(0.1, 0.1, 0.1, 0, "cm"))
             combined_plot <- forest + tab + plot_layout(ncol = 2, widths = c(3, 1))
           }else{combined_plot <-forest}
           
           plot_list[[outcome]] <- combined_plot
         }
         return(plot_list)
       }
       #Extra
       
       
       multi_pac_disease_plot <- multi_pac_disease(result_disease, 
                                                   outcome_list=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death",
                                                                  "copd","mi","stroke","ckd","ld","dm"),Pval = T)
       
       multi_extra_disease_plot <- multi_pac_disease(result_disease_extra, 
                                                     outcome_list=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death",
                                                                    "copd","mi","stroke","ckd","ld","dm"
                                                     ),Pval = T)
       
       multi_super_disease_plot <- multi_pac_disease(result_disease_super, 
                                                     outcome_list=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death",
                                                                    "copd","mi","stroke","ckd","ld","dm"
                                                     ),Pval = T)
       
       plot_grid(multi_super_disease_plot$death)
       
       #Supplementary figs
       pac_by_disease <- plot_grid(plotlist =multi_pac_disease_plot[1:17],ncol=3,labels = "auto",label_size = 7)
       super_by_disease <- plot_grid(plotlist =multi_super_disease_plot[1:17],labels = "auto",label_size = 7)
       extra_by_disease <- plot_grid(plotlist =multi_extra_disease_plot[1:17],labels = "auto",label_size = 7)
       
       
       
       #A----
       {
         col_fun = colorRamp2(c(-0.75,0,0.75), c("#4DBBD5B2","white","#DC0000B2"))
         col_fun2 = colorRamp2(c(0,1), c("white","#DC0000B2"))
         #Nature style color panel
         #show_col(pal_npg("nrc", alpha = 1)(9))
         color_panel <- c("Black", c( "#A52A2A", "#4B0082","#FF7F50", "#FFFF00", "#808080", "#808000",
                                               "#D2B48C", "#00FFFF", "#0000FF","#008000"))
                                               
         mat_super <- xtabs(Beta ~ Outcome_label + Exposure_label, result_disease_super)
         mat_p_super <- xtabs(FDR_P_Value ~ Outcome_label + Exposure_label, result_disease_super)
         # mat_super[mat_p_super > 0.05] <- 1            
         
         
         super <- Heatmap(mat_super, col = col_fun, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                          cluster_rows = F,cluster_columns = F,
                          width = unit(4, "cm"), height = unit(3.5, "cm"),
                          row_names_side = "left",row_names_gp = gpar(fontsize = 6),column_names_rot=45,
                          column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 6),
                          row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",5),rep("All-cause mortality",1),rep("Others",6)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality","Others")),
                          # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                          # column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                          column_title_gp = gpar(fontsize = 6),
                          column_title_side = "top",
                          cell_fun = function(j, i, x, y, w, h, fill) {
                            if((mat_p_super)[i, j] < 0.001) {
                              grid.text("***", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_p_super)[i, j] < 0.01) {
                              grid.text("**", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_p_super)[i, j] < 0.05) {
                              grid.text("*", x, y,gp = gpar(fontsize = 4))
                            }
                          },
                          heatmap_legend_param = list(
                            # at = c(-0.75,-0.5, 0, 0.5,0.75),
                            # labels = c(-0.75,-0.5, 0, 0.5,0.75),
                            title = NULL,
                            legend_height = unit(1, "cm"),
                            grid_width = unit(0.5, "mm"),
                            #legend_width =unit(1, "mm"),
                            #title_position = "lefttop-rot",
                            labels_gp = gpar(col = "black", fontsize = 5)
                            
                          ),
                          # right_annotation=rowAnnotation(
                          #   "Disease cases"= anno_barplot(mat_case, border=F, ylim=c(0,6000), gp = gpar(fill =c("#3c5488ff","#00a087ff"), col = c("#3c5488ff","#00a087ff")),
                          #                                 axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,3000,6000),labels = c("0","3000", "6000"),labels_rot=0),
                          #   ),
                          #   
                          #   "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                          #                                            pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                          #   
                          #   width = unit(8, "cm"), annotation_name_gp= gpar(fontsize = 6),
                          #   annotation_width=unit(c(1, 2), c("null", "null")),
                          #   gap = unit(3, "mm")
                          #   
                          # ),
                          
                          row_title =NULL,column_title = "Super Ageotype (-)"
         )
         
         
         mat_extra <- xtabs(Beta ~ Outcome_label + Exposure_label, result_disease_extra)
         mat_p_extra <- xtabs(FDR_P_Value ~ Outcome_label + Exposure_label, result_disease_extra)
         # mat_extra[mat_p_extra > 0.05] <- 1            
         
         extra <- Heatmap(mat_extra, col = col_fun, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                          cluster_rows = F,cluster_columns = F,
                          width = unit(4, "cm"), height = unit(3.5, "cm"),
                          row_names_side = "left",row_names_gp = gpar(fontsize = 6),column_names_rot=45,
                          column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 6),
                          row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",5),rep("All-cause mortality",1),rep("Others",6)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality","Others")),
                          # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                          # column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                          column_title_gp = gpar(fontsize = 6),
                          column_title_side = "top",
                          cell_fun = function(j, i, x, y, w, h, fill) {
                            if((mat_p_extra)[i, j] < 0.001) {
                              grid.text("***", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_p_extra)[i, j] < 0.01) {
                              grid.text("**", x, y,gp = gpar(fontsize = 4))
                            } else if((mat_p_extra)[i, j] < 0.05) {
                              grid.text("*", x, y,gp = gpar(fontsize = 4))
                            }
                          },
                          heatmap_legend_param = list(
                            # at = c(-0.75,-0.5, 0, 0.5,0.75),
                            # labels = c(-0.75,-0.5, 0, 0.5,0.75),
                            title = NULL,
                            legend_height = unit(1, "cm"),
                            grid_width = unit(0.5, "mm"),
                            #legend_width =unit(1, "mm"),
                            #title_position = "lefttop-rot",
                            labels_gp = gpar(col = "black", fontsize = 5)
                            
                          ),
                          # right_annotation=rowAnnotation(
                          #   "Disease cases"= anno_barplot(mat_case, border=F, ylim=c(0,6000), gp = gpar(fill =c("#3c5488ff","#00a087ff"), col = c("#3c5488ff","#00a087ff")),
                          #                                 axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,3000,6000),labels = c("0","3000", "6000"),labels_rot=0),
                          #   ),
                          #   
                          #   "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                          #                                            pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                          #   
                          #   width = unit(8, "cm"), annotation_name_gp= gpar(fontsize = 6),
                          #   annotation_width=unit(c(1, 2), c("null", "null")),
                          #   gap = unit(3, "mm")
                          #   
                          # ),
                          
                          row_title =NULL,column_title = "Accelerated Ageotype (+)"
         )
       }
       
       
       
       result_ageo <-  rbind(result_disease_extra %>% mutate(label="extra"),result_disease_super %>% mutate(label="super"))
       
       ageo_forest <- ggplot(result_ageo %>% 
                               filter(Outcome %in% c("death","acd","dep","ckd")) %>% 
                               mutate(Outcome=factor(Outcome,levels=c("death","acd","dep","ckd"))) %>% 
                               group_by((label),Outcome) %>% arrange(-HR) %>% mutate(Exposure_label=factor(Exposure_label,levels = rev(unique(Exposure_label))),
                                                                                     Outcome_label=factor(Outcome_label)
                               )%>% ungroup(), 
                             # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("All-cause dementia","Alzheimer's disease","Vascular dementia","Parkinson's disease","Multiple sclerosis","All-cause mortality"))),
                             #        Exposure_label=factor(Exposure_label,levels=c("Heart","Artery","Organismal","Brain"))
                             # ),
                             aes(x = HR, y=(Outcome_label),groups=(Exposure_label),shape=factor(label))) +
         geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
         geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
         geom_point(aes(color = Exposure_label),size = 1,position=position_dodge(0.7))+
         # scale_x_continuous(limits = c(0, 4.8), breaks = c(0.5,1,2,3,4)) +
         scale_y_discrete()+
         labs(x = "HR (95% CI)", y = NULL, title = NULL) +
         theme_classic(base_size = 6,base_line_size = 0.2)+
         theme(            legend.position = "right",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                           axis.text.x = element_text(color="black",size = 5),
                           axis.text.y = element_text(color="black",size = 6),
                           legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                           legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                           legend.key.spacing.y = unit(0, "pt"),
                           legend.box.spacing = margin(0),
                           # legend.margin = margin(t = 1, b =1, unit = "pt"),
                           # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                           
         )+
         # facet_grid(~label,scales = "free")+
         coord_flip()+
         scale_color_manual(values=c( "Organismal" = "black",
                                      "Artery" = "#A52A2A",
                                      "Brain" = "#4B0082",
                                      "Heart" = "#FF7F50",
                                      "Immune" = "#FFFF00",
                                      "Intestine" = "#808080",
                                      "Kidney" = "#808000",
                                      "Liver" = "#D2B48C",
                                      "Lung" = "#00FFFF",
                                      "Muscle" = "#0000FF",
                                      "Pancreas" = "#008000"))+
         scale_shape_manual(values = c(17, 16))+scale_x_log10(limits = c(0.25, 5),breaks=c(0.5,1,2,4)) 
       
       
       ukb_olink_inci2 <- ukb_olink_inci %>%  
         mutate(across(contains("_num"), ~ case_when(
           . == 0 ~ "0",
           . %in% 1:2 ~ "1-2",
           . %in% 3:4 ~ "3-4",
           . >= 5 ~ "5+",
           TRUE ~ as.character(.)
         ), .names = "{.col}_cat")) %>% 
         mutate(neuropsy=pmax(acd_inci,pd_inci,ms_inci,scz_inci,dep_inci,anx_inci,sd_inci,na.rm = TRUE),
                chronic=pmax(copd_inci,mi_inci,stroke_inci,ckd_inci,ld_inci,dm_inci,na.rm = TRUE),
                co_neuropsy=rowSums(.[c(paste0(c("acd","ad","vd","pd","ms","scz","dep","anx","sd"),"_inci"))]),
                co_neuropsy2=case_when(
                  co_neuropsy == 0 ~ "0",
                  co_neuropsy %in% 1 ~ "one",
                  co_neuropsy >= 2 ~ "two",
                  TRUE ~ as.character(co_neuropsy)
                ),
                co_chronic=rowSums(.[c(paste0(c( "copd","mi","stroke","ckd","ld","dm"),"_inci"))]),
                co_chronic2=case_when(
                  co_chronic == 0 ~ "0",
                  co_chronic %in% 1 ~ "one",
                  co_chronic >= 2 ~ "two",
                  TRUE ~ as.character(co_chronic)
                ),
                neuropsy_chronic=case_when(
                  co_neuropsy2 == 0 & co_chronic2 == 0 ~ "0",
                  co_neuropsy2 != 0 & co_chronic2 == 0 ~ "a_neuro_psy",
                  co_neuropsy2 == 0 & co_chronic2 != 0 ~ "b_chronic",
                  co_neuropsy2 != 0 & co_chronic2 != 0  ~ "comorbidity"
                )
                
                
         ) %>% mutate(num_all=rowSums(.[c("co_neuropsy","co_chronic")]))
       
       
       cox_format <- function(cox)
       {summary(cox)$coefficients %>%
           as.data.frame() %>%
           rownames_to_column("Variable") %>%
           mutate(
             HR = `exp(coef)`,
             CI_Lower = exp(confint(cox)[, 1]),
             CI_Upper = exp(confint(cox)[, 2])
           ) %>%
           arrange(desc(HR)) %>% filter(!str_detect(Variable, "edu")) %>% 
           filter(!str_detect(Variable, "age")) %>%
           filter(!str_detect(Variable, "sex"))
         # mutate(Variable= sub("_res", "", Variable) %>% str_to_title)
       }  
       
       
       cox_super_ad <- coxph(as.formula(paste("Surv(acd_fu, acd_inci) ~", paste(c("as.factor(super_num_cat)", "age", "sex","edu"), collapse = "+"))),
                             data = ukb_olink_inci2) %>% cox_format() %>% mutate(label="Super")
       cox_extra_ad <- coxph(as.formula(paste("Surv(acd_fu, acd_inci) ~", paste(c("as.factor(extra_num_cat)", "age", "sex","edu"), collapse = "+"))),
                             data = ukb_olink_inci2) %>% cox_format() %>% mutate(label="Extra")
       
       cox_super_death <- coxph(as.formula(paste("Surv(death_fu, death) ~", paste(c("as.factor(super_num_cat)", "age", "sex","edu"), collapse = "+"))),
                                data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Super")
       cox_extra_death<- coxph(as.formula(paste("Surv(death_fu, death) ~", paste(c("as.factor(extra_num_cat)", "age", "sex","edu"), collapse = "+"))),
                               data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Extra")
       
       cox_super_dep <- coxph(as.formula(paste("Surv(dep_fu, dep_inci) ~", paste(c("as.factor(super_num_cat)", "age", "sex","edu"), collapse = "+"))),
                              data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Super")
       cox_extra_dep<- coxph(as.formula(paste("Surv(dep_fu, dep_inci) ~", paste(c("as.factor(extra_num_cat)", "age", "sex","edu"), collapse = "+"))),
                             data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Extra")
       cox_super_ckd <- coxph(as.formula(paste("Surv(ckd_fu, ckd_inci) ~", paste(c("as.factor(super_num_cat)", "age", "sex","edu"), collapse = "+"))),
                              data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Super")
       cox_extra_ckd<- coxph(as.formula(paste("Surv(ckd_fu, ckd_inci) ~", paste(c("as.factor(extra_num_cat)", "age", "sex","edu"), collapse = "+"))),
                             data = ukb_olink_inci2)%>% cox_format() %>% mutate(label="Extra")
       
       
       num_forest <- 
         rbind(cox_super_ad,cox_extra_ad) %>% mutate(Variable=sub(".*\\)", "", Variable),label=factor(label),outcome="All-cause mortality") %>% 
         rbind(.,
               rbind(cox_super_death,cox_extra_death) %>% mutate(Variable=sub(".*\\)", "", Variable),label=factor(label),outcome="All-cause dementia")
         ) %>% 
         rbind(.,
               rbind(cox_super_dep,cox_extra_dep) %>% mutate(Variable=sub(".*\\)", "", Variable),label=factor(label),outcome="Depression")
         )%>% 
         rbind(.,
               rbind(cox_super_ckd,cox_extra_ckd) %>% mutate(Variable=sub(".*\\)", "", Variable),label=factor(label),outcome="Chronic kidney disease")
         )%>% mutate(outcome=factor(outcome,labels = c("All-cause mortality","All-cause dementia","Depression","Chronic kidney disease")),
                     HR=ifelse(HR>10,10,HR)
         ) %>% 
         ggplot(.
                # result_ageo %>% 
                #        filter(Outcome %in% c("death","acd","dep","ckd")) %>% 
                #        mutate(Outcome=factor(Outcome,levels=c("death","acd","dep","ckd"))) %>% 
                #        group_by((label),Outcome) %>% arrange(-HR) %>% mutate(Exposure_label=factor(Exposure_label,levels = rev(unique(Exposure_label))),
                #                                                              Outcome_label=factor(Outcome_label)
                #        )%>% ungroup()
                , 
                # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("All-cause dementia","Alzheimer's disease","Vascular dementia","Parkinson's disease","Multiple sclerosis","All-cause mortality"))),
                #        Exposure_label=factor(Exposure_label,levels=c("Heart","Artery","Organismal","Brain"))
                # ),
                aes(x = HR, y=(Variable),groups=(Variable),shape=factor(label))) +
         geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
         geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
         geom_point(aes(color = Variable),size = 1,position=position_dodge(0.7))+
         scale_x_continuous(limits = c(0.1, 10)) +
         scale_y_discrete()+
         labs(x = "HR (95% CI)", y = "No. of extreme organ aging", title = NULL) +
         theme_classic(base_size = 6,base_line_size = 0.2)+
         facet_grid(~outcome,scales = "free")+
         theme(            legend.position = "right",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                           axis.text.x = element_text(color="black",size = 5),
                           axis.text.y = element_text(color="black",size = 6),
                           legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                           legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                           legend.key.spacing.y = unit(0, "pt"),
                           legend.box.spacing = margin(0),
                           # legend.margin = margin(t = 1, b =1, unit = "pt"),
                           # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                           strip.background = element_blank(),
                           # panel.spacing = unit(0.01, "cm")
                           
         )+scale_color_manual(values = c("1-2"="#4daf4a","3-4"="#377eb8","5+"="#e41a1c"))+
         scale_shape_manual(values = c(17, 16))+scale_x_log10(limits = c(0.2, 10),breaks=c(0.5,1,2,4,8))+
         coord_flip()
       
       
       
       library(nnet)
       multi_logistic <- function(exposure_list, outcome_list, data) {
         # Initialize an empty list to store results
         results <- data.frame()
         
         # Loop over each combination of exposure and outcome
         for (outcome in outcome_list) {
           for (exposure in exposure_list) {
             
             # Fit multinomial logistic regression model
             model <- multinom(as.formula(paste0(outcome, "~", exposure,"+age+sex")), data = data)
             
             # Get model summary
             summary_model <- summary(model)
             
             # Extract odds ratios (OR) and 95% confidence intervals (CI)
             coefs <- summary_model$coefficients [,exposure]
             se <- summary_model$standard.errors [,exposure]
             
             
             
             OR <- exp(coefs)
             CI_lower <- exp(coefs - 1.96 * se)
             CI_upper <- exp(coefs + 1.96 * se)
             
             # Create a data frame to store OR and CI
             result_df <- data.frame(
               Exposure = exposure,
               Outcome = outcome,
               OR = OR,
               CI_Lower = CI_lower,
               CI_Upper = CI_upper
             )
             
             # Store the result in the list
             results<- rbind(results,result_df)
           }
         }
         
         # Return the list of results
         return(results)
       }
       
       multi_results <- multi_logistic(exposure_list =c(organ_pac) ,outcome_list =c("co_neuropsy2","co_chronic2","neuropsy_chronic") ,data = ukb_olink_inci2)
       
       
       table(ukb_olink_inci2$co_neuropsy2)
       table(ukb_olink_inci2$co_chronic2)
       table(ukb_olink_inci2$neuropsy_chronic)
       
       # linear(ukb_olink_inci2,
       #        a=c("num_all"),b=organ_pac,covar = c("age","sex","ethnic","edu","tdi","center"))
       
       
       
       
       multi_result <- multi_results %>% 
         mutate(label=row.names(.)) %>% mutate(label= (gsub("[0-9]", "", label))) %>% 
         mutate(Exposure=gsub("_res", "",Exposure),
                label=case_when(label=="one"~"1",
                                label=="two"~"2+",
                                label=="a_neuro_psy"~"Neuropsychiatric",
                                label=="b_chronic"~"Physical",
                                label=="comorbidity"~"Comorbidity"
                ),
                label=factor(label,levels  = c("1","2+","Neuropsychiatric","Physical","Comorbidity")),
                
                Outcome=case_when(Outcome=="co_neuropsy2"~"No. of neuropsychiatric diseases",
                                  Outcome=="co_chronic2"~"No. of chronic diseases",
                                  Outcome=="neuropsy_chronic"~"Multicomorbidity"
                ),
                
                Outcome=factor(Outcome,levels= c("No. of neuropsychiatric diseases",
                                                 "No. of chronic diseases",
                                                 "Multicomorbidity"
                ))
         )
       
       comorbid_forest <- ggplot(multi_result
                                 # result_ageo %>% 
                                 #        filter(Outcome %in% c("death","acd","dep","ckd")) %>% 
                                 #        mutate(Outcome=factor(Outcome,levels=c("death","acd","dep","ckd"))) %>% 
                                 #        group_by((label),Outcome) %>% arrange(-HR) %>% mutate(Exposure_label=factor(Exposure_label,levels = rev(unique(Exposure_label))),
                                 #                                                              Outcome_label=factor(Outcome_label)
                                 #        )%>% ungroup()
                                 , 
                                 # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("All-cause dementia","Alzheimer's disease","Vascular dementia","Parkinson's disease","Multiple sclerosis","All-cause mortality"))),
                                 #        Exposure_label=factor(Exposure_label,levels=c("Heart","Artery","Organismal","Brain"))
                                 # ),
                                 aes(x = OR, y=(Exposure),groups=(label))) +
         geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
         geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
         geom_point(aes(color=label),size = 1,shape=15,position=position_dodge(0.7))+
         # scale_x_continuous(limits = c(0.8, 2.2),breaks = c(1,1.5,2.0)) +
         scale_y_discrete()+
         labs(x = "OR (95% CI)", y = "Organ age gap", title = NULL) +
         theme_classic(base_size = 6,base_line_size = 0.2)+
         facet_grid(~factor(Outcome),scales = "free")+
         theme(            legend.position = NULL,legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                           axis.text.x = element_text(color="black",size = 5),
                           axis.text.y = element_text(color="black",size = 6),
                           legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
                           legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                           legend.key.spacing.y = unit(0, "pt"),
                           legend.box.spacing = margin(0),
                           # legend.margin = margin(t = 1, b =1, unit = "pt"),
                           # legend.box.margin = margin(t = 1, b = 1, unit = "pt"),
                           plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                           strip.background = element_blank(),
                           # panel.spacing = unit(0.01, "cm")
                           
         )+scale_color_manual(values = 
                                c("1"="#fee08b",
                                  "2+"="#d53e4f",
                                  "Neuropsychiatric"="#377eb8",
                                  "Physical"="#984ea3",
                                  "Comorbidity"="#cb181d"
                                ))
       
       
       library(nnet)
       multinom(co_neuropsy ~ Brain_res+age+sex, data = ukb_olink_inci2) %>% summary()
       multinom(co_neuropsy ~ Organismal_res+age+sex, data = ukb_olink_inci2) %>% summary()
       
       multinom(co_chronic ~ Brain_res, data = ukb_olink_inci2) %>% summary()
       multinom(co_chronic ~ Organismal_res, data = ukb_olink_inci2) %>% summary()
       
       multinom(neuropsy_chronic ~ Brain_res, data = ukb_olink_inci2) %>% summary()
       summary_model <- multinom(neuropsy_chronic ~ Organismal_res+age+sex, data = ukb_olink_inci2) %>% summary()
       
       
       
       
       
       #Fig 2 ageotype and disease, mulicomorbidity, and death
       
       fig2 <-grid.grabExpr(draw(super+extra,padding = unit(c(0.5, 0.5, 0.5, 0.5), "mm"))) %>% 
         plot_grid(.,ageo_forest,num_forest,comorbid_forest,nrow = 4,align="h",rel_heights = c(1.5,1,1,1.2),
                   labels = "auto",label_size = 6)
       
       ggsave(filename = ".\\Project\\Proteomic aging\\Submission\\Submission\\Nat Aging\\Revision\\2-16 R1\\Figure\\Fig.2.pdf",
              plot = fig2,
              device = cairo_pdf,width =5.5,height=6.2,units="in")
       