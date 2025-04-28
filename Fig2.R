
#Fig2------------------------------------------------------------------------------------------------


#Cor between PAC and traits at baseline 
library(ComplexHeatmap)
library(cowplot)
library(BioAge)
library(dplyr)
library(smplot2)
library(ggpubr)

calculate_beta_matrix <- function(data, var_list1, var_list2, covariates = NULL) {
  beta_matrix <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  fdr_p_values <- matrix(NA, nrow = length(var_list1), ncol = length(var_list2))
  
  for (i in seq_along(var_list1)) {
    for (j in seq_along(var_list2)) {
      formula <- as.formula(paste0(var_list1[i], "~", "`",var_list2[j],"`"))
      if (!is.null(covariates)) {
        covariate_formula <- paste(covariates, collapse = " + ")
        formula <- update(formula, paste(". ~ . +", covariate_formula))
      }
      lm_model <- lm(formula, data = data)
      beta_matrix[i, j] <- coef(lm_model)[2]  # Store the beta coefficient
      
      p_value <- summary(lm_model)$coefficients[2, "Pr(>|t|)"]
      fdr_p_values[i, j] <- p.adjust(p_value, method = "fdr")
      
    }
  }
  
  colnames(beta_matrix) <- var_list2
  rownames(beta_matrix) <- var_list1
  
  return(list(beta_matrix = beta_matrix, fdr_p_values = fdr_p_values))
}

beta_matrix_all <- calculate_beta_matrix(ukbb_trait, 
                                         var_list1=c(organ_pac,"phenoage_res","kdm_res"), 
                                         var_list2=ukbb_trait_list$fid,
                                         covariates = c("age","sex","ethnic","edu","tdi","center","pa","smok"))

beta_matrix_nmr <- calculate_beta_matrix(ukbb_trait_nmr %>% select(1:23,any_of(nmr_in$fid)), 
                                         var_list1=c(organ_pac,"phenoage_res","kdm_res"), 
                                         var_list2=nmr_in$fid,
                                         covariates = c("age","sex","ethnic","edu","tdi","center","pa"))

#Filter only meaningful NMR traits

nmr_in <- ukbb_trait_list_NMR %>% filter(!str_detect(Abbr, "-(FC|CE|PL|L)$"))

group_list <- list(
  
  ` `=c( "ApoA1", "ApoB", "ApoB/ApoA1","Total Cholesterol", "Total Cholesterol Minus HDL-C", "TG","LDL-P","IDL-P","HDL-P"),
  
  `VLDL`=c( "VLDL-C", "XL-VLDL-C", "L-VLDL-C", "M-VLDL-C", "S-VLDL-C", "XS-VLDL-C",
         "VLDL-TG", "XXL-VLDL-TG", "XL-VLDL-TG","L-VLDL-TG", "M-VLDL-TG", "S-VLDL-TG","XS-VLDL-TG"),
  
  `IDL`=c("IDL-C", "IDL-TG"),
  
  `LDL`=c("LDL-C", "L-LDL-C", "M-LDL-C", "S-LDL-C", "XL-LDL-C", 
          "LDL-TG","L-LDL-TG", "M-LDL-TG", "S-LDL-TG"),
  
  `HDL`=c("HDL-C", "L-HDL-C", "M-HDL-C", "S-HDL-C", "XL-HDL-C",
          "HDL-TG", "L-HDL-TG", "M-HDL-TG", "S-HDL-TG", "XL-HDL-TG"),
  
  `Fatty acids` = c("Omega-3 Fatty Acids", "Omega-6 Fatty Acids", "Linoleic Acid", 
                    "Docosahexaenoic Acid", "Saturated Fatty Acids", "Degree of Unsaturation",
                    "Monounsaturated Fatty Acids", "Polyunsaturated Fatty Acids", 
                    "Total Fatty Acids"),
  
  `Ketone bodies` = c( "Acetate", "Acetoacetate", "Acetone","3-Hydroxybutyrate"),
  
  `Amino acids` = c("Alanine", "Glutamine", "Glycine", "Histidine", "Isoleucine", 
                    "Leucine", "Phenylalanine", "Tyrosine", "Valine", "Branched-Chain Amino Acids"),
  
  `Glycolysis` = c("Glucose", "Glucose-lactate", "Lactate", "Pyruvate", "Citrate", 
                   "Spectrometer-corrected alanine"),
  
  `Other` = c("Albumin", "Creatinine","Glycoprotein Acetyls")
)

group_df <- stack(group_list) %>% rename(values="Abbr", ind="Group")

nmr_in$Abbr <- factor(nmr_in$Abbr, levels = unlist(group_list, use.names = FALSE))

nmr_in <- nmr_in %>% filter(Abbr %in% unlist(group_list, use.names = FALSE)) %>%
                arrange(match(Abbr, unlist(group_list, use.names = FALSE))) %>% 
                left_join(group_df, by = "Abbr") %>% 
                mutate(Abbr=case_when(Abbr=="Total Cholesterol Minus HDL-C"~"nonHDL-C",
                                      Abbr=="Omega-3 Fatty Acids"~"Omega-3",
                                      Abbr=="Omega-6 Fatty Acids"~"Omega-6",
                                      Abbr=="Docosahexaenoic Acid"~"DHA",
                                      Abbr=="Linoleic Acid"~"LA",
                                      Abbr=="Saturated Fatty Acids"~"SFA",
                                      Abbr=="Degree of Unsaturation"~"Unsaturation",
                                      Abbr=="Monounsaturated Fatty Acids"~"MUFA",
                                      Abbr=="Polyunsaturated Fatty Acids"~"PUFA",
                                      Abbr=="Total Fatty Acids"~"Total FA",
                                      Abbr=="Branched-Chain Amino Acids"~"BCAA",
                                      Abbr=="Spectrometer-corrected alanine"~"Alanine",
                                      T~Abbr
                                      ))








# Corr between PAC and other aging clock, and PAC and baseline traits
colnames(beta_matrix_all$beta_matrix) <- trait_base
rownames(beta_matrix_all$beta_matrix) <- c(organ_list,"Phenoage gap","KDMage gap")


colnames(beta_matrix_nmr$beta_matrix) <- nmr_in$Abbr
rownames(beta_matrix_nmr$beta_matrix) <- c(organ_list,"Phenoage gap","KDMage gap")



mat_beta <- (beta_matrix_all$beta_matrix) 
mat_p <- (beta_matrix_all$fdr_p_values)

mat_beta_nmr <- (beta_matrix_nmr$beta_matrix) 
mat_p_nmr <- (beta_matrix_nmr$fdr_p_values)

ht_opt$COLUMN_ANNO_PADDING = unit(0.5, "mm")
xf2_4 <-grid.grabExpr(draw( 
  Heatmap(mat_beta,
          name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
          col = colorRamp2(c(-0.1,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
          show_row_names = TRUE, 
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 5),
          cluster_rows = F,
          cluster_columns = F,
          column_names_rot=45,
          column_split=ukbb_trait_list$Categroy_short,
          column_title_gp = gpar(fontsize = 5),
          row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
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
            # at = c(-0.1, 0, 0.5),
            # labels = c(-0.1, 0, 0.5),
            title = "r",
            legend_height = unit(1, "cm"),
            grid_width = unit(0.5, "mm"),
            legend_width =unit(1, "mm"),
            # title_position = "lefttop-rot",
            labels_gp = gpar(col = "black", fontsize = 5),
            title_gp = gpar(col = "black", fontsize = 5)
            
          ),
          width = ncol(mat_beta_nmr)*unit(2, "mm"),
          height =nrow(mat_beta_nmr)*unit(2, "mm")
          
  )
),padding = unit(c(1, 1, 0, 1), "mm")
)

xf2_5 <-grid.grabExpr(draw( 
  Heatmap(mat_beta_nmr,
          name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
          col = colorRamp2(c(-0.1,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
          show_row_names = TRUE, 
          show_column_names = TRUE,
          row_names_gp = gpar(fontsize = 5),
          column_names_gp = gpar(fontsize = 4.5),
          cluster_rows = F,
          cluster_columns = F,
          column_names_rot=45,
          # column_split=ukbb_trait_list$Categroy_short,
          column_title_gp = gpar(fontsize = 5),
          column_split=nmr_in$Group,
          row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
          column_title_side = "top",
          row_names_side = "left",
          row_title = NULL,
          #column_title = NULL,
          cell_fun = function(j, i, x, y, w, h, fill) {
            if((mat_p_nmr)[i, j] < 0.001) {
              grid.text("***", x, y,gp = gpar(fontsize = 2))
            } else if((mat_p_nmr)[i, j] < 0.01) {
              grid.text("**", x, y,gp = gpar(fontsize = 2))
            } else if((mat_p_nmr)[i, j] < 0.05) {
              grid.text("*", x, y,gp = gpar(fontsize = 2))
            }
          },
          heatmap_legend_param = list(
            # at = c(-0.1, 0, 0.5),
            # labels = c(-0.1, 0, 0.5),
            title = "r",
            legend_height = unit(1, "cm"),
            grid_width = unit(0.5, "mm"),
            legend_width =unit(1, "mm"),
            # title_position = "lefttop-rot",
            labels_gp = gpar(col = "black", fontsize = 5),
            title_gp = gpar(col = "black", fontsize = 5)
          ),
          width = ncol(mat_beta_nmr)*unit(2, "mm"),
          height =nrow(mat_beta_nmr)*unit(2, "mm")
          
  ),padding = unit(c(0, 1, 0, 1), "mm")
)
)
plot_grid(xf2_4)
plot_grid(xf2_5)


library(circlize) 
mat_beta_nmr <- (beta_matrix_nmr$beta_matrix) %>% t()
mat_p_nmr <- (beta_matrix_nmr$fdr_p_values)%>% t()

circle <- function() {
  circos.clear()
  par(mar = c(0, 0, 0, 0))
  circos.par(start.degree = 90-15, gap.after = c(15), cell.padding = c(0.01, 0, 0, 0))
  

  circos.heatmap(
    mat_beta_nmr,  
    col = colorRamp2(c(-0.1, 0, 0.5), c("#4DBBD5B2", "white", "#DC0000B2")),
    track.height = 0.35,
    rownames.side = "inside",
    rownames.cex = 0.2,
    rownames.font = 0.5,
    cell.border = "black",
    cell.lwd = 0.2,
    ignore.white = FALSE,
    track.margin = c(0, 0)
  )
  
  circos.track(
    track.index = 1,
    panel.fun = function(x, y) {
      if (CELL_META$sector.numeric.index == 1) {
        cn <- colnames(mat_beta_nmr)
        n <- length(cn)
        circos.text(
          rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"),
          1:n - 0.5,
          cn,
          cex = 0.25,
          adj = c(0, 0.5),
          facing = "inside"
        )
      }
    },
    bg.border = NA
  )
  
  for (i in 1:ncol(mat_beta_nmr)) {
    for (j in 1:nrow(mat_beta_nmr)) {
      pval <- mat_p_nmr[j, i] 
      if (!is.na(pval)) {
        stars <- if (pval < 0.001) {
          "***"
        } else if (pval < 0.01) {
          "**"
        } else if (pval < 0.05) {
          "*"
        } else {
          ""
        }
        
        if (stars != "") {
          circos.text(
            x = j - 0.5,  
            y = i - 0.5,  
            labels = stars,
            cex = 0.1,    
            col = "black",
            facing = "inside"
          )
        }
      }
    }
  }
}


circlize_grob <- cowplot::as_grob(circle)




# circle <- function(){
#   circos.clear()
#   par(mar = c(0, 0, 0, 0) + 0.2)
#   circos.par(gap.after = c(15),cell.padding = c(0.01, 0, 0, 0))
#   circos.heatmap(mat_beta_nmr %>% t(), col = colorRamp2(c(-0.1,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
#                  track.height = 0.3,rownames.side = "inside",rownames.cex=0.2,rownames.font = 0.5,cell.border = "black",cell.lwd = 0.2
#                  ,ignore.white = F,track.margin = c(0,0)
#   )
#   
#   circos.track(track.index =1, panel.fun = function(x, y) {
#     if(CELL_META$sector.numeric.index == 1) { # the last sector
#       cn = colnames(mat_beta_nmr %>% t())
#       n = length(cn)
#       circos.text( rep(CELL_META$cell.xlim[2], n) + convert_x(1, "mm"), 1:n-0.5, cn, 
#                    cex = 0.25, adj = c(0, 0.5), facing = "inside")
#     }
#   }, bg.border = NA)
#   
# }
# circlize_grob <-  cowplot::as_grob(circle)
# plot_grid(circlize_grob)

show_col(pal_npg("nrc")(10))
xf2_1 <- ggplot(ukb_olink %>% filter(!is.na(phenoage)), aes(x = age, y = phenoage)) +
  geom_point(size = 1, alpha = 0.2,color="#8491B4FF") +  # Add points
  labs(title = "",x = "Chronological age",y = "KDMAge"
  )+
  theme(plot.caption = element_text(hjust = 1, vjust = 1, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
        )+
  sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2,size=0.2)+theme_classic(base_size = 6,base_line_size = 0.2)

xf2_2 <- ggplot(ukb_olink %>% filter(!is.na(kdm)), aes(x = age, y = kdm)) +
  geom_point(size = 1, alpha = 0.2,color="#8491B4FF") +  # Add points
  labs(title = "",x = "Chronological age",y = "PhenoAge"
  )+
  theme(plot.caption = element_text(hjust = 1, vjust = 1, color = "black"),
        plot.margin = unit(c(0, 0, 0, 0), "mm")
        )+
  sm_statCorr(color = 'black', corr_method = 'pearson',text_size = 2,size=0.2)+theme_classic(base_size = 6,base_line_size = 0.2)


age_corr <- ukb_olink %>%
  dplyr::select(c("phenoage_res","kdm_res",organ_pac)) %>% cor(use="na.or.complete")
mat_beta <-age_corr[c(1,2),c(-1,-2)]
colnames(mat_beta) <- c(organ_list)
rownames(mat_beta) <- c("PhenoAge gap","KDMAge gap")

xf2_3 <- Heatmap(mat_beta,show_heatmap_legend = FALSE,
                 name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                 col = colorRamp2(c(-0.2,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
                 show_row_names = TRUE, 
                 show_column_names = TRUE,
                 row_names_gp = gpar(fontsize = 5),
                 column_names_gp = gpar(fontsize = 5),
                 cluster_rows = F,
                 cluster_columns = F,
                 column_names_rot=45,
                 # row_split=c("","",rep("ProtAge gap",11)),
                 row_title_gp = gpar(fontsize = 5),
                 column_title_side = "bottom",
                 row_names_side = "left",
                 column_title = NULL,
                 cell_fun = function(j, i, x, y, width, height, fill) {if(i==1|i==2)
                   grid.text(sprintf("%.1f", mat_beta[i, j]), x, y, gp = gpar(fontsize = 5))
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
                 # width = ncol(mat_beta)*unit(3, "mm"), 
                 height = nrow(mat_beta)*unit(4, "mm")
                 
)

# plot_grid(circlize_grob,xf2_6,nrow = 1,rel_widths = c(4,1))

# xf2 <- plot_grid(xf2_1,xf2_2,grid.grabExpr(draw(xf2_3)),ncol=3,nrow=1,labels = "auto",label_size = 6,rel_widths = c(1,1,1.5)) %>% 
#         plot_grid(.,plot_grid(xf2_4, plot_grid(circlize_grob),nrow = 2,align = "v",rel_heights  = c(1,1.5),labels = c("d","e"),label_size = 6),nrow=2,rel_heights = c(0.6,2.2)) %>% 
#           plot_grid(.,xf2_6,ncol=1,rel_heights = c(1,0.2))

# 
# xf2 <- plot_grid(xf2_1,xf2_2,grid.grabExpr(draw(xf2_3)),ncol=,nrow=1,labels = "auto",label_size = 6,rel_widths = c(1,1,1.5)) %>% 
#   plot_grid(.,plot_grid(xf2_4, plot_grid(circlize_grob,plot_grid(xf2_6),nrow = 1,rel_widths = c(1,1),rel_heights = c(3,1)),nrow = 2,align = "v",rel_heights  = c(1,1.5),labels = c("d","e"),label_size = 6),nrow=2,rel_heights = c(0.6,2.2))


# xf2 <- plot_grid(xf2_4,ncol=1,rel_heights = c(1,1.7),labels = c("c",""),label_size = 8,
# plot_grid(plot_grid(plot_grid(plot_grid(xf2_1,xf2_2)+theme(plot.margin = margin(0, 0, 0, 0, "cm")),labels = c("a","b"),label_size = 8,
#                               grid.grabExpr(draw(xf2_3)),nrow=1),xf2_6,ncol=1,labels = c("","e"),label_size = 8
#                     )+theme(plot.margin = margin(0, 1, 0, 0, "mm")),circlize_grob,nrow=1,rel_widths = c(1,1),labels = c("","d"),label_size = 8)
# )


f2 <- plot_grid(plot_grid(plot_grid(xf2_1,xf2_2)+theme(plot.margin = margin(0, 0, 0, 0, "cm")),labels = c("a","b"),label_size = 8,grid.grabExpr(draw(xf2_3)),nrow=1),
          plot_grid(xf2_4),
          plot_grid(xf2_5),
          xf2_6,labels = c("","c","d","e"),
          ncol=1,rel_heights = c(0.5,1,1,0.8),
          label_size = 8)

ggsave(filename = ".\\Project\\Proteomic aging\\Submission\\Submission\\Nat Aging\\Revision\\2-16 R1\\Figure\\Fig.2.pdf",
       plot = f2,
       device = cairo_pdf,width =7,height=7.3,units="in")

###Association with healthy lifestyle


ukb_lifestyle <- ukb_olink %>% #mutate(pa=ifelse(is.na(pa),0,1)) %>% 
  dplyr::select(eid,age,sex,edu,tdi,center,ethnic,
                organ_pac,"phenoage_res","kdm_res",
                smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat,lifestyle) %>%
mutate(lifestyle=rowSums(dplyr::select(., smoke_cat,alcohol_cat,pa_cat,tv_cat,sleep_cat,fruit_cat,fish_cat,red_cat,process_cat)),
       lifestyle_cat=ifelse(lifestyle<=2,"Favorable",ifelse(lifestyle<=5,"Intermediate","Unfavorable")),
       lifestyle_cat=factor(lifestyle_cat, levels = c("Unfavorable", "Favorable","Intermediate"))
)

describe(ukb_lifestyle$lifestyle)
describe(ukb_lifestyle$lifestyle_cat)

library(broom)
lm(Brain_res~lifestyle_cat+age+sex,data=ukb_lifestyle) %>% tidy()
coxph(formula =Surv(death_fu, death) ~ lifestyle_cat + age + sex + ethnic + edu + tdi + center, data=ukb_lifestyle) %>% summary()




linear2 <- function(data, a, b,covar=NULL) {
  
  result_list <- data.frame()
  # Iterate over all combinations of variables in lists 'a' and 'b'
  for (var_a_name in a) {
    for (var_b_name in b) {
      
      formula_str <- paste(var_a_name, "~", var_b_name)
      if (!is.null(covar)) {
        formula_str <- paste(formula_str, paste(covar, collapse = " + "), sep = " + ")
      }
      formula <- as.formula(formula_str)
      model <- lm(formula, data = data) %>% tidy() %>% #filter(term %in% var_b_name) %>% 
        mutate(outcome=var_a_name)
      
      result_list <- rbind(result_list,model)
    }
    
  }
  
  
  
  return(result_list)
}

lifestyle_result <- linear2(ukb_lifestyle %>% filter(),
                            b=c("lifestyle_cat"),
                            a=organ_pac,
                            covar = c("age","sex","ethnic","edu","tdi","center"))


lifestyle_result2 <- lifestyle_result %>% filter(term %in% c("lifestyle_catFavorable","lifestyle_catIntermediate")) %>% 
  mutate(term=ifelse(term=="lifestyle_catFavorable","Favorable","Intermediate")) %>% 
  mutate(term=factor(term,levels=c("Intermediate","Favorable"))) %>% 
  mutate(outcome=sub("_res", "", outcome))





lifestyle_cat_plot <- ggplot(lifestyle_result2
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
                             aes(x = estimate, y=(outcome),groups=(term))) +
  geom_vline(xintercept = 0, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = estimate-1.96*std.error, xmax = estimate+1.96*std.error), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color=term),size = 1,shape=15,position=position_dodge(0.7))+
  scale_color_manual(
    values = c(
      "Intermediate" = "#FC8D59",
      "Favorable" = "#91BFDB"
    ),
    labels = c(
      "Intermediate" = "Intermediate (3-5 unhealthy lifestyle factors)",
      "Favorable" = "Favorable (0-2 unhealthy lifestyle factors)"
    )
  ) +
  # scale_x_continuous(limits = c(0.9, 2),breaks = c(1,1.5,2.0)) +
  scale_y_discrete()+
  labs(x = "Beta (95% CI)", y = "Organ age gap", title = "Healthy lifestyle category and organ age gaps") +
  theme_classic(base_size = 5,base_line_size = 0.2)+
  # facet_grid(~factor(Outcome),scales = "free")+
  theme(            plot.title = element_text(hjust = 0.5),
                    legend.position = "top",
                    legend.title =element_blank(),
                    legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 5,angle = 45,vjust = 0.65),
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
                    
  )+coord_flip()



#Association between pac and unhealthy lifestyle

beta_matrix_lifestyle <- calculate_beta_matrix(ukb_lifestyle, 
                                               var_list1=c(organ_pac,"phenoage_res","kdm_res"), 
                                               var_list2=c("smoke_cat","alcohol_cat","pa_cat","tv_cat","sleep_cat","fruit_cat","fish_cat","red_cat","process_cat","lifestyle"),
                                               covariates = c("age","sex"
                                                              # "age","sex","ethnic","edu","tdi","center"
                                               ))

colnames(beta_matrix_lifestyle$beta_matrix) <- c("Smoking","Alcohol","Physical activity","TV viewing",
                                                 "Sleep duration","Fruit and vegetable","Oily fish","Red meat",
                                                 "Processed meat","No. of unhealthy lifestyle"
)

rownames(beta_matrix_lifestyle$beta_matrix) <- c(organ_list,"Phenoage gap","KDMage gap")

mat_beta_life <- (beta_matrix_lifestyle$beta_matrix) 
mat_p_life  <- (beta_matrix_lifestyle$fdr_p_values)


xf2_6 <-  plot_grid(
  # ggmatrix_gtable(neuro_plot),
  grid.grabExpr(draw( 
    Heatmap(mat_beta_life,
            name = "beta",border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
            col = colorRamp2(c(-0.5,0,0.5), c("#4DBBD5B2","white","#DC0000B2")),
            show_row_names = TRUE,
            show_column_names = TRUE,
            row_names_gp = gpar(fontsize = 5),
            column_names_gp = gpar(fontsize = 5),
            cluster_rows = F,
            cluster_columns = F,
            column_names_rot=45,
            # column_split=ukbb_trait_list$Categroy_short,
            column_title_gp = gpar(fontsize = 5),
            row_split = factor(c(rep("Organ-specific aging", each = 11),rep("Phenoage", each = 2)),levels = c("Organ-specific aging","Phenoage")),
            column_title_side = "top",
            row_names_side = "left",
            row_title = NULL,
            #column_title = NULL,
            cell_fun = function(j, i, x, y, w, h, fill) {
              if((mat_p_life)[i, j] < 0.001) {
                grid.text("***", x, y,gp = gpar(fontsize = 4))
              } else if((mat_p_life)[i, j] < 0.01) {
                grid.text("**", x, y,gp = gpar(fontsize = 4))
              } else if((mat_p_life)[i, j] < 0.05) {
                grid.text("*", x, y,gp = gpar(fontsize = 4))
              }
            },
            heatmap_legend_param = list(
              # at = c(-0.1, 0, 0.5),
              # labels = c(-0.1, 0, 0.5),
              title = "r",
              legend_direction = "horizontal",
              legend_height = unit(0.5, "mm"),
              grid_width = unit(0.2, "mm"),
              grid_height = unit(0.5, "mm"),
              legend_width =unit(1, "cm"),
              title_position = "leftcenter",
              labels_gp = gpar(col = "black", fontsize = 5),
              title_gp  = gpar(col = "black",fontsize = 6)
            ),
            width = ncol(mat_beta_life)*unit(3, "mm"),
            height =nrow(mat_beta_life)*unit(2, "mm")
            
    ),padding = unit(c(1, 1, 1, 1), "mm"),heatmap_legend_side = "top"
  )) %>% plot_grid(.,lifestyle_cat_plot,nrow =1,rel_widths = c(1,1))
  # ,rel_widths = c(1,1),labels = c("f",""),label_size = 6
)




###Association with Neurobiomarkers
library(broom)
library(GGally)
ukb_neuro <- ukb_olink %>% filter(!is.na(ab40_1)) %>% 
  mutate(brain_ageo=ifelse(Brain_res_super==1,"Super",ifelse(Brain_res_extra==1,"Extra","Normal")),
         brain_ageo=factor(brain_ageo, levels = c("Super","Normal", "Extra")),
         ab42_40_1=ab42_1/ab40_1,
         ab42_40_2=ab42_2/ab40_2                
  ) %>% 
  mutate(across(ends_with(c("_1", "_2")), ~ as.numeric(scale(as.numeric(.))), .names = "{.col}")) %>% filter(acd_date==""&pd_date=="")




neuro_plot <- ggpairs(ukb_neuro %>% dplyr::select(ends_with(c("_1")),"ad_prs","Brain_res") %>% 
                        dplyr::rename(Aβ40=ab40_1,Aβ42=ab42_1,GFAP=gfap_1,NEFL=nefl_1,pTau181=pTau181_1,`Aβ42 / Aβ40`=ab42_40_1,`AD PRS`=ad_prs, Brain=Brain_res)
                      ,
                      upper = list(continuous = wrap("cor", size = 1,title=NULL)), 
                      diag = list(continuous = wrap("barDiag", bins = 15,color="black",fill="grey",lwd=0.1)),
                      lower = list(continuous = wrap("points", alpha = 0.3, size = 0.5)))+
  theme_bw(base_size = 6,base_line_size = 0.2)+
  theme(
    text = element_text(size = 5,color = "black"),
    strip.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # axis.text = element_text(size = 8),
    # axis.title = element_text(size = 12),
    # plot.title = element_text(size = 14, face = "bold")
  )






