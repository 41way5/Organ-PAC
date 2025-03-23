# xFig2
  #Organ aging model performance across cohorts---------------------

ckb_pred <-read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/CKB/validation/ckb_preds_df.csv")
nhs1_pred <- read.csv("E:\\Tutor\\UK Biobank\\Project\\Proteomic aging\\Submission\\Submission\\Nat Aging\\Revision\\2-16 R1\\Validation\\NHS\\nhs1\\nhs1\\ckb_preds_df.csv")

#Read in UKB model data (dot plot and SHAP plot)

color_panel <- c("Black", "#A52A2A", "#4B0082", "#FF7F50", "#FFFF00", "#808080", "#808000",
                 "#D2B48C", "#00FFFF", "#0000FF", "#008000")

organ_list

plot_organ_age <- function(data, organ_list, color_panel,label=T) {
  plot_list <- list()
  
  for (i in seq_along(organ_list)) {
    organ <- organ_list[i]
    color <- color_panel[i]  
    
    model <- lm(reformulate("age_at_recruitment", response = organ), data = data)
    r2 <- summary(model)$r.squared
    r <- sqrt(r2)
    mse <- mean(model$residuals^2)
    
    if (label) {y=paste0(organ)
      
    }else{y=NULL}
    

    x <- if (organ %in% c("Pancreas", "Intestine")) "Chronological age" else NULL
    
    
    plot_list[[organ]] <- ggplot(data, aes(x = age_at_recruitment, y = .data[[organ]])) +
      geom_point( color = color, size=0.3, pch=20) +  
      geom_smooth(method = "lm", se = FALSE, color = ifelse(organ=="Organismal","white","black"), size = 0.6) +  # Match color
     
      labs(x = x, 
           y = y) +
      
      # ggtitle(organ) +
      theme_bw(base_size = 6,base_line_size = 0.2)+
      annotate("text", 
               x = max(data$age_at_recruitment, na.rm = TRUE), 
               y = min(data[[organ]], na.rm = TRUE), 
               label = sprintf("R = %.2f\nR² = %.2f\nMSE = %.2f", r, r2, mse), 
               hjust = 1, vjust = 0, size = 1.5,lineheight = 0.8)+
      theme(legend.position = "none",
            panel.grid = element_blank(),
            axis.title.y = element_text(face = "bold")
            
            )
  }
  
  return(plot_list)
}

plot_organ_shape <- function(organ_list) {
  library(ggbeeswarm) 
  plot_list <- list()
  
  for (i in seq_along(organ_list)) {
    organ <- organ_list[i]
    
    shap_data <- fread(paste0("./Data/Proteomic/Data/shap/shap_",organ,".csv"))
    
    feature_importance <- shap_data[, lapply(.SD, function(x) mean(abs(x), na.rm = TRUE)), .SDcols = -"eid"]
    top_features <- names(sort(feature_importance, decreasing = TRUE))[1:20] 
    top_features <- top_features[!is.na(top_features)] 
    
    shap_long <- melt(shap_data, id.vars = "eid", variable.name = "Feature", value.name = "SHAP")
    shap_long <- shap_long[Feature %in% top_features]
    
    
    value <- olink_full %>% dplyr::select(all_of(top_features),eid) %>% filter(eid %in% (shap_long$eid)) %>% normalize_median_center_df()
    value_long <- melt(value, id.vars = "eid", variable.name = "Feature", value.name = "value")
    
    shap_long <- shap_long %>% left_join(.,value_long)
    
    value_median <- median(shap_long$value, na.rm = TRUE)
    value_low <- quantile(shap_long$value, 0.1, na.rm = TRUE)
    value_high <- quantile(shap_long$value, 0.9, na.rm = TRUE)
    
    x <- if (organ %in% c("Pancreas", "Intestine")) "SHAP" else NULL
    
    plot_list[[organ]] <- ggplot(shap_long, aes(x = SHAP, y = reorder(Feature, abs(SHAP), FUN = median), color = value)) +
      geom_violin(scale = "width", fill = NA, color = "white",alpha = 0.7) + 
      geom_quasirandom(size = 0.2) +  
      scale_color_gradient2(low = "#2b83ba", mid = "white", high = "#d7191c", 
                            midpoint = value_median, 
                            limits = c(value_low, value_high), oob = scales::squish,
                            na.value = "white"
      ) +
      theme_bw(base_size = 5,base_line_size = 0.2) +
      labs(x = x, y = NULL) +
      theme(legend.position = "none",
            panel.grid = element_blank()
      )
  }
  
  return(plot_list)
}

#Dot plot for UKB
plot_list_ukb <- plot_organ_age(data=ukb_olink,organ_list = organ_list,color_panel = color_panel)
#Dot plot for CKB
plot_list_ckb <- plot_organ_age(data=ckb_pred,organ_list = organ_list,color_panel = color_panel,label = F)
#Dot plot for NHS
plot_list_nhs <- plot_organ_age(data=nhs1_pred,organ_list = organ_list,color_panel = color_panel,label = F)
#SHAP plot
plot_list_shape <- plot_organ_shape(organ_list=organ_list)


xf2 <- plot_grid(plot_grid(
             plot_grid(plotlist = plot_list_ukb[1:6],ncol=1,labels = "auto",label_size = 7), 
             plot_grid(plotlist = plot_list_shape[1:6],ncol=1), 
             plot_grid(plotlist = plot_list_ckb[1:6],ncol=1), 
             plot_grid(plotlist = plot_list_nhs[1:6],ncol=1),
             align = "hv",labels = c("i","ii","iii","iv"),label_size = 6,
             ncol = 4
             ),
       plot_grid(
             plot_grid(NULL, plot_grid(plotlist = plot_list_ukb[7:11],ncol=1,labels = c("g","h","i","j","k"),label_size = 7),ncol = 1,rel_heights = c(1,5)),
             plot_grid(NULL, plot_grid(plotlist = plot_list_shape[7:11],ncol=1),ncol = 1,rel_heights = c(1,5)),
             plot_grid(NULL, plot_grid(plotlist = plot_list_ckb[7:11],ncol=1),ncol = 1,rel_heights = c(1,5)),
             plot_grid(NULL, plot_grid(plotlist = plot_list_nhs[7:11],ncol=1),ncol = 1,rel_heights = c(1,5)),
             align = "hv",
             ncol = 4
      ),align = "h"
      
      )+theme(plot.margin = unit(c(0.1, 0.3, 0.1, 0.1), "cm"))
       


ggsave(filename = "Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Figure/Final/xFig.2.pdf",
       plot = xf2,
       device = cairo_pdf,width =8,height=8.2,units="in")

ggsave(filename = "Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Figure/Final/xFig.2.jpg",
       plot = xf2,
       device = "jpeg", width = 8, height = 8.2, units = "in")








