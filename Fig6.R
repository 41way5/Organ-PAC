#Fig 6-----------------------------------------------
# Load necessary libraries
library(pROC)
library(ggplot2)
library(gridExtra)
library(cowplot)

#ROC curve with AUC
  ##Age+sex comparision

roc_brain2 <- function(data,outcome){
  data[[outcome]] <- factor(data[[outcome]])
  model_demo <- glm(paste(outcome, "~ age+sex+edu"), data = data, family = binomial)
  
  model_brain <- glm(paste(outcome, "~ Brain_res"), data = data, family = binomial)
  model_pheno <- glm(paste(outcome, "~ phenoage_res"), data = data, family = binomial)
  model_organ <- glm(paste(outcome, "~ Organismal_res"), data = data, family = binomial)
  
  
  
  model_demo_brain <- glm(paste(outcome, "~ age+sex+edu+Brain_res"), data = data, family = binomial)
  model_demo_pheno <- glm(paste(outcome, "~ age+sex+edu+phenoage_res"), data = data, family = binomial)
  model_demo_organ <- glm(paste(outcome, "~ age+sex+edu+Organismal_res"), data = data, family = binomial)
  
  
  
  
  # model_cog <- glm(paste(outcome, "~ age+sex+edu+rt"), data = data, family = binomial)
  # model_pheno <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res"), data = data, family = binomial)
  # model_pheno_brain <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+Brain_res"), data = data, family = binomial)
  
  model_pheno_organ <- glm(paste(outcome, "~ age+sex+edu+rt+Organismal_res+phenoage_res"), data = data, family = binomial)
  model_pheno2 <- glm(paste(outcome, "~ age+sex+edu+Brain_res+phenoage_res"), data = data, family = binomial)
  
  model_prs <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+ad_prs"), data = data, family = binomial)
  model_prs_brain <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+ad_prs+Brain_res"), data = data, family = binomial)
  
  
  model_organ <- glm(paste(outcome, "~ age+sex+Organismal_res"), data = data, family = binomial)
  
  
  
  
  model_all <- glm(paste(outcome, "~ age+sex+edu+",paste(organ_pac,collapse = "+")), data = data, family = binomial)
  model_pheno3 <- glm(paste(outcome, "~ age+sex+edu+phenoage_res"), data = data, family = binomial)
  
  
  get_predictions <- function(model, data) {
    predict(model, newdata = data, type = "response")
  }
  
  prediction <- data.frame(
    outcome = data[[outcome]],
    demo = get_predictions(model_demo, data),
    brain= get_predictions(model_brain, data),
    pheno = get_predictions(model_pheno, data),
    demo_brain = get_predictions(model_demo_brain, data),
    demo_pheno = get_predictions(model_demo_pheno, data)
    
    # cog = get_predictions(model_cog, data),
    # brain1 = get_predictions(model_brain1, data),
    # brain2 = get_predictions(model_brain2, data),
    # pheno = get_predictions(model_pheno, data),
    # pheno_brain = get_predictions(model_pheno_brain, data),
    # prs=get_predictions(model_prs, data),
    # prs_brain=get_predictions(model_prs_brain, data)
  )
  
  if(outcome %in% c("ad_inci","acd_inci","vd_inci")){
    
    prediction_roc <- list(
      `Age+sex+edu` = roc(outcome ~ demo, data=prediction),
      `Brain age gap`=roc(outcome ~ brain, data=prediction),
      `Phenoage gap`=roc(outcome ~ pheno, data=prediction),
      `Age+sex+edu+Brain age gap`=roc(outcome ~ demo_brain, data=prediction),
      `Age+sex+edu+Phenoage gap`=roc(outcome ~ demo_pheno, data=prediction)
      
      
      
      # `Cognition` = roc(outcome ~ cog, data=prediction),
      # `Demo+brain age gap` = roc(outcome ~ brain1, data=prediction),
      # `Demo+cog+brain age gap` = roc(outcome ~ brain2, data=prediction),
      # `Cognition+PhenoAge` = roc(outcome ~ pheno, data=prediction),
      # `Cognition+PhenoAge+Brain` = roc(outcome ~ pheno_brain, data=prediction),
      # `Cognition+PhenoAge+AD PRS` = roc(outcome ~ prs, data=prediction),
      # `Cognition+PhenoAge+AD PRS+Brain` = roc(outcome ~ prs_brain, data=prediction)
    )
    
  }else{
    if(outcome %in% c("death")){
      
      prediction <- data.frame(
        outcome = data[[outcome]],
        
        #Demo only
        demo = get_predictions(model_demo, data),
        #Aging only
        brain= get_predictions(model_brain, data),
        organ= get_predictions(model_organ, data),
        pheno = get_predictions(model_pheno, data),

        #Demo+Aging
        demo_brain = get_predictions(model_demo_brain, data),
        demo_organ = get_predictions(model_demo_organ, data),
        demo_pheno = get_predictions(model_demo_pheno, data)
        
        # # demo = get_predictions(model_demo, data),
        # brain= get_predictions(model_brain, data),
        # organ= get_predictions(model_organ, data),
        # all=get_predictions(model_all, data),
        # cog = get_predictions(model_cog, data),
        # # brain1 = get_predictions(model_brain1, data),
        # # brain2 = get_predictions(model_brain2, data),
        # pheno = get_predictions(model_pheno2, data),
        # pheno_organ=get_predictions(model_pheno_organ, data),
        # pheno2 = get_predictions(model_pheno3, data),
        # prs=get_predictions(model_prs, data)
      )
      
      prediction_roc <- list(
        `Age+sex+edu` = roc(outcome ~ demo, data=prediction),
        
        `Brain age gap`=roc(outcome ~ brain, data=prediction),
        `Organismal age gap`=roc(outcome ~ organ, data=prediction),
        `Phenoage gap`=roc(outcome ~ pheno, data=prediction),
        
        `Age+sex+edu+Brain age gap` = roc(outcome ~ demo_brain, data=prediction),
        `Age+sex+edu+Organismal age gap` = roc(outcome ~ demo_organ, data=prediction),
        `Age+sex+edu+Phenoage gap` = roc(outcome ~ demo_pheno, data=prediction)
        
        # `All organ gap`=roc(outcome ~ all, data=prediction),
        # PhenoAge=roc(outcome ~ pheno2, data=prediction),
        # `Brain+cognition` = roc(outcome ~ cog, data=prediction),
        # `Demo+brain age gap` = roc(outcome ~ brain1, data=prediction),
        # `Demo+cog+brain age gap` = roc(outcome ~ brain2, data=prediction),
        # `Organismal+PhenoAge` = roc(outcome ~ pheno_organ, data=prediction),
        # `Brain+PhenoAge` = roc(outcome ~ pheno, data=prediction)
        
        # `Brain age gap+cognition+PhenoAge+PRS` = roc(outcome ~ prs, data=prediction)
      )
    }else{
      prediction_roc <- list(
        # Basic = roc(outcome ~ demo, data=prediction),
        `Brain age gap`=roc(outcome ~ brain, data=prediction),
        # `Cognition` = roc(outcome ~ cog, data=prediction),
        `Cognition+PhenoAge` = roc(outcome ~ pheno, data=prediction),
        `Cognition+PhenoAge+Brain` = roc(outcome ~ pheno_brain, data=prediction)
        
        # `Brain age gap+cognition+PhenoAge+PRS` = roc(outcome ~ prs, data=prediction)
      )
    }
    
  }
  
  create_labels <- function(roc_curves) {
    sapply(names(roc_curves), function(name) {
      auc <- auc(roc_curves[[name]])
      p_value <- roc.test(roc_curves[[name]], roc_curves$`Age+sex+edu`,boot.n=2000)$p.value
      sprintf("%s (%.3f, p=%.3f)", name, auc, p_value)
    })
  }
  
  # create_labels <- function(roc_curves) {
  #   names_list <- names(roc_curves)
  #   n <- length(roc_curves)
  #   
  #   sapply(seq_along(roc_curves), function(i) {
  #     auc_value <- auc(roc_curves[[i]])
  #     
  #     if (i == 1) {
  #       # For the first ROC curve, there is no preceding curve to compare with
  #       sprintf("%s (AUC=%.3f)", names_list[i], auc_value)
  #     } else {
  #       p_value <- roc.test(roc_curves[[i]], roc_curves[[i - 1]], boot.n=2000)$p.value
  #       sprintf("%s (%.3f, p=%.3f)", names_list[i], auc_value, p_value)
  #     }
  #   })
  # }
  # 
  
  create_ggroc_plot <- function(roc_curves, labels, title) {
    ggroc(roc_curves,size=0.2,legacy.axes = TRUE) +
      scale_color_manual(values = 1:length(roc_curves), labels = labels) +
      ggtitle(NULL) +
      theme_classic(base_size = 6,base_line_size = 0.2) +
      theme(
        legend.title = element_blank(),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
        legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
        legend.key.spacing.y = unit(0, "pt"),
        legend.box.spacing = margin(0),
        legend.margin = margin(t = 1, b =-3, unit = "pt"),
        legend.box.margin = margin(t = 1, b = -3, unit = "pt"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # aspect.ratio=1
        # legend.margin=margin(0,0,0,0),
        # legend.box.margin=margin(-10,-10,-10,-10)
      ) +
      guides(color = guide_legend(nrow = 7, byrow = TRUE))+labs(title=NULL,subtitle = NULL,caption = NULL)+
      annotate("segment",x = 0, xend = 1, y = 0, yend = 1, color="grey", linetype="dashed",size=.2)
    
  }
  
  labels<- create_labels(prediction_roc)
  
  plot<- create_ggroc_plot(prediction_roc, labels, NULL)
  
  # 
  # get_legend <- function(myggplot) {
  #   tmp <- ggplotGrob(myggplot)
  #   leg <- gtable::gtable_filter(tmp, "guide-box")
  #   return(leg)
  # }
  # 
  # legend<- get_legend(plot)
  # 
  # plot <- plot + theme(legend.position = "none")
  # 
  # # Combine plots with legends
  # plot <- plot_grid(legend, plot, ncol = 1,nrow=2, rel_heights = c(0.2, 1))
  # 
  
  return(plot)
  
}

plot_brain_acd <- roc_brain2(ukb_olink_inci,outcome = "acd_inci")
plot_brain_ad <- roc_brain2(ukb_olink_inci,outcome = "ad_inci")
plot_brain_vd <- roc_brain2(ukb_olink_inci,outcome = "vd_inci")
plot_brain_death <- roc_brain2(ukb_olink_inci, outcome = "death")

#Supp fig
plot_grid(plot_brain_acd, plot_brain_ad, plot_brain_vd,plot_brain_death,
          ncol = 2,labels = "auto",label_size = 6)




#ROC curve with AUC main figure------------------------


roc_brain2 <- function(data,outcome){
  data[[outcome]] <- factor(data[[outcome]])
  model_demo <- glm(paste(outcome, "~ age+sex+edu"), data = data, family = binomial)
  model_brain <- glm(paste(outcome, "~ age+sex+edu+Brain_res"), data = data, family = binomial)
  model_cog <- glm(paste(outcome, "~ age+sex+edu+rt"), data = data, family = binomial)
  model_pheno <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res"), data = data, family = binomial)
  model_pheno_brain <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+Brain_res"), data = data, family = binomial)
  
  model_pheno_organ <- glm(paste(outcome, "~ age+sex+edu+rt+Organismal_res+phenoage_res"), data = data, family = binomial)
  model_pheno2 <- glm(paste(outcome, "~ age+sex+edu+Brain_res+phenoage_res"), data = data, family = binomial)
  
  model_prs <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+ad_prs"), data = data, family = binomial)
  model_prs_brain <- glm(paste(outcome, "~ age+sex+edu+rt+phenoage_res+ad_prs+Brain_res"), data = data, family = binomial)
  
  
  model_organ <- glm(paste(outcome, "~ age+sex+edu+Organismal_res"), data = data, family = binomial)
  model_all <- glm(paste(outcome, "~ age+sex+edu+",paste(organ_pac,collapse = "+")), data = data, family = binomial)
  model_pheno3 <- glm(paste(outcome, "~ age+sex+edu+phenoage_res"), data = data, family = binomial)
  
  
  get_predictions <- function(model, data) {
    predict(model, newdata = data, type = "response")
  }
  
  prediction <- data.frame(
    outcome = data[[outcome]],
    demo = get_predictions(model_demo, data),
    brain= get_predictions(model_brain, data),
    cog = get_predictions(model_cog, data),
    # brain1 = get_predictions(model_brain1, data),
    # brain2 = get_predictions(model_brain2, data),
    pheno = get_predictions(model_pheno, data),
    pheno_brain = get_predictions(model_pheno_brain, data),
    prs=get_predictions(model_prs, data),
    prs_brain=get_predictions(model_prs_brain, data)
  )
  
  if(outcome %in% c("ad_inci")){
    
    prediction_roc <- list(
      # Demo = roc(outcome ~ demo, data=prediction),
      `Brain age gap`=roc(outcome ~ brain, data=prediction),
      # `Cognition` = roc(outcome ~ cog, data=prediction),
      # `Demo+brain age gap` = roc(outcome ~ brain1, data=prediction),
      # `Demo+cog+brain age gap` = roc(outcome ~ brain2, data=prediction),
      `Cognition+PhenoAge` = roc(outcome ~ pheno, data=prediction),
      `Cognition+PhenoAge+Brain` = roc(outcome ~ pheno_brain, data=prediction),
      `Cognition+PhenoAge+AD PRS` = roc(outcome ~ prs, data=prediction),
      `Cognition+PhenoAge+AD PRS+Brain` = roc(outcome ~ prs_brain, data=prediction)
    )
    
  }else{
    if(outcome %in% c("death")){
      
      prediction <- data.frame(
        outcome = data[[outcome]],
        # demo = get_predictions(model_demo, data),
        brain= get_predictions(model_brain, data),
        organ= get_predictions(model_organ, data),
        all=get_predictions(model_all, data),
        cog = get_predictions(model_cog, data),
        # brain1 = get_predictions(model_brain1, data),
        # brain2 = get_predictions(model_brain2, data),
        pheno = get_predictions(model_pheno2, data),
        pheno_organ=get_predictions(model_pheno_organ, data),
        pheno2 = get_predictions(model_pheno3, data),
        prs=get_predictions(model_prs, data)
      )
      
      prediction_roc <- list(
        # Demo = roc(outcome ~ demo, data=prediction),
        `Brain age gap`=roc(outcome ~ brain, data=prediction),
        `Organismal age gap`=roc(outcome ~ organ, data=prediction),
        `All organ gap`=roc(outcome ~ all, data=prediction),
        # PhenoAge=roc(outcome ~ pheno2, data=prediction),
        # `Brain+cognition` = roc(outcome ~ cog, data=prediction),
        # `Demo+brain age gap` = roc(outcome ~ brain1, data=prediction),
        # `Demo+cog+brain age gap` = roc(outcome ~ brain2, data=prediction),
        `Organismal+PhenoAge` = roc(outcome ~ pheno_organ, data=prediction),
        `Brain+PhenoAge` = roc(outcome ~ pheno, data=prediction)
        
        # `Brain age gap+cognition+PhenoAge+PRS` = roc(outcome ~ prs, data=prediction)
      )
    }else{
      prediction_roc <- list(
        # Basic = roc(outcome ~ demo, data=prediction),
        `Brain age gap`=roc(outcome ~ brain, data=prediction),
        # `Cognition` = roc(outcome ~ cog, data=prediction),
        `Cognition+PhenoAge` = roc(outcome ~ pheno, data=prediction),
        `Cognition+PhenoAge+Brain` = roc(outcome ~ pheno_brain, data=prediction)
        
        # `Brain age gap+cognition+PhenoAge+PRS` = roc(outcome ~ prs, data=prediction)
      )
    }
    
  }
  
  create_labels <- function(roc_curves) {
    sapply(names(roc_curves), function(name) {
      auc <- auc(roc_curves[[name]])
      p_value <- roc.test(roc_curves[[name]], roc_curves$`Brain age gap`)$p.value
      sprintf("%s (%.3f, p=%.3f)", name, auc, p_value)
    })
  }
  
  # create_labels <- function(roc_curves) {
  #   names_list <- names(roc_curves)
  #   n <- length(roc_curves)
  #   
  #   sapply(seq_along(roc_curves), function(i) {
  #     auc_value <- auc(roc_curves[[i]])
  #     
  #     if (i == 1) {
  #       # For the first ROC curve, there is no preceding curve to compare with
  #       sprintf("%s (AUC=%.3f)", names_list[i], auc_value)
  #     } else {
  #       p_value <- roc.test(roc_curves[[i]], roc_curves[[i - 1]], boot.n=2000)$p.value
  #       sprintf("%s (%.3f, p=%.3f)", names_list[i], auc_value, p_value)
  #     }
  #   })
  # }
  # 
  
  create_ggroc_plot <- function(roc_curves, labels, title) {
    ggroc(roc_curves,size=0.2,legacy.axes = TRUE) +
      scale_color_manual(values = 1:length(roc_curves), labels = labels) +
      ggtitle(NULL) +
      theme_classic(base_size = 6,base_line_size = 0.2) +
      theme(
        legend.title = element_blank(),
        legend.position = "top",
        legend.box = "horizontal",
        legend.key.size = unit(0.1, "cm"),  # Adjust key size to reduce spacing
        legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
        legend.key.spacing.y = unit(0, "pt"),
        legend.box.spacing = margin(0),
        legend.margin = margin(t = 1, b =-3, unit = "pt"),
        legend.box.margin = margin(t = 1, b = -3, unit = "pt"),
        plot.margin = unit(c(0, 0, 0, 0), "cm"),
        # aspect.ratio=1
        # legend.margin=margin(0,0,0,0),
        # legend.box.margin=margin(-10,-10,-10,-10)
      ) +
      guides(color = guide_legend(nrow = 7, byrow = TRUE))+labs(title=NULL,subtitle = NULL,caption = NULL)+
      annotate("segment",x = 0, xend = 1, y = 0, yend = 1, color="grey", linetype="dashed",size=.2)
    
  }
  
  labels<- create_labels(prediction_roc)
  
  plot<- create_ggroc_plot(prediction_roc, labels, NULL)
  
  # 
  # get_legend <- function(myggplot) {
  #   tmp <- ggplotGrob(myggplot)
  #   leg <- gtable::gtable_filter(tmp, "guide-box")
  #   return(leg)
  # }
  # 
  # legend<- get_legend(plot)
  # 
  # plot <- plot + theme(legend.position = "none")
  # 
  # # Combine plots with legends
  # plot <- plot_grid(legend, plot, ncol = 1,nrow=2, rel_heights = c(0.2, 1))
  # 
  
  return(plot)
  
}

plot_brain_acd <- roc_brain2(ukb_olink_inci,outcome = "acd_inci")
plot_brain_ad <- roc_brain2(ukb_olink_inci,outcome = "ad_inci")
plot_brain_vd <- roc_brain2(ukb_olink_inci,outcome = "vd_inci")

plot_brain_acd10 <- roc_brain2(ukb_olink_inci %>% filter(!(acd_inci==1&acd_fu<10)),outcome = "acd_inci")
plot_brain_ad10 <- roc_brain2(ukb_olink_inci %>% filter(!(ad_inci==1&ad_fu<10)),outcome = "ad_inci")
plot_brain_vd10 <- roc_brain2(ukb_olink_inci %>% filter(!(vd_inci==1&vd_fu<10)),outcome = "vd_inci")

plot_brain_death <- roc_brain2(ukb_olink_inci, outcome = "death")

plot_grid(plot_brain_acd, plot_brain_ad, plot_brain_vd, 
          plot_brain_acd10,plot_brain_ad10,plot_brain_vd10,
          plot_brain_death,c_index_brain,plot_grid(acd_hr,death_hr,nrow=2),
          ncol = 3,labels = "auto",label_size = 6)+theme(plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"))




