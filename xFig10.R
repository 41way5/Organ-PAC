# xFig10
  #Refined aging clock---------------------
  
  #Performance of refined vs full brain/organismal clock 
    #Brain (full/refined)
      #UKB 0.777/0.737
      #CKB 0.760/0.710
      #NHS 0.614/0.563
  #organismal (full/refined)
    #UKB 0.939/0.888
    #CKB 0.919/0.859
    #NHS 0.850/0.741
  
  refined_df <- data.frame(
    Organ = rep(c("Brain", "Organismal"), each = 6),
    ClockType = rep(c("Full", "Refined"), times = 3, each = 1),
    Study = rep(c("UKB", "CKB", "NHS"), each = 2),
    R_Performance = c(0.777, 0.737, 0.760, 0.710, 0.614, 0.563,  # Brain
                      0.939, 0.888, 0.919, 0.859, 0.850, 0.741)  # Organismal
  )


  protein_df <- data.frame(
    Organ = c("Brain", "Organismal"),
    Full = c(70, 230),
    Refined = c(10, 20)
  ) %>%
    mutate(Ratio_Brain = 10 / 70,   # Brain: 10/70
           Ratio_Organismal = 20 / 230)  # Organismal: 20/230
  

  protein_expanded <- expand.grid(Organ = unique(refined_df$Organ), Study = "UKB") %>%
    left_join(protein_df, by = "Organ") %>%
    mutate(ClockType = "Protein",
           R_Performance_Brain = 1-Ratio_Brain,
           R_Performance_Organismal = 1-Ratio_Organismal) 
  
  refined_df$Study <- factor(refined_df$Study, levels = c("UKB", "CKB", "NHS"))
  protein_expanded$Study <- factor(protein_expanded$Study, levels = c("UKB", "CKB", "NHS"))
  

x10_b <-   ggplot() +
    geom_bar(data = refined_df, aes(x = Study, y = R_Performance, fill = ClockType),
             stat = "identity", position = "dodge", width = 0.4) +
    
    geom_bar(data = protein_expanded, aes(x = as.numeric(Study) - 0.4, 
                                          y = ifelse(Organ == "Brain", R_Performance_Brain, R_Performance_Organismal)),
             stat = "identity", fill = "#fdae61", width = 0.2, alpha = 0.6) +
    
    facet_wrap(~ Organ) +
    labs(y = "Performance (r) / Feature reduction (%)", x = NULL) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.2)) +
    coord_cartesian(ylim = c(0, 1)) +
    scale_fill_manual(values = c("Full" = "#a6dba0", "Refined" = "#c2a5cf")) +
    
    theme_bw(base_size = 7) +
    theme(legend.title = element_blank(),
          panel.grid = element_blank(),
          
          legend.position = c(0.88, 0.88),
          legend.key.size = unit(0.2, "cm"))
  
  
  
  
  
  
  #Short-panel based PAC for Brain and Organismal
  
  #UKB----
  result_list_rfe_brain <- 
    cox(
      data=ukb_olink,
      exposures=c("Brain_res","Brain_res_rfe"),
      outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death","copd","mi","stroke","ckd","ld","dm"), 
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok"
      ))
  
  
  result_disease_rfe_brain <- result_disease_fomrat(result_list_rfe_brain) %>% 
    mutate(Exposure=ifelse(Exposure=="Brain_res","Full","Refined"))
  
  result_list_rfe_organismal <- 
    cox(
      data=ukb_olink,
      exposures=c("Organismal_res","Organismal_res_rfe"),
      outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death","copd","mi","stroke","ckd","ld","dm"), 
      covar=c("age","sex","ethnic","edu","tdi","center","pa","smok"
      ))
  
  result_disease_rfe_organismal <- result_disease_fomrat(result_list_rfe_organismal) %>% 
    mutate(Exposure=ifelse(Exposure=="Organismal_res","Full","Refined"))
  
  
  refined_cox <- rbind(      result_disease_rfe_brain %>% mutate(organ="Brain"),
                            result_disease_rfe_organismal %>% mutate(organ="Organismal")
                            )
  
  
  
  x10_c <-   ggplot(refined_cox,
         # group_by(rev(Outcome)) %>%
         # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("Depression","Anxiety","Schizophrenia","Sleep disorder","Substance use disorder","All-cause mortality")))) %>% 
         # arrange(desc(Outcome_label),desc(Beta)) %>% 
         # mutate(Exposure_label = factor(Exposure_label, levels = rev(unique(Exposure_label)))),
         aes(x = HR, y=reorder(Outcome_label,HR),groups=Exposure)) +
    geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
    geom_point(aes(color = Exposure),size = 1,shape=16,position=position_dodge(0.7))+
    #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
    scale_y_discrete()+
    labs(x = "HR (95% CI)", y = NULL, title = NULL) +
    facet_wrap(~ organ, scales = "free_x")+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    theme(            plot.title = element_text(hjust = 0.5,face = "bold",size = 6),
                      legend.position = c(0.85,0.1),legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                      axis.text.x = element_text(color="black",size = 6),
                      axis.text.y = element_text(color="black",size = 6),
                      legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                      legend.key.spacing.y = unit(0, "pt"),
                      legend.box.spacing = margin(0),
                      legend.margin = margin(t = -2, b =-3, unit = "pt"),
                      legend.box.margin = margin(t = -2, b = -3, unit = "pt"),
                      plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm"),
                      
    )+scale_color_manual(values = c("Full" = "#a6dba0", "Refined" = "#c2a5cf")) +

    # scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
    guides(color = guide_legend(nrow = 1))
  
  #CKB
  ckb_cox_full <- read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/CKB/Cox_results 1.csv")
  ckb_cox_rfe  <- read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/CKB/Cox_results_pAD10.csv")

  
  ckb_cox <- rbind(ckb_cox_full %>% mutate(Exposure="Full") %>% filter(Score %in% c("Organismal_pAD","Brain_pAD"))%>% mutate(organ=case_when(Score=="Organismal_pAD"~"Organismal",
                                                                                                                                             Score=="Brain_pAD"~"Brain")),
                   ckb_cox_rfe %>% mutate(Exposure="Refined") %>% mutate(organ=case_when(Score=="Organismal_pAD_10"~"Organismal",
                                                                                         Score=="Brain_pAD_10"~"Brain"
                                                                                         ))
                   ) %>% mutate(HR=exp(coef),CI_Lower=exp(coef-1.96*se),CI_Upper=exp(coef+1.96*se)) %>% filter(N>50)
  
  ckb_cox <- ckb_cox %>% filter()
  
  
  x10_d <- ggplot(ckb_cox,
         # group_by(rev(Outcome)) %>%
         # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("Depression","Anxiety","Schizophrenia","Sleep disorder","Substance use disorder","All-cause mortality")))) %>% 
         # arrange(desc(Outcome_label),desc(Beta)) %>% 
         # mutate(Exposure_label = factor(Exposure_label, levels = rev(unique(Exposure_label)))),
         aes(x = HR, y=reorder(Outcome,HR),groups=Exposure)) +
    geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
    geom_point(aes(color = Exposure),size = 1,shape=16,position=position_dodge(0.7))+
    #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
    scale_y_discrete()+
    labs(x = "HR (95% CI)", y = NULL, title = NULL) +
    facet_wrap(~ organ, scales = "free_x")+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    theme(            plot.title = element_text(hjust = 0.5,face = "bold",size = 6),
                      legend.position = c(0.85,0.1),legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                      axis.text.x = element_text(color="black",size = 6),
                      axis.text.y = element_text(color="black",size = 6),
                      legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                      legend.key.spacing.y = unit(0, "pt"),
                      legend.box.spacing = margin(0),
                      legend.margin = margin(t = -2, b =-3, unit = "pt"),
                      legend.box.margin = margin(t = -2, b = -3, unit = "pt"),
                      plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"),
                      
    )+scale_color_manual(values = c("Full" = "#a6dba0", "Refined" = "#c2a5cf")) +
    
    # scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
    guides(color = guide_legend(nrow = 1))
  
  #NHS
  
  nhs_cox_full <- read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/NHS/nhs1/disease_nhs.csv") %>% 
    filter(group=="full") %>% filter(Organ %in% c("Brain","Organismal")) %>% mutate(Exposure="Full") 
  nhs_cox_rfe  <- read.csv("E:/Tutor/UK Biobank/Project/Proteomic aging/Submission/Submission/Nat Aging/Revision/2-16 R1/Validation/NHS/nhs1/disease_nhs_rfe.csv")%>% filter(group=="full") %>% mutate(Exposure="Refined") 
  
  
  nhs_cox <- rbind(nhs_cox_full,
                   nhs_cox_rfe
                   )%>% mutate(HR=exp(coef),CI_Lower=exp(coef-1.96*`se.coef.`),CI_Upper=exp(coef+1.96*`se.coef.`)) %>% 
    mutate(Outcome=case_when(Outcome=="db2"~"Type 2 diabetes",
                                     Outcome=="death"~"All-cause mortality",
                                     Outcome=="stroke"~"Stroke",
                                     Outcome=="mi"~"Myocardial infarction",
                                     Outcome=="adcase"~"Alzheimer's disease",
    ))
  x10_e <-  ggplot(nhs_cox,
         # group_by(rev(Outcome)) %>%
         # mutate(Outcome_label=factor(Outcome_label,levels=rev(c("Depression","Anxiety","Schizophrenia","Sleep disorder","Substance use disorder","All-cause mortality")))) %>% 
         # arrange(desc(Outcome_label),desc(Beta)) %>% 
         # mutate(Exposure_label = factor(Exposure_label, levels = rev(unique(Exposure_label)))),
         aes(x = HR, y=reorder(Outcome,HR),groups=Exposure)) +
    geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
    geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
    geom_point(aes(color = Exposure),size = 1,shape=16,position=position_dodge(0.7))+
    #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
    scale_y_discrete()+
    labs(x = "HR (95% CI)", y = NULL, title = NULL) +
    facet_wrap(~ Organ, scales = "free_x")+
    theme_classic(base_size = 6,base_line_size = 0.2)+
    theme(            plot.title = element_text(hjust = 0.5,face = "bold",size = 6),
                      legend.position = c(0.85,0.1),legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                      axis.text.x = element_text(color="black",size = 6),
                      axis.text.y = element_text(color="black",size = 6),
                      legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                      legend.key.spacing.y = unit(0, "pt"),
                      legend.box.spacing = margin(0),
                      legend.margin = margin(t = -2, b =-3, unit = "pt"),
                      legend.box.margin = margin(t = -2, b = -3, unit = "pt"),
                      plot.margin = unit(c(0.1, 0.2, 0.1, 0.1), "cm"),
                      
    )+scale_color_manual(values = c("Full" = "#a6dba0", "Refined" = "#c2a5cf")) +
    
    # scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
    guides(color = guide_legend(nrow = 1))
  
  
  
plot_grid(x10_b,x10_c,x10_d,x10_e,nrow=2,ncol=2,align = "v",labels = c("b","c","d","e"),label_size = 8)
  
  
  
  
  
  
  
  
  
  
 
  