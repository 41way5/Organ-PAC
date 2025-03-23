##Fig 5. APOE4 and Brain age gap-------------
ukb_olink_inci <- result_list$data

coxph(formula =Surv(acd_fu, acd_inci) ~ relevel(apoe4,ref="3/3") + Brain_res +age + sex + edu , data=ukb_olink_inci) %>% summary()

coxph(formula =Surv(acd_fu, acd_inci) ~ Brain_res*relevel(apoe4,ref="3/3") +age + sex + edu , data=ukb_olink_inci) %>% summary()
coxph(formula =Surv(acd_fu, acd_inci) ~ Organismal_res*relevel(apoe4,ref="3/3") +age + sex + edu , data=ukb_olink_inci) %>% summary()
coxph(formula =Surv(acd_fu, acd_inci) ~ Artery_res*relevel(apoe4,ref="3/3") +age + sex + edu , data=ukb_olink_inci) %>% summary()


coxph(formula =Surv(acd_fu, acd_inci) ~  Brain_res +age + sex + edu , data=ukb_olink_inci %>% filter(apoe4=="4/4")) %>% summary()
coxph(formula =Surv(acd_fu, acd_inci) ~  Brain_res_super+Brain_res_extra+age + sex + edu , data=ukb_olink_inci %>% filter(apoe4=="3/4")) %>% summary()
ukb_olink_inci$Brain_res_super


results_apoe <- lapply(split(ukb_olink_inci, ukb_olink_inci$apoe4),
                       function(x) coxph(Surv(acd_fu, acd_inci) ~ Brain_res+ age + sex + edu, data = x) %>% tidy(exponentiate = TRUE,conf.int = TRUE))

add_label <- function(df, label) {
  df$label <- label
  return(df)
}
results_apoe_labeled <- lapply(names(results_apoe), function(name) add_label(results_apoe[[name]], name))
results_apoe_ad <- do.call(rbind, results_apoe_labeled) %>% filter(term=="Brain_res") %>% 
  dplyr::rename(HR=estimate,CI_Lower=conf.low,CI_Upper=conf.high) %>% 
  mutate(label=ifelse(label=="4/3","3/4",label)) %>% 
  mutate(label=factor(label,levels = c("2/x", "3/3", "3/4", "4/4")))


f4b.2 <- ggplot(results_apoe_ad,
                aes(x = HR, y=(label))) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color = label),size = 1,shape=15,position=position_dodge(0.7))+
  #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = expression( italic("APOE")*" haplotype"),) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            plot.title = element_text(hjust = 0.5,face = "bold",size = 6),
                    legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6),
                    legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                    legend.key.spacing.y = unit(0, "pt"),
                    legend.box.spacing = margin(0),
                    legend.margin = margin(t = -2, b =-3, unit = "pt"),
                    legend.box.margin = margin(t = -2, b = -3, unit = "pt"),
                    plot.margin = unit(c(0, 0, 0, 0), "cm"),
                    
  )+
  # scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))+
  coord_flip()


f4b.1 <- coxph(formula =Surv(acd_fu, acd_inci) ~ Brain_res*relevel(apoe4,ref="3/3") +age + sex + edu , data=ukb_olink_inci) %>% tidy(exponentiate = TRUE,conf.int = TRUE) %>% 
  filter(grepl(c("Brain_res|apoe4"),term)) %>% 
  mutate(
    term = case_when(
      str_detect(term, "^relevel") ~ str_sub(term, -3),
      str_detect(term, ":") ~ paste0(str_sub(term, -3), "*"),
      TRUE ~ term
    )
  ) %>% 
  mutate(term=ifelse(term=="Brain_res","Brain",term)) %>% 
  dplyr::rename(HR=estimate,CI_Lower=conf.low,CI_Upper=conf.high) %>% 
  bind_rows(.,
            data.frame(term = "3/3",HR = 1,CI_Lower=1,CI_Upper=1)
  ) %>% 
  mutate(term=factor(term, levels=c("2/x","3/3","3/4","4/4","Brain","2/x*","3/4*","4/4*")),
         # label=ifelse(grep("*",term),"interaction","")
  )%>% 
  
  ggplot(.,
         aes(x = HR, y=(term),color=term)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(size = 1,shape=16,position=position_dodge(0.7))+
  scale_x_continuous( breaks = c(1,2,4,6)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = expression(italic("APOE")*" haplotype and Brain age gap"),) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            plot.title = element_text(hjust = 0.5,face = "bold",size = 6),
                    legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6),
                    legend.spacing = unit(0, "cm"),    # Adjust spacing between legend items
                    legend.key.spacing.y = unit(0, "pt"),
                    legend.box.spacing = margin(0),
                    legend.margin = margin(t = -2, b =-3, unit = "pt"),
                    legend.box.margin = margin(t = -2, b = -3, unit = "pt"),
                    plot.margin = unit(c(0, 0, 0, 0), "cm"),
                    
  )+
  # scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))+
  coord_flip()+
  scale_color_manual(values=c("#4daf4a","#377eb8","#ff7f00","#e41a1c","#4B0082","black","black","black"))




ukb_olink_inci <- result_list$data %>% mutate(ad_year=ifelse(ad_inci==1,age+ad_fu,NA),
                                              acd_year=ifelse(acd_inci==1,age+acd_fu,NA),
) %>% 
  left_join(olink_full %>% select(c(top_proteins_list$Brain,
                                    top_proteins_list$Organismal,
                                    top_proteins_list$Artery,
                                    top_proteins_list$Heart),"eid"),by="eid") %>% 
  mutate(across(c(top_proteins_list$Brain,top_proteins_list$Organismal,top_proteins_list$Artery,top_proteins_list$Heart),scale)) %>% 
  mutate(across(c(top_proteins_list$Brain,top_proteins_list$Organismal,top_proteins_list$Artery,top_proteins_list$Heart),as.numeric)) %>% 
  mutate(apoe4=ifelse(apoe4=="4/3","3/4",apoe4)) %>% 
  mutate(apoe4=factor(apoe4, levels = c("2/x", "3/3", "3/4", "4/4")))

describe(ukb_olink_inci$apoe4)

psych::describeBy(ukb_olink_inci$acd_year,group =ukb_olink_inci$apoe4)



create_plot_apoe <- function(data=ukb_olink_inci, y_var) {
  ggplot(data %>% 
           filter(!is.na(apoe4)) %>%
           filter(is.na(acd_base)) %>% 
           filter(apoe4 %in% c("2/x", "3/3", "4/4")), 
         aes(x = age, y = .data[[y_var]], group = apoe4, color = apoe4)) +
    geom_smooth(method = "loess", se = TRUE, size = 0.5, fill = "#bdbdbd") +
    labs(title = NULL,
         x = "Age",
         y = case_when(y_var=="Brain_res"~"Brain age gap",y_var=="Artery_res"~"Artery age gap", T~y_var),
         color = NULL) +
    theme_bw(base_size = 6, base_line_size = 0.2) +
    scale_color_manual(values = c("#4daf4a", "#377eb8", "#e41a1c"),
                       labels = c("2/X", "3/3", "4/4")) +
    guides(color = guide_legend(override.aes = list(fill = NA))) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"),
          legend.key.width = unit(0.3, "cm"),
          legend.key = element_rect(color = NA, fill = NA),
          legend.key.size = unit(0.2, "cm"),
          legend.position = if (y_var == "Brain_res") c(0.4, 0.96) else "none",
          legend.justification = c("right", "top"),
          legend.box.just = "right")
}


f4b.4 <- plot_grid(plotlist = lapply(c("Brain_res","NEFL","GFAP","SYT1"), function(var) create_plot_apoe(ukb_olink_inci, var)), 
                   nrow = 1,labels = c("d"),label_size = 6)
f4b.5 <- plot_grid(plotlist = lapply(c("Artery_res","ELN","LTBP2","TNFRSF11B"), function(var) create_plot_apoe(ukb_olink_inci, var)), 
                   nrow = 1,labels = c("e"),label_size = 6)




ggplot(ukb_olink_inci %>% filter(!is.na(apoe4)) %>% filter(is.na(acd_base)) %>% 
         filter(apoe4=="2/x"|apoe4=="3/3"|apoe4=="4/4"), 
       aes(x = age, y = Brain_res, group=apoe4, color = apoe4)) +
  # geom_point() + 
  geom_smooth(method = "loess",se = T, size = 0.5,fill="#bdbdbd") +
  labs(title = NULL,
       x = "Age",
       y = "Biomarker",
       color = NULL) +
  theme_bw(base_size = 6,base_line_size = 0.2)+
  scale_color_manual(values = c("#4daf4a", "#377eb8","#e41a1c"),labels=c("2/X","3/3","4/4"))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.width = unit(0.5, "cm"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.2, "cm"),
        legend.position  = c(0.2, 0.96), legend.justification = c("right", "top"), legend.box.just = "right",
  )

ggplot(ukb_olink_inci %>% filter(!is.na(apoe4)) %>% filter(is.na(acd_base)) %>% 
         filter(apoe4=="2/x"|apoe4=="3/3"|apoe4=="4/4"), aes(x = age, y = NEFL, group=apoe4, color = apoe4)) +
  # geom_point() + 
  geom_smooth(method = "loess",se = T, size = 0.5,fill="#bdbdbd") +
  labs(title = NULL,
       x = "Age",
       y = "Biomarker",
       color = NULL) +
  theme_bw(base_size = 6,base_line_size = 0.2)+
  scale_color_manual(values = c("#4daf4a", "#377eb8","#e41a1c"),labels=c("2/X","3/3","4/4"))+
  guides(color=guide_legend(override.aes=list(fill=NA)))+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"),
        legend.key.width = unit(0.5, "cm"),
        legend.key = element_rect(color = NA, fill = NA),
        legend.key.size = unit(0.2, "cm"),
        legend.position  = c(0.2, 0.96), legend.justification = c("right", "top"), legend.box.just = "right",
  )



##KM curve for stratification
describe(ukb_olink_inci_km$Brain_ageo %>% as.factor())
describe(ukb_olink_inci_km$brain_apoe %>% as.factor())
describe(ukb_olink_inci_km$apoe4_2 %>% as.factor())
describe(ukb_olink_inci_km$apoe4 %>% as.factor())
unique(ukb_olink_inci_km$brain_apoe)
ukb_olink_inci_km <- ukb_olink_inci %>% filter(!is.na(apoe4)&apoe4!="2/x") %>% 
  mutate(Brain_ageo=ifelse(Brain_res_extra==1,"↑",ifelse(Brain_res_super==1,"↓","-")),
         apoe4_2=case_when(apoe4=="3/4"|apoe4=="4/4"~"4/X",T~apoe4),
         brain_apoe=paste0(apoe4_2," ",Brain_ageo),
         brain_apoe=factor(brain_apoe,levels=c("3/3 ↓","4/X ↓","3/3 -","4/X -","3/3 ↑","4/X ↑"))
  )

coxph(formula =Surv(ad_fu, acd_inci) ~ relevel(factor(brain_apoe),ref="3/3 -")  +age + sex + edu , data=ukb_olink_inci_km) %>% summary()


brain_apoe_km <- survfit(Surv(acd_fu, acd_inci) ~ brain_apoe, ukb_olink_inci_km) %>% 
  survminer::ggsurvplot(fun = "event",
                        conf.int = T,
                        censor = FALSE,
                        size = 0.2,
                        ggtheme = theme_classic(base_size=6, base_line_size = 0.2, base_family = "Arial"),
                        legend.title=expression(italic("APOE")*"+Brain ageotypes"),
                        # palette = c("#abdda4","pink","#0571b0","","#f4a582","#ca0020"),
                        break.x.by = 5,
                        xlab = "Time (years)",
                        legend.labs=c("3/3 ↓","4/X ↓","3/3 -","4/X -","3/3 ↑","4/X ↑")
  )

brain_apoe_km$plot <-  brain_apoe_km$plot+
  ggplot2::annotate("text",x = 0.5, y = 0.09,hjust = 0,
                    label = "3/3 ↓: 0.4 (0.2-0.8)\n3/3 -: REF\n3/3 ↑: 3.0 (2.0-4.5)\n4/X ↓: 0.6 (0.3-1.5)\n4/X -: 3.6 (3.1-4.2)\n4/X ↑: 11.0 (7.9-15.4)", size = 2)+
  theme_classic(base_line_size = 0.2,base_size = 6)+
  theme(legend.text = element_text(size = 5, color = "black"),
        legend.title =element_text(size = 5),
        legend.key.width = unit(0.15, "cm"),
        legend.key.size = unit(0.1, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        legend.position = "top"
  )+
  guides(colour = guide_legend(nrow = 2))
f4b.3 <- brain_apoe_km$plot

#Fig 4b: APOE and BRAIN AGING
plot_grid(plot_grid(NULL,plot_grid(f4b.1,f4b.2,nrow=1,labels=c("a","b"),label_size = 6,rel_widths = c(1.5,1)),NULL,nrow=3,rel_heights = c(0.2,1,0.2)),f4b.3,nrow=1,rel_widths = c(1.5,1),labels=c("","c"),label_size = 6) %>% 
  plot_grid(.,f4b.4,f4b.5,nrow=3,rel_heights=c(1,0.5,0.5))+theme(plot.margin = unit(c(1,1,1,1), "mm"))
