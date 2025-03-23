#-----*****Figure 3*****---------
#Fig 3 PAC in cognitive decline and neurodegenerative disease--------------------------------------------------------------

disease_prefixes <- c("acd", "ad", "vd", "pd", "ms", "scz", "dep", "anx", "sd")

#Add mental health scale measures
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

# 10057 with VM cog test at base and fu
n_distinct(cognition %>% filter(!is.na(`2050-0.0`)&!is.na(`20518-0.0`))) 
describe(cognition %>% filter(!is.na(vm)&!is.na(vm_fu)) %>% select(vm_change))
n_distinct(cognition %>% filter(!is.na(vm)&!is.na(vm_fu))) 

#Add prefix to excluding those with baseline disease
for (prefix in disease_prefixes) {
  interval_col <- paste0(prefix, "_interval")
  cognition <- cognition %>% 
    mutate(!!interval_col := if_else(!is.na(!!sym(paste0(prefix, "_date"))),
                                     as.numeric(difftime(.$recruit_date, !!sym(paste0(prefix, "_date")), units = "days")),
                                     0))
}

cognition_nodisease <-cognition %>% filter(acd_interval<=0&pd_interval<=0&ms_interval<=0)

#Function (linear/logistic/cox) for the association test

cog_result <- linear(cognition_nodisease %>% filter(),a=c("rt","nm","fi","vm"),b=organ_pac,covar = c("age","sex","ethnic","edu","smok","pa","tdi","center")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="rt"~"Reaction time",
                         var_a=="nm"~"Numeric memory",
                         var_a=="fi"~"Fluid intelligence",
                         var_a=="vm"~"Visual memory",
                         T~var_a))


#2.1 Selected PAC [Brain, Organismal, Artery,Heart] and baseline Cog function
#in those wo Neruro disease [N=43141]
f2_1 <- ggplot(cog_result %>% filter(var_b %in% c("Organismal","Brain","Artery","Heart")) %>% mutate(var_b=factor(var_b,levels=c("Organismal","Brain","Artery","Heart"))), 
               aes(x = beta, y=reorder(var_a, beta), group=var_b)) +
  geom_vline(xintercept = 0, color = "grey",linewidth=0.1)+
  geom_errorbarh(aes(xmin = beta-1.95*se, xmax = beta+1.95*se),height = 0,linewidth=0.2,position=position_dodge(0.7)) +
  geom_point(aes(color = var_b), size = 1,shape=15,position=position_dodge(0.7))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #scale_y_discrete(expand = c(0,1))+
  labs(x = "Beta", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
        #aspect.ratio=1
  )+
  scale_color_manual(values=c("black","#4B0082","#A52A2A","#FF7F50"))

##2.2 PAC and FU cognition decline [39684]

n_distinct(cognition_nodisease %>% filter(rt<1.5))

cog_result_fu <- logistic(cognition_nodisease %>% filter(),a=c("vm_change"),b=organ_pac,covar = c("age","sex","ethnic","edu","smok","pa","tdi","center")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="rt"~"Reaction time",
                         var_a=="nm"~"Numeric memory",
                         var_a=="fi"~"Fluid intelligence",
                         var_a=="vm"~"Visual memory",
                         T~var_a))

f2_2 <-ggplot(cog_result_fu %>% filter(var_b %in% c("Organismal","Brain","Artery","Heart")) %>% mutate(var_b=factor(var_b,levels=c("Organismal","Brain","Artery","Heart"))), 
              aes(x = exp(beta), y=reorder(var_b, beta))) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1)+
  geom_errorbarh(aes(xmin = exp(beta-1.95*se), xmax = exp(beta+1.95*se)),linewidth=0.2,height = 0,position=position_dodge(0.7)) +
  geom_point(aes(color = var_b), size = 1,shape=15,position=position_dodge(0.7))+
  scale_x_continuous(limits=c(0.9,1.2), breaks = c(0.9,1,1.10,1.2)) +
  #scale_y_discrete(expand = c(0,1))+
  labs(x = "OR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
        #aspect.ratio=1
  )+
  scale_color_manual(values=c("black","#4B0082","#A52A2A","#FF7F50"))

##2.3 PAC and incident disease in those with MCI [N=3337]
cognition_mci <- ukb_olink %>% filter(rt>=1.5)
# describe(cognition_mci$rt)
result_list_mci <- cox(data=cognition_mci,exposures=c(organ_pac),
                       outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","death"),
                       covar=c("age","sex","ethnic","edu","tdi","center","pa","smok"))

result_disease_mci <- result_list_mci$result_list %>%   
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
    Outcome == "death" ~ "All-cause mortality",
    TRUE ~ as.character(Outcome)
  ),
  Exposure_label = sub("_res$", "", Exposure)
  ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                             "Schizophrenia","Depression", "Anxiety", "Sleep disorder","All-cause mortality")),
               
               Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                                                              "Lung","Muscle","Pancreas"))
  )


f2_3 <-
  ggplot(result_disease_mci %>% filter(Exposure_label %in% c("Organismal","Brain","Artery","Heart")) %>%
           filter(Outcome %in% c("acd","vd","ad","pd","ms","death")) %>% 
           mutate(Outcome_label=factor(Outcome_label,levels=rev(c("All-cause dementia","Alzheimer's disease","Vascular dementia","Parkinson's disease","Multiple sclerosis","All-cause mortality"))),
                  Exposure_label=factor(Exposure_label,levels=c("Heart","Artery","Organismal","Brain"))
           ),
         aes(x = HR, y=(Outcome_label),groups=Exposure_label)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color = Exposure_label),size = 1,shape=16,position=position_dodge(0.7))+
  #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6))+
  scale_color_manual(values=c("#FF7F50","#A52A2A","black","#4B0082"
  ))

##2.4 PAC and incident disease in all participants [N=43624]

f2_4 <-
  ggplot(result_disease %>% filter(Exposure_label %in% c("Organismal","Brain","Artery","Heart")) %>%
           filter(Outcome %in% c("acd","vd","ad","pd","ms","death")) %>% 
           mutate(Outcome_label=factor(Outcome_label,levels=rev(c("All-cause dementia","Alzheimer's disease","Vascular dementia","Parkinson's disease","Multiple sclerosis","All-cause mortality"))),
                  Exposure_label=factor(Exposure_label,levels=c("Heart","Artery","Organismal","Brain"))
           ),
         aes(x = HR, y=(Outcome_label),groups=Exposure_label)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color = Exposure_label),size = 1,shape=16,position=position_dodge(0.7))+
  #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6))+
  scale_color_manual(values=c("#FF7F50","#A52A2A","black","#4B0082"
  ))

##2.5 Mutivar cox regression

ukb_olink_inci_cox <- result_list$data

ukb_olink_inci <- ukb_olink_inci_cox %>% left_join(olink_full %>% select(eid,GDF15,MAPT,APOE),by="eid") %>% 
  mutate(GDF15=as.numeric(scale(GDF15)),
         MAPT=as.numeric(scale(MAPT)),
         APOE=as.numeric(scale(APOE))
  ) %>% 
  mutate(ad_prs=as.numeric(scale(ad_prs)),
         pd_prs=as.numeric(scale(pd_prs)),
         ms_prs=as.numeric(scale(ms_prs))
  ) %>% 
  mutate(brain_prs=ad_prs+Brain_res) %>% 
  mutate(brain_prs=as.numeric(scale(brain_prs)),
         brain_prs5=cut(brain_prs, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = c("-2", "-1", "0", "+1", "+2"))
  )

# describe(ukb_olink_inci$brain_prs5)

#Outcome with PRS
outcome="ad"
# outcome="pd"
# outcome="ms"

out <-  ukb_olink_inci %>% filter(is.na(!!sym(paste0(outcome, "_base")))) %>% 
  mutate(apoe4x=ifelse(apoe4=="2/x"|apoe4=="3/3"|is.na(apoe4),"0","1") %>% as.factor())

describe(out$apoe4x)
describe(out$edu)
summary(out$apoe4x)
exposure="Brain_res"
# exposure="Organismal_res"
# exposure="brain_prs"

covar = c("age","sex","ethnic","edu","tdi","smok","pa","center","rt",paste0(outcome,"_prs"),"GDF15","APOE","apoe4x")
# covar = c("age","sex","ethnic","edu","tdi","center","rt",paste0(outcome,"_prs"),"GDF15")
# covar = c("age","sex","ethnic","edu","tdi","center","rt","GDF15")

formula_input <-  as.formula(paste0(paste0( "survival::Surv(",outcome, "_fu,",outcome, "_inci", ") ~"), paste(c(exposure, covar), collapse="+")))

cox_summary <- survival::coxph(formula = formula_input, data=out) %>% summary()
cox_summary$coefficients

results <- cox_summary$coefficients[c(exposure,"age","rt","GDF15","APOE","ad_prs","apoe4x1"), ] %>% as.data.frame() %>% 
  mutate(label=row.names(.)) %>% 
  mutate(hr=`exp(coef)`,
         lci=exp(coef-1.96*`se(coef)`),
         hci=exp(coef+1.96*`se(coef)`),
  ) %>% 
  mutate(label=case_when(label=="Brain_res"~"Brain age gap",
                         label=="age"~"Age",
                         label=="rt"~"Cognitive function",
                         str_detect(label, "prs") ~ "PRS",
                         str_detect(label, "apoe4x1") ~ "APOE4/X",
                         T~label
  ))


f2_5 <-ggplot(results, aes(x = hr, y=reorder(label, coef))) +
  geom_point(size = 1,shape=15)+
  geom_errorbarh(aes(xmin = lci, xmax = hci), height = 0,linewidth=0.2) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  scale_x_continuous(limits = c(0.5, 3)) +
  scale_y_discrete(labels = function(x) {
    ifelse(x == "APOE4/X", expression(italic("APOE4/X")), x)
  }) +
  labs(x = "HR (95% CI)", y = NULL, title =NULL ) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
##2.6 KM curve for Brain age+ AD_PRS
library(survminer)

coxph(formula =Surv(ad_fu, ad_inci) ~ brain_prs + age + sex + ethnic + edu + pa + tdi + center, data=out) %>% summary()
coxph(formula =Surv(ad_fu, ad_inci) ~ relevel(brain_prs5,ref="0") + age + sex , data=out) %>% summary()

ad_prs_agegap <- survfit(Surv(ad_fu, ad_inci) ~ brain_prs5, out) %>% 
  ggsurvplot(fun = "event",
             conf.int = F,
             censor = FALSE,
             size = 0.2,
             ggtheme = theme_classic(base_size=6, base_line_size = 0.2, base_family = "Arial"),
             legend.title="PRS+Brain age gap",
             palette = c("#0571b0","#92c5de","pink","#f4a582","#ca0020"),
             break.x.by = 5,
             xlab = "Time (years)",
             legend.labs=c("-2", "-1","0","+1","+2")
  )

ad_prs_agegap$plot <-  ad_prs_agegap$plot+ ggplot2::annotate("text", 
                                                             x = 4, y = 0.08, 
                                                             label = "HR (1-SD) = 2.8 (2.6-3.1)", size = 2)+
  theme_classic(base_line_size = 0.2,base_size = 6)+
  theme(legend.text = element_text(size = 5, color = "black"),
        legend.title =element_text(size = 5),
        legend.key.width = unit(0.15, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        legend.position = "top"
  )+
  guides(colour = guide_legend(nrow = 1))
f2_6 <- ad_prs_agegap$plot


##2.7 Individual association of protein with disease    
##Extract top n proteins---------------------------------------------------------
extract_top_proteins <- function(organ, num_proteins) {
  file_name <- paste("Data/Proteomic/Data/shap/shap_", organ, ".csv", sep = "")
  data <- read.csv(file = file_name) %>% 
    dplyr::select(-eid) %>% 
    mutate_all(abs) %>%  
    summarise_all(mean, na.rm = TRUE) %>% 
    gather(key = "protein", value = "mean_value") %>%
    arrange(desc(mean_value))
  
  if (nrow(data) < num_proteins) {
    num_proteins <- nrow(data)
  }
  
  top_proteins <- data$protein[1:num_proteins]
  return(top_proteins)
}
#Extract the top 20 prot for each organ
top_proteins_list <- list()
for (organ in organ_list) {
  top_proteins <- extract_top_proteins(organ, 20)
  top_proteins_list[[organ]] <- top_proteins
}
top_proteins_list

#Extract the top 20 prot for selected organ
olink_subset <-ukb_olink %>% left_join(olink_full %>% select(c(top_proteins_list$Brain,
                                                               top_proteins_list$Organismal,
                                                               top_proteins_list$Artery,
                                                               top_proteins_list$Heart),"eid"),by="eid") %>% 
  mutate(across(c(top_proteins_list$Brain,
                  top_proteins_list$Organismal,
                  top_proteins_list$Artery,
                  top_proteins_list$Heart),scale)) %>% mutate(across(c(top_proteins_list$Brain,
                                                                       top_proteins_list$Organismal,
                                                                       top_proteins_list$Artery,
                                                                       top_proteins_list$Heart),as.numeric)) %>% 
  mutate(age2=round(age, 0))


cox_prot <- function(data,exposures,outcomes,covar){
  result_list <- data.frame()
  p_values_all <- list()
  
  data <- data %>% mutate(across(where(is.character), ~ifelse(. == "", NA, .))) %>% mutate(death=ifelse(is.na(death),0,1)) %>% 
    mutate(death_date=date_of_death,death_s=NA)
  #Calculated the incident case and follow-up years
  for (var in outcomes){
    data <- data  %>% 
      mutate(!!paste0(var, "_base") := ifelse(as.Date(!!sym(paste0(var, "_date"))) <= as.Date(recruit_date), 1, NA),
             !!paste0(var, "_inci") := ifelse(as.Date(!!sym(paste0(var, "_date"))) > as.Date(recruit_date)&!!paste0(var, "_s")!=0, 1, NA)) %>% 
      mutate(!!paste0(var, "_date") := ifelse(!is.na(!!sym(paste0(var, "_date"))),!!sym(paste0(var, "_date")),ifelse(death==1,date_of_death,max(!!sym(paste0(var, "_date")),na.rm = T))),
             !!paste0(var, "_fu") := as.numeric(as.Date(!!sym(paste0(var, "_date")))-as.Date(recruit_date))/365.25,
             !!paste0(var, "_inci") := ifelse(is.na(!!sym(paste0(var, "_inci"))),0,!!sym(paste0(var, "_inci"))),
             !!paste0(var, "_base2") := ifelse(is.na(!!sym(paste0(var, "_base"))),0,!!sym(paste0(var, "_base"))),
             !!paste0(var) := ifelse(!!sym(paste0(var, "_inci"))==1,1,0)
      )
  }
  
  for (outcome in outcomes) {
    for (exposure in exposures) {
      
      formula_input <-  as.formula( paste0( paste0( "survival::Surv(",outcome, "_fu,",outcome, "_inci", ") ~"), paste( c(exposure, covar), collapse="+")))
      data2 <-  data %>% filter(is.na(!!sym(paste0(outcome, "_base"))))
      cox_model <- survival::coxph(formula = formula_input, data = data2)
      p_values <- summary(cox_model)$coefficients[exposure, "Pr(>|z|)"]
      p_values_all[[paste(outcome, exposure)]] <- p_values
    }
  }
  all_p_values <- unlist(p_values_all)
  
  # Perform FDR correction
  fdr_p_values <- p.adjust(all_p_values, method = "fdr")
  
  for (outcome in outcomes) {
    for (exposure in exposures) {
      formula_input <-  as.formula( paste0( paste0( "survival::Surv(",outcome, "_fu,",outcome, "_inci", ") ~"), paste( c(exposure, covar), collapse="+")))
      data2 <-  data %>% filter(is.na(!!sym(paste0(outcome, "_base"))))
      cox_model <- survival::coxph(formula = formula_input, data =data2)
      cox_summary <- summary(cox_model)
      coef_data <- as.data.frame(cox_summary$coefficients[exposure, "coef"])
      se_data <- as.data.frame(cox_summary$coefficients[exposure, "se(coef)"])
      z_data <- as.data.frame(cox_summary$coefficients[exposure, "z"])
      
      n_cases_inci <- sum(data2[, paste0(outcome, "_inci")]==1)
      n_cases_base <- sum(data[, paste0(outcome, "_base2")]==1)
      
      mean_diff <- data %>%
        summarise(mean_diff = mean(.data[[exposure]][.data[[outcome]] == 1]) - 
                    mean(.data[[exposure]][.data[[outcome]] == 0])) %>%
        pull(mean_diff)
      
      median_follow_up <- median(data2[, paste0(outcome, "_fu")], na.rm = TRUE)
      hr_data <- exp(coef_data)
      hr_ci_data <- exp(confint(cox_model))
      
      # Store results in a list
      result <- data.frame(
        Outcome = outcome,
        Exposure = exposure,
        N_Cases_inci = n_cases_inci,
        N_Cases_base = n_cases_base,
        Mean_diff=mean_diff,
        Median_Follow_Up = median_follow_up,
        
        Beta = coef_data$`cox_summary$coefficients[exposure, "coef"]`,
        SE = se_data$`cox_summary$coefficients[exposure, "se(coef)"]`,
        HR = hr_data$`cox_summary$coefficients[exposure, "coef"]`,
        CI_Lower = hr_ci_data[exposure, "2.5 %"],
        CI_Upper = hr_ci_data[exposure, "97.5 %"],
        z_score = z_data$`cox_summary$coefficients[exposure, "z"]`,
        P_Value = as.numeric(p_values_all[[paste(outcome, exposure)]]),
        FDR_P_Value = as.numeric(fdr_p_values[paste(outcome, exposure)]
                                 
        )
      )
      
      result_list <- rbind(result_list,result)
    }
  }
  return(list(result_list = result_list, data = data))
  
}

result_list_brain<-cox_prot(data=olink_subset,
                            exposures=top_proteins_list$Brain,
                            outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","death"),
                            covar=c("age","sex","ethnic","tdi","center","pa","smok",top_proteins_list$Brain))

#Plotting
result_disease_brain_pac <- result_list_brain$result_list %>%   
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
    Outcome == "death" ~ "All-cause mortality",
    TRUE ~ as.character(Outcome)
  )
  ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                             "Schizophrenia","Depression", "Anxiety", "Sleep disorder","All-cause mortality")),
               
               # Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
               #                                                "Lung","Muscle","Pancreas"))
  ) %>% mutate(z_score_trans=ifelse(P_Value>0.05,0,z_score))

mat_prot <- xtabs(z_score ~ Outcome_label + Exposure, result_disease_brain_pac) 
mat_prot_p <- xtabs(P_Value ~ Outcome_label + Exposure, result_disease_brain_pac) 

mat_prot_trans <- xtabs(z_score_trans ~ Outcome_label + Exposure, result_disease_brain_pac) %>% as.data.frame() %>% 
  group_by(Outcome_label) %>% 
  mutate(Freq=ifelse(Freq<0,-Freq,Freq)) %>% 
  mutate(Freq_pct=Freq/sum(Freq))

mat_prot_trans <- xtabs(Freq_pct ~ Outcome_label + Exposure, mat_prot_trans)
mat_prot_trans_bar <- mat_prot_trans[]
# mat_prot <- mat_prot[order(rowMeans(mat_prot)), ]
# mat_prot_p <- mat_prot_p[order(rowMeans(mat_prot)), ]

col_fun = colorRamp2(c(-10,0,20), c("#3288bd","white","#d53e4f"))
col_panel_prot=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
                          "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#1f78b4",
                          "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
                          "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
                          
f2_7 <- grid.grabExpr(draw(Heatmap(mat_prot,
                                   col = col_fun, na_col="white", show_column_dend = T, show_row_dend = T, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                                   cluster_rows = T,cluster_columns = T,
                                   row_dend_width = unit(0.3, "cm"),
                                   row_dend_gp = gpar(lwd = 0.2),
                                   column_dend_height  = unit(0.3, "cm"),
                                   column_dend_gp = gpar(lwd = 0.2),
                                   width = unit(5, "cm"), height = unit(2, "cm"),
                                   row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                                   column_names_rot=45,
                                   column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                                   #row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",4),rep("All-cause mortality",1)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality")),
                                   # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                                   #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                                   column_title_gp = gpar(fontsize = 6),
                                   column_title_side = "bottom",
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     
                                     if((mat_prot_p)[i, j] < 0.001) {
                                       grid.text("***", x, y,gp = gpar(fontsize = 4))
                                     } else if((mat_prot_p)[i, j] < 0.01) {
                                       grid.text("**", x, y,gp = gpar(fontsize = 4))
                                     } else if((mat_prot_p)[i, j] < 0.05) {
                                       grid.text("*", x, y,gp = gpar(fontsize = 4))
                                     }
                                     
                                     
                                   },
                                   # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                                   #                                  show_legend = FALSE,
                                   #                                  width= unit(1, "mm")
                                   # ),
                                   heatmap_legend_param = list(
                                     #at = c(-2, 0, 2),
                                     #labels = c("0", "10%", "20%", "30%", "40%"),
                                     title =  "z-score",
                                     legend_height = unit(1, "cm"),
                                     grid_width = unit(0.2, "mm"),
                                     #legend_width =unit(1, "mm"),
                                     #title_position = "lefttop-rot",
                                     labels_gp = gpar(col = "black", fontsize = 5),
                                     title_gp=gpar(col = "black", fontsize = 5)
                                     
                                   ),
                                   
                                   right_annotation=rowAnnotation(
                                     " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                                                       gp = gpar(fill =col_panel_prot, 
                                                                 col = col_panel_prot,lwd=0.2),bar_width = 0.8,
                                                       axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                                     ),
                                     
                                     # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                                     #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                                     # 
                                     width = unit(4, "cm"), annotation_name_gp= gpar(fontsize = 6)
                                     # annotation_width=unit(c(1, 2), c("null", "null")),
                                     # gap = unit(3, "mm")
                                     
                                   ),
                                   bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                                                                       col= list(" "=c(
                                                                         "BCAN"="#1b9e77",    "C1QL2"="#d95f02",   "CA14"="#7570b3",    "CDH15"="#e7298a",   "CNTN1"="#66a61e",   
                                                                         "CNTN2"="#e6ab02","COL28A1"="#a6761d", "FGFBP2"="#666666",  "FLT3"="#a6cee3",    "GFAP"="#1f78b4",    
                                                                         "IGSF21"="#b2df8a",  "KLK6"="#33a02c", "MEPE"="#fb9a99",    "MOG"="#e31a1c",     "NEFL"="#fdbf6f",    
                                                                         "OMG"="#ff7f00",     "PTPRR"="#cab2d6",   "RBFOX3"="#6a3d9a","SPOCK1"="#ffff99",  "SYT1"="#b15928")),
                                                                       simple_anno_size = unit(2, "mm"),
                                                                       show_legend = FALSE), 
                                   
                                   row_title =NULL
)))

##2.8 Brain rnq-seq for the top 20 brain prot

###BP functional enrichment
library(clusterProfiler)
library(org.Hs.eg.db)

select_prot_list <- read.csv(file = "Data/Proteomic/Data/boruta_protein_df.csv") %>% filter_all(any_vars(. != ""))
# x <- lapply(prot_list[, c("Organismal", "Brain", "Artery", "Heart")], function(x) x[x != ""])
x <- lapply(prot_list, function(x) x[x != ""])

go_enrichment <- enrichGO(gene = x$Brain, 
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",   
                          ont = "BP",
                          pvalueCutoff = 0.05
)


head(go_enrichment,n=20)


dotplot(go_enrichment, showCategory =10)+theme_bw()
barplot(go_enrichment, showCategory =10)+theme_bw()


enrichment_bp_brain <- as.data.frame(go_enrichment)

GO_brain <- enrichment_bp_brain %>% filter(ID=="GO:0099177"|ID=="GO:0061564") %>% dplyr::select(ID,geneID) %>% 
  mutate(geneID=(strsplit(geneID, "/")))


brain_20_seq <-  read.csv(file = "Data/Proteomic/Data/rna_seq_brain.csv",skip=9) %>% 
  mutate(Cell.Type=case_when(Cell.Type=="aggregated"~"Brain",
                             Cell.Type=="microglial cell"~"Microglia",
                             Cell.Type=="endothelial cell"~"Endothelial cell",
                             T~str_to_title(Cell.Type))) %>% 
  mutate(Cell.Type=factor(Cell.Type,levels=c("Brain","Neuron","Oligodendrocyte","Astrocyte","Microglia","Endothelial cell"))) %>%
  mutate(Tissue.Composition=gsub("%", "", Tissue.Composition)) %>% mutate(Tissue.Composition=as.numeric(Tissue.Composition) / 100) %>% 
  
  mutate(Tissue.Composition=ifelse(Tissue.Composition==0,0.001,Tissue.Composition)) %>%
  mutate(expression_all=ifelse(!is.na(Tissue.Composition)&!is.na(Expression),Expression*(Tissue.Composition),NA)) %>% 
  group_by(Gene.Symbol) %>% mutate(expression_nom=scale(expression_all)) %>% ungroup() %>% 
  mutate(`Synaptic signaling`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[1]] %>% unlist(),1,0),
         `Axon development`=ifelse(Gene.Symbol %in% (GO_brain$geneID)[[2]] %>% unlist(),1,0)
  ) %>%  arrange(Cell.Type,Gene.Symbol)


mat_brain_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, brain_20_seq)
mat_brain_seq[-1,]

col_fun_3 = colorRamp2(c(-1,0,1), c("#2166ac","white","#b2182b"))

f2_8 <- grid.grabExpr(draw(Heatmap(mat_brain_seq[-1,],
                                   col = col_fun_3, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                                   cluster_rows = F,cluster_columns = T,
                                   row_dend_width = unit(0.3, "cm"),
                                   row_dend_gp = gpar(lwd = 0.2),
                                   column_dend_height  = unit(0.3, "cm"),
                                   column_dend_gp = gpar(lwd = 0.2),
                                   width = unit(5, "cm"),
                                   height = unit(1, "cm"),
                                   row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                                   column_names_rot=45,
                                   column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                                   #row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",4),rep("All-cause mortality",1)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality")),
                                   # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                                   #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                                   column_title_gp = gpar(fontsize = 6),
                                   column_title_side = "bottom",
                                   # cell_fun = function(j, i, x, y, width, height, fill) {
                                   #   
                                   #   if((mat_prot_p)[i, j] < 0.001) {
                                   #     grid.text("***", x, y,gp = gpar(fontsize = 4))
                                   #   } else if((mat_prot_p)[i, j] < 0.01) {
                                   #     grid.text("**", x, y,gp = gpar(fontsize = 4))
                                   #   } else if((mat_prot_p)[i, j] < 0.05) {
                                   #     grid.text("*", x, y,gp = gpar(fontsize = 4))
                                   #   }
                                   #   
                                   #   
                                   # },
                                   # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                                   #                                  show_legend = FALSE,
                                   #                                  width= unit(1, "mm")
                                   # ),
                                   heatmap_legend_param = list(
                                     at = c(-1, 0, 1),
                                     #labels = c("0", "10%", "20%", "30%", "40%"),
                                     title =  "Mean expression",
                                     legend_height = unit(1, "cm"),
                                     grid_width = unit(0.1, "mm"),
                                     legend_width =unit(0.1, "mm"),
                                     title_position = "lefttop-rot",
                                     labels_gp = gpar(col = "black", fontsize = 5),
                                     title_gp=gpar(col = "black", fontsize = 5)
                                     
                                   ),
                                   top_annotation = columnAnnotation(`Synaptic signaling`=brain_20_seq %>% filter(Cell.Type=="Brain") %>% dplyr::select(`Synaptic signaling`) %>% as.matrix(),
                                                                     `Axon development`=brain_20_seq %>% filter(Cell.Type=="Brain") %>% dplyr::select(`Axon development`)%>% as.matrix(),
                                                                     "Top proteins"=rep(1,20),
                                                                     col=list(`Synaptic signaling`=colorRamp2(c(0,1), c("white","black")),
                                                                              `Axon development`=colorRamp2(c(0,1), c("white","black")),
                                                                              "Top proteins"=colorRamp2(c(0,1), c("white","darkgreen"))
                                                                     ),
                                                                     gp = gpar(col = "black",lwd=0.2),
                                                                     border=F,
                                                                     simple_anno_size = unit(1.5, "mm"),
                                                                     show_legend = FALSE,
                                                                     annotation_name_gp= gpar(fontsize = 5)
                                   ),
                                   
                                   # right_annotation=rowAnnotation(
                                   #   " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                                   #                     gp = gpar(fill =col_panel_prot, 
                                   #                               col = col_panel_prot,lwd=0.2),bar_width = 0.8,
                                   #                     axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                                   #   ),
                                   #   
                                   #   # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                                   #   #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                                   #   # 
                                   #   width = unit(4, "cm"), annotation_name_gp= gpar(fontsize = 6)
                                   #   # annotation_width=unit(c(1, 2), c("null", "null")),
                                   #   # gap = unit(3, "mm")
                                   #   
                                   # ),
                                   # bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                                   #                                     col= list(" "=c(
                                   #                                       "BCAN"="#1b9e77",    "C1QL2"="#d95f02",   "CA14"="#7570b3",    "CDH15"="#e7298a",   "CNTN1"="#66a61e",   
                                   #                                       "CNTN2"="#e6ab02","COL28A1"="#a6761d", "FGFBP2"="#666666",  "FLT3"="#a6cee3",    "GFAP"="#1f78b4",    
                                   #                                       "IGSF21"="#b2df8a",  "KLK6"="#33a02c", "MEPE"="#fb9a99",    "MOG"="#e31a1c",     "NEFL"="#fdbf6f",    
                                   #                                       "OMG"="#ff7f00",     "PTPRR"="#cab2d6",   "RBFOX3"="#6a3d9a","SPOCK1"="#ffff99",  "SYT1"="#b15928")),
                                   #                                     simple_anno_size = unit(2, "mm"),
                                   #                                     show_legend = FALSE), 
                                   
                                   row_title =NULL
)))
plot_grid(f2_8) 

##2.9 Changes with aging and AD in featured proteins of brain aging or dementia 
# across plasma and brain (levels of protein, bulk RNA and sc RNA expression
brain_seq <- readRDS(file = "Data/Proteomic/Data/brain_seq.rds")
library(Seurat)
seurat_brain = UpdateSeuratObject(brain_seq)

brain_rna <- brain_seq@assays[["RNA"]]@data
brain_summary <- brain_seq@meta.data

unique(brain_summary$cell_type)

brain_gene_list <- brain_seq@assays[["RNA"]]@meta.features
brain_gene_list <- brain_gene_list %>% filter(feature_name %in% c("NEFL","GFAP","KLK6","SYT1","BCAN","CNTN1",
                                                                  "NRXN3","APP","APPL2","MAPT","GDF15"
))

rename_info <- as.list(setNames(row.names(brain_gene_list),brain_gene_list$feature_name))


brain_summary_filtered <- brain_summary %>%  
  filter(grepl(paste(c("oligodendrocyte",
                       "endothelial cell",
                       "mature astrocyte",
                       "neuron",
                       "cell",
                       "mature microglial cell"), collapse = "|"), cell_type))



brain_rna_filtered <- brain_rna[row.names(brain_gene_list), row.names(brain_summary_filtered)] %>% t() %>% as.data.frame()

brain_filtered_all <- cbind(brain_summary_filtered,brain_rna_filtered) %>% 
  rename_with(~ as.character(brain_gene_list$feature_name), as.character(row.names(brain_gene_list))) %>% 
  group_by(SampleID,clusterCellType) %>% 
  mutate(across(as.character(brain_gene_list$feature_name), ~ log(. + 1))) %>% 
  mutate(across(as.character(brain_gene_list$feature_name), ~ (.x), .names = "{.col}_sum"))  %>% ungroup() %>% 
  group_by(donor_id) %>% dplyr::select(donor_id,
                                       clusterCellType,
                                       cell_type,
                                       disease,
                                       ends_with("_sum"),
                                       development_stage,
                                       nGene,
                                       clusterAssignment
  ) %>% group_by(donor_id,clusterCellType) %>% 
  mutate(disease_ad=ifelse(disease=="normal","0","1"),
         #disease_ad=as.numeric(disease_ad)
  ) %>% 
  mutate(age=substr(development_stage, 1, 2))
# group_by(cell_type) %>% 
# mutate_at(vars(ends_with("_sum")), scale) %>% ungroup() %>% 
# mutate_at(vars(ends_with("_sum")), as.numeric)


gene_list <- as.character(brain_gene_list$feature_name)
cell_list <- c("mature astrocyte",
               "endothelial cell",
               "mature microglial cell",
               "oligodendrocyte",
               "glutamatergic neuron",
               "GABAergic neuron",
               "all")

run_by_cell <- function(data, cell_list, genelist) {
  results <- data.frame()
  for (cell in cell_list) {
    #data <- brain_filtered_all
    if (cell=="all") {subset_data <- data 
    
    }else{subset_data <- data %>% filter(cell_type==cell)}
    
    for (gene in genelist) {
      x <- subset_data$disease_ad
      y <- subset_data[[paste0(gene, "_sum")]]
      
      # model <- glm(y ~ x, family = binomial(link = "logit"))
      model <- lm(y ~ x+development_stage,data = subset_data) %>% summary()
      
      coef <- model$coefficients["x1","Estimate"]
      se <- model$coefficients["x1","Std. Error"]
      p_value <- model$coefficients["x1","Pr(>|t|)"]
      
      result_row <- data.frame(
        cell=cell,
        gene=gene,
        cell_gene = paste0(cell,"_",gene),
        coef = coef,
        se = se,
        p_value = p_value
      )
      results <- rbind(results, result_row)
    }
  }
  
  results$q <- p.adjust(results$p_value, method = "fdr")
  return(results)
}

brain_coef_list <- run_by_cell(brain_filtered_all,cell_list, gene_list)

brain_coef_list_filtered <- brain_coef_list %>% 
  filter(cell %in% c("mature astrocyte",
                     "oligodendrocyte",
                     "glutamatergic neuron",
                     "GABAergic neuron",
                     "mature microglial cell",
                     #"endothelial cell",
                     "all"))
mat_braincell_seq <- xtabs(coef ~ cell + gene, brain_coef_list_filtered) 
mat_braincell_seq
rownames(mat_braincell_seq) <-            c("Brain RNA (AD/HC)",
                                            "GABAergic neuron scRNA (AD/HC)",
                                            "Glutamatergic neuron scRNA (AD/HC)",
                                            "Astrocyte scRNA (AD/HC)",
                                            "Microglial scRNA (AD/HC)",
                                            "Oligodendrocyte scRNA (AD/HC)"
)   

mat_braincell_seq_q <- xtabs(q ~ cell + gene, brain_coef_list_filtered)

###Plasma protein and AD/HC or aging
ukb_olink_prot_brain <- ukb_olink %>% left_join(olink_full %>% dplyr::select("eid",as.character(brain_gene_list$feature_name)),
                                                by="eid") %>% 
  mutate(across(as.character(brain_gene_list$feature_name),scale)) %>% 
  mutate(across(as.character(brain_gene_list$feature_name),as.numeric)) %>% 
  mutate(age_scale=scale(age) %>% as.numeric())


result_list_brain_ad<-cox_prot(data=ukb_olink_prot_brain,
                               exposures=as.character(brain_gene_list$feature_name),
                               outcomes=c("ad"), #"acd","ad","vd","pd","ms","scz","dep","anx","sd","death"
                               covar=c("age","sex","ethnic","tdi","center","pa"))


result_disease_brain_ad <- result_list_brain_ad$result_list %>%   
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
    Outcome == "death" ~ "All-cause mortality",
    TRUE ~ as.character(Outcome)
  )
  )


result_brain_ad <- result_disease_brain_ad %>% dplyr::select(Exposure,Outcome,Beta,FDR_P_Value) %>% rbind(.,
                                                                                                          
                                                                                                          linear(ukb_olink_prot_brain,
                                                                                                                 a=as.character(brain_gene_list$feature_name),
                                                                                                                 b="age",covar = c("age","sex","ethnic","edu","tdi","center","pa")) %>% 
                                                                                                            mutate(q=p.adjust(p, method = "fdr")) %>% dplyr::select(var_a,var_b,beta,q) %>% 
                                                                                                            rename("Exposure"="var_a","Outcome"="var_b","Beta"="beta","FDR_P_Value"="q")
)

mat_brain_prot_ad <- xtabs(Beta ~ Outcome+Exposure, result_brain_ad) 
mat_brain_prot_ad_q <- xtabs(FDR_P_Value ~ Outcome+Exposure, result_brain_ad) 


mat_brain_ad <- rbind(mat_brain_prot_ad,mat_braincell_seq)
mat_brain_ad_q <- rbind(mat_brain_prot_ad_q,mat_braincell_seq_q)


rownames(mat_brain_ad) <-            c(    "Plasma protein (aging)",
                                           "Plasma protein (AD/HC)",
                                           "Brain RNA (AD/HC)",
                                           "GABAergic neuron scRNA (AD/HC)",
                                           "Glutamatergic neuron scRNA (AD/HC)",
                                           "Astrocyte scRNA (AD/HC)",
                                           "Microglial scRNA (AD/HC)",
                                           "Oligodendrocyte scRNA (AD/HC)"
)   
###Ploting
col_fun_4 = colorRamp2(c(-0.25,0,0.5), c("#2166ac","white","#b2182b"))
f2_9 <- grid.grabExpr(draw(Heatmap(mat_brain_ad,
                                   col = col_fun_4, na_col="white", show_column_dend = F, show_row_dend = F, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                                   cluster_rows = F,cluster_columns = T,
                                   row_dend_width = unit(0.3, "cm"),
                                   row_dend_gp = gpar(lwd = 0.2),
                                   column_dend_height  = unit(0.3, "cm"),
                                   column_dend_gp = gpar(lwd = 0.2),
                                   width = unit(2.5, "cm"),
                                   height = unit(1.5, "cm"),
                                   row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                                   column_names_rot=45,
                                   column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                                   row_split=factor(c("Plasma","Plasma",rep("Brain",6)),levels = c("Plasma", "Brain")),
                                   row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                                   column_split =factor(c(0,0,1,1,0,1,1,0,1,0,1)),
                                   column_title = NULL,
                                   column_title_gp = gpar(fontsize = 6),
                                   column_title_side = "bottom",
                                   cell_fun = function(j, i, x, y, width, height, fill) {
                                     
                                     if((mat_brain_ad_q)[i, j] < 0.001) {
                                       grid.text("***", x, y,gp = gpar(fontsize = 4))
                                     } else if((mat_brain_ad_q)[i, j] < 0.01) {
                                       grid.text("**", x, y,gp = gpar(fontsize = 4))
                                     } else if((mat_brain_ad_q)[i, j] < 0.05) {
                                       grid.text("*", x, y,gp = gpar(fontsize = 4))
                                     }
                                     
                                     
                                   },
                                   # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                                   #                                  show_legend = FALSE,
                                   #                                  width= unit(1, "mm")
                                   # ),
                                   heatmap_legend_param = list(
                                     at = c(-0.25, 0, 0.25),
                                     #labels = c("0", "10%", "20%", "30%", "40%"),
                                     title =  "beta",
                                     legend_height = unit(1, "cm"),
                                     grid_width = unit(0.1, "mm"),
                                     #legend_width =unit(1, "mm"),
                                     #title_position = "lefttop-rot",
                                     labels_gp = gpar(col = "black", fontsize = 5),
                                     title_gp=gpar(col = "black", fontsize = 5)
                                     
                                   ),
                                   row_title =NULL
)))

plot_grid(f2_9)     


plot_grid(f2_1,f2_2,f2_3,f2_4,f2_5,f2_6,rel_widths = c(1,0.7,1,1,0.7,1),
                nrow = 2, ncol= 3,align = "h",labels = "auto",label_size = 7,hjust =-1,vjust = 1.5) %>% 
  plot_grid(.,plot_grid(f2_7,plot_grid(f2_8,f2_9,nrow=1,ncol=2,labels = c("h","i"),label_size = 7,
                                       hjust =-1,vjust = 1.5,rel_widths = (5:2)),
                        nrow = 2,labels = c("g",""),label_size = 7,rel_heights = c(2.5:1)), 
            nrow = 2,hjust =-1,vjust = 1.5,rel_heights = c(4.6:2.4))


#Fig 4 Body and brain PAC in cognitive decline and neurodegenerative disease--------------------------------------------------------------  

olink_subset <-ukb_olink %>% left_join(olink_full %>% dplyr::select(c(top_proteins_list$Brain,
                                                                      top_proteins_list$Organismal,
                                                                      top_proteins_list$Artery,
                                                                      top_proteins_list$Heart),"eid"),by="eid") %>% 
  mutate(across(c(top_proteins_list$Brain,
                  top_proteins_list$Organismal,
                  top_proteins_list$Artery,
                  top_proteins_list$Heart),scale)) %>% mutate(across(c(top_proteins_list$Brain,
                                                                       top_proteins_list$Organismal,
                                                                       top_proteins_list$Artery,
                                                                       top_proteins_list$Heart),as.numeric)) %>% 
  mutate(age2=round(age, 0))


data_test <- olink_subset %>% filter(age2>=40)


prot_diease_by_organ <- 
  function(organ_list,disease){
    heatmap_list <- list()
    sig_prot <- data.frame()
    for (organ in organ_list){
      
      
      # ukb_olink_prot <- ukb_olink %>% left_join(olink_full %>% select(brain_20$gene,"eid"),by="eid") %>% 
      #   mutate(across(brain_20$gene,scale)) %>% mutate(across(brain_20$gene,as.numeric))
      # 
      result_list_organ<-cox_prot(data=olink_subset,
                                  exposures=c(top_proteins_list[[organ]]),
                                  outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","death"),
                                  covar=c("age","sex","ethnic","tdi","center","smok","pa",top_proteins_list[[organ]]))
      
      
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
          Outcome == "death" ~ "All-cause mortality",
          TRUE ~ as.character(Outcome)
        )
        ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                                   "Schizophrenia","Depression", "Anxiety", "Sleep disorder","All-cause mortality")),
                     
                     # Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                     #                                                "Lung","Muscle","Pancreas"))
        ) %>% mutate(z_score_trans=ifelse(P_Value>0.05,0,z_score))
      
      mat_prot <- xtabs(z_score ~ Outcome_label + Exposure, result_disease_organ_pac) 
      mat_prot_p <- xtabs(P_Value ~ Outcome_label + Exposure, result_disease_organ_pac) 
      
      mat_prot_trans_pre <- xtabs(z_score_trans ~ Outcome_label + Exposure, result_disease_organ_pac) %>% as.data.frame() %>% 
        group_by(Outcome_label) %>% 
        mutate(Freq=ifelse(Freq<0,-Freq,Freq)) %>% 
        mutate(Freq_pct=Freq/sum(Freq))
      
      mat_prot_trans <- xtabs(Freq_pct ~ Outcome_label + Exposure, mat_prot_trans_pre)
      mat_prot_trans_bar <- mat_prot_trans[]
      
      mat_prot_sig <- mat_prot_trans_pre %>% ungroup() %>% 
        filter(Outcome_label %in% disease & Freq_pct > 0) %>%
        dplyr::select(Outcome_label,Exposure)  %>% mutate(organ=organ)
      
      sig_prot <- rbind(sig_prot, mat_prot_sig)
      # mat_prot <- mat_prot[order(rowMeans(mat_prot)), ]
      # mat_prot_p <- mat_prot_p[order(rowMeans(mat_prot)), ]
      
      col_fun = colorRamp2(c(-10,0,20), c("#3288bd","white","#d53e4f"))
      
      col_panel_prot=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
                                "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#1f78b4",
                                "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
                                "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
                                
      
      col_panel_prot_subset <- col_panel_prot[1:length(top_proteins_list[[organ]])]
      
      # paste0(colnames(mat_prot),"'","=","'",col_panel_prot)
      
      heatmap_list[[organ]] <- grid.grabExpr(draw(Heatmap(mat_prot,
                                                          col = col_fun, na_col="white", show_column_dend = T, show_row_dend = T, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                                                          cluster_rows = T,cluster_columns = T,
                                                          row_dend_width = unit(0.3, "cm"),
                                                          row_dend_gp = gpar(lwd = 0.2),
                                                          column_dend_height  = unit(0.3, "cm"),
                                                          column_dend_gp = gpar(lwd = 0.2),
                                                          width = unit(5, "cm"), height = unit(2, "cm"),
                                                          row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                                                          column_names_rot=45,
                                                          column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                                                          #row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",4),rep("All-cause mortality",1)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality")),
                                                          # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                                                          #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                                                          column_title_gp = gpar(fontsize = 6),
                                                          column_title_side = "top",
                                                          column_title =paste0(organ," age gap"),
                                                          cell_fun = function(j, i, x, y, width, height, fill) {
                                                            
                                                            if((mat_prot_p)[i, j] < 0.001) {
                                                              grid.text("***", x, y,gp = gpar(fontsize = 4))
                                                            } else if((mat_prot_p)[i, j] < 0.01) {
                                                              grid.text("**", x, y,gp = gpar(fontsize = 4))
                                                            } else if((mat_prot_p)[i, j] < 0.05) {
                                                              grid.text("*", x, y,gp = gpar(fontsize = 4))
                                                            }
                                                            
                                                            
                                                          },
                                                          # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                                                          #                                  show_legend = FALSE,
                                                          #                                  width= unit(1, "mm")
                                                          # ),
                                                          heatmap_legend_param = list(
                                                            #at = c(-2, 0, 2),
                                                            #labels = c("0", "10%", "20%", "30%", "40%"),
                                                            title =  "z-score",
                                                            legend_height = unit(1, "cm"),
                                                            grid_width = unit(0.2, "mm"),
                                                            #legend_width =unit(1, "mm"),
                                                            #title_position = "lefttop-rot",
                                                            labels_gp = gpar(col = "black", fontsize = 5),
                                                            title_gp=gpar(col = "black", fontsize = 5)
                                                            
                                                          ),
                                                          
                                                          right_annotation=rowAnnotation(
                                                            " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                                                                              gp = gpar(fill =col_panel_prot,
                                                                                        col = col_panel_prot,
                                                                                        lwd=0.2),bar_width = 0.8,
                                                                              axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                                                            ),
                                                            
                                                            # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                                                            #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                                                            # 
                                                            width = unit(4, "cm"), annotation_name_gp= gpar(fontsize = 6)
                                                            # annotation_width=unit(c(1, 2), c("null", "null")),
                                                            # gap = unit(3, "mm")
                                                            
                                                          ),
                                                          bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                                                                                              col= list(" "=setNames(col_panel_prot_subset, levels(mat_prot_trans_pre$Exposure))),
                                                                                              simple_anno_size = unit(2, "mm"),
                                                                                              show_legend = FALSE)
      )))
    }
    return(list(heatmap_list=heatmap_list,sig_prot=sig_prot))
  }

#Selected the aging-related prot [brain, organ, artery, heart] that also associated with AD
heatmap <- prot_diease_by_organ(c(organ_list="Brain","Organismal","Artery","Heart"),disease = c("All-cause dementia",
                                                                                                "Alzheimer's disease",
                                                                                                "Vascular dementia"))
#Individual prot and disease by organ [For AD]
plot_grid(plotlist = heatmap$heatmap_list, ncol = 1,nrow=length(heatmap$heatmap_list))

#Prioritized proteins related to AD
heatmap$sig_prot %>% distinct(Exposure) %>% pull()

#LOWESS Plot for Change with Age 
library(ggrepel) 

generate_loess_geoms <- function(data,var_list) {
  plot <- ggplot()+geom_hline(yintercept = 0, linetype = "dashed", color = "gray", lwd=0.2)+
    theme_classic(base_line_size = 0.2,base_size = 6)+labs(title = NULL, x = "Age", y = "Z-scored protein level")+
    ylim(-2,3)+xlim(39,75)
  
  # var_list <- data %>% dplyr::select(-eid, -age, -age2,-sex) %>% colnames()
  # var_list <-   heatmap$sig_prot %>% distinct(Exposure) %>% pull()
  # var_list <- ad_sig_prot$prot
  var_list <- var_list
  
  
  top_protein_endpoints <- data.frame(var = character(), x = numeric(), y = numeric(), color= character(), stringsAsFactors = FALSE)
  
  for (var in var_list) {
    
    #start_value <- colMeans(data[data$age2 == 40, var], na.rm = TRUE)
    
    mean_40 <-  data %>% filter(age2 == 40) %>% dplyr::select(var) %>%
      colMeans(na.rm = TRUE)
    
    sd <- sd(data[[var]], na.rm = TRUE)
    
    data[[var]] <- ifelse(data$age2 == 40, 0, (data[[var]] - mean_40) / sd)
    
    #Brain
    # values=c("black","#4B0082","#A52A2A","#FF7F50")
    color <-if (var %in% top_proteins_list$Brain) "#4B0082"
      else if (var %in% top_proteins_list$Artery) "#A52A2A"
        else if (var %in% top_proteins_list$Heart) "#FF7F50"
          else "black"
            
    
    loess_fit <- loess(as.formula(paste0(var, "~age2")), data = data)
    predicted_values <- predict(loess_fit, newdata = data.frame(age2 = data$age2))
    offset <- predicted_values[data$age2 == 40][1]
    adjusted_predicted_values <- predicted_values - offset
    adjusted_data <- data.frame(age2 = data$age2, adjusted_values = adjusted_predicted_values)
    
    endpoint <- data.frame(var = var, x = 71, y = adjusted_predicted_values[data$age2 == 71], color = color)   
    
    top_protein_endpoints <- rbind(top_protein_endpoints, endpoint[1,])
    
    
    plot <- plot +
      geom_line(aes(x = age2, y = adjusted_values), color = color, size = 0.3, data = adjusted_data)
    
  }
  
  top_protein_endpoints <- top_protein_endpoints %>%
    arrange(desc(abs(y))) %>%
    head(10)
  

  
  # Add labels for the top 10 proteins
  plot <- plot + 
    geom_text_repel(data = top_protein_endpoints, aes(x = x, y = y, label = var,color = color), direction = "both",
                    size = 1.5, segment.size = 0.2, nudge_x=2, max.overlaps = 10,box.padding=0.3, segment.color = "black")+scale_color_identity() 
  
  return(plot)
}

traj_plot <- generate_loess_geoms(data=data_test,
                                  var_list=heatmap$sig_prot$Exposure %>% as.character() %>% unique())
f4.1 <- traj_plot

# Plot for Change with Age 
{
  data <- olink_subset
  
  # scale_variable <- function(x) {
  #   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  # }
  # data <- data %>%
  #   mutate(across(-c(eid,age,sex,age2), scale_variable))
  
  names(data)
  
  lm_results <- lapply(names(data %>% dplyr::select(heatmap$sig_prot$Exposure)), function(var) {
    lm_model <- lm(data[[var]] ~ age+sex, data = data)
    p_value <- summary(lm_model)$coefficients["age", "Pr(>|t|)"]
    return(p_value)
  })
  
  p_values <- unlist(lm_results)
  fdr_corrected_p_values <- p.adjust(p_values, method = "fdr")
  neg_log10_fdr_p_values <- -log10(fdr_corrected_p_values)
  
  effect_of_age <- sapply(names(data %>% dplyr::select(heatmap$sig_prot$Exposure)), function(var) {
    lm_model <- lm(data[[var]] ~ age+sex, data = data)
    effect <- coef(lm_model)["age"]
    return(effect)
  })
  volcano_data <- data.frame(effect_of_age = effect_of_age, neg_log10_fdr_p_values = neg_log10_fdr_p_values)
  volcano_data$neg_log10_fdr_p_values[is.infinite(volcano_data$neg_log10_fdr_p_values)] <- 291
  volcano_data$protein <- sub("\\.age$", "", row.names(volcano_data)) 
  
  protein_colors <- c(
    Brain = "#4B0082",
    Artery = "#A52A2A",
    Heart = "#FF7F50",
    Other = "black"  # Default color for proteins not falling into the specified categories
  )
  
  # Assign colors based on protein categories
  volcano_data$protein_category <- ifelse(volcano_data$protein %in% top_proteins_list$Brain, "Brain",
                                          ifelse(volcano_data$protein %in% top_proteins_list$Artery, "Artery",
                                                 ifelse(volcano_data$protein %in% top_proteins_list$Heart, "Heart", "Other")))
  
  f4.2 <- ggplot(volcano_data, aes(x = effect_of_age, y = neg_log10_fdr_p_values)) +
    geom_point(aes(size = neg_log10_fdr_p_values,color=protein_category),alpha = 0.7) + 
    scale_color_manual(values = protein_colors) +
    scale_size_continuous(range = c(0.2, 3))+
    geom_text_repel(aes(label = protein), size = 1.5, max.overlaps = 30, direction= "both",segment.size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted",size=0.2)+
    labs(title = NULL,
         x = "Effect of Age",
         y = "-log10(FDR-corrected p-value)")+theme_classic(base_size = 5,base_line_size = 0.2)+ylim(0,320)+
    theme(text = element_text(size = 5),
          line = element_line(linewidth = 0.2),
          legend.position = "none"
    )
  
  
  targeted_proteins_neuro <- volcano_data %>% arrange(desc(abs(effect_of_age)),desc(neg_log10_fdr_p_values)) %>%
    # filter(protein_category %in% c("Other","Artery")) %>%
    head(12) %>% dplyr::select(protein) %>% pull()
  
  
  
}


#Proteins interaction based on STING database
library(STRINGdb)
library(igraph)

#Conduct interaction analyses on: https://string-db.org/
raw_data <- fread(input = "Data/string_interactions_ad.tsv")
# Check the structure of the raw_data data frame
str(raw_data)

# Create an igraph object from the data frame
g <- graph_from_data_frame(d = raw_data[, c("node1", "node2")], directed = FALSE)
E(g)$width <- raw_data$combined_score *1.1


node_color_mapping <- list(
  "ELN" = "#A52A2A",
  "LTBP2" = "#A52A2A",
  "NEFL" = "#4B0082",
  "GDF15" = "black",
  "GFAP" = "#4B0082",
  "TNFRSF11B" = "#A52A2A",
  "IGDCC4" = "black",
  "COMP" = "#A52A2A"
    
)

# Create a vector of colors for nodes based on the mapping
node_colors <- sapply(V(g)$name, function(node) {
  if (node %in% names(node_color_mapping)) {
    return(node_color_mapping[[node]])
  } else {
    return("gray")  # Assign gray color if node is not found in the mapping
  }
})


vertex_label_colors <- rep("black", length(V(g)))  # Set all labels to dark red

# Plot the network

par(mar = c(0, 0, 0, 0) + 0.1)
f4.3 <- plot(g,
             vertex.label = V(g)$name,  # Remove automatic labeling
             vertex.color = node_colors,
             vertex.label.color = vertex_label_colors,
             vertex.label.cex = 0.5, 
             vertex.label.dist = 1.5,
             vertex.size = 8,
             edge.width = E(g)$width,
             vertex.label.family = "sans"
)

c(raw_data$node1,raw_data$node2) %>% unique()

#Enrichment of Proteins interaction

library(clusterProfiler)
library(org.Hs.eg.db)

go_enrichment <- enrichGO(gene = c(raw_data$node1,raw_data$node2) %>% unique(), 
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",   
                          ont = "BP",
                          # pvalueCutoff = 0.05
)

go_enrichment_bg <- enrichGO(gene = c(raw_data$node1,raw_data$node2) %>% unique(), 
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",   
                          ont = "BP",
                          universe = olink_list$GeneName
                          #pvalueCutoff = 0.05
)


head(go_enrichment,n=20)
head(go_enrichment_bg,n=20)

f4.4 <- dotplot(go_enrichment, showCategory =6)+
  guides(size = FALSE,color = guide_colorbar(barwidth = 0.5,barheight = 2))+
  geom_point(aes(size =1.5))+
  scale_size_continuous(range = c(0.5, 2))+ 
  scale_color_gradient(high = "#d0d1e6", low ="#023858")+
  labs(color = "q value",x="Gene ratio")+
  theme_bw(base_size = 6,base_line_size = 0.1)+
  theme(panel.grid = element_blank(),  
        axis.text.x = element_text(size = 6,colour ="black"),
        axis.text.y = element_text(size = 6,colour ="black"),
  )+
  scale_y_discrete(labels = function(x) (str_wrap(str_to_sentence(x), width = 30)))+
  scale_x_continuous(limits = c(0.15,0.45), breaks = c(0.2, 0.4))+
  theme(legend.text = element_text(size = 5))



library(ggplotify)

par(mar = c(0, 0, 0, 0) + 0.1)
plot_grid(f4.1,f4.2,ncol=1,align ="v",labels = "auto",label_size = 6) %>% 
  plot_grid(.,
            plot_grid(
              plot_grid(base2grob(~plot(g,
                                        vertex.label = V(g)$name,  # Remove automatic labeling
                                        vertex.color = node_colors,
                                        vertex.label.color = vertex_label_colors,
                                        vertex.label.cex = 0.4, 
                                        vertex.label.dist = 2,
                                        vertex.size = 3,
                                        edge.width = E(g)$width,
                                        vertex.label.family = "sans")),f4.4,
                        
                        labels = c("c","d"),label_size = 6,nrow=1,rel_widths = c(1.2,1),rel_heights = c(1.5,1)),
              plot_grid(f4.5,f4.6,align = "hv",rel_widths = c(1.5,1),labels=c("e","f"),label_size = 6),labels=c("",""),label_size = 6,
              nrow = 2,ncol = 1,rel_heights = c(1.2,1)
            )
            
            ,nrow = 1,ncol = 2,rel_widths = c(1,2)
            
  ) %>% plot_grid(.,NULL,labels=c("","g"),label_size = 6, ncol=1,nrow=2,rel_heights = c(1,0.3))


#Expression of selected protein in vascular [of body or brain]
#Whole body
library(stringr)
vas_seq <-  read.csv(file = "Data/Proteomic/Data/rna_seq_vascular (2).csv",skip=9) %>% 
  mutate(Cell.Type=case_when(
    Cell.Type=="smooth muscle cell"~"SMC",
    T~str_to_title(Cell.Type))) %>%
  filter(Tissue=="vasculature") %>% 
  filter(Cell.Type %in% c("Fibroblast","Endothelial Cell","SMC","Pericyte","Epithelial Cell")) %>% 
  mutate(Cell.Type=factor(Cell.Type,levels=c("Epithelial Cell","SMC","Endothelial Cell","Fibroblast","Pericyte"))) %>%
  mutate(Tissue.Composition=gsub("%", "", Tissue.Composition)) %>% mutate(Tissue.Composition=as.numeric(Tissue.Composition)) %>% 
  mutate(Tissue.Composition=ifelse(Tissue.Composition==0,0.001,Tissue.Composition)) %>%
  mutate(fraction=Number.of.Cells.Expressing.Genes/Cell.Count) %>% 
  mutate(expression_all=ifelse(!is.na(Tissue.Composition)&!is.na(Expression),Expression*(Tissue.Composition),NA)) %>% 
  group_by(Gene.Symbol) %>% mutate(expression_nom=scale(expression_all)) %>% ungroup() %>%  arrange(Cell.Type,Gene.Symbol)


mat_brain_seq <- xtabs(`expression_nom` ~ Cell.Type + Gene.Symbol, vas_seq)

col_fun_3 = colorRamp2(c(-1,0,1), c("#2166ac","white","#b2182b"))



library(Seurat)
library(ggplot2)
library(Matrix)

brainvas_expression_matrix <- readMM("Data/Proteomic/Data/matrix.mtx")
brainvas_genes <- read.delim("Data/Proteomic/Data/features.tsv", header = FALSE)
barcodes <- read.delim("Data/Proteomic/Data/barcodes.tsv", header = FALSE)
metadata <- read.delim("Data/Proteomic/Data/meta.tsv", header = TRUE)

rownames(brainvas_expression_matrix) <- brainvas_genes$V1
colnames(brainvas_expression_matrix) <- barcodes$V1

seurat_object <- CreateSeuratObject(counts = brainvas_expression_matrix)
metadata <- metadata[match(colnames(seurat_object), metadata$Cell),]
seurat_object <- AddMetaData(seurat_object, metadata = metadata$Cell_Type, col.name = "cell_type")
seurat_object <- AddMetaData(seurat_object, metadata = metadata$Treat, col.name = "treat")
seurat_object <- AddMetaData(seurat_object, metadata = metadata$Sample, col.name = "sample")


seurat_object_control <- subset(seurat_object, subset = treat == "Control")

#Target GENE
selected_genes <- c(targeted_proteins,"FBLN2","MFAP2","MFAP5","FBLN5","HSPG2")
selected_expression <-  GetAssayData(seurat_object_control, slot = "counts")[selected_genes, ]


# Get cell type information
cell_type_info <- seurat_object_control@meta.data$cell_type

# Calculate mean expression for each gene in each cell type
mean_expression <- aggregate(t(as.matrix(selected_expression)), by = list(cell_type_info), FUN = mean)
colnames(mean_expression) <- c("CellType", selected_genes)

# Calculate fraction of cells expressing each gene in each cell type
fraction_expression <- aggregate(t(as.matrix(selected_expression)) > 0, by = list(cell_type_info), FUN = function(x) mean(x) * 100)
colnames(fraction_expression) <- c("CellType", selected_genes)

# Merge mean expression and fraction expression data
expression_data <- merge(mean_expression, fraction_expression, by = "CellType", suffixes = c("_mean", "_frac")) %>% 
  filter(CellType!="Ependymal"&CellType!="Astrocyte"&CellType!="Neuron"&CellType!="Microglia/Mφ"&CellType!="Oligo"&CellType!="T cell")

# Convert to long format for plotting
expression_long <- melt(expression_data, id.vars = "CellType")
expression_long$Type <- ifelse(grepl("_mean$", expression_long$variable), "MeanExpression", "Fraction")
expression_long$Gene <- gsub("_mean$|_frac$", "", expression_long$variable)
expression_long <- dcast(expression_long, CellType + Gene ~ Type, value.var = "value")

# Brain+Body seq


vas_seq_all <- expression_long %>% filter(CellType %in% c("Arterial","P. Fibro","Pericyte","SMC")) %>% mutate(group="Brain") %>% 
  rbind(.,
        vas_seq %>% dplyr::select(Cell.Type,Gene.Symbol,Expression..Scaled,fraction) %>% filter(Cell.Type !="Epithelial Cell") %>% 
          rename("CellType"="Cell.Type",
                 "Gene"="Gene.Symbol",
                 "MeanExpression"="Expression..Scaled",
                 "Fraction"="fraction") %>% mutate(group="Body",Fraction=Fraction*100)
        
  ) %>% mutate(CellType=case_when(CellType=="P. Fibro"~"Fibroblast",
                                  CellType=="Arterial"~"Endothelial Cell",
                                  T~CellType
  )
  
  ) %>% group_by(group) %>% mutate(CellType=factor(CellType,c("SMC","Fibroblast","Endothelial Cell","Pericyte"))) %>%
  filter(!is.na(MeanExpression))

f4.5 <-  ggplot(vas_seq_all, aes(x = Gene, y = CellType)) +
  geom_point(aes(size = Fraction, color = MeanExpression)) +
  scale_size_continuous(range = c(0, 3),  breaks = seq(0, 60, by = 20)) +
  scale_color_gradient(breaks = c(0, 1),  low = "lightgrey", high = "black") +
  facet_wrap(~group)+
  theme_classic(base_size = 6,base_line_size = 0.2) +
  labs(size = "Fraction of cells", color = "Mean expression") +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle = 45,hjust = 1,color="black"),
        axis.text.y = element_text(color="black"),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),  
        legend.text = element_text(angle = 0, hjust = 0.5, size = 5),
        legend.position = "bottom",
        legend.key.size = unit(0.3, "cm")
  )+
  theme(legend.direction = "horizontal") + 
  guides(size = guide_legend(label.position = "bottom", title.position = "top", title.hjust = 0.5, keywidth = 0.5, keyheight = 0.1),
         color=guide_colorbar(label.position = "bottom", title.position = "top", title.hjust = 0.5, keywidth = 2.5, keyheight = 0.2)
  )

#Change of cell type between AD and normal control

library(ggplot2)
library(ggsignif)

head(seurat_object@meta.data,n=50
)

cell_proportions <- seurat_object@meta.data %>%
  group_by(sample, treat, cell_type) %>%
  summarise(count = n()) %>%
  group_by(treat) %>%
  mutate(total_count = sum(count),
         percentage = (count / total_count) * 100) %>%
  ungroup() 

cell_proportions <- seurat_object@meta.data %>%
  group_by(sample, treat, cell_type) %>%
  summarise(total_count = sum(nCount_RNA)) %>%
  group_by(sample) %>%
  mutate(percentage = (total_count / sum(total_count))) %>%
  ungroup()



p_cell <- cell_proportions %>%
  filter(cell_type %in% c("Pericyte", "SMC", "P. Fibro","Arterial")) %>%
  group_by(cell_type) %>%
  summarise(p_value = t.test(percentage ~ treat)$p.value)
p_cell


f4.6 <- ggplot(cell_proportions %>% filter(cell_type %in% c("Pericyte","SMC","P. Fibro","Arterial")),
               aes(x = factor(cell_type,c("Pericyte","SMC","P. Fibro","Arterial")), y = percentage, fill=factor(treat,c("Control","AD")))) +
  geom_boxplot(outlier.shape=NA,position=position_dodge(width=0.9),lwd = 0.2) +
  geom_signif(y_position=c(0.34, 0.15, 0.05, 0.12), xmin=c(0.8, 1.8, 2.8, 3.8), xmax=c(1.2, 2.2, 3.2, 4.2),
              annotation=c("**", "0.052","*","**"), tip_length=0,textsize=1.5,size=0.2) +
  geom_point(position=position_jitterdodge(dodge.width=0.9),size=0.5)+
  labs(x = NULL, y = "Nuclei (%)",title = NULL) +
  scale_fill_manual(values = c("Control" = "#66c2a5", "AD" = "#f46d43"))+
  theme_classic(base_size = 6,base_line_size = 0.2)+
  guides(fill = guide_legend(title = NULL,override.aes = list(linetype = "blank", shape = NA)))+
  theme(legend.position = c(0.95, 0.95), legend.justification = c("right", "top"), legend.box.just = "right",
        axis.text.x = element_text(color="black"),
        axis.text.y = element_text(color="black"),
        
  ) 


# # Body vas seq
# ggplot(vas_seq_all, aes(x = Gene.Symbol, y = Cell.Type)) +
#   geom_point(aes(size = fraction*100, color = Expression..Scaled)) +
#   scale_size_continuous(range = c(0, 3),  breaks = seq(0, 60, by = 10)) +
#   scale_color_gradient(limits = c(0, 1),  breaks = seq(0, 1), low = "lightgrey", high = "black") +
#   theme_classic(base_size = 6,base_line_size = 0.2) +
#   labs(size = "Fraction of cells", color = "Mean expression") +
#   theme_classic(base_size = 6,base_line_size = 0.2)+
#   theme(axis.text.x = element_text(angle = 45,hjust = 1),
#         axis.title.x = element_blank(),
#         axis.title.y = element_blank(),  
#         legend.text = element_text(angle = 0, hjust = 0.5, size = 6))+
#   theme(legend.direction = "horizontal") + 
#   guides(size = guide_legend(label.position = "bottom", title.position = "top", title.hjust = 0.5, keywidth = 0.5, keyheight = 0.1),
#          color=guide_colorbar(label.position = "bottom", title.position = "top", title.hjust = 0.5, keywidth = 2.5, keyheight = 0.2)
#   )

##plot the distribution across organ in GTEx for give protein----------------------------------------------------

targeted_proteins <- volcano_data %>% arrange(desc(abs(effect_of_age)),desc(neg_log10_fdr_p_values)) %>%
  # filter(protein_category %in% c("Other","Artery")) %>%
  head(8) %>% dplyr::select(protein) %>% pull()

# targeted_proteins <-x$Artery 
targeted_proteins


#Get the ORGAN level GTEx data
{
  library(TissueEnrich)
  library(cowplot)
  library(SummarizedExperiment)
  library(RTCGAToolbox)
  library(Hmisc)
  library(DESeq2)
  gtex<- read.table("Data/Proteomic/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip = 2, header = TRUE,sep = "\t") %>% 
    distinct(Description,.keep_all = TRUE)#Remove duplicated genes (that encoding epitopes of proteins)
  gtex_name <- gtex %>% dplyr::select(1,2)
  row.names(gtex) <- gtex$Description
  
  #Mapping tissues to 21 organs
  gtex_df <-gtex %>% dplyr::select(-1,-2) %>% plyr::rename(c(
    "Adipose...Subcutaneous" = "Adipose_1",
    "Adipose...Visceral..Omentum." = "Adipose_2",
    "Adrenal.Gland" = "Adrenal",
    "Artery...Aorta" = "Artery_1",
    "Artery...Tibial" = "Artery_2",
    "Artery...Coronary" = "Artery_3",
    "Brain...Anterior.cingulate.cortex..BA24." = "Brain_1",
    "Brain...Hippocampus" = "Brain_2",
    "Brain...Cerebellum" = "Brain_3",
    "Brain...Hypothalamus" = "Brain_4",
    "Brain...Substantia.nigra" = "Brain_5",
    "Brain...Putamen..basal.ganglia." = "Brain_6",
    "Brain...Nucleus.accumbens..basal.ganglia." = "Brain_7",
    "Brain...Caudate..basal.ganglia." = "Brain_8",
    "Brain...Frontal.Cortex..BA9." = "Brain_9",
    "Brain...Amygdala" = "Brain_10",
    "Brain...Cerebellar.Hemisphere" = "Brain_11",
    "Brain...Spinal.cord..cervical.c.1." = "Brain_12",
    "Nerve...Tibial" = "Brain_13",
    "Brain...Cortex" = "Brain_14",
    "Esophagus...Mucosa" = "Esophagus_1",
    "Esophagus...Muscularis" = "Esophagus_2",
    "Esophagus...Gastroesophageal.Junction" = "Esophagus_3",
    "Ovary" = "Female_1",
    "Uterus" = "Female_2",
    "Vagina" = "Female_3",
    "Fallopian.Tube" = "Female_4",
    "Cervix...Endocervix" = "Female_5",
    "Cervix...Ectocervix" = "Female_6",
    "Heart...Left.Ventricle" = "Heart_1",
    "Heart...Atrial.Appendage" = "Heart_2",
    "Whole.Blood" = "Immune_1",
    "Spleen" = "Immune_2",
    "Colon...Transverse" = "Intestine_1",
    "Small.Intestine...Terminal.Ileum" = "Intestine_2",
    "Colon...Sigmoid" = "Intestine_3",
    "Kidney...Cortex" = "Kidney_1",
    "Kidney...Medulla" = "Kidney_2",
    "Testis" = "Male_1",
    "Prostate" = "Male_2",
    "Muscle...Skeletal" = "Muscle",
    "Minor.Salivary.Gland" = "Salivary",
    "Skin...Sun.Exposed..Lower.leg." = "Skin_1",
    "Skin...Not.Sun.Exposed..Suprapubic." = "Skin_2")
  ) %>% dplyr::mutate(Adipose = pmax(!!!dplyr::select(., contains("Adipose")), na.rm = TRUE),
                      Artery = pmax(!!!dplyr::select(., contains("Artery")), na.rm = TRUE),
                      Brain= pmax(!!!dplyr::select(., contains("Brain")), na.rm = TRUE),
                      Esophagus= pmax(!!!dplyr::select(., contains("Esophagus")), na.rm = TRUE),
                      Female= pmax(!!!dplyr::select(., contains("Female")), na.rm = TRUE),
                      Heart= pmax(!!!dplyr::select(., contains("Heart")), na.rm = TRUE),
                      Immune= pmax(!!!dplyr::select(., contains("Immune")), na.rm = TRUE),
                      Intestine= pmax(!!!dplyr::select(., contains("Intestine")), na.rm = TRUE),
                      Kidney= pmax(!!!dplyr::select(., contains("Kidney")), na.rm = TRUE),
                      Male= pmax(!!!dplyr::select(., contains("Male")), na.rm = TRUE),
                      Skin= pmax(!!!dplyr::select(., contains("Skin")), na.rm = TRUE),
  ) %>% dplyr::select(Adipose,Adrenal,Artery,Bladder,Brain,Esophagus,Female,Heart,Immune,Intestine,Kidney,Liver,Lung,Male,
                      Muscle,Pancreas,Pituitary,Salivary,Skin,Stomach,Thyroid)
  
}

plots <- lapply(targeted_proteins, function(targeted_protein) {
  # Filter data for the current protein
  protein_data <- gtex_df %>%
    filter(row.names(gtex_df) == targeted_protein) %>%
    pivot_longer(cols = everything(), names_to = "Tissue", values_to = "Expression")
  
  # Create ggplot for the current protein
  ggplot(protein_data, aes(x = reorder(Tissue, -Expression), y = Expression, fill = Tissue)) +
    geom_bar(stat = 'identity') +
    labs(x = '', y = paste0(targeted_protein, " TPM"),title =targeted_protein ) +
    theme_classic(base_size = 4,base_line_size = 0.2,base_rect_size = 6) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 5),
          axis.title = element_text(size = 5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 4),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = c(0, 0))
})

gridExtra::grid.arrange(grobs = plots, ncol = 3) 

#xFig 9_1 PAC in mental health and psychiatric disease--------------------------------------------------------------   

mental_nodisease <-cognition %>% filter(acd_interval<=0&pd_interval<=0&ms_interval<=0&dep_interval<=0&anx_interval<=0&scz_interval<=0&sd_interval<=0) %>% 
  mutate(phq4=scale(phq4) %>% as.numeric)


mental_result1 <- linear(mental_nodisease %>% filter(),a=c("phq4"),b=organ_pac,covar = c("age","sex","ethnic","edu","tdi","smok","pa","center")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="phq4"~"PHQ-4",
                         var_a=="dep_bas"~"Depressive symptoms",
                         var_a=="anx_bas"~"Anxiety symptoms",
                         T~var_a))



mental_result2 <- logistic(mental_nodisease %>% filter(),a=c("dep_bas","anx_bas"),b=organ_pac,covar = c("age","sex","ethnic","edu","tdi","smok","pa","center")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="phq4"~"PHQ-4",
                         var_a=="dep_bas"~"Depressive symptoms",
                         var_a=="anx_bas"~"Anxiety symptoms",
                         T~var_a))

mental_result <- rbind(mental_result1,mental_result2)

f3_1 <- ggplot(mental_result %>% 
                 filter(var_b %in% c("Brain","Intestine","Pancreas","Kidney")) %>% 
                 group_by(rev(var_a)) %>%  
                 arrange(desc(var_a),desc(beta)) %>% 
                 mutate(var_b = factor(var_b, levels = rev(unique(var_b)))), 
               aes(x = beta, y=var_a, group=var_b)) +
  geom_vline(xintercept = 0, color = "grey",linewidth=0.1)+
  geom_errorbarh(aes(xmin = beta-1.95*se, xmax = beta+1.95*se),height = 0,linewidth=0.2,position=position_dodge(0.7)) +
  geom_point(aes(color = var_b), size = 1,shape=15,position=position_dodge(0.7))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #scale_y_discrete(expand = c(0,1))+
  labs(x = "Beta", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
        #aspect.ratio=1
  )+
  scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))

f3_1


# %>% filter(dep_bas!=1&anx_bas!=1)
mental_result_fu1 <- linear(mental_nodisease, a=c("phq9_fu","gad7_fu"),b=organ_pac,covar = c("age","sex","ethnic","edu","tdi","smok","pa","center","dep_bas")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="phq9_fu"~"PHQ-9",
                         var_a=="gad7_fu"~"GAD-7",
                         var_a=="dep_bas"~"Depressive symptoms",
                         var_a=="anx_bas"~"Anxiety symptoms",
                         T~var_a))


mental_result_fu2 <- logistic(mental_nodisease,a=c("dep_fu","anx_fu","self_harm"),b=organ_pac,covar = c("age","sex","ethnic","edu","tdi","smok","pa","center")) %>% 
  mutate(var_b=gsub("_res","",var_b)) %>% 
  mutate(var_a=case_when(var_a=="phq9_fu"~"PHQ-9",
                         var_a=="dep_fu"~"Depressive symptoms",
                         var_a=="anx_fu"~"Anxiety symptoms",
                         var_a=="distress"~"Mental distress",
                         var_a=="self_harm"~"Self harm",
                         var_a=="Suicide"~"Suicide attempt",
                         T~var_a))

mental_result_fu <- rbind(mental_result_fu1,mental_result_fu2)

describe(mental_nodisease$dep_fu %>% as.factor())
#n=30187
f3_2 <- ggplot(mental_result_fu %>% filter(var_b %in% c("Brain","Intestine","Pancreas","Kidney")) %>%
                 mutate(var_a =factor(var_a,levels=rev(c("PHQ-9","GAD-7","Depressive symptoms","Anxiety symptoms","Mental distress","Self harm")))) %>% 
                 group_by(rev(var_a)) %>%
                 arrange(desc(var_a),desc(beta)) %>% 
                 mutate(var_b = factor(var_b, levels = rev(unique(var_b)))),
               aes(x = beta, y=var_a, group=var_b)) +
  geom_vline(xintercept = 0, color = "grey",linewidth=0.1)+
  geom_errorbarh(aes(xmin = beta-1.95*se, xmax = beta+1.95*se),height = 0,linewidth=0.2,position=position_dodge(0.7)) +
  geom_point(aes(color = var_b), size = 1,shape=15,position=position_dodge(0.7))+
  scale_x_continuous(labels = scales::number_format(accuracy = 0.01)) +
  #scale_y_discrete(expand = c(0,1))+
  labs(x = "Beta", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
        #aspect.ratio=1
  )+
  scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))
f3_2

#n=3954
# mental_stress <- mental_nodisease %>% filter(dep_bas==1|anx_bas==1)
# 

mental_stress <- ukb_olink %>% right_join(mental_nodisease %>% dplyr::select(eid,dep_bas,anx_bas),by="eid") %>% filter(dep_bas==1|anx_bas==1)


result_list_stress <- cox(data=mental_stress,exposures=c(organ_pac),
                          outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death"), #"acd","ad","vd","pd","ms","scz","dep","anx","sd","death"
                          covar=c("age","sex","ethnic","edu","tdi","center","smok","pa"))

result_disease_stress <- result_list_stress$result_list %>%   
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
    TRUE ~ as.character(Outcome)
  ),
  Exposure_label = sub("_res$", "", Exposure)
  ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                             "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","All-cause mortality")),
               
               Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                                                              "Lung","Muscle","Pancreas"))
  )


f3_3 <-
  ggplot(result_disease_stress %>% filter(Exposure_label %in% c("Brain","Intestine","Pancreas","Kidney")) %>%
           filter(Outcome %in% c("dep","anx","sd","sud","death")) %>% 
           mutate(Outcome_label=factor(Outcome_label,levels=rev(c("Depression","Anxiety","Sleep disorder","Substance use disorder","All-cause mortality")))) %>% 
           group_by(rev(Outcome_label)) %>%  
           arrange(desc(Outcome_label),desc(Beta)) %>% 
           mutate(Exposure_label = factor(Exposure_label, levels = rev(unique(Exposure_label)))),
         aes(x = HR, y=(Outcome_label),groups=Exposure_label)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color = Exposure_label),size = 1,shape=16,position=position_dodge(0.7))+
  # scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6))+
  scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))

f3_3


result_list <- cox(data=ukb_olink,exposures=c(organ_pac),
                   outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death"), #"acd","ad","vd","pd","ms","scz","dep","anx","sd","death"
                   covar=c("age","sex","ethnic","edu","tdi","center","smok","pa"))

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
    TRUE ~ as.character(Outcome)
  ),
  Exposure_label = sub("_res$", "", Exposure)
  ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                             "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","All-cause mortality")),
               
               Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                                                              "Lung","Muscle","Pancreas"))
  )

f3_4 <-
  ggplot(result_disease %>% filter(Exposure_label %in% c("Brain","Intestine","Pancreas","Kidney")) %>%
           filter(Outcome %in% c("dep","anx","sd","scz","sud","death")) %>% 
           group_by(rev(Outcome)) %>%
           mutate(Outcome_label=factor(Outcome_label,levels=rev(c("Depression","Anxiety","Schizophrenia","Sleep disorder","Substance use disorder","All-cause mortality")))) %>% 
           arrange(desc(Outcome_label),desc(Beta)) %>% 
           mutate(Exposure_label = factor(Exposure_label, levels = rev(unique(Exposure_label)))),
         aes(x = HR, y=(Outcome_label),groups=Exposure_label)) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  geom_errorbarh(aes(xmin = CI_Lower, xmax = CI_Upper), height = 0,linewidth=0.2,position=position_dodge(0.7))+
  geom_point(aes(color = Exposure_label),size = 1,shape=16,position=position_dodge(0.7))+
  #scale_x_continuous(limits = c(0.8, 3), breaks = c(1,2,3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title = NULL) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(            legend.position = "top",legend.title =element_blank(),legend.key.width = unit(0.1, "cm"),
                    axis.text.x = element_text(color="black",size = 6),
                    axis.text.y = element_text(color="black",size = 6))+
  scale_color_manual(values=c("Organismal"="black","Pancreas"="#008000","Intestine"="#808080","Brain"="#4B0082","Kidney"="#808000","Liver"="#D2B48C"))+
  guides(color = guide_legend(nrow = 1))

f3_4



#Mutivar cox regression
ukb_olink_inci <- result_list$data
ukb_olink_inci <- ukb_olink_inci %>% left_join(cognition %>% dplyr::select(eid,bd_prs,phq4),by="eid") %>% 
  # mutate(GDF15=as.numeric(scale(GDF15)),
  #        MAPT=as.numeric(scale(MAPT))
  # ) %>% 
  # mutate(ad_prs=as.numeric(scale(ad_prs)),
  #        pd_prs=as.numeric(scale(pd_prs)),
  #        ms_prs=as.numeric(scale(ms_prs))
  # ) %>% 
  mutate(brain_prs=bd_prs+Brain_res) %>% 
  mutate(brain_prs=as.numeric(scale(brain_prs)),
         brain_prs5=cut(brain_prs, breaks = c(-Inf, -1.5, -0.5, 0.5, 1.5, Inf), labels = c("-2", "-1", "0", "+1", "+2"))
  )

describe(ukb_olink_inci$brain_prs5)
#Outcome with PRS
# outcome="ad"
# outcome="pd"
# outcome="ms"
outcome="dep"
out2 <-  ukb_olink_inci %>% filter(is.na(!!sym(paste0(outcome, "_base"))))

exposure="Brain_res"
# exposure="Organismal_res"

# covar = c("age","sex","ethnic","edu","tdi","center","rt",paste0(outcome,"_prs"),"GDF15","MAPT")
# covar = c("age","sex","ethnic","edu","tdi","center","rt",paste0(outcome,"_prs"),"GDF15")
covar = c("age","sex","ethnic","edu","tdi","center",paste0("bd","_prs"),"phq4")
# covar = c("age","sex","ethnic","edu","tdi","center",paste0("brain","_prs"))



formula_input <-  as.formula(paste0(paste0( "survival::Surv(",outcome, "_fu,",outcome, "_inci", ") ~"), paste( c(exposure, covar), collapse="+")))

cox_summary <- survival::coxph(formula = formula_input, data =out2) %>% summary()
cox_summary$coefficients

results <- cox_summary$coefficients[c(exposure,"age",paste0("bd","_prs"),"phq4"), ] %>% as.data.frame() %>% 
  mutate(label=row.names(.)) %>% 
  mutate(hr=`exp(coef)`,
         lci=exp(coef-1.96*`se(coef)`),
         hci=exp(coef+1.96*`se(coef)`),
  ) %>% 
  mutate(label=case_when(label=="Brain_res"~"Brain age gap",
                         label=="age"~"Age",
                         label=="rt"~"Cognitive function",
                         label=="phq4"~"PHQ-4",
                         str_detect(label, "prs") ~ "PRS",
                         T~label
  ))

# 
# 
# results <- cox_summary$coefficients[c("age",paste0("brain","_prs")), ] %>% as.data.frame() %>% 
#   mutate(label=row.names(.)) %>% 
#   mutate(hr=`exp(coef)`,
#          lci=exp(coef-1.96*`se(coef)`),
#          hci=exp(coef+1.96*`se(coef)`),
#   ) %>% 
#   mutate(label=case_when(label=="Brain_res"~"Brain age gap",
#                          label=="age"~"Age",
#                          label=="rt"~"Cognitive function",
#                          label=="phq4"~"PHQ-4",
#                          str_detect(label, "prs") ~ "PRS",
#                          T~label
#   ))


f3_5 <-ggplot(results, aes(x = hr, y=reorder(label, coef))) +
  geom_point(size = 1,shape=15)+
  geom_errorbarh(aes(xmin = lci, xmax = hci), height = 0,linewidth=0.2) +
  geom_vline(xintercept = 1, color = "grey",linewidth=0.1) +
  #scale_x_continuous(limits = c(0.5, 3)) +
  scale_y_discrete()+
  labs(x = "HR (95% CI)", y = NULL, title =NULL ) +
  theme_classic(base_size = 6,base_line_size = 0.2)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        #plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
f3_5

library(survminer)
dep_prs_agegap <- survfit(Surv(dep_fu, dep_inci) ~ brain_prs5, out2) %>% 
  ggsurvplot(fun = "event",
             conf.int = F,
             censor = FALSE,
             size = 0.2,
             ggtheme = theme_classic(base_size=6, base_line_size = 0.2, base_family = "Arial"),
             legend.title="PRS+Brain age gap",
             palette = c("#0571b0","#92c5de","pink","#f4a582","#ca0020"),
             break.x.by = 5,
             xlab = "Time (years)",
             legend.labs=c("-2", "-1","0","+1","+2")
  )

dep_prs_agegap$plot <-  dep_prs_agegap$plot+ ggplot2::annotate("text", 
                                                               x = 4, y = 0.06, 
                                                               label = "HR = 1.2 (1.1-1.2)", size = 2)+
  theme_classic(base_line_size = 0.2,base_size = 6)+
  theme(legend.text = element_text(size = 5, color = "black"),
        legend.title =element_text(size = 5),
        legend.key.width = unit(0.15, "cm"),
        axis.text.x = element_text(color="black",size = 6),
        axis.text.y = element_text(color="black",size = 6),
        legend.position = "top"
  )+
  guides(colour = guide_legend(nrow = 1))
f3_6 <- dep_prs_agegap$plot

xf9 <- plot_grid(f3_1,f3_2,f3_3,f3_4,f3_5,f3_6,rel_widths = c(1,0.7,1,1,0.7,1),
                 nrow = 2, ncol= 3,align = "h",labels = "auto",label_size = 7,hjust =-1,vjust = 1.5) 


#xFig.9_2 Body and brain PAC in psychiatric disease--------------------------------------------------------------     
olink_subset2 <-ukb_olink %>% left_join(olink_full %>% dplyr::select(c(
  top_proteins_list$Brain,
  top_proteins_list$Organismal,
  top_proteins_list$Artery,
  top_proteins_list$Heart,
  
  top_proteins_list$Kidney,
  top_proteins_list$Pancreas,
  top_proteins_list$Intestine
  
  
),"eid"),by="eid") %>% 
  mutate(across(c( top_proteins_list$Brain,
                   top_proteins_list$Organismal,
                   top_proteins_list$Artery,
                   top_proteins_list$Heart,
                   
                   top_proteins_list$Kidney,
                   top_proteins_list$Pancreas,
                   top_proteins_list$Intestine),scale)) %>% mutate(across(c( top_proteins_list$Brain,
                                                                             top_proteins_list$Organismal,
                                                                             top_proteins_list$Artery,
                                                                             top_proteins_list$Heart,
                                                                             
                                                                             top_proteins_list$Kidney,
                                                                             top_proteins_list$Pancreas,
                                                                             top_proteins_list$Intestine),as.numeric)) %>% 
  mutate(age2=round(age, 0))


data_test2 <- olink_subset2 %>% filter(age2>=40)


prot_diease_by_organ2 <- 
  function(organ_list,disease){
    heatmap_list <- list()
    sig_prot <- data.frame()
    for (organ in organ_list){
      
      # ukb_olink_prot <- ukb_olink %>% left_join(olink_full %>% select(brain_20$gene,"eid"),by="eid") %>% 
      #   mutate(across(brain_20$gene,scale)) %>% mutate(across(brain_20$gene,as.numeric))
      # 
      result_list_organ<-cox_prot(data=olink_subset2,
                                  exposures=c(top_proteins_list[[organ]]),
                                  outcomes=c("acd","ad","vd","pd","ms","scz","dep","anx","sd","sud","death"),
                                  covar=c("age","sex","ethnic","tdi","center","pa",top_proteins_list[[organ]]))
      
      
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
          TRUE ~ as.character(Outcome)
        )
        ) %>% mutate(Outcome_label =factor(Outcome_label, levels=c("All-cause dementia", "Alzheimer's disease","Vascular dementia", "Parkinson's disease","Multiple sclerosis",
                                                                   "Schizophrenia","Depression", "Anxiety", "Sleep disorder","Substance use disorder","All-cause mortality")),
                     
                     # Exposure_label=factor(Exposure_label, levels=c("Organismal","Artery","Brain","Heart","Immune","Intestine","Kidney","Liver",
                     #                                                "Lung","Muscle","Pancreas"))
        ) %>% mutate(z_score_trans=ifelse(P_Value>0.05,0,z_score))
      
      mat_prot <- xtabs(z_score ~ Outcome_label + Exposure, result_disease_organ_pac) 
      mat_prot_p <- xtabs(P_Value ~ Outcome_label + Exposure, result_disease_organ_pac) 
      
      mat_prot_trans_pre <- xtabs(z_score_trans ~ Outcome_label + Exposure, result_disease_organ_pac) %>% as.data.frame() %>% 
        group_by(Outcome_label) %>% 
        mutate(Freq=ifelse(Freq<0,-Freq,Freq)) %>% 
        mutate(Freq_pct=Freq/sum(Freq))
      
      mat_prot_trans <- xtabs(Freq_pct ~ Outcome_label + Exposure, mat_prot_trans_pre)
      mat_prot_trans_bar <- mat_prot_trans[]
      
      mat_prot_sig <- mat_prot_trans_pre %>% ungroup() %>% 
        filter(Outcome_label %in% disease & Freq_pct > 0) %>%
        dplyr::select(Outcome_label,Exposure)  %>% mutate(organ=organ)
      
      sig_prot <- rbind(sig_prot, mat_prot_sig)
      # mat_prot <- mat_prot[order(rowMeans(mat_prot)), ]
      # mat_prot_p <- mat_prot_p[order(rowMeans(mat_prot)), ]
      
      col_fun = colorRamp2(c(-10,0,20), c("#3288bd","white","#d53e4f"))
      
      col_panel_prot=c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e",
                                "#e6ab02", "#a6761d", "#666666", "#a6cee3", "#1f78b4",
                                "#b2df8a", "#33a02c", "#fb9a99", "#e31a1c", "#fdbf6f",
                                "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")
                                
      
      col_panel_prot_subset <- col_panel_prot[1:length(top_proteins_list[[organ]])]
      
      # paste0(colnames(mat_prot),"'","=","'",col_panel_prot)
      
      heatmap_list[[organ]] <- grid.grabExpr(draw(Heatmap(mat_prot,
                                                          col = col_fun, na_col="white", show_column_dend = T, show_row_dend = T, border_gp = gpar(col = "black",lwd = 0.2),rect_gp = gpar(col = "black",lwd = 0.2),
                                                          cluster_rows = T,cluster_columns = T,
                                                          row_dend_width = unit(0.3, "cm"),
                                                          row_dend_gp = gpar(lwd = 0.2),
                                                          column_dend_height  = unit(0.3, "cm"),
                                                          column_dend_gp = gpar(lwd = 0.2),
                                                          width = unit(3, "cm"), height = unit(2, "cm"),
                                                          row_names_side = "right",row_names_gp = gpar(fontsize = 5),
                                                          column_names_rot=45,
                                                          column_names_side = "bottom",show_column_names = T,column_names_gp = gpar(fontsize = 5),
                                                          #row_split=factor(c(rep("Neurodegenerative",5),rep("Psychiatric",4),rep("All-cause mortality",1)),levels = c("Neurodegenerative", "Psychiatric","All-cause mortality")),
                                                          # row_title_gp = gpar(fontsize = 5,fontface = "bold"),
                                                          #column_split =factor(c(rep("Overall", each = 1),rep("Organ-specific aging", each = 10)),levels = c("Overall", "Organ-specific aging")),
                                                          column_title_gp = gpar(fontsize = 4),
                                                          column_title_side = "top",
                                                          column_title =paste0(organ," age gap"),
                                                          cell_fun = function(j, i, x, y, width, height, fill) {
                                                            
                                                            if((mat_prot_p)[i, j] < 0.001) {
                                                              grid.text("***", x, y,gp = gpar(fontsize = 4))
                                                            } else if((mat_prot_p)[i, j] < 0.01) {
                                                              grid.text("**", x, y,gp = gpar(fontsize = 4))
                                                            } else if((mat_prot_p)[i, j] < 0.05) {
                                                              grid.text("*", x, y,gp = gpar(fontsize = 4))
                                                            }
                                                            
                                                            
                                                          },
                                                          # right_annotation = rowAnnotation(" "=row.names(mat_prot),
                                                          #                                  show_legend = FALSE,
                                                          #                                  width= unit(1, "mm")
                                                          # ),
                                                          heatmap_legend_param = list(
                                                            #at = c(-2, 0, 2),
                                                            #labels = c("0", "10%", "20%", "30%", "40%"),
                                                            title =  "z-score",
                                                            legend_height = unit(1, "cm"),
                                                            grid_width = unit(0.2, "mm"),
                                                            #legend_width =unit(1, "mm"),
                                                            #title_position = "lefttop-rot",
                                                            labels_gp = gpar(col = "black", fontsize = 5),
                                                            title_gp=gpar(col = "black", fontsize = 5)
                                                            
                                                          ),
                                                          
                                                          right_annotation=rowAnnotation(
                                                            " "= anno_barplot(mat_prot_trans_bar, border=F, ylim=c(0,1),
                                                                              gp = gpar(fill =col_panel_prot,
                                                                                        col = col_panel_prot,
                                                                                        lwd=0.2),bar_width = 0.8,
                                                                              axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(0,0.5,1),labels = c("0","0.5", "1"),labels_rot=0),
                                                            ),
                                                            
                                                            # "Age gap for each disease" = anno_points(mat_disease_gap, ylim = c(-0.3,0.5),border = F, axis_param = list(gp=gpar(fontsize = 5,col="black",lwd=0.2),at=c(-0.25,0,0.25,0.5),labels = c("-0.25","0","0.25", "0.5"),labels_rot=0),
                                                            #                                          pch = 16, size = unit(0.15, "cm"), gp = gpar(col = color_panel,lwd = 0.2)),
                                                            # 
                                                            width = unit(3, "cm"), annotation_name_gp= gpar(fontsize = 6)
                                                            # annotation_width=unit(c(1, 2), c("null", "null")),
                                                            # gap = unit(3, "mm")
                                                            
                                                          ),
                                                          bottom_annotation =columnAnnotation(" "=as.character(colnames(mat_prot)),
                                                                                              col= list(" "=setNames(col_panel_prot_subset, levels(mat_prot_trans_pre$Exposure))),
                                                                                              simple_anno_size = unit(2, "mm"),
                                                                                              show_legend = FALSE)
      )))
    }
    return(list(heatmap_list=heatmap_list,sig_prot=sig_prot))
  }

#Selected the aging-related prot [brain, organ, artery, heart] that also associated with AD
#eFig.8-----
heatmap2 <- prot_diease_by_organ2(c(organ_list="Brain","Organismal","Artery","Heart","Kidney","Pancreas","Intestine"),
                                  disease = c("Depression","Anxiety","All-cause dementia", "Alzheimer's disease","Vascular dementia"))
#Individual prot and disease by organ [For dep/anx]
plot_grid(plotlist = heatmap2$heatmap_list,labels = "auto",ncol = 2, nrow=4, label_size = 6)

#Prioritized proteins related to AD
heatmap2$sig_prot %>% filter(organ %in% c("Brain","Kidney","Pancreas","Intestine")) %>% 
  filter(Outcome_label %in% c("Depression","Anxiety")) %>% 
  distinct(Exposure) %>% pull()

#LOWESS Plot for Change with Age 
library(ggrepel) 

generate_loess_geoms2 <- function(data,var_list) {
  plot <- ggplot()+geom_hline(yintercept = 0, linetype = "dashed", color = "gray", lwd=0.2)+
    theme_classic(base_line_size = 0.2,base_size = 6)+labs(title = NULL, x = "Age", y = "Z-scored protein level")+
    ylim(-2,3)+xlim(39,75)
  
  # var_list <- data %>% dplyr::select(-eid, -age, -age2,-sex) %>% colnames()
  # var_list <-   heatmap$sig_prot %>% distinct(Exposure) %>% pull()
  # var_list <- ad_sig_prot$prot
  var_list <- var_list
  
  
  top_protein_endpoints <- data.frame(var = character(), x = numeric(), y = numeric(), color= character(), stringsAsFactors = FALSE)
  
  for (var in var_list) {
    
    #start_value <- colMeans(data[data$age2 == 40, var], na.rm = TRUE)
    
    mean_40 <-  data %>% filter(age2 == 40) %>% dplyr::select(var) %>%
      colMeans(na.rm = TRUE)
    
    sd <- sd(data[[var]], na.rm = TRUE)
    
    data[[var]] <- ifelse(data$age2 == 40, 0, (data[[var]] - mean_40) / sd)
    
    #Brain
    # values=c("black","#4B0082","#A52A2A","#FF7F50")
    color <-if (var %in% top_proteins_list$Brain) "#4B0082"
      else if (var %in% top_proteins_list$Kidney) "#808000"
        else if (var %in% top_proteins_list$Intestine) "#808080"
          else if (var %in% top_proteins_list$Pancreas) "#008000"
            else if (var %in% top_proteins_list$Liver) "#D2B48C"
              else "black"
                
    
    loess_fit <- loess(as.formula(paste0(var, "~age2")), data = data)
    predicted_values <- predict(loess_fit, newdata = data.frame(age2 = data$age2))
    offset <- predicted_values[data$age2 == 40][1]
    adjusted_predicted_values <- predicted_values - offset
    adjusted_data <- data.frame(age2 = data$age2, adjusted_values = adjusted_predicted_values)
    
    endpoint <- data.frame(var = var, x = 71, y = adjusted_predicted_values[data$age2 == 71], color = color)   
    
    top_protein_endpoints <- rbind(top_protein_endpoints, endpoint[1,])
    
    
    plot <- plot +
      geom_line(aes(x = age2, y = adjusted_values), color = color, size = 0.3, data = adjusted_data)
    
  }
  
  top_protein_endpoints <- top_protein_endpoints %>%
    arrange(desc(abs(y))) %>%
    head(10)
  
  # Add labels for the top 10 proteins
  plot <- plot + 
    geom_text_repel(data = top_protein_endpoints, aes(x = x, y = y, label = var,color = color), direction = "both",
                    size = 1.5, segment.size = 0.2, nudge_x=2, max.overlaps = 20,box.padding=0.3, segment.color = "black")+scale_color_identity() 
  
  return(plot)
}

traj_plot <- generate_loess_geoms2(data=data_test2,
                                   var_list=heatmap2$sig_prot %>% filter(organ %in% c("Brain","Kidney","Pancreas","Intestine")) %>% 
                                     filter(Outcome_label %in% c("Depression","Anxiety")) %>% 
                                     distinct(Exposure) %>% pull()
)
f5.1 <- traj_plot


# Plot for Change with Age 
{
  data <- olink_subset2
  
  # scale_variable <- function(x) {
  #   (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
  # }
  # data <- data %>%
  #   mutate(across(-c(eid,age,sex,age2), scale_variable))
  
  names(data)
  
  lm_results <- lapply(names(data %>% dplyr::select(heatmap2$sig_prot %>% filter(organ %in% c("Brain","Kidney","Pancreas","Intestine")) %>% 
                                                      filter(Outcome_label %in% c("Depression","Anxiety")) %>% dplyr::select(Exposure) %>% pull() )), function(var) {
                                                        lm_model <- lm(data[[var]] ~ age+sex, data = data)
                                                        p_value <- summary(lm_model)$coefficients["age", "Pr(>|t|)"]
                                                        return(p_value)
                                                      })
  
  p_values <- unlist(lm_results)
  fdr_corrected_p_values <- p.adjust(p_values, method = "fdr")
  neg_log10_fdr_p_values <- -log10(fdr_corrected_p_values)
  
  effect_of_age <- sapply(names(data %>% dplyr::select(heatmap2$sig_prot %>% filter(organ %in% c("Brain","Kidney","Pancreas","Intestine")) %>% 
                                                         filter(Outcome_label %in% c("Depression","Anxiety")) %>% dplyr::select(Exposure) %>% pull())), function(var) {
                                                           lm_model <- lm(data[[var]] ~ age+sex, data = data)
                                                           effect <- coef(lm_model)["age"]
                                                           return(effect)
                                                         })
  volcano_data <- data.frame(effect_of_age = effect_of_age, neg_log10_fdr_p_values = neg_log10_fdr_p_values)
  volcano_data$neg_log10_fdr_p_values[is.infinite(volcano_data$neg_log10_fdr_p_values)] <- 291
  volcano_data$protein <- sub("\\.age$", "", row.names(volcano_data)) 
  # organ_list="Brain","Organismal","Kidney","Liver","Pancreas","Intestine"
  protein_colors <- c(
    Brain = "#4B0082",
    Kidney = "#808000",
    # Liver = "#D2B48C",
    Pancreas="#008000",
    Intestine="#808080"
      # Other = "black"  # Default color for proteins not falling into the specified categories
  )
  
  # Assign colors based on protein categories
  volcano_data <- volcano_data %>%
    mutate(protein_category = case_when(
      protein %in% top_proteins_list$Brain ~ "Brain",
      protein %in% top_proteins_list$Kidney ~ "Kidney",
      protein %in% top_proteins_list$Liver ~ "Liver",
      protein %in% top_proteins_list$Pancreas ~ "Pancreas",
      protein %in% top_proteins_list$Intestine ~ "Intestine",
      TRUE ~ "Other"
    ))
  
  f5.2 <- ggplot(volcano_data, aes(x = effect_of_age, y = neg_log10_fdr_p_values)) +
    geom_point(aes(size = neg_log10_fdr_p_values,color=protein_category),alpha = 0.7) + 
    scale_color_manual(values = protein_colors) +
    scale_size_continuous(range = c(0.2, 2.5))+
    geom_text_repel(aes(label = protein), size = 1, max.overlaps = 20, direction= "both",segment.size = 0.2) +
    geom_hline(yintercept = 0, linetype = "dotted",size=0.2)+
    labs(title = NULL,
         x = "Effect of Age",
         y = "-log10(FDR-corrected p-value)")+theme_classic(base_size = 5,base_line_size = 0.2)+ylim(0,320)+
    theme(text = element_text(size = 5),
          line = element_line(linewidth = 0.2),
          legend.position = "none"
    )
  
  
}


targeted_proteins <- volcano_data %>% arrange(desc(abs(effect_of_age)),desc(neg_log10_fdr_p_values)) %>%
  # filter(protein_category %in% c("Other","Artery")) %>%
  head(12) %>% dplyr::select(protein) %>% pull()

# targeted_proteins <-x$Artery 
targeted_proteins_psy <- targeted_proteins
targeted_proteins_psy

#Proteins interaction based on STING database
library(STRINGdb)
library(igraph)
library(ggplotify)
#Conduct interaction analyses on: https://string-db.org/
raw_data <- fread(input = "Data/string_interactions_dep.tsv")
# Check the structure of the raw_data data frame
str(raw_data)

# Create an igraph object from the data frame
g <- graph_from_data_frame(d = raw_data[, c("node1", "node2")], directed = FALSE)
E(g)$width <- raw_data$combined_score *1.1

# protein_colors <- c(
#   Brain = "#4B0082",
#   Kidney = "#808000",
#   Liver = "#D2B48C",
#   Pancreas="#008000",
#   Intestine="#808080",
#   Other = "black"  # Default color for proteins not falling into the specified categories
# )


node_color_mapping <- list(
  
  "NEFL" = "#4B0082",
  "RBFOX3" = "#4B0082",
  "LGALS4" = "#808080",
  "MOG" = "#4B0082",
  "SPINK1" = "#008000",
  "SPINK4" = "#808080",
  "FCAMR"="#808000",
  "CNTN2"="#4B0082",
  "CA14"="#4B0082",
  "SPOCK1" ="#4B0082"
    
)

# Create a vector of colors for nodes based on the mapping
node_colors <- sapply(V(g)$name, function(node) {
  if (node %in% names(node_color_mapping)) {
    return(node_color_mapping[[node]])
  } else {
    return("gray")  # Assign gray color if node is not found in the mapping
  }
})


vertex_label_colors <- rep("black", length(V(g)))  # Set all labels to dark red

# Plot the network

par(mar = c(0, 0, 0, 0) + 0.1)
plot(g,
     vertex.label = V(g)$name,  # Remove automatic labeling
     vertex.color = node_colors,
     vertex.label.color = vertex_label_colors,
     vertex.label.cex = 0.5, 
     vertex.label.dist = 1.5,
     vertex.size = 5,
     edge.width = E(g)$width,
     vertex.label.family = "sans"
)

f5.3 <- base2grob(~plot(g,vertex.label = V(g)$name,
                        vertex.color = node_colors,
                        vertex.label.color = vertex_label_colors,
                        vertex.label.cex = 0.4, 
                        vertex.label.dist = 1.5,
                        vertex.size = 3,
                        edge.width = E(g)$width,
                        vertex.label.family = "sans"))




#Enrichment of Proteins interaction

library(clusterProfiler)
library(org.Hs.eg.db)

go_enrichment <- enrichGO(gene = 
                            c(targeted_proteins,"SPINK2","SPINK13","SPINK6","REG4","LGALS4","LBH"),
                          # heatmap2$sig_prot %>% distinct(Exposure) %>% pull() %>% unique(),
                          OrgDb = org.Hs.eg.db,
                          keyType = "SYMBOL",   
                          ont = "BP",
                          pvalueCutoff = 0.05
)


head(go_enrichment,n=20)

f5.4 <- dotplot(go_enrichment, showCategory =6)+
  guides(size = FALSE,color = guide_colorbar(barwidth = 0.5,barheight = 2))+
  geom_point(aes(size =1.5))+
  scale_size_continuous(range = c(0.5, 2))+ 
  scale_color_gradient(high = "#d0d1e6", low ="#023858")+
  labs(color = "q value",x="Gene ratio")+
  theme_bw(base_size = 6,base_line_size = 0.1)+
  theme(panel.grid = element_blank(),  
        axis.text.x = element_text(size = 6))+
  scale_y_discrete(labels = function(x) (str_wrap(str_to_sentence(x), width = 30)))+
  scale_x_continuous(limits = c(0.15,0.45), breaks = c(0.2, 0.4))+
  theme(legend.text = element_text(size = 5))


#Get the ORGAN level GTEx data
{
  library(TissueEnrich)
  library(cowplot)
  library(SummarizedExperiment)
  library(RTCGAToolbox)
  library(Hmisc)
  library(DESeq2)
  gtex<- read.table("Data/Proteomic/GTEx/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct", skip = 2, header = TRUE,sep = "\t") %>% 
    distinct(Description,.keep_all = TRUE)#Remove duplicated genes (that encoding epitopes of proteins)
  gtex_name <- gtex %>% dplyr::select(1,2)
  row.names(gtex) <- gtex$Description
  
  #Mapping tissues to 21 organs
  gtex_df <-gtex %>% dplyr::select(-1,-2) %>% plyr::rename(c(
    "Adipose...Subcutaneous" = "Adipose_1",
    "Adipose...Visceral..Omentum." = "Adipose_2",
    "Adrenal.Gland" = "Adrenal",
    "Artery...Aorta" = "Artery_1",
    "Artery...Tibial" = "Artery_2",
    "Artery...Coronary" = "Artery_3",
    "Brain...Anterior.cingulate.cortex..BA24." = "Brain_1",
    "Brain...Hippocampus" = "Brain_2",
    "Brain...Cerebellum" = "Brain_3",
    "Brain...Hypothalamus" = "Brain_4",
    "Brain...Substantia.nigra" = "Brain_5",
    "Brain...Putamen..basal.ganglia." = "Brain_6",
    "Brain...Nucleus.accumbens..basal.ganglia." = "Brain_7",
    "Brain...Caudate..basal.ganglia." = "Brain_8",
    "Brain...Frontal.Cortex..BA9." = "Brain_9",
    "Brain...Amygdala" = "Brain_10",
    "Brain...Cerebellar.Hemisphere" = "Brain_11",
    "Brain...Spinal.cord..cervical.c.1." = "Brain_12",
    "Nerve...Tibial" = "Brain_13",
    "Brain...Cortex" = "Brain_14",
    "Esophagus...Mucosa" = "Esophagus_1",
    "Esophagus...Muscularis" = "Esophagus_2",
    "Esophagus...Gastroesophageal.Junction" = "Esophagus_3",
    "Ovary" = "Female_1",
    "Uterus" = "Female_2",
    "Vagina" = "Female_3",
    "Fallopian.Tube" = "Female_4",
    "Cervix...Endocervix" = "Female_5",
    "Cervix...Ectocervix" = "Female_6",
    "Heart...Left.Ventricle" = "Heart_1",
    "Heart...Atrial.Appendage" = "Heart_2",
    "Whole.Blood" = "Immune_1",
    "Spleen" = "Immune_2",
    "Colon...Transverse" = "Intestine_1",
    "Small.Intestine...Terminal.Ileum" = "Intestine_2",
    "Colon...Sigmoid" = "Intestine_3",
    "Kidney...Cortex" = "Kidney_1",
    "Kidney...Medulla" = "Kidney_2",
    "Testis" = "Male_1",
    "Prostate" = "Male_2",
    "Muscle...Skeletal" = "Muscle",
    "Minor.Salivary.Gland" = "Salivary",
    "Skin...Sun.Exposed..Lower.leg." = "Skin_1",
    "Skin...Not.Sun.Exposed..Suprapubic." = "Skin_2")
  ) %>% dplyr::mutate(Adipose = pmax(!!!dplyr::select(., contains("Adipose")), na.rm = TRUE),
                      Artery = pmax(!!!dplyr::select(., contains("Artery")), na.rm = TRUE),
                      Brain= pmax(!!!dplyr::select(., contains("Brain")), na.rm = TRUE),
                      Esophagus= pmax(!!!dplyr::select(., contains("Esophagus")), na.rm = TRUE),
                      Female= pmax(!!!dplyr::select(., contains("Female")), na.rm = TRUE),
                      Heart= pmax(!!!dplyr::select(., contains("Heart")), na.rm = TRUE),
                      Immune= pmax(!!!dplyr::select(., contains("Immune")), na.rm = TRUE),
                      Intestine= pmax(!!!dplyr::select(., contains("Intestine")), na.rm = TRUE),
                      Kidney= pmax(!!!dplyr::select(., contains("Kidney")), na.rm = TRUE),
                      Male= pmax(!!!dplyr::select(., contains("Male")), na.rm = TRUE),
                      Skin= pmax(!!!dplyr::select(., contains("Skin")), na.rm = TRUE),
  ) %>% dplyr::select(Adipose,Adrenal,Artery,Bladder,Brain,Esophagus,Female,Heart,Immune,Intestine,Kidney,Liver,Lung,Male,
                      Muscle,Pancreas,Pituitary,Salivary,Skin,Stomach,Thyroid)
  
  }
#Featured prot for psychiatric disease
targeted_proteins_psy <- c("NEFL","RBFOX3","LGALS4","MOG",
                           "SPINK1","SPINK4","FCAMR","CNTN2",
                           "CA14","SPOCK1","MEP1A","FLT3"
)

targeted_proteins_neuro <- c("ELN","LTBP2","NEFL","GDF15",
                             "GFAP","TNFRSF11B","IGDCC4","COMP",
                             "FBLN2","NPPB","CA14","C1QL2"
)

plots_psy <- lapply(targeted_proteins_psy, function(targeted_protein) {
  # Filter data for the current protein
  protein_data <- gtex_df %>%
    filter(row.names(gtex_df) == targeted_protein) %>%
    pivot_longer(cols = everything(), names_to = "Tissue", values_to = "Expression")
  
  # Create ggplot for the current protein
  ggplot(protein_data, aes(x = reorder(Tissue, -Expression), y = Expression, fill = Tissue)) +
    geom_bar(stat = 'identity') +
    labs(x = '', y = paste0(targeted_protein, " TPM"),title =targeted_protein ) +
    theme_classic(base_size = 4,base_line_size = 0.2,base_rect_size = 6) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 5),
          axis.title = element_text(size = 5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 4),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = c(0, 0))
})
gridExtra::grid.arrange(grobs = plots_psy, ncol = 3) 

plots_neuro <- lapply(targeted_proteins_neuro, function(targeted_protein) {
  # Filter data for the current protein
  protein_data <- gtex_df %>%
    filter(row.names(gtex_df) == targeted_protein) %>%
    pivot_longer(cols = everything(), names_to = "Tissue", values_to = "Expression")
  
  # Create ggplot for the current protein
  ggplot(protein_data, aes(x = reorder(Tissue, -Expression), y = Expression, fill = Tissue)) +
    geom_bar(stat = 'identity') +
    labs(x = '', y = paste0(targeted_protein, " TPM"),title =targeted_protein ) +
    theme_classic(base_size = 4,base_line_size = 0.2,base_rect_size = 6) +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, size = 5),
          axis.title = element_text(size = 5),
          axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),
          axis.title.y = element_text(size = 4),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    scale_y_continuous(expand = c(0, 0))
})

gridExtra::grid.arrange(grobs = plots_neuro, ncol = 4) %>% 
  plot_grid(.,plot_grid(plotlist = plots_psy,ncol =4),nrow=2,labels = c("h","i"),label_size = 6)


#GWAS on PAC (Brain and organismal)------
#Brain
brain_gwas <- fread(file = "Data/GWAS/brain/assoc.linear.plink.merged.txt")
brain_gwas_freq <- fread(file = "Data/GWAS/brain/freq.merged.txt") %>% dplyr::select(ID,A1) %>% 
  dplyr::rename(freq=A1)

brain_gwas <- brain_gwas %>% left_join(brain_gwas_freq,by="ID") %>% filter(ERRCODE!="CONST_OMITTED_ALLELE")
write.table(brain_gwas, quote = F,file = "Data/GWAS/brain/brain.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

brain_gwas2 <- brain_gwas %>% dplyr::select(CHROM,POS,ID,REF,ALT,freq,BETA,SE,P,N) %>% 
  dplyr::rename(FREQ=freq,
                CHR=CHROM,
                BP=POS,
                A1=REF,
                A2=ALT,
                SNP=ID
  )
write.table(brain_gwas2, quote = F,file = "Data/GWAS/brain/brain2.txt", sep = "\t", row.names = FALSE, col.names = TRUE)

#Organismal
organ_gwas <- fread(file = "Data/GWAS/organismal/assoc.linear.plink.merged.txt")
organ_gwas_freq <- fread(file = "Data/GWAS/organismal/freq.merged.txt") %>% dplyr::select(ID,A1) %>% 
  dplyr::rename(freq=A1)

organ_gwas <- organ_gwas %>% left_join(organ_gwas_freq,by="ID") %>% filter(ERRCODE!="CONST_OMITTED_ALLELE")
write.table(organ_gwas, quote = F,file = "Data/GWAS/organismal/organismal.txt", sep = "\t", row.names = FALSE, col.names = TRUE)  



library(qqman)
library(dplyr)
library(ggrepel)

#FUMA functional annotations
#Brain
brain_gene <- fread(file="Data/GWAS/brain/FUMA_brain/genes.txt")
brain_snp <- fread(file="Data/GWAS/brain/FUMA_brain/leadSNPs.txt")
brain_gene_set <- fread(file = "Data/GWAS/brain/FUMA_job599346/magma.genes.out")
#Organismal
organ_gene <- fread(file="Data/GWAS/organismal/FUMA_organ/genes.txt")
organ_gene_set <- fread(file = "Data/GWAS/organismal/FUMA_organ/magma.genes.out")


#Brain
# Prepare the dataset
don <- brain_gene_set %>% mutate(CHR=ifelse(CHR=="X",23,CHR) %>% as.numeric()) %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(START)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(brain_gene_set %>% mutate(CHR=ifelse(CHR=="X",23,CHR) %>% as.numeric()), ., by=c("CHR"="CHR")) %>%
  arrange(CHR, START) %>%
  mutate( BPcum=START+tot)
# # Add highlight and annotation information
# mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
# mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 
# 
# Prepare X axis
axisdf <- don %>% group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
brain_geneset_plot <- ggplot(don, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#c6dbef", "#08519c"), 22 )) +
  scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # 
  # Add label using ggrepel to avoid overlapping
  geom_text_repel( data=subset(don,P<=2.520e-6), aes(label=SYMBOL), size=2) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())+xlab("Chromosome")



don2 <- brain_gwas %>% mutate(CHROM=ifelse(CHROM=="X",23,CHROM) %>% as.numeric()) %>% 
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(brain_gwas %>% mutate(CHROM=ifelse(CHROM=="X",23,CHROM) %>% as.numeric()), ., by=c("CHROM"="CHROM")) %>%
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot)

axisdf2 <- don2 %>% group_by(CHROM) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# brain_gene <- brain_gene %>% dplyr::rename(POS=start) 

brain_gwas_plot <- ggplot(don2, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#c6dbef", "#08519c"), 22 )) +
  scale_x_continuous(label = axisdf2$CHROM, breaks= axisdf2$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # 
  # Add label using ggrepel to avoid overlapping
  # geom_text_repel(data=subset(brain_gene,minGwasP<=5e-8), aes(x = BPcum,y=-log10(minGwasP), label=symbol), size=2,max.overlaps=1000) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())



#Organisaml
# Prepare the dataset
don_organ <- organ_gene_set %>% mutate(CHR=ifelse(CHR=="X",23,CHR) %>% as.numeric()) %>% 
  group_by(CHR) %>% 
  summarise(chr_len=max(START)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(brain_gene_set %>% mutate(CHR=ifelse(CHR=="X",23,CHR) %>% as.numeric()), ., by=c("CHR"="CHR")) %>%
  arrange(CHR, START) %>%
  mutate( BPcum=START+tot)
# # Add highlight and annotation information
# mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
# mutate( is_annotate=ifelse(-log10(P)>4, "yes", "no")) 
# 
# Prepare X axis
axisdf3 <- don_organ %>% group_by(CHR) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# Make the plot
organ_geneset_plot <- ggplot(don_organ, aes(x=BPcum, y=-log10(P))) +
  geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#c6dbef", "#08519c"), 22 )) +
  scale_x_continuous( label = axisdf3$CHR, breaks= axisdf3$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # 
  # Add label using ggrepel to avoid overlapping
  geom_text_repel( data=subset(don_organ,P<=2.520e-6), aes(label=SYMBOL), size=2) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  theme(legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

don2_organ <- organ_gwas %>% mutate(CHROM=ifelse(CHROM=="X",23,CHROM) %>% as.numeric()) %>% 
  group_by(CHROM) %>% 
  summarise(chr_len=max(POS)) %>% 
  mutate(tot=cumsum(as.numeric(chr_len))-chr_len) %>%
  select(-chr_len) %>%
  left_join(brain_gwas %>% mutate(CHROM=ifelse(CHROM=="X",23,CHROM) %>% as.numeric()), ., by=c("CHROM"="CHROM")) %>%
  arrange(CHROM, POS) %>%
  mutate( BPcum=POS+tot)

axisdf4 <- don2_organ %>% group_by(CHROM) %>% dplyr::summarize(center=( max(BPcum) + min(BPcum) ) / 2 )

# brain_gene <- brain_gene %>% dplyr::rename(POS=start) 

organ_gwas_plot <- ggplot(don2_organ, aes(x=BPcum, y=-log10(P))) +
  geom_point(aes(color=as.factor(CHROM)), alpha=0.8, size=1.3) +
  scale_color_manual(values = rep(c("#c6dbef", "#08519c"), 22 )) +
  scale_x_continuous(label = axisdf4$CHROM, breaks= axisdf4$center ) +
  scale_y_continuous(expand = c(0, 0) ) + 
  # # Add highlighted points
  # geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  # 
  # Add label using ggrepel to avoid overlapping
  # geom_text_repel(data=subset(organ_gene,minGwasP<=5e-8), aes(x = BPcum,y=-log10(minGwasP), label=symbol), size=2,max.overlaps=1000) +
  theme_classic(base_size = 6,base_line_size = 0.2) +
  theme( 
    legend.position="none",
    panel.border = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank())



plot_grid(organ_gwas_plot,organ_geneset_plot,brain_gwas_plot,brain_geneset_plot,nrow = 4,labels = "auto", label_size = 6)

