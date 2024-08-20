library(data.table)
library(tidyfst)
library(dplyr)
library(plyr)
library(survival)
library(survminer)
library(stringr)
library(rms)
library(vroom)
library(dietaryindex)
library(writexl)
library(lavaan)
# function annotation -----------------------------------------------------

my_merge <- function(df1, df2){                                # Create own merging function
  merge(df1, df2, by = "eid")
}

EER_function <- function(sex,Age,PA,weight,height){
  if(sex=="2"){
    EER=354.1-6.91*Age+PA*(9.36*weight+7.26*height)
  } else {
    EER=662-9.53*Age+PA*(15.91*weight+5.396*height)
  }
  return(EER)
}

Quantile<-function(x){
  ifelse(x>quantile(x,.80),"Q1",ifelse(x>quantile(x,.60),"Q2",ifelse(x>quantile(x,.40),"Q3",ifelse(x>quantile(x,.20),"Q4","Q5"))))
}

Third.quartile<-function(x){
  ifelse(x>quantile(x,0.8),"Q3",ifelse(x>quantile(x,0.2),"Q2","Q1"))
}


read_and_select <- function(file_name) {
  read_csv(file_name) %>%
    select(eid, everything())
}
# DII score calculation ---------------------------------------------------


## primary diet data import ------------------------------------------------
Field_description <- readxl::read_xlsx("F:/paper/paper 6/UKB_dietary/primary data/Dietary Field ID.xlsx",sheet = "Sheet2")
diet_primary_data <- fread("F:/paper/paper 6/UKB_dietary/primary data/UKB_dietary_weight.csv")

food_Field_ID_split <- split(Field_description,Field_description$`Field ID`)
for (i in 1:length(food_Field_ID_split)) {
  food_Field_ID_split[[i]] <- diet_primary_data %>% 
    select(1,contains(c(as.character(food_Field_ID_split[[i]][["Field ID"]]))))
}

for (i in 1:length(food_Field_ID_split)) {
  N = apply(food_Field_ID_split[[i]], 2, is.na)
  food_Field_ID_split[[i]] <- food_Field_ID_split[[i]] %>% mutate(
    mean = rowMeans(food_Field_ID_split[[i]][,-1], na.rm=TRUE),
    count = length(colnames(food_Field_ID_split[[i]]))-1-apply(N, 1, sum)) %>% 
    filter(count>=1)
  setnames(food_Field_ID_split[[i]],
           c("count","mean"),
           c(paste0("count_sum_",names(food_Field_ID_split)[i],"_",Field_description$Description[Field_description$`Field ID`==names(food_Field_ID_split)[i]]),
             paste0("mean_sum_",names(food_Field_ID_split)[i],"_",Field_description$Description[Field_description$`Field ID`==names(food_Field_ID_split)[i]])))
}

food_Field_ID_summry <-Reduce(my_merge, food_Field_ID_split) %>% select(1,contains("mean_sum"))
####data prepare###
tea.coffee <- food_Field_ID_summry %>% select(1,contains(c("26141","26081")))
tea.coffee$caffeine <- tea.coffee$`mean_sum_26081_Coffee, caffeinated`/100*0.08
DII_data <- vroom("E:/paper/paper 6/UKB_dietary/primary data/UKB_estimated_nutrients63_data1.csv.zip")
DII_data <- merge(DII_data,tea.coffee,by="eid")
DII_data$Energy <- DII_data$Energy/4.184

## exclude unreasonable data -----------------------------------------------


physical <- c("22038","22039","22040","22037","904","914","884","894","864","874")
covariate <- fread("E:/paper/paper 6/UKB_dietary/primary data/UKB_covariate13.csv")
MET_score_primary_data <- fread("E:/paper/paper 6/UKB_dietary/primary data/physical-All-ins0.csv")
MET_score <- MET_score_primary_data %>% select(1,contains(physical)) %>% filter(eid %in% DII_score$eid)

MET_score_EER <- MET_score %>% mutate(
  time_22037 = MET_score$`864`*MET_score$`874`,
  time_22038 = MET_score$`884`*MET_score$`894`,
  time_22039 = MET_score$`904`*MET_score$`914`,
  `22039` = time_22039*8,
  `22037` = time_22037*3.3,
  `22038` = time_22038*4
)
MET_score_EER <- MET_score_EER %>% mutate(
  calculate_22040 = rowSums(MET_score_EER[,c(2:3,5)],na.rm = T),
  BMI_group+PAL = ifelse(calculate_22040<600,"1",ifelse(calculate_22040>=3000,"3","2"))
) %>% filter(!(is.na(`22040`)&calculate_22040==0&!(`904`==0&`884`==0&`864`==0)))

EER <- merge(covariate,MET_score_EER[,c(1,16)],by="eid")
EER$PA <- ifelse(EER$sex==2&EER$BMI_group+PAL=="1",1.12,ifelse(EER$sex==2&EER$BMI_group+PAL=="2",1.27,ifelse(EER$sex==2&EER$BMI_group+PAL=="3",1.45,
                                                                                                             ifelse(EER$sex==1&EER$BMI_group+PAL=="1",1.11,ifelse(EER$sex==1&EER$BMI_group+PAL=="2",1.25,1.48)))))
EER$EER_score <- ifelse(EER$sex=="1",EER_function("1",EER$age,EER$PA,EER$weight,EER$height),
                        EER_function("2",EER$age,EER$PA,EER$weight,EER$height))
EER <- merge(EER,DII_data[,c(1,4)],by="eid") %>% filter(!is.na(EER_score))
EER$EI_EER <- EER$Energy/EER$EER_score

EER_filter <- EER[EER$EI_EER>=mean(EER$EI_EER)-1.96*sd(EER$EI_EER)&EER$EI_EER<=mean(EER$EI_EER)+1.96*sd(EER$EI_EER)]
DII_data_EER.filter <- DII_data %>% filter(eid %in% EER_filter$eid)
DII_score_EER.filter <- DII(DII_data_EER.filter, DII_data_EER.filter$eid, 1, ALCOHOL_DII=DII_data_EER.filter$Alcohol, VITB12_DII=DII_data_EER.filter$`Vitamin B12`, 
                            VITB6_DII= DII_data_EER.filter$`Vitamin B6`, BCAROTENE_DII=DII_data_EER.filter$`Beta-carotene`, 
                            CAFFEINE_DII=NULL, CARB_DII=DII_data_EER.filter$Carbohydrate, CHOLES_DII=DII_data_EER.filter$Cholesterol,
                            KCAL_DII=DII_data_EER.filter$Energy, EUGENOL_DII=NULL,
                            TOTALFAT_DII=DII_data_EER.filter$Fat, FIBER_DII=DII_data_EER.filter$`Englyst fibre`, FOLICACID_DII=DII_data_EER.filter$Folate,
                            GARLIC_DII=NULL, GINGER_DII=NULL,IRON_DII=DII_data_EER.filter$Iron, MG_DII=DII_data_EER.filter$Magnesium, 
                            MUFA_DII=DII_data_EER.filter$`Monounsaturated fatty acids`, NIACIN_DII=DII_data_EER.filter$`Niacin equivalent`,
                            N3FAT_DII=DII_data_EER.filter$`n-3 fatty acids`, N6FAT_DII=DII_data_EER.filter$`n-6 fatty acids`,
                            ONION_DII=NULL, PROTEIN_DII=DII_data_EER.filter$Protein, PUFA_DII=NULL, 
                            RIBOFLAVIN_DII=DII_data_EER.filter$Riboflavin,SAFFRON_DII=NULL, SATFAT_DII=DII_data_EER.filter$`Saturated fatty acids`,
                            SE_DII=DII_data_EER.filter$Selenium, THIAMIN_DII=DII_data_EER.filter$Thiamin, TRANSFAT_DII=DII_data_EER.filter$`Trans fatty acids`,
                            TURMERIC_DII=NULL, VITA_DII=DII_data_EER.filter$`Vitamin A retinol equivalents`,
                            VITC_DII=DII_data_EER.filter$`Vitamin C`, VITD_DII=DII_data_EER.filter$`Vitamin D`,
                            VITE_DII=DII_data_EER.filter$`Vitamin E`, ZN_DII=DII_data_EER.filter$Zinc, TEA_DII=DII_data_EER.filter$mean_sum_26141_Tea,
                            FLA3OL_DII=NULL,FLAVONES_DII=NULL,FLAVONOLS_DII=NULL,FLAVONONES_DII=NULL,ANTHOC_DII=NULL,ISOFLAVONES_DII=NULL,
                            PEPPER_DII=NULL,THYME_DII=NULL,ROSEMARY_DII=NULL)
colnames(DII_score_EER.filter)[1] <- "eid"
DII_score_EER.filter$DII_GROUP <- Quantile(DII_score_EER.filter$DII_ALL)
DII_score_EER.filter$DII_GROUP <- factor(DII_score_EER.filter$DII_GROUP,levels=c("Q3","Q1","Q2","Q4","Q5"))
DII_score_EER.filter$DII_NOETOH_GROUP <- Quantile(DII_score_EER.filter$DII_NOETOH)
DII_score_EER.filter$DII_NOETOH_GROUP <- factor(DII_score_EER.filter$DII_NOETOH_GROUP,levels=c("Q3","Q1","Q2","Q4","Q5"))
DII_score_EER.filter <- merge(DII_score_EER.filter,DII_data,by="eid")
DII_score_EER.filter$EDII_ALL <- DII_score_EER.filter$DII_ALL*1000/DII_score_EER.filter$Energy
DII_score_EER.filter$EDII_GROUP <- Quantile(DII_score_EER.filter$EDII_ALL)
DII_score_EER.filter$EDII_GROUP <- factor(DII_score_EER.filter$EDII_GROUP,levels=c("Q3","Q1","Q2","Q4","Q5"))


# covariate data import ---------------------------------------------------
covariate_new <- fread("E:/paper/paper 6/UKB_dietary/primary data/UKB_covariate.csv")
covariate <- fread("E:/paper/paper 6/UKB_dietary/primary data/UKB_covariate13.csv")
###MET###
covariate <- covariate %>% mutate(
  MET_group=ifelse(activity<600,"1",ifelse(activity<3000,"2","3")))
###BMI###
covariate <- covariate %>% mutate(
  BMI_group=ifelse(BMI<18.5,"1",ifelse(BMI<25,"2",ifelse(BMI<30,"3","4"))))
###family history###
covariate_family_history <- covariate_new %>%  select(1,contains(c("20107","20110"))) %>% dplyr::filter(eid %in% DII_score_EER.filter$eid)
covariate_family_history$history <- ifelse(str_detect(covariate_family_history$`20107-0.0`,"10")=="TRUE"|str_detect(covariate_family_history$`20110-0.0`,"10")=='TRUE',"1","0")
###Menopause###
covariate_Menopause <- covariate_new %>%  select(1,contains(c("2724"))) %>% dplyr::filter(eid %in% DII_score_EER.filter$eid)
###Hypertension###
covariate_Hypertension <- covariate_new %>%  select(1,contains(c("6150","4080","93","4079","94","6177","6153"))) %>% dplyr::filter(eid %in% DII_score_EER.filter$eid)
covariate_Hypertension$Hypertension<- ifelse(covariate_Hypertension$`4080-0.0`>140|covariate_Hypertension$`93-0.0`>140|covariate_Hypertension$`4079-0.0`>90|covariate_Hypertension$`94-0.0`>90,"1","0")
covariate_Hypertension$HyperMed_history <- ifelse(str_detect(covariate_Hypertension$`6150-0.0`,"4")=="TRUE"|str_detect(covariate_Hypertension$`6153-0.0`,"2")=='TRUE'|str_detect(covariate_Hypertension$`6177-0.0`,"2")=='TRUE',"1","0")
covariate_Hypertension$Hypertension_history <- ifelse(covariate_Hypertension$Hypertension=="1"|covariate_Hypertension$HyperMed_history=="1","1","0")
covariate_Hypertension <- covariate_Hypertension %>% replace_na_dt(Hypertension_history,to="0")
###Diabetes###
covariate_Diabetes <- covariate_new %>% select(1,contains(c("2443","6177","6153","30750"))) %>% dplyr::filter(eid %in% DII_score_EER.filter$eid)
covariate_Diabetes$history6177 <- ifelse(str_detect(covariate_Diabetes$`6177-0.0`,"-3")=="TRUE","-9",covariate_Diabetes$`6177-0.0`)
covariate_Diabetes$history6153 <- ifelse(str_detect(covariate_Diabetes$`6153-0.0`,"-3")=="TRUE","-9",covariate_Diabetes$`6153-0.0`)
covariate_Diabetes$diabetes_history <- ifelse(str_detect(covariate_Diabetes$history6177,"3")=="TRUE"|str_detect(covariate_Diabetes$history6153,"3")=='TRUE',"1","0")


###High cholesterol###
covariate_cholesterol <- covariate_new %>% select(1,contains(c("6177","6153","30690","30760"))) %>% dplyr::filter(eid %in% DII_score_EER.filter$eid)
covariate_cholesterol$history6177 <- ifelse(str_detect(covariate_cholesterol$`6177-0.0`,"-1")=="TRUE","-9",covariate_cholesterol$`6177-0.0`)
covariate_cholesterol$history6153 <- ifelse(str_detect(covariate_cholesterol$`6153-0.0`,"-1")=="TRUE","-9",covariate_cholesterol$`6153-0.0`)
covariate_cholesterol$cholesterol_history <- ifelse(str_detect(covariate_cholesterol$history6177,"1")=="TRUE"|str_detect(covariate_cholesterol$history6153,"1")=='TRUE',"1","0")
covariate_cholesterol$cholesterol_lipids<- ifelse(covariate_cholesterol$`30690-0.0`>5|covariate_cholesterol$`30760-0.0`>3,"1","0")
covariate_cholesterol$cholesterol <- ifelse(covariate_cholesterol$cholesterol_history=="1"|covariate_cholesterol$cholesterol_lipids=="1","1","0")


covariate <- list(covariate_Diabetes[,c(1,13)],covariate,covariate_Hypertension[,c(1,24)],covariate_family_history[,c(1,4)],covariate_cholesterol[,c(1,13)])
covariate <- Reduce(my_merge, covariate)
#covariate <- covariate %>% dplyr::filter(age>=50)
rm(covariate_cholesterol,covariate_Diabetes,covariate_family_history,covariate_Hypertension,covariate_Menopause,covariate_new)

# brain disorder analysis -------------------------------------------------

## brain disorders data ----------------------------------------------------
# braindisorder_files <- list.files(path = paste0(primary_data.path,"BrainDisorder"),pattern = "csv",full.names = T)
# dt_braindisorder <- braindisorder_files %>%
#   map(~read_and_select(.)) %>%
#   reduce(full_join, by = "eid")
dt_braindisorder <- fread("E:/paper/paper-Connectomes WES/primary_data/BrainDisorder/BrainDisorder.csv")
dt_braindisorder <- dt_braindisorder %>% select("eid",contains(c("Anxiety", "MDD", "PD", "Sleep","AD", "Dementia","Stroke")))
dt_re_date <- read_table(paste0(primary_data.path,"field53.tsv.gz")) #%>% dplyr::filter(eid %in% dt_connectomes_eid$V1)
dt_braindisorder <- merge(dt_braindisorder,dt_re_date,by="eid")
dt_disorder_day <- grep("day",colnames(dt_braindisorder),value = T)
dt_braindisorder <- setDT(dt_braindisorder)
for (i in 1:length(dt_disorder_day)) {
  dt_braindisorder[,paste0(dt_disorder_day[i],"_time") := dt_braindisorder[,"53-0.0"]+dt_braindisorder[,dt_disorder_day[i],with = F]]
}
dt_diet_time <- read_csv(paste0(primary_data.path,"diet_date.csv"))
dt_diet_time <- dt_diet_time %>%
  dplyr::filter(rowSums(is.na(select(., 2:6))) < 5)
dt_diet_time <- dt_diet_time %>%
  dplyr::mutate(diet_time = pmax(`X105010-0.0`, `X105010-1.0`, `X105010-2.0`,
                                           `X105010-3.0`,`X105010-4.0`, na.rm = TRUE))
dt_braindisorder <- merge(dt_braindisorder,dt_diet_time[,c(1,7)],by="eid")
dt_braindisorder <- setDT(dt_braindisorder)
for (i in 1:length(dt_disorder_day)) {
  dt_braindisorder[,dt_disorder_day[i] := time_length(interval(as_date(dt_braindisorder[["diet_time"]]),
                                                               as_date(dt_braindisorder[[paste0(dt_disorder_day[i],"_time")]])),unit = "years")]
}

dt_braindisorder <- dt_braindisorder_data
dt_braindisorder_control <- dt_braindisorder %>% 
  dplyr::filter(if_all(ends_with("status"), ~ .x == 0))
braindisorder_cols <- c("Anxiety_status", "MDD_status", "PD_status", "Sleep_status", 
                        "AD_status", "Dementia_status", "Stroke_status")
dt_braindisorder_list <- lapply(braindisorder_cols, function(col) {
  rbind(dt_braindisorder[get(col) == 1],dt_braindisorder_control)
})
names(dt_braindisorder_list) <- braindisorder_cols




## Divide into three group -------------------------------------------------
DII_score_EER.filter$DII_GROUP_QQ3 <- Third.quartile(DII_score_EER.filter$DII_ALL)
DII_score_EER.filter$DII_GROUP_QQ3 <- factor(DII_score_EER.filter$DII_GROUP_QQ3,levels=c("Q1","Q2","Q3"))
DII_score_EER.filter$EDII_GROUP_QQ3 <- Third.quartile(DII_score_EER.filter$EDII_ALL)
DII_score_EER.filter$EDII_GROUP_QQ3 <- factor(DII_score_EER.filter$EDII_GROUP_QQ3,levels=c("Q1","Q2","Q3"))
DII_score_EER.filter <- DII_score_EER.filter %>% 
  dplyr::group_by(DII_GROUP_QQ3) %>% 
  dplyr::mutate(DII_GROUP_QQ3_trend = mean(DII_ALL)) %>% 
  dplyr::ungroup()
DII_score_EER.filter <- DII_score_EER.filter %>% 
  dplyr::group_by(EDII_GROUP_QQ3) %>% 
  dplyr::mutate(EDII_GROUP_QQ3_trend = mean(EDII_ALL)) %>% 
  dplyr::ungroup()
for (i in 1:length(dt_braindisorder_list)) {
  dt_braindisorder_list[[i]] <- merge(dt_braindisorder_list[[i]],DII_score_EER.filter[,c(1,124:127)],by="eid")
}

### cox analysis -----------------------------------------------------------
#### brain disorder analysis -------------------------------------------------
variable.names.y <- c("Anxiety", "MDD","PD","Sleep")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_LF_result_EER.filter_Q3 <- as.data.frame(Result)
#### cox dementia (APOE) -----------------------------------------------------
variable.names.y <- c("Dementia","AD")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}

Cox_LF_result_EER.filter_Q3 <- rbind(Cox_LF_result_EER.filter_Q3,as.data.frame(Result))
#### cox Stroke (risk factor) ------------------------------------------------------------
variable.names.y <- c("Stroke")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_LF_result_EER.filter_Q3 <- rbind(Cox_LF_result_EER.filter_Q3,as.data.frame(Result))
Cox_LF_result_EER.filter_continue <- Cox_LF_result_EER.filter_Q3 %>% dplyr::filter(!Var2 %like% 'trend')
Cox_LF_result_EER.filter_continue$p.adj <- p.adjust(Cox_LF_result_EER.filter_continue$p.value,method = "fdr")
Cox_LF_result_EER.filter_trend <-  Cox_LF_result_EER.filter_Q3 %>% dplyr::filter(Var2 %like% 'trend')
Cox_LF_result_EER.filter_trend$p.adj <- p.adjust(Cox_LF_result_EER.filter_trend$p.value,method = "fdr")
Cox_LF_result_EER.filter_Q3 <- rbind(Cox_LF_result_EER.filter_continue,Cox_LF_result_EER.filter_trend)
fwrite(Cox_LF_result_EER.filter_Q3,"Cox_LF_result_EER.filter_Q3.csv",row.names = FALSE)
### group analysis -----------------------------------------------------------------
#### brain disorder analysis -------------------------------------------------
variable.names.y <- c("Anxiety","MDD", "PD","Sleep")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_QQ3 <- factor(dt_braindisorder$DII_GROUP_QQ3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_QQ3 <- factor(dt_braindisorder$EDII_GROUP_QQ3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter <- as.data.frame(Result)
Cox_FL_result_EER.filter
#### cox dementia (APOE) -----------------------------------------------------
variable.names.y <- c("Dementia","AD")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_Q3 <- factor(dt_braindisorder$DII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_Q3 <- factor(dt_braindisorder$EDII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter <- rbind(Cox_FL_result_EER.filter,as.data.frame(Result))
#### cox Stroke (risk factor) ------------------------------------------------------------
variable.names.y <- c("Stroke")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_Q3 <- factor(dt_braindisorder$DII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_Q3 <- factor(dt_braindisorder$EDII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter <- rbind(Cox_FL_result_EER.filter,as.data.frame(Result))
Cox_FL_result_EER.filter
fwrite(Cox_FL_result_EER.filter,"Cox_FL_result_EER.filter_Q3.csv",row.names = FALSE)




### PH analysis -------------------------------------------------------------
variable.names.y <- c("Anxiety", "MDD", "PD", "Sleep")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend',"DII_GROUP_QQ3","EDII_GROUP_QQ3")
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    test.ph <- cox.zph(cox_model[[i]])
    ggcoxzph(test.ph)
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    row.var <- rownames(test.ph[[1]])
    p.value <- test.ph[[1]][,3]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'row.var'=row.var,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
PH_LF_result_EER.filter <- as.data.frame(Result)
PH_LF_result_EER.filter



variable.names.y <- c("AD", "Dementia")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend',"DII_GROUP_QQ3","EDII_GROUP_QQ3")
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    test.ph <- cox.zph(cox_model[[i]])
    ggcoxzph(test.ph)
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    row.var <- rownames(test.ph[[1]])
    p.value <- test.ph[[1]][,3]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'row.var'=row.var,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
PH_LF_result_EER.filter <- rbind(PH_LF_result_EER.filter,as.data.frame(Result))


variable.names.y <- c("Stroke")
variable.names.x <- c('DII_ALL','EDII_ALL','EDII_GROUP_QQ3_trend','DII_GROUP_QQ3_trend',"DII_GROUP_QQ3","EDII_GROUP_QQ3")
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    test.ph <- cox.zph(cox_model[[i]])
    ggcoxzph(test.ph)
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    row.var <- rownames(test.ph[[1]])
    p.value <- test.ph[[1]][,3]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'row.var'=row.var,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
PH_LF_result_EER.filter <- rbind(PH_LF_result_EER.filter,as.data.frame(Result))
fwrite(PH_LF_result_EER.filter,"PH_LF_result_EER.filter_Q3.csv",row.names = FALSE)
## Sensitive analysis ------------------------------------------------------
### twice follow available --------------------------------------------------
#### twice follow available data--------------------------------------------------
dt_braindisorder_twice <- dt_braindisorder_list
DII_score_EER_twice <- DII_score_EER.filter %>% dplyr::filter(eid %in% marked_Field_ID_summry$eid)
# eid_M50 <- covariate$eid[which(covariate$age>=55)]
# DII_score_EER.filter <- DII_score_EER.filter %>% dplyr::filter(eid %in% eid_M50)
DII_score_EER_twice$DII_GROUP_QQ3_twice <- Third.quartile(DII_score_EER_twice$DII_ALL)
DII_score_EER_twice$DII_GROUP_QQ3_twice <- factor(DII_score_EER_twice$DII_GROUP_QQ3_twice,levels=c("Q1","Q2","Q3"))
DII_score_EER_twice$EDII_GROUP_QQ3_twice <- Third.quartile(DII_score_EER_twice$EDII_ALL)
DII_score_EER_twice$EDII_GROUP_QQ3_twice <- factor(DII_score_EER_twice$EDII_GROUP_QQ3_twice,levels=c("Q1","Q2","Q3"))
DII_score_EER_twice <- DII_score_EER_twice %>% 
  dplyr::group_by(DII_GROUP_QQ3_twice) %>% 
  dplyr::mutate(DII_GROUP_QQ3_trend_twice = mean(DII_ALL)) %>% 
  dplyr::ungroup()
DII_score_EER_twice <- DII_score_EER_twice %>% 
  dplyr::group_by(EDII_GROUP_QQ3_twice) %>% 
  dplyr::mutate(EDII_GROUP_QQ3_trend_twice = mean(EDII_ALL)) %>% 
  dplyr::ungroup()
for (i in 1:length(dt_braindisorder_twice)) {
  dt_braindisorder_twice[[i]] <- merge(dt_braindisorder_twice[[i]],DII_score_EER_twice[,c(1,128:131)],by="eid")
}

#### twice follow available analysis cox model -------------------------------
variable.names.y <- c("Anxiety", "MDD", "PD", "Sleep")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_LF_result_twice <- as.data.frame(Result)

variable.names.y <- c("AD", "Dementia")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}

Cox_LF_result_twice <- rbind(Cox_LF_result_twice,as.data.frame(Result))


variable.names.y <- c("Stroke")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history')))
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}

Cox_LF_result_twice <- rbind(Cox_LF_result_twice,as.data.frame(Result))
Cox_LF_result_twice$p.adj <- p.adjust(Cox_LF_result_twice$p.value,method = "fdr")

#### twice follow available analysis cox model -------------------------------

variable.names.y <- c("Anxiety", "MDD", "PD", "Sleep")
# variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
variable.names.x <- c('DII_GROUP_QQ3_twice','EDII_GROUP_QQ3_twice')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP <- factor(dt_braindisorder$DII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  dt_braindisorder$EDII_GROUP <- factor(dt_braindisorder$EDII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_twice <- as.data.frame(Result)
Cox_FL_result_twice


variable.names.y <- c("AD", "Dementia")
# variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
variable.names.x <- c('DII_GROUP_QQ3_twice','EDII_GROUP_QQ3_twice')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP <- factor(dt_braindisorder$DII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  dt_braindisorder$EDII_GROUP <- factor(dt_braindisorder$EDII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_twice <- rbind(Cox_FL_result_twice,as.data.frame(Result))


variable.names.y <- c("Stroke")
# variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
variable.names.x <- c('DII_GROUP_QQ3_twice','EDII_GROUP_QQ3_twice')

Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_twice[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP <- factor(dt_braindisorder$DII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  dt_braindisorder$EDII_GROUP <- factor(dt_braindisorder$EDII_GROUP,levels = c("Q1","Q2","Q3","Q4","Q5"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_twice <- rbind(Cox_FL_result_twice,as.data.frame(Result))
Cox_FL_result_twice
fwrite(Cox_FL_result_twice,"Cox_FL_result_twice.csv",row.names = FALSE)
fwrite(Cox_LF_result_twice,"Cox_LF_result_twice.csv",row.names = FALSE)



### more covariate -----------------------------------------------------------
#### cox brain disorders -----------------------------------------------------
variable.names.y <- c("Anxiety", "MDD","PD","Sleep")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+diabetes_history+Hypertension_history+cholesterol')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_LF_result_EER.filter_M2 <- as.data.frame(Result)
#### cox dementia (APOE) -----------------------------------------------------
variable.names.y <- c("Dementia","AD")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4+diabetes_history+Hypertension_history+cholesterol')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}

Cox_LF_result_EER.filter_M2 <- rbind(Cox_LF_result_EER.filter_M2,as.data.frame(Result))
#### cox Stroke (risk factor) ------------------------------------------------------------
variable.names.y <- c("Stroke")
variable.names.x <- c('DII_ALL','EDII_ALL')
Result <- list()
for (j in 1:length(variable.names.y)) {
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history+diabetes_history+cholesterol')))
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    number.data=length(dt_braindisorder$eid)
    case.data=length(dt_braindisorder$eid[which(dt_braindisorder[[paste0(variable.names.y[[j]],'_status')]]==1)])
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'number'=number.data,
                    'case'=case.data,
                    'beta'=beta, 
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_LF_result_EER.filter_M2 <- rbind(Cox_LF_result_EER.filter_M2,as.data.frame(Result))
Cox_LF_result_EER.filter_continue <- Cox_LF_result_EER.filter_M2 %>% dplyr::filter(!Var2 %like% 'trend')
Cox_LF_result_EER.filter_continue$p.adj <- p.adjust(Cox_LF_result_EER.filter_continue$p.value,method = "fdr")
Cox_LF_result_EER.filter_trend <-  Cox_LF_result_EER.filter_M2 %>% dplyr::filter(Var2 %like% 'trend')
Cox_LF_result_EER.filter_trend$p.adj <- p.adjust(Cox_LF_result_EER.filter_trend$p.value,method = "fdr")
Cox_LF_result_EER.filter_M2 <- rbind(Cox_LF_result_EER.filter_continue,Cox_LF_result_EER.filter_trend)
fwrite(Cox_LF_result_EER.filter_M2,"Cox_LF_result_EER.filter_M2.csv",row.names = FALSE)
#### group analysis brain disorder analysis -------------------------------------------------
variable.names.y <- c("Anxiety","MDD", "PD","Sleep")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_QQ3 <- factor(dt_braindisorder$DII_GROUP_QQ3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_QQ3 <- factor(dt_braindisorder$EDII_GROUP_QQ3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+diabetes_history+Hypertension_history+cholesterol')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter_M2 <- as.data.frame(Result)
Cox_FL_result_EER.filter_M2
#### group analysis cox dementia (APOE) -----------------------------------------------------
variable.names.y <- c("Dementia","AD")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_Q3 <- factor(dt_braindisorder$DII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_Q3 <- factor(dt_braindisorder$EDII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4+diabetes_history+Hypertension_history+cholesterol')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter_M2 <- rbind(Cox_FL_result_EER.filter_M2,as.data.frame(Result))
#### group analysis cox Stroke (risk factor) ------------------------------------------------------------
variable.names.y <- c("Stroke")
# variable.names.x <- c('DII_GROUP_Q3','EDII_GROUP_Q3')
variable.names.x <- c('DII_GROUP_QQ3','EDII_GROUP_QQ3')
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_braindisorder <- dt_braindisorder_list[[variable.names.y[j]]]
  dt_braindisorder$DII_GROUP_Q3 <- factor(dt_braindisorder$DII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  dt_braindisorder$EDII_GROUP_Q3 <- factor(dt_braindisorder$EDII_GROUP_Q3,levels = c("Q1","Q2","Q3"))
  cox_formula <- sapply(variable.names.x,
                        function(x) as.formula(
                          paste0('Surv(',variable.names.y[[j]],'_days,',variable.names.y[[j]],'_status==1)~',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+Hypertension_history+diabetes_history+cholesterol')))
  cox_model <- lapply(cox_formula, function(x){coxph(x, data = dt_braindisorder)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(cox_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- rownames(x[["coefficients"]])[1:2]
    p.value <- sprintf('%0.3f',x[["coefficients"]][1:2,5])
    beta <- sprintf('%0.3f',x[["coefficients"]][1:2,1])
    HR <- sprintf('%0.3f',x[["coefficients"]][1:2,2])
    HR.confint.lower <- sprintf('%0.3f',x[["conf.int"]][,"lower .95"][1:2])
    HR.confint.upper <- sprintf('%0.3f',x[["conf.int"]][,"upper .95"][1:2])
    HR <- paste0(HR, " (", 
                 HR.confint.lower, "-", HR.confint.upper, ")")
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'beta'=beta, 
                    'dataset'=variable.names.y[j],
                    'HR (95% CI)'=HR,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
Cox_FL_result_EER.filter_M2 <- rbind(Cox_FL_result_EER.filter_M2,as.data.frame(Result))
Cox_FL_result_EER.filter_M2
fwrite(Cox_FL_result_EER.filter_M2,"Cox_FL_result_EER.filter_M2.csv",row.names = FALSE)
## SEM function annotation -------------------------------------------------

### SEM function annotation data --------------------------------------------

dt_sem <- merge(dt_inflamation,dt_braindisorder_data,by="eid") %>% dplyr::filter(eid %in% eid_M50)
dt_sem <- dt_sem %>% dplyr::select(c(2:7,9:14,129,135:140,147:148,contains('status')))
dt_sem <- na.omit(dt_sem)
dt_sem$BMI_group <- as.numeric(dt_sem$BMI_group)
dt_sem$PAL <- as.numeric(dt_sem$PAL)

### SEM function annotation analysis ----------------------------------------

model=
  "
  brain.disorders =~ MDD_status+Anxiety_status+Dementia_status
  inflammation_marker =~ Z_NLR+Z_PLR+Z_LMR+Z_SII+Z_CRP
  inflammation_marker ~ age+sex+townsend+education+PAL+BMI_group+b1*DII_ALL
  brain.disorders ~ age+sex+townsend+education+PAL+smoking+BMI_group+b2*DII_ALL+c1*inflammation_marker
  ie1 := b1*c1
  de := b2
"
fit_dep <- sem(model,data = dt_sem, se="boot", bootstrap = 1000)
fit_summary_dep<-summary(fit_dep,standardized=T)
fit_summary_dep
fwrite(fit_summary_dep$pe,"fit_summary_DII_dep.csv")

model=
  "
  brain.disorders =~ MDD_status+Anxiety_status+Dementia_status
  inflammation_marker =~ Z_NLR+Z_PLR+Z_LMR+Z_SII+Z_CRP
  inflammation_marker ~ age+sex+townsend+education+PAL+BMI_group+b1*EDII_ALL
  brain.disorders ~ age+sex+townsend+education+PAL+smoking+BMI_group+b2*EDII_ALL+c1*inflammation_marker
  ie1 := b1*c1
  de := b2
"
fit_dep <- sem(model,data = dt_sem, se="boot", bootstrap = 1000)
fit_summary_dep<-summary(fit_dep,standardized=T)
fit_summary_dep
fwrite(fit_summary_dep$pe,"fit_summary_EDII_dep.csv")



## inflammation ------------------------------------------------------------
dt_inflamation <- fread("F:/paper/paper 6/UKB_dietary/primary data/Blood_UKB_bl_Zscore.csv") %>% dplyr::select(2,contains("Z"))
dt_inflamation <- merge(dt_inflamation,DII_score_EER.filter,by="eid")
dt_inflamation <- merge(dt_inflamation,mydata[,c(1,5:21)],by="eid")
dt_inflamation <- dt_inflamation %>% dplyr::filter(eid %in% eid_final)

### blood inflammation analysis with DII score and EDII score ---------------------------------------------

variable.names.y <- colnames(dt_inflamation)[2:13]
variable.names.x <- c('DII_ALL','EDII_ALL')

###lm###
Result <- list()
for (j in 1:length(variable.names.y)) {
  lm_formula <- sapply(variable.names.x,
                       function(x) as.formula(
                         paste0(variable.names.y[[j]],' ~ ',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+PAL')))
  lm_model <- lapply(lm_formula, function(x){lm(x, data = dt_inflamation)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(lm_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- x[["coefficients"]][2,4]
    t_value <- x[["coefficients"]][2,3]
    beta_value <- x[["coefficients"]][2,2] 
    regression_value <- rownames(x[["coefficients"]])[2]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'regression'=regression_value,
                    'beta'=beta_value,
                    't'=t_value,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
lm_LF_result_inflammation <- as.data.frame(Result)
lm_LF_result_inflammation$P.adj <- p.adjust(lm_LF_result_inflammation$p.value,method = "fdr")
lm_LF_result_inflammation

fwrite(lm_LF_result_inflammation,"lm_LF_result_inflammation.csv",row.names = FALSE)
### blood inflammation analysis with DII score and EDII score ---------------------------------------------
dt_brain_inflammation <- dt_braindisorder_list
for (i in 1:length(dt_brain_inflammation)){
  dt_brain_inflammation[[i]] <- merge(dt_brain_inflammation[[i]],dt_inflamation[,1:13],by="eid")
}


variable.names.x <- colnames(dt_inflamation)[2:13]
variable.names.y <- c("Anxiety","MDD", "PD","Sleep")

###lm###
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_brain_inflammation1 <- dt_brain_inflammation[[variable.names.y[j]]]
  lm_formula <- sapply(variable.names.x,
                       function(x) as.formula(
                         paste0(paste0(variable.names.y[[j]],"_status"),' ~ ',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group')))
  lm_model <- lapply(lm_formula, function(x){lm(x, data = dt_brain_inflammation1)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(lm_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- x[["coefficients"]][2,4]
    t_value <- x[["coefficients"]][2,3]
    beta_value <- x[["coefficients"]][2,2] 
    regression_value <- rownames(x[["coefficients"]])[2]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'regression'=regression_value,
                    'beta'=beta_value,
                    't'=t_value,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
lm_LF_result_braininflammation <- as.data.frame(Result)


variable.names.x <- colnames(dt_inflamation)[2:13]
variable.names.y <- c("Dementia","AD")

###lm###
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_brain_inflammation1 <- dt_brain_inflammation[[variable.names.y[j]]]
  lm_formula <- sapply(variable.names.x,
                       function(x) as.formula(
                         paste0(paste0(variable.names.y[[j]],"_status"),' ~ ',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  lm_model <- lapply(lm_formula, function(x){lm(x, data = dt_brain_inflammation1)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(lm_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- x[["coefficients"]][2,4]
    t_value <- x[["coefficients"]][2,3]
    beta_value <- x[["coefficients"]][2,2] 
    regression_value <- rownames(x[["coefficients"]])[2]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'regression'=regression_value,
                    'beta'=beta_value,
                    't'=t_value,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
lm_LF_result_braininflammation <- rbind(lm_LF_result_braininflammation,as.data.frame(Result))


variable.names.x <- colnames(dt_inflamation)[2:13]
variable.names.y <- c("Stroke")

###lm###
Result <- list()
for (j in 1:length(variable.names.y)) {
  dt_brain_inflammation1 <- dt_brain_inflammation[[variable.names.y[j]]]
  lm_formula <- sapply(variable.names.x,
                       function(x) as.formula(
                         paste0(paste0(variable.names.y[[j]],"_status"),' ~ ',x,'+age+sex+ethnic+townsend+education+smoking+BMI_group+MET_group+APOE4')))
  lm_model <- lapply(lm_formula, function(x){lm(x, data = dt_brain_inflammation1)})
  result <- list()
  for (i in 1:length(variable.names.x)) {
    x <- summary(lm_model[[i]])
    Var1 <- variable.names.y[[j]]
    Var2 <- variable.names.x[[i]]
    p.value <- x[["coefficients"]][2,4]
    t_value <- x[["coefficients"]][2,3]
    beta_value <- x[["coefficients"]][2,2] 
    regression_value <- rownames(x[["coefficients"]])[2]
    res<-data.frame('Var1'=Var1,
                    'Var2'=Var2,
                    'regression'=regression_value,
                    'beta'=beta_value,
                    't'=t_value,
                    'p.value'=p.value)
    result <- rbind(result,res)
    rm(res)
  }
  Result <- rbind(Result,result)
}
lm_LF_result_braininflammation <- rbind(lm_LF_result_braininflammation,as.data.frame(Result))


lm_LF_result_braininflammation$P.adj <- p.adjust(lm_LF_result_braininflammation$p.value,method = "fdr")
lm_LF_result_braininflammation
fwrite(lm_LF_result_braininflammation,"lm_LF_result_braininflammation.csv",row.names = FALSE)
