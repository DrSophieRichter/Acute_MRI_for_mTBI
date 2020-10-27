#load packages
library(dplyr)
library(tidyr)
library(kableExtra)
library(broom)
library(grid)
library(gridExtra)

###LOAD DATA

##load diffusion data
  #this data is in the long format with up to three rows per patient (MR1, MR2, MR3)
  #The columns are the white matter tracts (identified by number) 
  #plus some other columns (subject ID, scan ID, timepoint etc)
  #the column "timepoint" is a factor variable with 4 levels:
    # "control", "MR1", "MR2", "MR3" denoting whether a scan is 
    # a healthy control scan or the first, second or third scan of a patient
   
fa <- read.csv("fa_allsubjects_alltracts_nooutliers.csv")
md <- read.csv("md_allsubjects_alltracts_nooutliers.csv")

#load outcome data
  #There is one row per patient
  #the column Outcome contains the dichotomised 3-months outcome
  # i.e "Favourable" if GOSE = 8, "Unfavourable" if GOSE < 8
out <- read.csv("Paper1_demog.csv")
out <- out %>% select(subject, 
                      Outcome)

#load volume data
  #this table is organised the same way as the dti table
  #with multiple rows per patient
  #and one column for each ROI, plus some identifier columns
  #I will only use the ROI Cerebral white matter
vol <- read.csv("normalised_volumes.csv")
vol <- vol %>% select(subject, 
                      study, 
                      scanner, 
                      timepoint, 
                      nCWM = Cerebralwhitematter)

#view number of patients per timepoint
summary(fa$timepoint)
temp <- fa %>% 
        group_by(timepoint, scanner) %>% 
        summarise(n())
#view number of scanners per timepoint
summary(temp$timepoint) 


###DATA PREPARATION

##Calculating mean and sd for controls for each scanner and tract

#getting reference data for fa
c_fa    <- fa %>% filter(timepoint == "control")
temp    <- c_fa %>% gather(key = "tract", value = "raw_fa", 9:80)
ref_fa  <- temp %>% 
            group_by(scanner, tract) %>% 
            summarise(mean_fa = mean(raw_fa, na.rm = TRUE), 
                      sd_fa = sd(raw_fa, na.rm = TRUE))

#getting reference data for md
c_md    <- md   %>% filter(timepoint == "control")
temp    <- c_md %>% gather(key = "tract", value = "raw_md", 9:80)
ref_md  <- temp %>% 
            group_by(scanner, tract) %>% 
            summarise(mean_md = mean(raw_md, na.rm = TRUE), 
                      sd_md = sd(raw_md, na.rm = TRUE))

##Calculating the standardised deviation of FA

#Ordering the table with patients and tracts so that it can be combined with reference data
p_fa    <- fa   %>% filter(timepoint != "control")
temp    <- p_fa %>% gather(key = "tract", value = "raw_fa", 9:80)

#merging the patient fa values with those of the respective control population
all_fa  <- merge(temp, ref_fa, 
                 by=c("tract", "scanner"), 
                 all.x = TRUE, 
                 all.y = FALSE)

#calculating by how much each patient fa value deviates from the control population
all_fa$delta_fa <- all_fa$raw_fa - all_fa$mean_fa
all_fa$dev_fa   <- all_fa$delta_fa/all_fa$sd_fa
all_fa$dev_fa   <- round(all_fa$dev_fa, 2) #round to two decimal places
all_fa          <- all_fa %>% select(-X)

##Calculating the standardised deviation of MD

#Ordering the table with patients and tracts so that it can be combined with reference data
p_md    <- md    %>% filter(timepoint != "control")
temp    <- p_md  %>% gather(key = "tract", value = "raw_md", 9:80)

#merging the patient fa values with those of the respective control population
all_md  <- merge(temp, ref_md, 
                 by=c("tract", "scanner"), 
                 all.x = TRUE, 
                 all.y = FALSE)

#calculating by how much each patient fa value deviates from the control population
all_md$delta_md <- all_md$raw_md - all_md$mean_md
all_md$dev_md   <- all_md$delta_md/all_md$sd_md
all_md$dev_md   <- round(all_md$dev_md, 2) #round to two decimal places
all_md          <- all_md %>% select(-X)

##Combine data for FA, MD and Outcome into one table

dti         <- merge(all_fa, all_md, 
                     by = c("tract", 
                            "scanner", 
                            "subject", 
                            "scan", 
                            "study", 
                            "age", 
                            "sex", 
                            "timepoint", 
                            "group"), 
                     all.x = TRUE, 
                     all.y = TRUE)

data        <- merge(dti, out, 
                     by = "subject", 
                     all.x = TRUE, 
                     all.y = FALSE)

#Make sex a factor variable
data$sex    <- factor(data$sex,
                      levels = c("M", "F"),
                      labels = c("M", "F"))

#make "tract" a factor variable
dti$tract   <- factor(dti$tract,
                      levels = unique(dti$tract),
                      labels = unique(dti$tract))

#Make "outcome" a factor variable
data$gose3  <- factor(data$Outcome,
                      levels = c("Unfavourable", "Favourable"),
                      labels = c("bad", "good"))
levels(data$gose3)
summary(data$gose3)

#The reference level is listed first. So the models will use "bad" as the refernece level and predict the odds of a GOOD outcome


## Adding a normalised white matter volume to the data table

#In the volume data I treated CENTER Cambridge and DETECT Cambridge as the same site/scanner but for DTI they were different scanners. So to merge the data I have to separate CENTER and DETECT volumes from Cambridge into two separate "scanners"
vol$scanner <- ifelse(vol$study == "DETECT", "DETECT_Verio", as.character(vol$scanner))
vol$scanner <- factor(vol$scanner,
                      levels = unique(vol$scanner),
                      labels = unique(vol$scanner))

###Calculating mean and sd of VOLUME for controls for each scanner 
c_vol     <- vol %>% filter(timepoint == "control")
ref_vol   <- c_vol %>% 
              group_by(scanner) %>% 
              summarise(mean_vol = mean(nCWM, na.rm = TRUE), 
                        sd_vol = sd(nCWM, na.rm = TRUE))
##Calculating the standardised deviation of VOLUME

#Make table with patients so that it can be combined with referenence data
p_vol     <- vol %>% filter(timepoint != "control")

#merging the patient vol values with those of the respective control population
all_vol   <- merge(p_vol, ref_vol, 
                   by=c("scanner"), 
                   all.x = TRUE, 
                   all.y = FALSE)

#calculating by how much each patient vol value deviates from the control population
all_vol$delta_vol <- all_vol$nCWM - all_vol$mean_vol
all_vol$dev_vol   <- all_vol$delta_vol/all_vol$sd_vol
all_vol$dev_vol   <- round(all_vol$dev_vol, 2) #round to two decimal places
data      <- merge(data, all_vol, 
                   by = c("subject", "timepoint", "scanner", "study"), 
                   all.x = TRUE, 
                   all.y = FALSE)

### COUNTING THE NUMBER OF ABNORMAL TRACTS PER PATIENT

# I will define as abnormal a   
# *dev_fa smaller or equal to -2  
# *dev_md larger or equal to +2
# 
# Abnormal will be coded as 1, normal as 0  
# If a tract was not measured, I will not count it as abnormal.
# 
# I will make three variables  
# abn_fa - ONLY FA abnormal, MD is normal  
# abn_md - ONLY MD abnormal, FA is normal  
# abn_both - BOTH MD and FA abnormal  
# 
# The "ONLY" term is important to avoid co-linearity between abn_fa (or abn_md) and abn_both

data$abn_fa   <- ifelse(data$dev_fa <= (-2) & data$dev_md <= 2, 1, 0)
data$abn_md   <- ifelse(data$dev_md >= 2 & data$dev_fa >= (-2), 1, 0)
data$abn_both <- ifelse(data$dev_fa <= (-2) & data$dev_md >= 2, 1, 0)

### EXCLUSION

#There is no control data for scanner 9109c8_Ingenia.
#This affects 2 patients (one of which does not have GOSE anyway):
data %>% filter(scanner == "9109c8_Ingenia") %>% distinct(subject) %>% dim()
#Thus the 2 patients scanned on this scanner have to be excluded. 
data  <- data %>% filter(scanner != "9109c8_Ingenia")

###MAKE DATAFRAME FOR MR1

#summarise abnormal tracts by patient

mr1.df <- data %>% 
          filter(timepoint == "MR1") %>%
          group_by(subject) %>% 
          summarise(sum_abn_fa   = sum(abn_fa, na.rm = TRUE),
                    sum_abn_md   = sum(abn_md, na.rm = TRUE),
                    sum_abn_both = sum(abn_both, na.rm = TRUE),
                    dev_vol      = mean(dev_vol),
                    gose3        = unique(gose3), 
                    age          = mean(age),
                    sex          = unique(sex))

#Patients with GOSE available at 3 months
mr1_gose3.df <- mr1.df %>% 
                filter(is.na(mr1.df$gose3)==FALSE)
dim(mr1_gose3.df) 
summary(mr1_gose3.df$gose3)

###MAKE DATAFRAME FOR MR2

#summarise abnormal tracts by patient
mr2.df <- data %>% 
          filter(subject %in% mr1.df$subject)%>% #ensure same sets of patients used in all models
          filter(timepoint == "MR2") %>%
          group_by(subject) %>% 
          summarise(sum_abn_fa   = sum(abn_fa, na.rm = TRUE),
                    sum_abn_md   = sum(abn_md, na.rm = TRUE),
                    sum_abn_both = sum(abn_both, na.rm = TRUE),
                    gose3        = unique(gose3), 
                    age          = mean(age),
                    sex          = unique(sex),
                    dev_vol      = mean(dev_vol))

summary(mr2.df)

#Patients with GOSE available at 3 months
mr2_gose3.df <- mr2.df %>% 
                filter(is.na(mr2.df$gose3)==FALSE)
dim(mr2_gose3.df)
summary(mr2_gose3.df$gose3)


###PREPARE DATA WITHOUT IMAGING INFORMATION

#Select patients
zero.df       <- data %>% 
                  group_by(subject) %>% 
                  summarise(gose3 = unique(gose3), 
                            age = mean(age),
                            sex = unique(sex))

#Patients with GOSE available at 3 months
zero_gose3.df <- zero.df %>% 
                  filter(is.na(zero.df$gose3)==FALSE)
dim(zero_gose3.df)
summary(zero_gose3.df$gose3)

###PREPARE FUNCTIONS TO FIT MODELS

#Make a function to fit the model
mymodel <- function(mydata)
  {
  mymodel   <- glm(gose3 ~.,
               data = mydata, 
               family = binomial(link='logit'))
  return(mymodel)
}


#Make a function that returns model coefficients as odds ratio
#with 95% CI and p-value
mycoef  <- function(mydata, Model_name)
  {
  mymodel       <- mymodel(mydata)#calls previous function
  res           <- summary(mymodel)$coefficients %>% data.frame()
  res$OR        <- exp(res$Estimate)
  res           <- cbind(res, exp(confint(mymodel))) %>% round(2)
  colnames(res) <- c("Estimate", "SE", "z", "P", "OR", "LLCI", "HLCI")
  res$CI        <- paste(res$LLCI, "-", res$HLCI)
  res$Model     <- Model_name
  res$Sig       <- ifelse(res$P<0.05, "*", "")
  res$Variable  <- rownames(res)
  rownames(res) <- NULL
  res           <- res %>% select(Model, Variable, OR, CI, P, Sig)
  return(res)
  }

#Make a function to get AUC
my_auc <- function(mydata)
  {
  #fit the model
  fit               <- mymodel(mydata)
  #generate predicted outcome
  gose3_pred        <- predict(fit, type=c("response"))
  mydata$gose3_pred <- gose3_pred
  
  #calculate AUC and its 95% confidence interval
  library(pROC)
  g                 <- roc(gose3 ~ gose3_pred, data = mydata)
  auc               <- g$auc[1] %>% round(2)
  LL <- ci.auc(g)[1] %>% round(2)
  UL <- ci.auc(g)[3] %>% round(2)
  
  #put AUC and CI into a single cell
  myauc <- paste0(auc, " (", LL, "-", UL, ")")
  myauc
  
  return(myauc)
  }



#Make a function for K-fold cross validation of model any model with gose3 as the outcome
my_cv <- function(mydata)
  {
  set.seed(123)
  
  #specify the type of cross-validation
  library(caret)
  ctrl    <- trainControl(method = "cv", 
                          number = 10, 
                          savePredictions = TRUE)
  #specify the model
  mod_fit <- train(gose3 ~., 
                   data = mydata, 
                   method="glm", 
                   family="binomial",
                   trControl = ctrl)
  
  #collect accuracy and 95% confidence interval
  cv <- mod_fit$results 
  LL <-(cv[1,2] - 1.96*cv[1,4]) %>% round(2) # calculate LL 95% CI from SD 
  UL <-(cv[1,2] + 1.96*cv[1,4]) %>% round(2) # calculate UL 95% CI from SD
  cv <- paste0(cv[1,2] %>% round(2), #mean
                " (", 
                ifelse(LL<0, 0, LL),
                "-", 
                ifelse(UL>1, 1, UL),
                ")")
  
  #return accuracy (95% CI)
  return(cv)
}

#Make a function to get the ROC curve for any model with gose 3 as the outcome
my_roc <- function(mydata)
  {
  #fit the model
  fit               <- mymodel(mydata)
  #generate predicted outcome
  gose3_pred        <- predict(fit,type=c("response"))
  mydata$gose3_pred <- gose3_pred
  
  #calculate AUC and its 95% confidence interval
  library(pROC)
  g                 <- roc(gose3 ~ gose3_pred, data = mydata)
  
  return(g)
  }

###COMPLETE CASE ANALYSIS

##Zeromodel
mydata        <- zero_gose3.df[,c("gose3", 
                                  "age", 
                                  "sex")]
zero_fit      <- mymodel(mydata)
zero_coef     <- mycoef(mydata, "zero")
zero_auc      <- my_auc(mydata)
zero_cv       <- my_cv(mydata)
zero_roc      <- my_roc(mydata)

##MR1
#T1
mydata        <- mr1_gose3.df[,c("gose3", 
                                 "age", 
                                 "sex", 
                                 "dev_vol")]
mr1_t1_fit    <- mymodel(mydata)
mr1_t1_coef   <- mycoef(mydata, "mr1_t1")
mr1_t1_auc    <- my_auc(mydata)
mr1_t1_cv     <- my_cv(mydata)
mr1_t1_roc    <- my_roc(mydata)

#DTI
mydata        <- mr1_gose3.df[,c("gose3", 
                                 "age", 
                                 "sex", 
                                 "sum_abn_fa", 
                                 "sum_abn_md", 
                                 "sum_abn_both")]
mr1_dti_fit   <- mymodel(mydata)
mr1_dti_coef  <- mycoef(mydata, "mr1_dti")
mr1_dti_auc   <- my_auc(mydata)
mr1_dti_cv    <- my_cv(mydata)
mr1_dti_roc   <- my_roc(mydata)

#T1&DTI
mydata        <- mr1_gose3.df[,c("gose3", 
                                 "age", 
                                 "sex", 
                                 "sum_abn_fa", 
                                 "sum_abn_md", 
                                 "sum_abn_both", 
                                 "dev_vol")]
mr1_t1dti_fit  <- mymodel(mydata)
mr1_t1dti_coef <- mycoef(mydata, "mr1_t1dti")
mr1_t1dti_auc  <- my_auc(mydata)
mr1_t1dti_cv   <- my_cv(mydata)
mr1_t1dti_roc  <- my_roc(mydata)

##MR2
#T1
mydata        <- mr2_gose3.df[,c("gose3", 
                                 "age", 
                                 "sex", 
                                 "dev_vol")]
mr2_t1_fit    <- mymodel(mydata)
mr2_t1_coef   <- mycoef(mydata, "mr2_t1")
mr2_t1_auc    <- my_auc(mydata)
mr2_t1_cv     <- my_cv(mydata)
mr2_t1_roc    <- my_roc(mydata)

#DTI
mydata        <- mr2_gose3.df[,c("gose3", 
                                 "age", 
                                 "sex", 
                                 "sum_abn_fa", 
                                 "sum_abn_md", 
                                 "sum_abn_both")]
mr2_dti_fit  <- mymodel(mydata)
mr2_dti_coef <- mycoef(mydata, "mr2_dti")
mr2_dti_auc  <- my_auc(mydata)
mr2_dti_cv   <- my_cv(mydata)
mr2_dti_roc  <- my_roc(mydata)

#T1&DTI
mydata      <- mr2_gose3.df[,c("gose3", 
                               "age", 
                               "sex", 
                               "sum_abn_fa", 
                               "sum_abn_md", 
                               "sum_abn_both", 
                               "dev_vol")]
mr2_t1dti_fit  <- mymodel(mydata)
mr2_t1dti_coef <- mycoef(mydata, "mr2_t1dti")
mr2_t1dti_auc  <- my_auc(mydata)
mr2_t1dti_cv   <- my_cv(mydata)
mr2_t1dti_roc  <- my_roc(mydata)
###TEST IF AUCs ARE STATISTICALLY SIGNIFICANTLY DIFFERENT
#using paired DeLong's test ("paired" as each subject had multiple tests)

##CHOICE OF TIMEPOINT: MR1 vs MR2 (T1 & DTI)
#This requires the same patients in both ROC curves i.e. the 65 who had both MR1 and MR2
newdat1   <- mr1_gose3.df %>% filter(subject %in% mr2_gose3.df$subject)
newdat2   <- mr2_gose3.df %>% filter(subject %in% mr1_gose3.df$subject)
mydata    <- newdat1[,c("gose3", 
                        "age", 
                        "sex", 
                        "sum_abn_fa", 
                        "sum_abn_md", 
                        "sum_abn_both", 
                        "dev_vol")]
roc1      <- my_roc(mydata)
mydata    <- newdat2[,c("gose3", 
                        "age", 
                        "sex", 
                        "sum_abn_fa", 
                        "sum_abn_md", 
                        "sum_abn_both", 
                        "dev_vol")]
roc2      <- my_roc(mydata)
p_time    <- roc.test(roc1, roc2, paired = TRUE, method = "delong")$p.value %>% round(3)

#CHOICE OF SEQUENCE AT MR1: T1 or DTI vs both combined
p_t1vsboth  <- roc.test(mr1_t1_roc, mr1_t1dti_roc, paired = TRUE)$p.value %>% round(3)
p_dtivsboth <- roc.test(mr1_dti_roc, mr1_t1dti_roc, paired = TRUE)$p.value %>% round(3)

###COMBINE ALL ROC CURVES

roclist1 <- list(zero_roc,
                mr1_t1dti_roc,
                mr2_t1dti_roc)

p1    <- ggroc(roclist1, legacy.axes= T, size=0.6) +
           theme_classic(base_size = 6) + 
           labs(x="Specificity", y="Sensitivity")+
           theme(legend.justification=c(1,0), 
                 legend.position=c(0.98,0.08), 
                 legend.text=element_text(size=6),
                 legend.title=element_text(size=6),
                 legend.key.size = unit(0.6,"line"),
            plot.title = element_text(size = 9, hjust = 0.5),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
           geom_abline(intercept = 0, 
                       slope = 1, 
                       color = "darkgrey", 
                       linetype = "dashed")+
           scale_colour_manual(name  ="Time of imaging",
                                 values = c("snow3", "red3", "palegreen3"),
                                 breaks=c("1", "2", "3"),
                                 labels=c(paste0("no MRI", ": AUC ", zero_auc),
                                          paste0("within 72h", ": AUC " , mr1_t1dti_auc),
                                          paste0("at 2-3 weeks", ": AUC ", mr2_t1dti_auc))) +
           annotate("text", x = 0.78, y = 0.01, size = 2.1,
                    label = paste0("within 72h vs. at 2-3 weeks: p = ", p_time)) +
           ggtitle("A. Choice of timing of imaging")
    

roclist2 <- list(zero_roc,
                mr1_t1_roc,
                mr1_dti_roc,
                mr1_t1dti_roc)
              
p2      <- ggroc(roclist2, legacy.axes= T, size = 0.6) +
           theme_classic(base_size = 6) + 
           labs(x="Specificity", y="Sensitivity")+
           theme(legend.justification=c(1,0), 
                 legend.position=c(0.93,0.12), 
                 legend.text=element_text(size=6),
                 legend.title=element_text(size=6),
                 legend.key.size = unit(0.6,"line"),
                 plot.title = element_text(size = 9, hjust = 0.5),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
           geom_abline(intercept = 0, 
                       slope = 1, 
                       color = "darkgrey", 
                       linetype = "dashed")+
           scale_colour_manual(name  ="Sequences used",
                               values = c("snow3", "mediumpurple", "turquoise3", "red3"),
                                 breaks=c("1", "2", "3", "4"),
                                 labels=c(paste0("no MRI",   ": AUC ", zero_auc), 
                                          paste0("T1 only", ":  AUC ", mr1_t1_auc),
                                          paste0("DTI only", ": AUC ", mr1_dti_auc), 
                                          paste0("T1 & DTI", ": AUC ", mr1_t1dti_auc))) +
           annotate("text", x = 0.73, y = 0.03, size = 2.1,
                    label = paste0("T1 alone vs. DTI & T1: p = ", p_t1vsboth, 
                                   "\nDTI alone vs. T1& DTI: p = ", p_dtivsboth)) +
           ggtitle("B. Choice of imaging sequences")


myplot <- grid.arrange(p1, p2, ncol = 1, nrow = 2)

ggsave(filename = "Fig3.tiff", 
       plot = myplot,
       dpi = 300,
       width = 9,
       height = 17,
       units = "cm")
#Make thumbnail figure

p1    <- ggroc(roclist1, legacy.axes= T, size=0.6) +
           theme_classic(base_size = 6) + 
           labs(x="Specificity", y="Sensitivity")+
           theme(legend.justification=c(1,0), 
                 legend.position=c(0.98,0.08), 
                 legend.text=element_text(size=6),
                 legend.title=element_text(size=6),
                 legend.key.size = unit(0.6,"line"),
            plot.title = element_text(size = 9, hjust = 0.5),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
           geom_abline(intercept = 0, 
                       slope = 1, 
                       color = "darkgrey", 
                       linetype = "dashed")+
           scale_colour_manual(name  ="Time of imaging",
                                 values = c("snow3", "red3", "palegreen3"),
                                 breaks=c("1", "2", "3"),
                                 labels=c(paste0("no MRI", ": AUC ", zero_auc),
                                          paste0("within 72h", ": AUC " , mr1_t1dti_auc),
                                          paste0("at 2-3 weeks", ": AUC ", mr2_t1dti_auc))) +
           annotate("text", x = 0.72, y = 0.01, size = 2.1,
                    label = paste0("within 72h vs. at 2-3 weeks: p = ", p_time)) +
           ggtitle("A. Choice of timing of imaging")
    

roclist2 <- list(zero_roc,
                mr1_t1_roc,
                mr1_dti_roc,
                mr1_t1dti_roc)
              
p2      <- ggroc(roclist2, legacy.axes= T, size = 0.6) +
           theme_classic(base_size = 6) + 
           labs(x="Specificity", y="Sensitivity")+
           theme(legend.justification=c(1,0), 
                 legend.position=c(0.93,0.12), 
                 legend.text=element_text(size=6),
                 legend.title=element_text(size=6),
                 legend.key.size = unit(0.6,"line"),
                 plot.title = element_text(size = 9, hjust = 0.5),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
           geom_abline(intercept = 0, 
                       slope = 1, 
                       color = "darkgrey", 
                       linetype = "dashed")+
           scale_colour_manual(name  ="Sequences used",
                               values = c("snow3", "mediumpurple", "turquoise3", "red3"),
                                 breaks=c("1", "2", "3", "4"),
                                 labels=c(paste0("no MRI",   ": AUC ", zero_auc), 
                                          paste0("T1 only", ":  AUC ", mr1_t1_auc),
                                          paste0("DTI only", ": AUC ", mr1_dti_auc), 
                                          paste0("T1 & DTI", ": AUC ", mr1_t1dti_auc))) +
           annotate("text", x = 0.68, y = 0.03, size = 2.1,
                    label = paste0("T1 alone vs. DTI & T1: p = ", p_t1vsboth, 
                                   "\nDTI alone vs. T1& DTI: p = ", p_dtivsboth)) +
           ggtitle("B. Choice of imaging sequences")








myplot <- grid.arrange(p1, p2, ncol = 2, nrow = 1)

ggsave(filename = "Thumbnail.jpeg", 
       plot = myplot,
       dpi = 300,
       width = 15,
       height = 8,
       units = "cm")

###GET COEFFICIENTS FOR PREDICTION MODELS

all_coef <- rbind(zero_coef,
                  mr1_t1_coef,
                  mr1_dti_coef,
                  mr1_t1dti_coef,
                  mr2_t1_coef,
                  mr2_dti_coef,
                  mr2_t1dti_coef)
all_coef
write.csv(all_coef, row.names = FALSE, file = "GOSEpred_coef.csv")

###SUMMARISE PREDICTIVE PERFORMANCE OF MR1 AND MR2 MODELS

Timepoint   <-  c(" ",
                  "MR1",
                  "MR1",
                  "MR1",
                  "MR2",
                  "MR2",
                  "MR2")

Sequences   <-  c("No imaging",
                  "T1 only",
                  "DTI only",
                  "T1 & DTI",
                  "T1 only",
                  "DTI only",
                  "T1 & DTI")

N           <-  c(dim(zero_gose3.df)[1],
                  dim(mr1_gose3.df)[1],
                  dim(mr1_gose3.df)[1],
                  dim(mr1_gose3.df)[1],
                  dim(mr2_gose3.df)[1],
                  dim(mr2_gose3.df)[1],
                  dim(mr2_gose3.df)[1])

AUC         <-  c(zero_auc,
                  mr1_t1_auc,
                  mr1_dti_auc,
                  mr1_t1dti_auc,
                  mr2_t1_auc,
                  mr2_dti_auc,
                  mr2_t1dti_auc)

CV          <-  c(zero_cv,
                  mr1_t1_cv,
                  mr1_dti_cv,
                  mr1_t1dti_cv,
                  mr2_t1_cv,
                  mr2_dti_cv,
                  mr2_t1dti_cv)

AIC         <- c(zero_fit$aic,
                  mr1_t1_fit$aic,
                  mr1_dti_fit$aic,
                  mr1_t1dti_fit$aic,
                  mr2_t1_fit$aic,
                  mr2_dti_fit$aic,
                  mr2_t1dti_fit$aic)%>%
                round(0)

res         <-  data.frame(Timepoint, N, Sequences, AUC, CV, AIC)

#make a pretty version of table "res"
res %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "bordered"), full_width = F) %>%
  add_header_above(c("Predictive ability of MR1 and MR2"=6))

###ASSESS MODEL ASSUMPTIONS (LOGISTIC REGRESSION)

##LINEARITY
#using the Box-Tidwell test i.e. a predictor x is considered to violate the linearity assumption if the interaction term x:log(x) is significant

#note I am adding a very small number to avoid taking the log of zero
x = 0.000000001

fit_test1 <- glm(gose3 ~  sum_abn_fa * log(sum_abn_fa + x) + 
                          sum_abn_md * log(sum_abn_fa + x) + 
                          sum_abn_both * log(sum_abn_fa + x) + 
                          dev_vol * log(dev_vol + x)+
                          age*log(age) + 
                          sex, 
                          data = mr1_gose3.df, 
                          family = binomial(link='logit'))
summary(fit_test1)

fit_test2 <- glm(gose3 ~  sum_abn_fa * log(sum_abn_fa + x) + 
                          sum_abn_md * log(sum_abn_fa + x) + 
                          sum_abn_both * log(sum_abn_fa + x) + 
                          dev_vol * log(dev_vol + x)+
                          age*log(age) + 
                          sex, 
                          data = mr2_gose3.df, 
                          family = binomial(link='logit'))
summary(fit_test2)

#The model summaries show that linearity assumption is not violated



##MULTICOLINARITY
#using the variance inflation factor from function vif() in the package car 3.0-3 with a vif > 5 suggestive of collinerity 

fit_test1 <- glm(gose3 ~  sum_abn_fa  + 
                          sum_abn_md  + 
                          sum_abn_both  + 
                          dev_vol +
                          age + 
                          sex, 
                          data = mr1_gose3.df, 
                          family = binomial(link='logit'))
car::vif(fit_test1)

fit_test2 <- glm(gose3 ~  sum_abn_fa  + 
                          sum_abn_md  + 
                          sum_abn_both  + 
                          dev_vol +
                          age + 
                          sex, 
                          data = mr2_gose3.df, 
                          family = binomial(link='logit'))
car::vif(fit_test2)

#all vif < 5 (indeed <2) so there is no multicollinearity


##OUTLIERS
#using the standardised residual, where a datapoint with a s.r. greater than 3 is considered an outlier

model.data <- augment(fit_test1) %>% mutate(index = 1:n()) 
model.data %>% filter(abs(.std.resid) > 3)

model.data <- augment(fit_test2) %>% mutate(index = 1:n()) 
model.data %>% filter(abs(.std.resid) > 3)

#there are no overly influential datapoints i.e. outliers


######SENSITIVITY ANALYSIS

##CREATE BEST CASE AND WORST CASE DATA FOR MR1

#Single imputation for best case analysis
mr1_worst       <- mr1.df
mr1_worst$gose3 <- ifelse(is.na(mr1_worst$gose3) == TRUE, 
                          "bad", 
                          as.character(mr1_worst$gose3))
mr1_worst$gose3 <- factor(mr1_worst$gose3,
                          labels = unique(mr1_worst$gose3),
                          levels = unique(mr1_worst$gose3))
summary(mr1_worst$gose3)

#Single imputation for best case analysis
mr1_best        <- mr1.df
mr1_best$gose3  <- ifelse(is.na(mr1_best$gose3) == TRUE, 
                          "good", 
                          as.character(mr1_best$gose3))
mr1_best$gose3  <- factor(mr1_best$gose3,
                          labels = unique(mr1_best$gose3),
                          levels = unique(mr1_best$gose3))
summary(mr1_best$gose3)



##CREATE BEST CASE AND WORST CASE DATA FOR MR2

#Single imputation for worst case analysis
mr2_worst       <- mr2.df
mr2_worst$gose3 <- ifelse(is.na(mr2_worst$gose3) == TRUE, 
                          "bad", 
                          as.character(mr2_worst$gose3))
mr2_worst$gose3 <- factor(mr2_worst$gose3,
                          labels = unique(mr2_worst$gose3),
                          levels = unique(mr2_worst$gose3))
summary(mr2_worst$gose3)

#Single imputation for best case analysis
mr2_best        <- mr2.df
mr2_best$gose3  <- ifelse(is.na(mr2_best$gose3) == TRUE, 
                          "good", 
                          as.character(mr2_best$gose3))

mr2_best$gose3  <- factor(mr2_best$gose3,
                          labels = unique(mr2_best$gose3),
                          levels = unique(mr2_best$gose3))
summary(mr2_best$gose3)


##CREATE BEST CASE AND WORST CASE DATA FOR MODEL WITHOUT IMAGING INFO

#Single imputation for worst case analysis
zero_worst        <- zero.df
zero_worst$gose3  <- ifelse(is.na(zero_worst$gose3) == TRUE, 
                            "bad", 
                            as.character(zero_worst$gose3))
zero_worst$gose3  <- factor(zero_worst$gose3,
                            labels = unique(zero_worst$gose3),
                            levels = unique(zero_worst$gose3))
summary(zero_worst$gose3)

#Single imputation for best case analysis
zero_best         <- zero.df
zero_best$gose3   <- ifelse(is.na(zero_best$gose3) == TRUE, 
                            "good", 
                            as.character(zero_best$gose3))
zero_best$gose3   <- factor(zero_best$gose3,
                            labels = unique(zero_best$gose3),
                            levels = unique(zero_best$gose3))
summary(zero_best$gose3)


###WORST CASE ANALYSIS

#Zeromodel
mydata          <- zero_worst[,c("gose3", 
                                 "age", 
                                 "sex")]
wzero_fit       <- mymodel(mydata)
wzero_auc       <- my_auc(mydata)
wzero_cv        <- my_cv(mydata)

##MR1
#T1
mydata          <- mr1_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "dev_vol")]
wmr1_t1_fit     <- mymodel(mydata)
wmr1_t1_auc     <- my_auc(mydata)
wmr1_t1_cv      <- my_cv(mydata)

#DTI
mydata          <- mr1_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "sum_abn_fa", 
                                "sum_abn_md", 
                                "sum_abn_both")]
wmr1_dti_fit     <- mymodel(mydata)
wmr1_dti_auc    <- my_auc(mydata)
wmr1_dti_cv     <- my_cv(mydata)

#T1&DTI
mydata          <- mr1_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "sum_abn_fa", 
                                "sum_abn_md", 
                                "sum_abn_both", 
                                "dev_vol")]
wmr1_t1dti_fit  <- mymodel(mydata)
wmr1_t1dti_auc  <- my_auc(mydata)
wmr1_t1dti_cv   <- my_cv(mydata)


##MR2
#T1
mydata          <- mr2_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "dev_vol")]
wmr2_t1_fit     <- mymodel(mydata)
wmr2_t1_auc     <- my_auc(mydata)
wmr2_t1_cv      <- my_cv(mydata)

#DTI
mydata          <- mr2_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "sum_abn_fa", 
                                "sum_abn_md", 
                                "sum_abn_both")]
wmr2_dti_fit    <- mymodel(mydata)
wmr2_dti_auc    <- my_auc(mydata)
wmr2_dti_cv     <- my_cv(mydata)

#T1&DTI
mydata          <- mr2_worst[,c("gose3", 
                                "age", 
                                "sex", 
                                "sum_abn_fa", 
                                "sum_abn_md", 
                                "sum_abn_both", 
                                "dev_vol")]
wmr2_t1dti_fit  <- mymodel(mydata)
wmr2_t1dti_auc  <- my_auc(mydata)
wmr2_t1dti_cv   <- my_cv(mydata)


###BEST CASE ANALYSIS

#Zeromodel
mydata          <- zero_best[,c("gose3", 
                                "age", 
                                "sex")]
bzero_fit       <- mymodel(mydata)
bzero_auc       <- my_auc(mydata)
bzero_cv        <- my_cv(mydata)

##MR1
#T1
mydata          <- mr1_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "dev_vol")]
bmr1_t1_fit     <- mymodel(mydata)
bmr1_t1_auc     <- my_auc(mydata)
bmr1_t1_cv      <- my_cv(mydata)

#DTI
mydata          <- mr1_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "sum_abn_fa", 
                               "sum_abn_md", 
                               "sum_abn_both")]
bmr1_dti_fit    <- mymodel(mydata)
bmr1_dti_auc    <- my_auc(mydata)
bmr1_dti_cv     <- my_cv(mydata)

#T1&DTI
mydata          <- mr1_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "sum_abn_fa", 
                               "sum_abn_md", 
                               "sum_abn_both", 
                               "dev_vol")]
bmr1_t1dti_fit  <- mymodel(mydata)
bmr1_t1dti_auc  <- my_auc(mydata)
bmr1_t1dti_cv   <- my_cv(mydata)


##MR2
#T1
mydata          <- mr2_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "dev_vol")]
bmr2_t1_fit     <- mymodel(mydata)
bmr2_t1_auc     <- my_auc(mydata)
bmr2_t1_cv      <- my_cv(mydata)

#DTI
mydata          <- mr2_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "sum_abn_fa", 
                               "sum_abn_md", 
                               "sum_abn_both")]
bmr2_dti_fit    <- mymodel(mydata)
bmr2_dti_auc    <- my_auc(mydata)
bmr2_dti_cv     <- my_cv(mydata)

#T1&DTI
mydata          <- mr2_best[,c("gose3", 
                               "age", 
                               "sex", 
                               "sum_abn_fa", 
                               "sum_abn_md", 
                               "sum_abn_both", 
                               "dev_vol")]
bmr2_t1dti_fit  <- mymodel(mydata)
bmr2_t1dti_auc  <- my_auc(mydata)
bmr2_t1dti_cv   <- my_cv(mydata)


###SUMMARISE PREDICTIVE PERFORMANCE OF MODELS FROM SENSITIVITY ANALYSIS

Timepoint     <-  c(" ",
                    "MR1",
                    "MR1",
                    "MR1",
                    "MR2",
                    "MR2",
                    "MR2")

Sequences     <-  c("No imaging",
                    "T1 only",
                    "DTI only",
                    "T1 & DTI",
                    "T1 only",
                    "DTI only",
                    "T1 & DTI")

AUC           <-  c(zero_auc,
                    mr1_t1_auc,
                    mr1_dti_auc,
                    mr1_t1dti_auc,
                    mr2_t1_auc,
                    mr2_dti_auc,
                    mr2_t1dti_auc)

worst_AUC     <-  c(wzero_auc,
                    wmr1_t1_auc,
                    wmr1_dti_auc,
                    wmr1_t1dti_auc,
                    wmr2_t1_auc,
                    wmr2_dti_auc,
                    wmr2_t1dti_auc)

best_AUC      <-  c(bzero_auc,
                    bmr1_t1_auc,
                    bmr1_dti_auc,
                    bmr1_t1dti_auc,
                    bmr2_t1_auc,
                    bmr2_dti_auc,
                    bmr2_t1dti_auc)

CV            <-  c(zero_cv,
                    mr1_t1_cv,
                    mr1_dti_cv,
                    mr1_t1dti_cv,
                    mr2_t1_cv,
                    mr2_dti_cv,
                    mr2_t1dti_cv)

worst_CV      <-  c(wzero_cv,
                    wmr1_t1_cv,
                    wmr1_dti_cv,
                    wmr1_t1dti_cv,
                    wmr2_t1_cv,
                    wmr2_dti_cv,
                    wmr2_t1dti_cv)

best_CV       <-  c(bzero_cv,
                    bmr1_t1_cv,
                    bmr1_dti_cv,
                    bmr1_t1dti_cv,
                    bmr2_t1_cv,
                    bmr2_dti_cv,
                    bmr2_t1dti_cv)

AIC           <-  c(zero_fit$aic,
                    mr1_t1_fit$aic,
                    mr1_dti_fit$aic,
                    mr1_t1dti_fit$aic,
                    mr2_t1_fit$aic,
                    mr2_dti_fit$aic,
                    mr2_t1dti_fit$aic)%>%
                  round(0)

worst_AIC     <-  c(wzero_fit$aic,
                    wmr1_t1_fit$aic,
                    wmr1_dti_fit$aic,
                    wmr1_t1dti_fit$aic,
                    wmr2_t1_fit$aic,
                    wmr2_dti_fit$aic,
                    wmr2_t1dti_fit$aic)%>%
                    round(0)

best_AIC      <-  c(bzero_fit$aic,
                    bmr1_t1_fit$aic,
                    bmr1_dti_fit$aic,
                    bmr1_t1dti_fit$aic,
                    bmr2_t1_fit$aic,
                    bmr2_dti_fit$aic,
                    bmr2_t1dti_fit$aic)%>%
                    round(0)

res           <-  data.frame(Timepoint, 
                             Sequences, 
                             AUC, 
                             CV, 
                             AIC,
                             worst_AUC, 
                             worst_CV,
                             worst_AIC,
                             best_AUC, 
                             best_CV,
                             best_AIC)

colnames(res) <- c("Timepoint", 
                   "Sequences", 
                   "AUC (95% CI)", 
                   "CV (95% CI)", 
                   "AIC",
                   "AUC (95% CI)", 
                   "CV (95% CI)", 
                   "AIC",
                   "AUC (95% CI)", 
                   "CV (95% CI)",
                   "AIC")

#Make a pretty version of table "res"
res %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "bordered"), full_width = F) %>%
  add_header_above(c(" " = 2, 
                     "Complete case analysis" =3, 
                     "Worst-case scenario" = 3, 
                     "Best-case scenario" = 3))
