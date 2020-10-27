# load packages and source function
library(dplyr)
source("fdr_function.R")
library(ggplot2)

###1. DATA PREPARATION

# The following files are in the long format
# the rows are the participants
  # there are up to three rows per participants (for MR1, MR2, MR3)
# the columns are the 72 white matter tracts
  # depending on the dataframe, the column contain either the FA or MD
  # of that particular tract
  # in addition there are also some identifier columns
  # e.g. patient ID, scan ID, scanner on which the scan was conducted
# Note the column "timepoint" has 4 levels: control, MR1, MR2, MR3
  # thereby denoting whether a participant is a control or a patient
  # and if they are a patient whether this is their first, second or third scan

md    <- read.csv("md_allsubjects_alltracts_nooutliers.csv")
fa    <- read.csv("fa_allsubjects_alltracts_nooutliers.csv")

#The following file contains the list of tracts that changed significantly 
  #within patients between MR1 and MR2
  #with respect to diffusion parameters
  #plus the corpus callosum
#there is one row for each tract
#there is a column called FDR
  #which says whether the tract had changed significantly
  #which, by definition, is the case for all tracts except the corpus callosum
mytracts <- read.csv("mytracts.csv")
changedtracts   <- mytracts %>% filter(FDR == "significant")

# In the following dataframe the rows are the patients
# there are various columns
# we will only use the column "phenotype" which denotes the imaging phenotype
  # that each patient belongs to
  # according to the k-means clustering done earlier
  # there are three phenotypes:
    # "Progressive injury" alias "prog"
    # "Minimal change" alias "No DTI change" alias "static"
    # "Pseudonormalisation" alias "pseud"

pheno <- read.csv("DTI_phenotypes.csv")

##1.1 PREPARE FA DATA

# Adding phenotype to FA data

pheno         <- pheno %>% select(subject, phenotype)

fa_pheno      <- merge(fa, pheno, 
                     by = "subject", 
                     all.x = TRUE)

# make a new column "pc" that says whether a participant s a control or a patient
# and if the latter, which phenotype they belong to
# this will make "control" the reference level which each of the three
# phenotypes can be compared to 
fa_pheno$pc   <- ifelse(fa_pheno$group == "control", 
                        "control", 
                        as.character(fa_pheno$phenotype))
fa_pheno$pc   <- factor(fa_pheno$pc, 
                        levels = c("control",
                                   "No DTI change",
                                   "Progressive injury" , 
                                   "Pseudonormalisation"),
                        labels = c("control",
                                   "static",
                                   "prog", 
                                   "pseud"))

summary(fa_pheno$pc)

fa_pheno <- fa_pheno %>% filter(is.na(fa_pheno$pc) == FALSE) 

#Select MR1 data
fa_pheno_1 <- fa_pheno %>% 
              filter(timepoint %in% c("control", "MR1")) %>%
              select(subject, 
                    pc, 
                    age, 
                    sex, 
                    scan, 
                    scanner, 
                    timepoint, 
                    as.character(changedtracts$Number))

#Select MR2 data
fa_pheno_2 <- fa_pheno %>% 
              filter(timepoint %in% c("control", "MR2")) %>%
              select(subject, 
                    pc, 
                    age, 
                    sex, 
                    scan, 
                    scanner, 
                    timepoint, 
                    as.character(changedtracts$Number))



##1.2 PREPARE MD DATA

##Prepare MD data
md_pheno      <- merge(md, pheno, 
                       by = "subject", 
                       all.x = TRUE)

md_pheno$pc   <- ifelse(md_pheno$group == "control", 
                        "control", 
                        as.character(md_pheno$phenotype))

md_pheno$pc   <- factor(md_pheno$pc, 
                        levels = c("control",
                                   "No DTI change",
                                   "Progressive injury" , 
                                   "Pseudonormalisation"),
                        labels = c("control",
                                   "static",
                                   "prog", 
                                   "pseud"))

summary(md_pheno$pc)

md_pheno      <-  md_pheno %>% filter(is.na(md_pheno$pc) == FALSE)

#select MR1 data
md_pheno_1    <-  md_pheno %>% 
                    filter(timepoint %in% c("control", "MR1")) %>%
                    select(subject, 
                          pc, 
                          age, 
                          sex, 
                          scan, 
                          scanner, 
                          timepoint, 
                          as.character(changedtracts$Number))



#select MR2 data
md_pheno_2      <-  md_pheno %>% 
                    filter(timepoint %in% c("control", "MR2")) %>%
                    select(subject, 
                          pc, 
                          age, 
                          sex, 
                          scan, 
                          scanner, 
                          timepoint, 
                          as.character(changedtracts$Number))

##2. CHECK IF PHENOTYPE (VARIABLE PC) IS RELATED TO FA (OR MD RESPECTIVELY)

#write a function to do this, so it can be reused

relate_pheno_to_dti <- function (Data, Metric)
  {
  #initialise dataframe to collect results in
  df              <- data.frame(matrix(ncol = 3, nrow = 0))
  colnames(df)    <- c("Tract", "Data", "a")
  df              <- data.frame(Tract   = character(), 
                                Metric  = character(),
                                a       = numeric())
  
  
  # This code will loop through all tracts
  # For each tract compare the model with and without the variable "group"
  # collect results in a vector
  # and bind all vectors together into a results dataframe
  # Only "group" is significant (according to anova p-value)
  # will I look at the p-vaue for the individual levels of group (see later)
  for (i in 8:20) 
    {
    Tract <- colnames(Data)[[i]]
    
    library(lmerTest)

    fit   <- lmer(Data[[i]] ~ pc + age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
    fit0  <- lmer(Data[[i]] ~  age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
    a   <- anova(fit0, fit)$`Pr(>Chisq)`[2] %>% round(3)

    vec <- c(Tract, Metric, a)

    df <- berryFunctions::addRows(df, 1, values = vec)
    }
  return(df)
  }

##3. SINCE 2. IS THE TRUE, COMPARE EACH PHENOTYPE's FA WITH CONTROLS

#write a function to do this so it can be reused

compare_phenos_with_ctrl <- function(Data, Metric, x = 1)
  {
  #initialise dataframe to collect results in
  df            <- data.frame(matrix(ncol = 11, nrow = 0))
  
  colnames(df)  <- c("Tract", 
                     "Metric", 
                     "Coef_static", "SE_static", "P_static", 
                     "Coef_prog",   "SE_prog",   "P_prog", 
                     "Coef_pseud",  "SE_pseud",  "P_pseud")

  df            <- data.frame(Tract       = character(), 
                              Metric      = character(),
                
                              Coef_static = numeric(),
                              SE_static   = numeric(),
                              P_static    = numeric(),
                
                              Coef_prog   = numeric(),
                              SE_prog     = numeric(),
                              P_prog      = numeric(),
                
                              Coef_pseud  = numeric(),
                              SE_pseud    = numeric(),
                              P_pseud     = numeric())
  

  ##Loop through all 13 tracts
  #For each tract extract the coefficients, standard error and p-value
  #of each phenotype as compared with the reference level (i.e. control participants)
  for (i in 8:20) 
    {
    Tract       <- colnames(Data)[[i]]

    library(lmerTest)

    fit         <- lmer(Data[[i]]*x ~ pc + age + sex + (1| scanner), 
                        data = Data, 
                        REML = FALSE, 
                        control = lmerControl(optimizer ="Nelder_Mead"))

    Coef_static <- summary(fit)$coefficients[2,1] %>% round(4) 
    SE_static   <- summary(fit)$coefficients[2,2] %>% round(4)
    P_static    <- summary(fit)$coefficients[2,5] %>% round(4)

    Coef_prog   <- summary(fit)$coefficients[3,1] %>% round(4)
    SE_prog     <- summary(fit)$coefficients[3,2] %>% round(4)
    P_prog      <- summary(fit)$coefficients[3,5] %>% round(4)

    Coef_pseud  <- summary(fit)$coefficients[4,1] %>% round(4)
    SE_pseud    <- summary(fit)$coefficients[4,2] %>% round(4)
    P_pseud     <- summary(fit)$coefficients[4,5] %>% round(4)
 
 
    vec         <- c(Tract, 
                    Metric, 
                    Coef_static, SE_static, P_static, 
                    Coef_prog, SE_prog, P_prog, 
                    Coef_pseud, SE_pseud, P_pseud)

    df          <- berryFunctions::addRows(df, 1, values = vec)
    }

  df[,3:11] <- apply(df[,3:11], 2, as.numeric)

  return(df)
  }
  ?berryFunctions::addRows
##4. FOR EACH PHENOTYPE HOW MANY TRACTS DIFFER **SIGNIFICANTLY** FROM CONTROLS

# The following function generates a results table with one row per phenotype
# For each phenotype is says what metric the table refers to (FA or MD)
# How many of the 13 tracts are abnormal compared to controls
# The mean coefficients
  # if the coefficient is negative, then the phenotype has lower FA (or MD) values
    # compared to controls
  # if the coefficient is positive, then the phenotype has higher FA (or MD) values
    # compared to controls
# The standard error of the mean coefficient

count_abn_tracts <- function(mydata = df, Metric)
  {
  dat               <- data.frame(matrix(ncol = 5, nrow = 3))
  colnames(dat)     <- c("Pheno", "Metric", "Num_abnormal", "Mean_Coef", "Mean_SE")
  dat$Pheno         <- c("static", "prog", "pseud")
  dat$Metric        <- rep(Metric, times = 3)
  
  maxp              <- fdr_maxp(c(mydata$P_static,
                                  mydata$P_prog, 
                                  mydata$P_pseud), 
                                0.05)
  
  dat$Num_abnormal  <- c( length(which(mydata$P_static < maxp)),
                          length(which(mydata$P_prog < maxp)),
                          length(which(mydata$P_pseud < maxp)))
  
  dat$Mean_Coef     <- c( mean(mydata$Coef_static),
                          mean(mydata$Coef_prog),
                          mean(mydata$Coef_pseud)) %>% 
                        round(4)
  
  dat$Mean_SE       <- c( mean(mydata$SE_static),
                          mean(mydata$SE_prog),
                          mean(mydata$SE_pseud)) %>% 
                        round(4)
  
  return(dat)
  }


###5. RUN STEPS 2-4 ON ALL four datasets

#MD at MR1
  df        <- relate_pheno_to_dti(Data = md_pheno_1, Metric = "md_pheno_1")
  df 
  #group is significant so I can proceed
  df        <- compare_phenos_with_ctrl(Data = md_pheno_1, Metric = "md_pheno_1", x = 1000) 
  df        <- count_abn_tracts(Metric = "md")
  dat_md_1  <- df


#MD at MR2
  df        <- relate_pheno_to_dti(Data = md_pheno_2, Metric = "md_pheno_2")
  df 
  #group is significant so I can proceed
  df        <- compare_phenos_with_ctrl(Data = md_pheno_2, Metric = "md_pheno_2", x = 1000)
  df        <- count_abn_tracts(Metric = "md")
  dat_md_2  <- df


#FA at MR1
  df        <- relate_pheno_to_dti(Data = fa_pheno_1, Metric = "fa_pheno_1")
  df 
  #group is significant so I can proceed
  df        <- compare_phenos_with_ctrl(Data = fa_pheno_1, Metric = "fa_pheno_1", x = 1)
  df        <- count_abn_tracts(Metric = "fa")
  dat_fa_1  <- df

#FA at MR2
  df        <- relate_pheno_to_dti(Data = fa_pheno_2, Metric = "fa_pheno_2")
  df 
  #group is significant so I can proceed
  df        <- compare_phenos_with_ctrl(Data = fa_pheno_2, Metric = "fa_pheno_2", x = 1)
  df        <- count_abn_tracts(Metric = "fa")
  dat_fa_2  <- df

###6. SAVE AND VISUALISE RESULTS FOR MR1

dat_combi         <- rbind(dat_fa_1, dat_md_1) %>% as.data.frame()
dat_combi
write.csv(dat_combi, file = "DTIPheno_MR1.csv", row.names = FALSE)

#for the bar plot, I still want to see a slim line of colour, 
#even if the number of abnormal tracts is zero
dat_combi$Num_abnormal <- ifelse(dat_combi$Num_abnormal < 1, 
                                 0.2, 
                                 dat_combi$Num_abnormal)

#Barplot of Number of abnormal tracts
  #grouped by phenotype (as in published paper)
ggplot(dat_combi, aes(x = Pheno, y = Num_abnormal)) +
        geom_bar(aes(fill = Metric), position = "dodge", stat="identity") +
        ggtitle("DTI abnormalities at MR1") +
        theme_bw()
  
  #grouped by metric
ggplot(dat_combi, aes(x = Metric, y = Num_abnormal)) +
        geom_bar(aes(fill = Pheno), position = "dodge", stat="identity") +
        ggtitle("DTI abnormalities at MR1") +
        theme_bw()

###7. SAVE AND VISUALISE RESULTS FOR MR2

dat_combi         <- rbind(dat_fa_2, dat_md_2) %>% as.data.frame()
dat_combi
write.csv(dat_combi, file = "DTIPheno_MR2.csv", row.names = FALSE)

#for the bar plot, I still want to see a slim line of colour, 
#even if the number of abnormal tracts is zero
dat_combi$Num_abnormal <- ifelse(dat_combi$Num_abnormal < 1, 
                                 0.2, 
                                 dat_combi$Num_abnormal)

#Barplot of Number of abnormal tracts
  #grouped by phenotype (as in published paper)
ggplot(dat_combi, aes(x = Pheno, y = Num_abnormal)) +
        geom_bar(aes(fill = Metric), position = "dodge", stat="identity") +
        ggtitle("DTI abnormalities at MR2") +
        theme_bw()
  
  #grouped by metric
ggplot(dat_combi, aes(x = Metric, y = Num_abnormal)) +
        geom_bar(aes(fill = Pheno), position = "dodge", stat="identity") +
        ggtitle("DTI abnormalities at MR2") +
        theme_bw()
###8. CHECK MODEL ASSUMPTIONS

#8.1 Normality MD
par(mfrow=c(3,5))
for (i in 8:20) 
  {
    Data   <-  md_pheno_1 
    Tract  <-  colnames(Data)[[i]]
    x      <-  1000
    
    library(lmerTest)

    fit     <- lmer(Data[[i]]*x ~ pc + age + sex + (1| scanner), 
                    data = Data, 
                    REML = FALSE, 
                    control = lmerControl(optimizer ="Nelder_Mead"))
    
    #Normality assumption
    qqnorm(residuals(fit))
    qqline(residuals(fit))
  }


#8.2 homoscedasticity MD
for (i in 8:20) 
  {
    Data    <- md_pheno_1 
    Tract   <- colnames(Data)[[i]]
    x       <- 1000
    
    library(lmerTest)

    fit     <- lmer(Data[[i]]*x ~ pc + age + sex + (1| scanner), 
                    data = Data, 
                    REML = FALSE, 
                    control = lmerControl(optimizer ="Nelder_Mead"))
    
    plot(residuals(fit))
  }


#8.3 Normality FA
par(mfrow=c(3,5))
for (i in 8:20) 
  {
    Data    <- fa_pheno_1 
    Tract   <- colnames(Data)[[i]]
    x       <- 1
    
    library(lmerTest)

    fit     <- lmer(Data[[i]]*x ~ pc + age + sex + (1| scanner), 
                    data = Data, 
                    REML = FALSE, 
                    control = lmerControl(optimizer ="Nelder_Mead"))
    
    #Normality assumption
    qqnorm(residuals(fit))
    qqline(residuals(fit))
  }


#8.4 homoscedasticity FA
for (i in 8:20) 
  {
    Data      <-  fa_pheno_1 
    Tract     <-  colnames(Data)[[i]]
    x         <-  1
    
    library(lmerTest)

    fit       <-  lmer(Data[[i]]*x ~ pc + age + sex + (1| scanner), 
                       data = Data, 
                       REML = FALSE, 
                       control = lmerControl(optimizer ="Nelder_Mead"))
    
    plot(residuals(fit))
  }


