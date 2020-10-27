###INTRODUCTION

#This document contains the code to make tables of patient characteristics
#I will create an overall Table 1, that includes all 81 patients 
#i.e. the total cohort used in our paper

#Our paper contains multiple analyses
#Due to missing data only a subgroup of patients could be used in each analysis
#I compare characteristics of patients included and excluded in each analysis
#to better understand why data is missing
#in order to perform an informed sensitivity analysis

###DATA PREPARATION

#load required libraries
library(dplyr)

#load demographics data
data <- read.csv("Paper1_demog.csv")

#select characteristics I want to compare across subjects
data <- data %>% select(subject, 
                        Age, 
                        Sex, 
                        Edu, 
                        GCS, 
                        ISS, 
                        Mech = InjMechshort, # injury mechanism
                        CWM_ratio, #Cerebral white matter ratio MR2/MR1
                        phenotype, #DTI phenotype based on DTI data on MR1 and MR2
                        deltaRPQ,  #change in concussion symptoms between MR1 and MR2
                        GOSE,
                        Outcome,
                        Marshall,
                        Stratum = stratum,
                        MR1, #time to MR1 in hours
                        MR2, #time to MR2 in days
                        MR3) #time to MR3 in days

#I also need the DTI data, so that I know which patients are missing DTI data and have 
#load DTI data
dti <- read.csv("DTIvsoutcome_alltimepoints.csv")
dti_mr1 <- dti %>% filter(timepoint == "MR1") %>% select(subject, sum_fa_mr1 = sum_fa)
dti_mr2 <- dti %>% filter(timepoint == "MR2") %>% select(subject, sum_fa_mr2 = sum_fa)

#I kept the variable sum_fa to act as an indicator of
#whether of not DTI data was missing for this scan 
data <- merge(data, dti_mr1, by = "subject", all.x = TRUE)
data <- merge(data, dti_mr2, by = "subject", all.x = TRUE)


###MAKING LISTS OF PATIENTS INCLUDED AND EXCLUDED FOR EACH ANALYSIS

#I need to compare included and excluded patients for each analysis.
#Refer to the flowchart in eFigure 1 to see why patients were included in various anaylses.

#I will name the analysis A to E
  #A - analysis of radiologically visible lesions and ROI volumes (all 81 patients)
  #B - analysis of DTI evolution (63 patients)
  #C - analysis of symptom evolution (53 patients)
  #D - outcome prediction MR1 (65 patients)
  #E - outcome prediction MR2 (68 patients)


#define which patients belong to which analysis

#structural analysis
 A <- data$subject %>% #keep all subjects
  unlist
length(A)

#DTI evolution
B <- data %>% 
  filter(
    is.na(data$phenotype) == FALSE) %>% #has diffusion data for both timepoints
  select(subject) %>% 
  unlist
length(B)

#symptom evolution
C <- data %>% 
  filter(
    is.na(data$deltaRPQ) == FALSE  & #has RPQ data for both timepoints
      is.na(data$phenotype) == FALSE) %>% #has diffusion data for both timepoints
  select(subject) %>% 
  unlist
length(C)

#Outcome MR1
D <- data %>% 
  filter(
      is.na(data$GOSE) == FALSE  & #has outcome data
      is.na(data$sum_fa_mr1) == FALSE & #has diffusion data
      !subject %in% c("1yhA225", "3EIi342")#has control data
      ) %>% 
  select(subject) %>% 
  unlist
length(D)

#Outcome MR2
E <- data %>% 
  filter(
      is.na(data$GOSE) == FALSE  & #has outcome data
      is.na(data$sum_fa_mr2) == FALSE & #has diffusion data
      !subject %in% c("1yhA225", "3EIi342")#has control data
      ) %>% 
  select(subject) %>% 
  unlist
length(E)
#Add columns to data that indicate whether or not a subject belongs to each of the analyses
data$A <- ifelse(data$subject %in% A, "0", "1") %>% as.factor
data$B <- ifelse(data$subject %in% B, "0", "1") %>% as.factor
data$C <- ifelse(data$subject %in% C, "0", "1") %>% as.factor
data$D <- ifelse(data$subject %in% D, "0", "1") %>% as.factor
data$E <- ifelse(data$subject %in% E, "0", "1") %>% as.factor


### MAKING TABLES COMPARING INCLUDED AND EXCLUDED PATIENTS FOR EACH ANALYSIS

#load package
library(table1)

##define how table should be displayed using render functions

#set up display of continuous variables
my.render.cont <- function(x) 
  {
  #display 2 sig figs:
    with(stats.apply.rounding(stats.default(x), digits=2), 
  #in first column display variable name and "Median (IQR)":
    c("", "Median (Q1-Q3)" = 
  #in second column calculate and display:
        sprintf("%s (%s %s %s)", 
                MEDIAN,  
                round(quantile(x, 0.25, type = 7, na.rm = T),1) , 
                "-", 
                round(quantile(x, 0.75, type = 7, na.rm = T),1))))   #calculations requested
  }


#set up display of categorical variables
my.render.cat <- function(x) 
  {
  #in first column display variable name:
    c("", 
  #in second column calculate and display:
      sapply(stats.default(x), function(y) with(y, sprintf("%d (%0.0f %%)", 
                                                           FREQ, 
                                                           round(PCT, 0)))))
  }

#Set up statistical tests used to compare included with excluded patients
  #note that Z in data$Z refers to the respective analysis (A, B, C, D or E)
  #column data$Z will include three levels: 0, 1, 2. 
  #0 and 1 include will refer to whether the patient was included or excluded from analysis Z
  #the third level ("2") is a hack to display the p-value in the third column 
  #Please see eTable 2 to understand what layout we are aiming for
rndr <- function(x, name, ...) 
  {
    if (length(x) == 0) 
      {
        y <- data[[name]]
        s <- rep("", length(render.default(x=y, name=name, ...)))
        #use wilcox test for numeric variables
        if (is.numeric(y)) 
          {
            p <- wilcox.test(y ~ data$Z)$p.value
        #use Fisher exact test for categorical variables
          } 
        else 
          {
          p <- fisher.test(table(y, droplevels(data$Z)))$p.value
          }
        s[2] <- sub("<", "&lt;", format.pval(p, digits=3, eps=0.001))
        s
        } 
    else 
        {
        render.default(x=x, name=name, ...)
        }
  }
        
#set up display of the three
rndr.strat <- function(label, n, ...) 
  {
  ifelse(n==0, label, render.strat.default(label, n, ...))
  }

#Summary Table analysis A

tableA <- table1(~ Age + 
         Sex + 
         Edu +  
         Mech + 
         GCS + 
         ISS + 
         Stratum + 
         Marshall + 
         CWM_ratio + 
         phenotype + 
         deltaRPQ + 
         GOSE + 
         MR1 + 
         MR2 + 
         MR3, 
         data = data, 
         droplevels=F,
         render.continuous = my.render.cont, 
         render.categorical = my.render.cat)

#Create functions to make tables for analyses B-E

make_my_column <- function(column)
  {
  factor(column, levels=c(0, 1, 2), labels=c("Included", "Excluded", "P-value"))
  }

make_my_table <- function() 
  {
    tab <- table1(~ Age + 
                    Sex + 
                    Edu +  
                    Mech + 
                    GCS + 
                    ISS + 
                    Stratum + 
                    Marshall + 
                    CWM_ratio + 
                    phenotype + 
                    deltaRPQ + 
                    GOSE + 
                    MR1 + 
                    MR2 + 
                    MR3|Z, 
                    data = data, 
                    droplevels=F, 
                    render=rndr, 
                    render.continuous = my.render.cont, 
                    render.categorical = my.render.cat,
                    render.strat=rndr.strat, 
                    overall=F)
    
    return(tab)
  }
data$Z <- make_my_column(data$B)
tableB <- make_my_table()
data$Z <- make_my_column(data$B)
tableC <- make_my_table()
data$Z <- make_my_column(data$B)
tableD <- make_my_table()
data$Z <- make_my_column(data$B)
tableE <- make_my_table()

#Print the comparison tables of included and excluded patients
#for each analysis A-E
tableA
tableB
tableC
tableD
tableE
#Note I cannot cbind these tables as they are written directly as html outputs
#see package description
#So they have to manually combined through copy and paste

##ADJUSTING FOR FALSE DISCOVERY RATE

#load p-values
pdata <- read.csv("comparemissingpt_20200608.csv")
#this table is the result of cbind-ing tables A-E together manually in excel.

#combine the 4 columns of p-values into a vector
pdata <- pdata %>% select(PB, PC, PD, PE) %>% na.omit 
pdata <- sapply(pdata, as.numeric) %>% as.vector
my_p_values <- pdata

#make a function that calculates the maximal p-value which is still
#statistically significant according to the False Discovery Rate
source("fdr_function.R")

#Identify the maximum p-value among my p-values that is still significant
#accordingto a 5% false discovery rate threshold
maxp <- fdr_maxp(pdata, 0.05)
maxp

###Making Table 1

#Update the render function for continous variables
#to also display the range (min-max)
my.render.cont <- function(x) 
  {
  with(stats.apply.rounding(stats.default(x), digits=2), 
  c("", "Median [Q1-Q3] (Range)" =
    sprintf("%s [%s %s %s] (%s %s %s)", 
            MEDIAN,  
            round(quantile(x, 0.25, type = 7, na.rm = T),1) , 
            "-", 
            round(quantile(x, 0.75, type = 7, na.rm = T),1),
            MIN, 
            "-", 
            MAX)))      
}

#How variable names should be displayed
label(data$Age)       <- "Age in years"
label(data$Sex)       <- "Male sex"
label(data$GCS)       <- "Glasgow Coma Score"
label(data$Mech)      <- "Mechanism of injury"
label(data$ISS)       <- "Injury Severity Score"
label(data$Stratum)   <- "Stratum"
label(data$Outcome)   <- "Outcome at 3 months"
label(data$Marshall)  <- "Marshall score (pre MR1)"
label(data$MR1)       <- "Time to MR1 in hours"
label(data$MR2)       <- "Time to MR2 in days"
label(data$MR3)       <- "Time to MR3 in days"

#actually make the table
table1(~ Age + 
         Sex +  
         Mech + 
         GCS + 
         ISS + 
         Stratum + 
         Outcome + 
         Marshall + 
         MR1 + 
         MR2 + 
         MR3, 
        data = data, 
        droplevels=F,
        render.continuous = my.render.cont, 
        render.categorical = my.render.cat)

#Get test statistic and degrees of freedom for selected comparisons

hasRPQ <- data %>% filter(C == 1) %>% select(ISS) %>% na.omit() %>% unlist()
noRPQ  <- data %>% filter(C == 0) %>% select(ISS) %>% na.omit() %>% unlist()
wilcox.test(hasRPQ, noRPQ, alternative = "two.sided")
