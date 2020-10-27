#load required packages
library(dplyr)
library(tidyr)
library(knitr)
library(kableExtra)
source("fdr_function.R")

###DATA PREPARATION

#Loading list of radiology reports
reports <- read.csv("Paper1_allreportsfor81pts_20190701_vfjn2.csv")

#The dataframe "reports" is in the long format.
#The rows are the patients.
#For each patient there are three rows containing data for the CT, MR1 and MR2 respectively
#The columns are the types of lesions e.g. Subdural heamatoma, Intraventrucular haemorraghe etc
#they are reported as either "present" or "absent"
#There are also some columns containing the identifier, time and site at which the scan was done


##Cleaning the levels of categorical variables

#Combining all types of SDH into one
reports$SubduralHematomaCombined <- 
  ifelse(reports$SubduralHematomaAcute == "present", "present",
    ifelse(reports$SubduralCollectionMixedDensity == "present", "present",
      ifelse(reports$SubduralHematomaSubacuteChronic == "present", "present",
             "absent")))

reports <- reports %>% 
  select(-SubduralHematomaAcute, 
         -SubduralCollectionMixedDensity, 
         -SubduralHematomaSubacuteChronic)


#Combing all types of CT, as I know they all occurred within 24h
reports$ImagingTimepoint <- 
  ifelse(reports$ImagingTimepoint %in% c("CT Early", 
                                         "CT Followup",  
                                         "CT Post-Op"), 
    "CT Early", 
    as.character(reports$ImagingTimepoint))

reports$timepoint <- factor(reports$ImagingTimepoint, 
                            levels=c("MR Early", "CT Early", "MR 2 weeks"),
                            labels=c("MR1", "CT", "MR2"))

#One patient has been included for two CTs, I will remove the second CT
reports <- reports %>% 
  filter(ExperimentId != "CTBI_E28504")
summary(reports$timepoint)

#reorder columns
reports <- reports %>% select(ExperimentId, SiteCode, entity_id, 
                              ExperimentDate, ExperimentTime, timepoint,
                              AnyIntracranTraumaticAbnormality,
                              MassLesion,
                              MidlineShift,
                              CisternalCompression,
                              Contusion,
                              TAI,
                              EpiduralHematoma,
                              SubduralHematomaCombined,
                              TraumaticSubarachnoidHemorrhage,
                              SkullFracture,
                              IntraventricularHemorrhage
                              )


#for manually entered data (lesion reports)
#remove trailing spaces
reports[, 7:17] <- apply(reports[, 7:17], 2, gdata::trim)
#turn into factor variables
reports <- mutate_at(reports, 
                     vars(c(AnyIntracranTraumaticAbnormality:IntraventricularHemorrhage)),
                     as.factor)
#drop unused levels
reports <- droplevels(reports)


#Renaming columns to something that looks good in a table
colnames(reports) <- c("ExperimentId", 
                       "SiteCode", 
                       "entity_id", 
                       "ExperimentDate",
                       "ExperimentTime",
                       "timepoint",
                       "Any Abnormality",
                              "Mass > 25cc",
                              "Midline shift",
                              "Cisternal compression",
                              "Contusion",
                              "Traumatic axonal injury",
                              "Epidural haemorrhage",
                              "Subdural haemorrhage",
                              "Subarachnoid haemorrhage",
                              "Skull fracture",
                              "Intraventricular haemorrhage"
                              )

lesions <- colnames(reports[,7:17])
#initialise a dataframe to store results in
df <- data.frame(CTposMRneg = numeric(),
                    CTnegMRpos = numeric(),
                    CTposMRpos = numeric(),
                    CTnegMRneg = numeric(),
                    p.value    = numeric(),
                    lesion     = character())

i=1
  #for each iteration label the relevant column "lesion" 
  #so that it can be plugged into the subsequent code
  data <- reports %>% rename("lesion" = lesions[i])
  
  #This code initially creates a table for a particular type of lesion
  #it has only one row per patient
  #because the timepoint is now spread across three columns (CT, MR1, MR2)
  #each column says whether a lesion is "present" or "absent"
  #I then select the columns for CT and MR1 and turn the data into a 2x2 frequency table
  tab <- data %>% 
    select(lesion, timepoint, entity_id) %>% 
    spread(timepoint, lesion) %>% 
    select(CT, MR1) %>%
    table()
  mcnemar.test(tab)$statistic
  mcnemar.test(tab)$parameter

###COMPARING CT WITH MR1 

#Due to space limits we had to omit results for this analysis from the published paper, 
#but an analogous analysis comparing MR1 and MR2 is included in the paper, 
#the code for MR1 vs MR2 is identical to the one for CT1 vs MR1 
#and is displayed further down in this document

#Remember the patients are the rows
#there are 3 rows per patient (CT, MR1 and MR2 data)
#the three types of data per patient are identifiable through the column "timepoint"
#The columns are the different types of lesions (e.g. SDH, EDH, skull fracture etc)
#they are recorded either as "present" or "absent"
#For each lesion I want to know whether it differs between timepoints within the same patient
#i.e. whether the lesion is present in one timepoint and not the other or vice versa

#First I will compare CT with MR1

#make a vector of lesion names
lesions <- colnames(reports[,7:17])
#initialise a dataframe to store results in
df <- data.frame(CTposMRneg = numeric(),
                    CTnegMRpos = numeric(),
                    CTposMRpos = numeric(),
                    CTnegMRneg = numeric(),
                    p.value    = numeric(),
                    lesion     = character())

for(i in seq_along(lesions))
  {
  #for each iteration label the relevant column "lesion" 
  #so that it can be plugged into the subsequent code
  data <- reports %>% rename("lesion" = lesions[i])
  
  #This code initially creates a table for a particular type of lesion
  #it has only one row per patient
  #because the timepoint is now spread across three columns (CT, MR1, MR2)
  #each column says whether a lesion is "present" or "absent"
  #I then select the columns for CT and MR1 and turn the data into a 2x2 frequency table
  tab <- data %>% 
    select(lesion, timepoint, entity_id) %>% 
    spread(timepoint, lesion) %>% 
    select(CT, MR1) %>%
    table()
  
  #Make a row of data for each lesion type
  #This row contains the observed frequenciy of a lesion being present or absent
  #at each of the two timepoints
  #as well as the p-value from McNemar's test for paired categorical data
  mynum <- c(
    tab[2,1],
    tab[1,2],
    tab[2,2],
    tab[1,1],
    mcnemar.test(tab)$parameter,
    mcnemar.test(tab)$statistic %>% round(3),
    mcnemar.test(tab)$p.value %>% round(3)
    ) %>%
      round(4) %>%
      matrix(nrow=1, ncol=7) %>%
      as.data.frame() %>%
      cbind(lesions[i])
  
  #store the row of data inside the previously initialised dataframe
  df <- rbind(df, mynum)
  }


#df now contains the comparison data for each lesion type at CT vs MR1
#In the code below I am merely beautifying the table for publication

#rename columns
df <- df %>% select(Abnormality = "lesions[i]",
                    CTposMRneg = "V1",
                    CTnegMRpos = "V2",
                    CTposMRpos = "V3",
                    CTnegMRneg = "V4",
                    df         = "V5",
                    chi        = "V6",
                    p.value    = "V7")

#add columns denoting the total number of CTs positive and negative for that lesion
#irrespective of whether that lesions was seen on MR1
df$CTpos <- rowSums(df[,c("CTposMRneg", "CTposMRpos")])
df$CTneg <- rowSums(df[,c("CTnegMRneg", "CTnegMRpos")])
df <- df%>% select(Abnormality, CTpos, CTposMRpos, CTposMRneg, CTneg, CTnegMRneg, CTnegMRpos, df, chi, p.value)

#calculate percentage data in addition to the absolute numbers already calculated
#I will use the total number of CTs positive (or negative) for the respective lesion type as 100%
df$CTposP <- "(100%)"
df$CTposMRposP <- (df$CTposMRpos/df$CTpos*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$CTposMRnegP <- (df$CTposMRneg/df$CTpos*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$CTnegP <- "(100%)"
df$CTnegMRnegP <- (df$CTnegMRneg/df$CTneg*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$CTnegMRposP <- (df$CTnegMRpos/df$CTneg*100) %>% round(0) %>% paste0("(", . , "%", ")")

#display absolute numbers and percentage data in the same cell of the table
df <- df %>% unite("CTpos", c(CTpos, CTposP), sep = " ")
df <- df %>% unite("CTneg", c(CTneg, CTnegP), sep = " ")
df <- df %>% unite("CTposMRpos", c(CTposMRpos, CTposMRposP), sep = " ")
df <- df %>% unite("CTposMRneg", c(CTposMRneg, CTposMRnegP), sep = " ")
df <- df %>% unite("CTnegMRneg", c(CTnegMRneg, CTnegMRnegP), sep = " ")
df <- df %>% unite("CTnegMRpos", c(CTnegMRpos, CTnegMRposP), sep = " ")
df

#rename columns
colnames(df) <- c("Abnormality",
                  "CT positive",
                  "CT abnormality also seen on MR1",
                  "CT abnormality not seen on MR1",
                  "CT negative",
                  "MR1 also negative",
                  "MR1 shows abnormality not seen on CT",
                  "Degrees of freedom",
                  "McNemar's chi-squared",
                  "Unadjusted p-value")

##Renaming so I can reuse this table for a different way of adjusting for multiple comparisons later
df1 <- df

#Make a pretty version of the table
df1 %>% 
  kable(escape = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1, "CT positive" = 3, "CT negative" = 3, "CT vs MR1" = 3)) %>%
  pack_rows("Mass effect", 2, 4) %>%
  pack_rows("Intra-axial", 5, 6) %>%
  pack_rows("Extra-axial", 7, 9) %>%
  pack_rows("Other", 10, 11)


### COMPARING MR1 with MR2

#Please see comments in the section above to understand the meaning of the code

lesions <- colnames(reports[,7:17])
df <- data.frame(MR1posMR2neg = numeric(),
                    MR1negMR2pos = numeric(),
                    MR1posMR2pos = numeric(),
                    MR1negMR2neg = numeric(),
                    p.value    = numeric(),
                    lesion     = character())

for(i in seq_along(lesions))
{
  data <- reports %>% rename("lesion" = lesions[i])
  tab <- data %>% 
    select(lesion, timepoint, entity_id) %>% 
    spread(timepoint, lesion) %>% 
    select(MR1, MR2) %>%
    table()

  mynum <- c(
    tab[2,1],
    tab[1,2],
    tab[2,2],
    tab[1,1],
    mcnemar.test(tab)$parameter,
    mcnemar.test(tab)$statistic %>% round(3),
    mcnemar.test(tab)$p.value %>% round(3)
    ) %>%
      round(4) %>%
      matrix(nrow=1, ncol=7) %>%
      as.data.frame() %>%
      cbind(lesions[i])

  df <- rbind(df, mynum)
}

df <- df %>% select(Abnormality = "lesions[i]",
                    MR1posMR2neg = "V1",
                    MR1negMR2pos = "V2",
                    MR1posMR2pos = "V3",
                    MR1negMR2neg = "V4",
                    df           = "V5",
                    chi          = "V6",
                    p.value      = "V7")
df$MR1pos <- rowSums(df[,c("MR1posMR2neg", "MR1posMR2pos")])
df$MR1neg <- rowSums(df[,c("MR1negMR2neg", "MR1negMR2pos")])

df <- df %>% 
  select(Abnormality, MR1pos, MR1posMR2pos, MR1posMR2neg, MR1neg, MR1negMR2neg, MR1negMR2pos, df, chi, p.value)

df$MR1posP <- "(100%)"
df$MR1posMR2posP <- (df$MR1posMR2pos/df$MR1pos*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$MR1posMR2negP <- (df$MR1posMR2neg/df$MR1pos*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$MR1negP <- "(100%)"
df$MR1negMR2negP <- (df$MR1negMR2neg/df$MR1neg*100) %>% round(0) %>% paste0("(", . , "%", ")")
df$MR1negMR2posP <- (df$MR1negMR2pos/df$MR1neg*100) %>% round(0) %>% paste0("(", . , "%", ")")

df <- df %>% unite("MR1pos", c(MR1pos, MR1posP), sep = " ")
df <- df %>% unite("MR1neg", c(MR1neg, MR1negP), sep = " ")
df <- df %>% unite("MR1posMR2pos", c(MR1posMR2pos, MR1posMR2posP), sep = " ")
df <- df %>% unite("MR1posMR2neg", c(MR1posMR2neg, MR1posMR2negP), sep = " ")
df <- df %>% unite("MR1negMR2neg", c(MR1negMR2neg, MR1negMR2negP), sep = " ")
df <- df %>% unite("MR1negMR2pos", c(MR1negMR2pos, MR1negMR2posP), sep = " ")

colnames(df) <- c("Abnormality",
                  "MR1 positive",
                  "MR1 abnormality persists on MR2",
                  "MR1 abnormality resolved on MR2",
                  "MR1 negative",
                  "MR2 remains negative",
                  "MR2 shows new abnormality",
                  "Degrees of freedom",
                  "McNemar's chi-squared",
                  "Unadjusted p-value")

##Assess for statistical significance using the False Discovery Rate
#Calculate which p-value in df would be the maximum that is 
#still statistically significant when adopting a 5% false discovery rate threshold
maxp <- fdr_maxp(df$`Unadjusted p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
df$FDR <- ifelse(df$`Unadjusted p-value` <= maxp, "significant", 
                       "not sig.")


df2 <- df

df2 %>% 
  kable(escape = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1, "MR1 positive" = 3, "MR1 negative" = 3, "MR1 vs MR2" = 4)) %>%
  pack_rows("Mass effect", 2, 4) %>%
  pack_rows("Intra-axial", 5, 6) %>%
  pack_rows("Extra-axial", 7, 9) %>%
  pack_rows("Other", 10, 11)
