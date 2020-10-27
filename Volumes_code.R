#load libraries
library(dplyr)
library(tidyr)
library(kableExtra)
library(knitr)
library(ggplot2)
library(lme4)
source("fdr_function.R")

### 1. DATA PREPARATION

#load data
malpem <- read.csv("MALPEM_PtsVsCtr_20190711.csv")
malpem <- malpem %>% filter(study != "CAMCAN")

#This data is in the long format
#The rows are the participants (both patients and controls)
#For each participant there are up to three rows
#for MR1, MR2 and MR3 respectively

#The columns are the regions of brain regions as parcellated by MALPEM
#Each column contains the volume of that brain region
#In addition there are some identifier columns (participant ID, site, scan ID etc)


#collapse MALPEM regions into 15 ROIs

malpem$ConvexityCSF               <- malpem$CSF

malpem$Ventricles                 <- malpem$Ventricle

malpem$Cerebralwhitematter        <- malpem$CerebralWhiteMatter

malpem$Cerebellarwhitematter      <- malpem$CerebellarWhite

malpem$Brainstem                  <- malpem$BrainStem

malpem$Cerebellargreymatter       <- rowSums(malpem[,c("CerebellarVermis", 
                                                       "CerebellumGrey")])

malpem$Frontallobe                <- rowSums(malpem[,c("BasalForebrain", 
                                                       "DorsolateralFrontal", 
                                                       "FrontalPole", 
                                                       "MedialFrontal", 
                                                       "Orbitofrontal", 
                                                       "Precentral", 
                                                       "Subcallosal")])

malpem$Basalganglia               <- rowSums(malpem[,c("Acumbens", 
                                                       "Caudate", 
                                                       "Pallidum", 
                                                       "Putamen")])

malpem$Temporallobe               <- rowSums(malpem[,c("Amygdala", 
                                                       "FusiformGyrus", 
                                                       "TemporalLobe", "TemporalPole")])

malpem$Hippocampalcomplex         <- malpem$HippocampalComplex

malpem$Other                      <- malpem$Other

malpem$Insula                     <- malpem$Insula

malpem$Parietallobe               <- rowSums(malpem[,c("MedialParietalLobe",
                                                       "LateralParietalLobe")])

malpem$Occipitallobe              <- malpem$OccipitalLobe

malpem$Thalamus                   <- malpem$Thalamus


#select relevant column
malpem <- malpem %>% select(#identifier columns
                            study, site, subject, age, sex, timepoint, scan, model, 
                            #15 collapsed ROIs
                            ConvexityCSF,
                            Ventricles,
                            Cerebellarwhitematter,
                            Cerebralwhitematter,
                            Brainstem,
                            Cerebellargreymatter,
                            Frontallobe,
                            Temporallobe,
                            Parietallobe,
                            Occipitallobe,
                            Basalganglia,
                            Hippocampalcomplex,
                            Insula,
                            Thalamus,
                            Other
                            )


###2. COMPOSITIONAL DATA ANALYSIS

#I will follow the approach described by Aitchinson in his Concise Guide to data analysis

##2.1 Calculating perturbations

#Selecting relevant columns and rows from malpem dataframe
  #in particular only look at patients (not controls)
  #only look at MR1 and MR2 (not MR3)
pt <- malpem %>% 
  filter(timepoint != "control") %>%
  filter(timepoint != "MR3") %>%
  select(-site, -age, -sex, -scan, -model, -study) %>%
  select(subject, timepoint, everything())
pt <- pt %>% mutate(Total_per_Subject = rowSums(.[3:17]))
pt <- pt %>% select(subject, timepoint, Total_per_Subject, everything())

#Dividing each ROI volume by the patients total intracranial volume to give "bytot" i.e. the by_total_volume
bytot = lapply(pt[,-(1:3)], 
               function(x) {
                            x / pt$Total_per_Subject
                            }
               )
bytot = as.data.frame(bytot)
bytot$Subject = pt$subject
bytot$Timepoint =pt$timepoint
bytot <- bytot %>% select(Subject, Timepoint, everything())

#Calculate the ratio between the bytotal volumes of MR1 and MR2
pert <- bytot %>% gather(key = ROI, value = Volume_bytot, 3:17)
pert <- pert %>% spread(key = Timepoint, value = Volume_bytot)
pert <- pert %>% mutate(Ratio_MR1_MR2 = MR1/MR2)

##2.2 ALR transformation

#spreading the table again with ROIs as columns
pert <- pert %>% select(-MR1, -MR2)
pert <- pert %>% spread(key = ROI, value = Ratio_MR1_MR2)

#taking column Brainstem to the front as it will be the reference region
alr <- pert %>% select(Subject, Brainstem, everything())
alr <- lapply(alr[,-(1:2)], 
               function(x) {
                            log(x / alr$Brainstem)
                            }
               )
alr <-  as.data.frame(alr)
alr$Subject = pert$Subject
alr <- alr %>% select(Subject, everything())
dim(alr)
#Now I have a dataframe with subject (i.e. participant ID) as the first column
#The subsequent columns contain the additive log ratio for each Brain region
#The alr is a measure of the change in volume from MR1 to MR2
#an alr of zero would suggest that volume has not changed for that brain region

##2.3 Visualising data - does it look different from zero?

#Plot all data except the ID column
plotdata <- alr[,-1]

#Only plot a few graphs at a time, so they are displayed in a decent size
par(mfrow=c(1,5))
for (i in 1:5) {
        boxplot(plotdata[,i], main=names(plotdata[i]), type="l")
  abline(h=0, col="red")

}

par(mfrow=c(1,5))
for (i in 6:10) {
        boxplot(plotdata[,i], main=names(plotdata[i]), type="l")
  abline(h=0, col="red")
}

par(mfrow=c(1,5))
for (i in 11:14) {
        boxplot(plotdata[,i], main=names(plotdata[i]), type="l")
  abline(h=0, col="red")

}

##2.4 Formally testing if the alr's different from zero?

#I am using the Hotelling's t-test for multivariate data
#I used the same test from two different packages and they generate the same result

ICSNP::HotellingsT2(X=plotdata, mu = NULL)

rrcov::T2.test(plotdata)

### 3. UNIVARIATE ANALYSIS - WITHIN PATIENT MR1 vs MR2

##3.1 Calculate changes between scans for each patient

#Making a dataframe that shows the differences between MR1 and MR2 as
  #absolute difference
  #relative difference (ratio MR2/MR1)
  #ln of relative difference
diff <- pt %>% select(-Total_per_Subject) %>% gather(key = ROI, value = Volume, 3:17)
diff <- diff %>% spread(key = timepoint, value = Volume)
diff <- diff %>% mutate(Absdiff = MR2-MR1)
diff <- diff %>% mutate(Reldiff = MR2/MR1)
diff <- diff %>% mutate(Logrel = log(MR2/MR1))

#Make a table where
  #the rows are the subjects
  #there is one column for each ROI
  #the columns contain the log(MR2/MR1) of the volume for that ROI
  #if the log(MR2/MR1) is zero, the volume has not changed between scans
  #This is the data that I will feed into the one sample t-test (one column at a time)
logrel.df <- diff %>% select(subject, ROI, Logrel)
logrel.df <- logrel.df %>% spread(key = ROI, value = Logrel) %>% select(-subject)



##3.1 Summarise the within-patient changes for each ROI

#initialise an empty dataframe to collect results in
malp_sum <- data.frame(matrix(ncol = 6, nrow = 0))
colnames(malp_sum) <- c("ROI", "MR1_median", "MR1_iqr", "MR2_median", "MR2_iqr", "P-value")

malp_sum <- data.frame(ROI = character(), 
                       MR1_median = numeric(),
                       MR1_Q1 = numeric(),
                       MR1_Q3 = numeric(),
                       MR2_median = numeric(),
                       MR2_Q1 = numeric(),
                       MR2_Q3 = numeric(),
                      `P-value` = numeric(),
                       t = numeric(),
                       df = numeric())

#order logrel.df columns in same order as pt (dataframe with patient data see line 113)
logrel.df <- logrel.df %>% select(colnames(pt)[4:18])


#The following loop will loop through all ROIs 
  #and calculate summary statistics for MR1 and for MR2
  #and also a one-sample t-test for the log(MR1/MR2)
  #results are collected in a vector
  #vectors are being collected as rows in the empty dataframe malp_sum

for(i in 4:18)
  {
  name <- colnames(pt)[[i]]
  x <- pt[[i]]

  ttest <- t.test(logrel.df[,i-3], mu=0, alternative = "two.sided")

  mr1_median  <- subset(pt,  timepoint == "MR1", name) %>% unlist %>% median %>% round(2)
  mr1_Q1      <- subset(pt,  timepoint == "MR1", name) %>% unlist %>% quantile(.,0.25) %>% round(2)
  mr1_Q3      <- subset(pt,  timepoint == "MR1", name) %>% unlist %>% quantile(.,0.75) %>% round(2)
  mr2_median  <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% median %>% round(2)
  mr2_Q1      <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% quantile(.,0.25) %>% round(2)
  mr2_Q3      <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% quantile(.,0.75) %>% round(2) 

  vec <- c(name, 
           mr1_median, mr1_Q1, mr1_Q3, 
           mr2_median, mr2_Q1, mr2_Q3, 
           round(ttest$p.value, 3), 
           round(ttest$statistic, 3),
           ttest$parameter)

  malp_sum <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

#The first column contains the ROI name
#all other columns should be numeric as they contain medians, quartiles or p-values
malp_sum[,2:8] <- lapply(malp_sum[,2:8], function(x) as.numeric(as.character(x)))

#rename columns
colnames(malp_sum) <- c("ROI",
                        "MR1_median", 
                        "MR1_Q1",
                        "MR1_Q3", 
                        "MR2_median", 
                        "MR2_Q1",
                        "MR2_Q3", 
                        "Unadj. p-value",
                        "t",
                        "df") 

#Merge columns so that median (Q1-Q3) appears in the same cell
malp_sum$MR1_iqr  <- paste0("(",malp_sum$MR1_Q1, "-",malp_sum$MR1_Q3, ")")
malp_sum$MR2_iqr  <- paste0("(",malp_sum$MR2_Q1, "-",malp_sum$MR2_Q3, ")")
malp_sum          <- malp_sum %>% 
                      unite(col = "MR1", MR1_median, MR1_iqr, sep = " ")
malp_sum          <- malp_sum %>% 
                      unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")

#Calculate the summary statistics for the absolute and relative changes between scans
temp <- diff %>% group_by(ROI) %>% summarise(Absdiff_median = round(median(Absdiff),2), 
                                             Absdiff_Q1     = round(quantile(Absdiff, 0.25),2),
                                             Absdiff_Q3     = round(quantile(Absdiff, 0.75),2),
                                             Reldiff_median = round(median(Reldiff),2), 
                                             Reldiff_Q1     = round(quantile(Reldiff, 0.25),2),
                                             Reldiff_Q3     = round(quantile(Reldiff, 0.75),2))

#Add this extra information to the summary table malp_sum
malp_sum <- merge(malp_sum, temp, by = "ROI", all = TRUE)

#Format new columns including merging columns so that median (Q1-Q3) appears in the same cell
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+", malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)
malp_sum$Absdiff_iqr    <- paste0("(",malp_sum$Absdiff_Q1, "-",malp_sum$Absdiff_Q3, ")")
malp_sum$Reldiff_iqr    <- paste0("(",malp_sum$Reldiff_Q1, "-",malp_sum$Reldiff_Q3, ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            select(ROI, MR1, MR2, Absdiff, Reldiff, 'Unadj. p-value', t, df)

malp_sum 

## 3.3 Assess for statistical significance using the False Discovery Rate

#Calculate which p-value in malp_sum would be the maximum that is 
#still statistically significant when adopting a 5% false discovery rate threshold
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, "significant", 
                       "not sig.")

#display p < 0.001 as such
malp_sum$`Unadj. p-value` <- ifelse(malp_sum$`Unadj. p-value` < 0.001, "<0.001", 
                                    malp_sum$`Unadj. p-value`)
malp_sum$ROI <- factor(malp_sum$ROI,
                   ordered = TRUE,
                   levels = c("ConvexityCSF",
                              "Ventricles",
                              "Cerebellarwhitematter",
                              "Cerebralwhitematter",
                              "Brainstem",
                              "Cerebellargreymatter",
                              "Frontallobe",
                              "Temporallobe",
                              "Parietallobe",
                              "Occipitallobe",
                              "Basalganglia",
                              "Hippocampalcomplex",
                              "Insula",
                              "Thalamus",
                              "Other"),
                   labels = c("Convexity CSF", 
                              "Ventricles",
                              "Cerebellar white matter",
                              "Cerebral white matter",
                              "Brainstem",
                              "Cerebellar grey matter",
                              "Frontal lobe",
                              "Temporal lobe",
                              "Parietal lobe",
                              "Occipital lobe",
                              "Basal ganglia",
                              "Hippocampal complex",
                              "Insula",
                              "Thalamus",
                              "Other"))
malp_sum <- malp_sum[order(malp_sum$ROI),]

colnames(malp_sum)

##3.4 Make a pretty table for Volumes of MR1 vs MR2

malp_sum <-  malp_sum %>% select(ROI, 
                                 MR1, 
                                 MR2, 
                                 "Absolute difference" = Absdiff, 
                                 "Ratio MR2/MR1" = Reldiff,
                                 "Degrees of freedom" = df,
                                 "t-test statistic" = t,
                                 "Unadjusted p-value" = `Unadj. p-value`, 
                                FDR)

malp_sum %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  "ROI volume in cm^3, Median (IQR)" = 4, "MR1 vs MR2" = 4)) %>%
  pack_rows("Cerebrospinal fluid", 1, 2) %>%
  pack_rows("White matter", 3, 4) %>%
  pack_rows("Infratentorial grey matter", 5, 6) %>%
  pack_rows("Supratentorial grey matter - lobes", 7, 10) %>%
  pack_rows("Supratentorial grey matter - special regions", 11, 15) 

###4. CHECK IF CHANGE IN VOLUME IS REFLECTIVE OF THE OVERALL TREND

# i.e. not just driven by an extreme subpopulation of patients
#I will pick the three regions (CSF, Ventricles, CWM) that changed significantly
#and plot the change of each patient compared the median change of all patients

#Plotting change in volume between MR1 and MR2 irrespective of time between scans
delta_data <- diff %>% 
  filter(ROI %in% c("Ventricles", "ConvexityCSF", "Cerebralwhitematter")) %>%
  select(subject, ROI, Logrel)


ggplot(data = delta_data, aes(x = ROI, y = Logrel)) +
    geom_jitter(alpha = 0.4, width = 0.3, height = 0) +
    geom_boxplot(alpha = 0) +
    ylim(-1, +1) +
    geom_abline(intercept = 0, slope = 0, colour = "tomato") +
    theme_bw()  

###5. COMPARE PATIENTS WITH CONTROLS

#I will pick the three regions (CSF, Ventricles, CWM) that changed significantly within-patients from MR1 to MR2. 
#I will compare the volume of these lesions between patients and controls at MR1 and then at MR2

##5.1 Data preparation

#Normalising for total intracranial volume in each patient
temp          <- malpem
temp$cortgrey <- rowSums(temp[,15:22])
temp$tiv      <- rowSums(temp[,9:23])
bytot         <-  lapply(temp[,-c(1:8,25)], function(x) {x / temp$tiv})
bytot         <- as.data.frame(bytot)
bytot         <- cbind(temp[,c(1:8, 25)], bytot) %>% as.data.frame()
bytot         <- bytot %>% 
                  select(study, 
                         site, 
                         subject, 
                         age, 
                         sex, 
                         timepoint, 
                         scan, 
                         model, 
                         tiv, 
                         everything())

#Combining column for site and model into one variable ("scanner") which identifies individual machines
bytot <- bytot %>% unite("scanner", site, model, sep = "_")
bytot$scanner <- factor(bytot$scanner,
                           levels = unique(bytot$scanner),
                           labels = unique(bytot$scanner))

#removing camcan data (a dataset that I am not using for the current paper)
bytot       <- bytot %>% 
                filter(study != "CAMCAN") 

#adding column for group (i.e. whether a participant is a patient or a control)
bytot$group <- ifelse(bytot$timepoint %in% c("MR1", "MR2", "MR3"), "patient", "control")
bytot$group <- factor(bytot$group, 
                      levels = c("control", "patient"), 
                      labels = c("control", "patient"))

#The sex for one control subject was missing but could be retrieved and is male
bytot$sex <-  ifelse(bytot$subject == "5mDK752", "M", bytot$sex)

#save this table of normalised ROI volumes for future use
#write.csv(bytot, "normalised_volumes.csv", row.names = FALSE)

#selecting relevant ROIs and make two separate dataframes, one for MR1 and one for MR2
data_MR1    <- bytot %>%
                filter(timepoint %in% c("MR1", "control")) %>%
                select(cortgrey, 
                       Ventricles, 
                       ConvexityCSF, 
                       Cerebralwhitematter, 
                       group, 
                       age, 
                       sex, 
                       scanner, 
                       timepoint)

data_MR2    <- bytot %>%
                filter(timepoint %in% c("MR2", "control")) %>%
                select(cortgrey, 
                       Ventricles, 
                       ConvexityCSF, 
                       Cerebralwhitematter, 
                       group, 
                       age, 
                       sex, 
                       scanner, 
                       timepoint)


##5.2 Checking assumptions

#MR1
  #Ventricles
    fit <- lmer(log(Ventricles) ~ group + age + sex + (1| scanner), 
            data = data_MR1, 
            REML = FALSE, 
            control = lmerControl(optimizer ="Nelder_Mead"))

    #homogeneity of variance
    plot(fit)
    #linearity
    qqnorm(resid(fit))
    qqline(resid(fit))
#ConvexityCSF
    fit <- lmer(log(ConvexityCSF) ~ group + age + sex + (1| scanner), 
            data = data_MR1, 
            REML = FALSE, 
            control = lmerControl(optimizer ="Nelder_Mead"))

    #homogeneity of variance
    plot(fit)
    #linearity
    qqnorm(resid(fit))
    qqline(resid(fit))
#Cerebralwhitematter
    fit <- lmer(log(Cerebralwhitematter) ~ group + age + sex + (1| scanner), 
            data = data_MR1, 
            REML = FALSE, 
            control = lmerControl(optimizer ="Nelder_Mead"))

    #homogeneity of variance
    plot(fit)
    #linearity
    qqnorm(resid(fit))
    qqline(resid(fit))

##5.3 Modeling the volume of CSF/Ventricles/WM using mixed models

#Making a dataframe to collect results of pt vs ctr comparison
#Volume at MR1
df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("ROI", "Data", "Coef_pt", "SE", "DF", "t", "P_value")

df <- data.frame(ROI        = character(), 
                Metric      = character(),
                Coef_pt     = numeric(),
                SE          = numeric(),
                DF          = numeric(),
                t           = numeric(),
                P_value     = numeric())

for (i in 2:4) {

Data <- data_MR1
Metric <- "data_MR1"
ROI <- colnames(Data)[[i]]

library(lmerTest)

fit <- lmer(log(Data[[i]]) ~ group + age + sex + (1| scanner), data = Data, REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
Coef_pt <- summary(fit)$coefficients[2,1] %>% round(3)
SE      <- summary(fit)$coefficients[2,2] %>% round(3)
DF      <- summary(fit)$coefficients[2,3] %>% round(0)
t       <- summary(fit)$coefficients[2,4] %>% round(3)
P_value <- summary(fit)$coefficients[2,5] %>% round(3)

 
vec <- c(ROI, Metric, Coef_pt, SE, DF, t, P_value)

df <- berryFunctions::addRows(df, 1, values = vec)

}

df[,3:7]    <- lapply(df[,3:7], function(x) as.numeric(as.character(x)))
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")
df          <- df %>% separate(Metric, c("Metric", "Timepoint")) %>% select(-Metric)
dfMR1       <- df
dfMR1



#Making a dataframe to collect results of pt vs ctr comparison
#Volume at MR2
df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("ROI", "Data", "Coef_pt", "SE", "DF", "t", "P_value")

df <- data.frame(ROI        = character(), 
                Metric      = character(),
                Coef_pt     = numeric(),
                SE          = numeric(),
                DF          = numeric(),
                t           = numeric(),
                P_value     = numeric())

#Run the models for each of the 3 ROIs
for (i in 2:4) {

Data <- data_MR2
Metric <- "data_MR2"
ROI <- colnames(Data)[[i]]

library(lmerTest)

fit <- lmer(log(Data[[i]]) ~ group + age + sex + (1| scanner), data = Data, REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
Coef_pt <- summary(fit)$coefficients[2,1] %>% round(3)
SE      <- summary(fit)$coefficients[2,2] %>% round(3)
DF      <- summary(fit)$coefficients[2,3] %>% round(0)
t       <- summary(fit)$coefficients[2,4] %>% round(3)
P_value <- summary(fit)$coefficients[2,5] %>% round(3)

 
vec <- c(ROI, Metric, Coef_pt, SE, DF, t, P_value)

df <- berryFunctions::addRows(df, 1, values = vec)

}

df[,3:7]    <- lapply(df[,3:7],function(x) as.numeric(as.character(x)))

df$Coef_pt  <- ifelse(df$Coef_pt > 0, 
                      paste0("+", df$Coef_pt), 
                      paste0(df$Coef_pt))

maxp        <- fdr_maxp(df$P_value, 0.05)

df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

df          <- df %>% 
                separate(Metric, c("Metric", "Timepoint")) %>% 
                select(-Metric)

dfMR2       <- df




#Make combined table of pt vs controls for MR1 and MR2

dfboth          <- rbind(dfMR1, dfMR2) %>% as.data.frame()
dfboth$Coef_pt  <- as.numeric(dfboth$Coef_pt)
dfboth$Relvol   <- exp(dfboth$Coef_pt) %>% round(2)
dfboth$Coef_pt  <- round(dfboth$Coef_pt, 3)
dfboth$P_value  <- round(dfboth$P_value, 3)
dfboth          <- dfboth %>% select(
                            Timepoint,
                            ROI,
                            "Raw Coefficient" = Coef_pt,
                            "Standard Error" = SE,
                            "Relative volume compared to control" = Relvol,
                            "Degrees of freedom" = DF,
                            "t-test statistic" = t,
                            "Unadjusted p-value" = P_value,
                            "FDR-adjusted significance" = FDR)


#Make a pretty version of the above table
dfboth %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

###6. DO VOLUMES NORMALISE AT MR3?

##6.1 Patients vs controls at MR3

#select relevant data
data_MR3 <- bytot %>%
  filter(timepoint %in% c("MR3", "control")) %>%
  select(cortgrey, 
         Ventricles, 
         ConvexityCSF, 
         Cerebralwhitematter, 
         group, 
         age, 
         sex, 
         scanner, 
         timepoint)

#Making a dataframe to collect results of pt vs ctr comparison
#Volume at MR3
df <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(df) <- c("ROI", "Data", "Coef_pt", "SE", "DF", "t", "P_value")

df <- data.frame(ROI        = character(), 
                Metric      = character(),
                Coef_pt     = numeric(),
                SE          = numeric(),
                DF          = numeric(),
                t           = numeric(),
                P_value     = numeric())

for (i in 2:4) {

Data <- data_MR3
Metric <- "data_MR3"
ROI <- colnames(Data)[[i]]

library(lmerTest)

fit <- lmer(log(Data[[i]]) ~ group + age + sex + (1| scanner), data = Data, REML = FALSE, control = lmerControl(optimizer ="Nelder_Mead"))
  
Coef_pt <- summary(fit)$coefficients[2,1] %>% round(3)
SE      <- summary(fit)$coefficients[2,2] %>% round(3)
DF      <- summary(fit)$coefficients[2,3] %>% round(0)
t       <- summary(fit)$coefficients[2,4] %>% round(3)
P_value <- summary(fit)$coefficients[2,5] %>% round(3)

 
vec <- c(ROI, Metric, Coef_pt, SE, DF, t, P_value)

df <- berryFunctions::addRows(df, 1, values = vec)

}

df[,3:7]      <- lapply(df[,3:7], function(x) as.numeric(as.character(x)))
df$Coef_pt    <- ifelse(df$Coef_pt > 0, paste0("+", df$Coef_pt), paste0(df$Coef_pt))
maxp          <- fdr_maxp(df$P_value, 0.05)
df$FDR        <- ifelse(df$P_value <= maxp, "significant", "not sig.")
df            <- df %>% separate(Metric, c("Metric", "Timepoint")) %>% select(-Metric)
dfMR3         <- df

dfMR3$Coef_pt <- as.numeric(dfMR3$Coef_pt)
dfMR3$Relvol  <- exp(dfMR3$Coef_pt) %>% round(2)
dfMR3$Coef_pt <- round(dfMR3$Coef_pt, 5)
dfMR3$P_value <- round(dfMR3$P_value, 5)
dfMR3         <- dfMR3 %>% select(
                            Timepoint,
                            ROI,
                            "Raw Coefficient" = Coef_pt,
                            "Standard Error" = SE,
                            "Relative volume compared to control" = Relvol,
                            "Degrees of freedom" = DF,
                            "t-test statistic" = t,
                            "Unadjusted p-value" = P_value,
                            "FDR-adjusted significance" = FDR)

#Make a pretty version of the above table
dfMR3 %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

##6.2 univariate analysis, within patient MR2 vs MR3

#this code is analogous to the one used for MR1 vs MR2
#Please refer to the more extensive comments in that previous code

#select relevant data
pt <- malpem %>% 
      filter(timepoint != "control") %>%
      filter(timepoint != "MR1") %>%
      select(-site, -age, -sex, -scan, -model, -study) %>%
      select(subject, timepoint, everything())
pt <- pt %>% mutate(Total_per_Subject = rowSums(.[3:17]))
pt <- pt %>% select(subject, timepoint, Total_per_Subject, everything())

#Making a dataframe that shows the differences between MR2 and MR3 as
  #absolute difference
  #relative difference
  #ln of relative difference
  diff <- pt %>% select(-Total_per_Subject) %>% gather(key = ROI, value = Volume, 3:17)
  diff <- diff %>% spread(key = timepoint, value = Volume)
  diff <- diff %>% mutate(Absdiff = MR3-MR2)
  diff <- diff %>% mutate(Reldiff = MR3/MR2)
  diff <- diff %>% mutate(Logrel = log(MR3/MR2))
  # remove patients without MR3 
  diff <- na.omit(diff) 

#spreading the table again with ROIs as columns
logrel.df <- diff %>% select(subject, ROI, Logrel)
logrel.df <- logrel.df %>% spread(key = ROI, value = Logrel) %>% select(-subject)
logrel.df <- na.omit(logrel.df)# remove patients without MR3
malp_sum  <- data.frame(matrix(ncol = 10, nrow = 0))
colnames(malp_sum) <- c("ROI", 
                        "MR2_median", 
                        "MR2_Q1", 
                        "MR2_Q3",
                        "MR3_median", 
                        "MR3_Q1", 
                        "MR3_Q3",
                        "P-value",
                        "t",
                        "df")

malp_sum  <- data.frame(ROI        = character(), 
                        MR2_median = numeric(),
                        MR2_Q1     = numeric(),
                        MR2_Q3     = numeric(),
                        MR3_median = numeric(),
                        MR3_Q1     = numeric(),
                        MR3_Q3     = numeric(),
                        `P-value` = numeric(),
                        t = numeric(),
                        df = numeric())
#order logrel.df columns in same order as lr, since I will use t-test on rel_raw
logrel.df <- logrel.df %>% select(colnames(pt)[4:18])

#run baseline statistics on each ROI
for(i in 4:18)
  {
  name <- colnames(pt)[[i]]
  x <- pt[[i]]

  ttest <- t.test(logrel.df[,i-3], mu=0, alternative = "two.sided")

  mr2_median  <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% median %>% round(2)
  mr2_Q1      <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% quantile(.,0.25) %>% round(2)
  mr2_Q3      <- subset(pt,  timepoint == "MR2", name) %>% unlist %>% quantile(.,0.75) %>% round(2)
  mr3_median  <- subset(pt,  timepoint == "MR3", name) %>% unlist %>% median %>% round(2)
  mr3_Q1      <- subset(pt,  timepoint == "MR3", name) %>% unlist %>% quantile(.,0.25) %>% round(2)
  mr3_Q3      <- subset(pt,  timepoint == "MR3", name) %>% unlist %>% quantile(.,0.75) %>% round(2)

  vec <- c(name, 
           mr1_median, mr1_Q1, mr1_Q3, 
           mr2_median, mr2_Q1, mr2_Q3, 
           round(ttest$p.value, 3), 
           round(ttest$statistic, 3),
           ttest$parameter)

  malp_sum <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

malp_sum[,2:10] <- lapply(malp_sum[,2:10], function(x) as.numeric(as.character(x)))
colnames(malp_sum) <- c("ROI",
                        "MR2_median", 
                        "MR2_Q1",
                        "MR2_Q3", 
                        "MR3_median", 
                        "MR3_Q1",
                        "MR3_Q3", 
                        "Unadj. p-value",
                        "t",
                        "df") 

malp_sum$MR2_iqr  <- paste0("(",malp_sum$MR2_Q1, "-",malp_sum$MR2_Q3, ")")
malp_sum$MR3_iqr  <- paste0("(",malp_sum$MR3_Q1, "-",malp_sum$MR3_Q3, ")")
malp_sum          <- malp_sum %>% unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")
malp_sum          <- malp_sum %>% unite(col = "MR3", MR3_median, MR3_iqr, sep = " ")

temp <- diff %>% group_by(ROI) %>% summarise(Absdiff_median = round(median(Absdiff),2), 
                                             Absdiff_Q1     = round(quantile(Absdiff, 0.25),2),
                                             Absdiff_Q3     = round(quantile(Absdiff, 0.75),2),
                                             Reldiff_median = round(median(Reldiff),2), 
                                             Reldiff_Q1     = round(quantile(Reldiff, 0.25),2),
                                             Reldiff_Q3     = round(quantile(Reldiff, 0.75),2))

malp_sum                <- merge(malp_sum, temp, by = "ROI", all = TRUE)

#display medians and quartiles in the same cell of the results table
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+",malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)
malp_sum$Absdiff_iqr    <- paste0("(",malp_sum$Absdiff_Q1, "-",malp_sum$Absdiff_Q3, ")")
malp_sum$Reldiff_iqr    <- paste0("(",malp_sum$Reldiff_Q1, "-",malp_sum$Reldiff_Q3, ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")

#Select relevant columns and rows
malp_sum  <- malp_sum %>% select(ROI, MR2, MR3, Absdiff, Reldiff, 'Unadj. p-value', t, df)
malp_sum  <- malp_sum %>% 
              filter(ROI %in% c("Ventricles", "ConvexityCSF", "Cerebralwhitematter")) 
malp_sum  <- malp_sum %>%
              mutate(ROI =  factor(ROI, 
                                   levels=c("Ventricles", 
                                            "ConvexityCSF", 
                                            "Cerebralwhitematter"),
                                   labels = c("Ventricles", 
                                              "ConvexityCSF",
                                              "Cerebralwhitematter")))%>%
              arrange(ROI) 

##Assess for statistical significance using the False Discovery Rate
#Calculate which p-value in malp_sum would be the maximum that is 
#still statistically significant when adopting a 5% false discovery rate threshold
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, "significant", 
                       "not sig.")

#display p < 0.001 as such
malp_sum$`Unadj. p-value` <- ifelse(malp_sum$`Unadj. p-value` < 0.001, "<0.001", 
                                    malp_sum$`Unadj. p-value`)

#make a pretty version of the table

malp_sum <-  malp_sum %>% select(ROI, 
                                 MR2, 
                                 MR3, 
                                 "Absolute difference" = Absdiff, 
                                 "Ratio MR2/MR1" = Reldiff,
                                 "Degrees of freedom" = df,
                                 "t-test statistic" = t,
                                 "Unadjusted p-value" = `Unadj. p-value`, 
                                FDR)
malp_sum %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  "ROI volume in cm^3, Median (IQR)" = 4, "MR2 vs MR3" = 4)) 
