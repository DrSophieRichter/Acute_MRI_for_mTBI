#load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(kableExtra)
library(lme4)
source("fdr_function.R")

### 1. DATA PREPARATION

#The following files are in the long format
#the rows are the participants
  #there are up to three rows per participants (for MR1, MR2, MR3)
  #(actually for CENTER-data there are up to 6 rows per patient, 
  # DTI sequences have been obtained with and without b0 at each timepoint)
#the columns are the 72 white matter tracts
  #depending on the dataframe, the column contain either the FA, MD or volume
  #of that particular tract
  #in addition there are also some identifier columns
  #e.g. patient ID, scan ID, whether a participant is a patient or a control

c_md <- read.csv("CENTER_rawMeanMD.csv")
c_md <- c_md %>% select(Study:X71)

c_fa <- read.csv("CENTER_rawMeanFA.csv")
c_fa <- c_fa %>% select(Study:X71)

c_vo <- read.csv("CENTER_TractSeg_volumes.csv")
c_vo <- c_vo %>% select(Study:X71)

d_md <- read.csv("DETECT_rawMeanMD.csv")
d_md <- d_md %>% filter(Study == "DETECT", Category %in% c("healthy", "patients"))

d_fa <- read.csv("DETECT_rawMeanFA.csv")
d_fa <- d_fa %>% filter(Study == "DETECT", Category %in% c("healthy", "patients"))

d_vo <- read.csv("DETECT_TractSeg_volumes.csv")
d_vo <- d_vo %>% filter(Study == "DETECT", Category %in% c("healthy", "patients"))

#The following dataframe is also in the long format
#while it contains MALPEM volumes, these will not be used
#the dataframe is being loaded as it contains useful identifier columns such as 
#which scanner the scan was conducted on
pt <- read.csv("MALPEM_PtsVsCtr_20190711.csv")
pt <- pt %>% select(study:model)


#remove CamCan data from dataframe "pt", as I will not be using it for this paper 
pt        <- pt %>% filter(study != "CAMCAN")
pt$study  <- factor(pt$study,
                   levels = unique(pt$study),
                   labels = unique(pt$study))
pt$study  <- pt$study %>% droplevels()

#For diffusion, I cannot consider Cambridge CENTER and Cambridge DETECT the same site
pt$site   <- ifelse(pt$study == "DETECT", "DETECT", as.character(pt$site))
pt$site   <- factor(pt$site,
                  levels = unique(pt$site),
                  labels = unique(pt$site))
pt$site   <- pt$site %>% droplevels()

#checking all patients and controls for paper 1 are now listed in dataframe pt
str(pt)
pt %>% filter(timepoint == "control") %>% distinct(subject) %>% dim
pt %>% filter(timepoint != "control") %>% distinct(subject) %>% dim


#Making subject and scan IDs compatible for merging
c_md$SubjectID <- gsub(pattern = "sub-", replacement = "", x=c_md$SubjectID)
c_fa$SubjectID <- gsub(pattern = "sub-", replacement = "", x=c_fa$SubjectID)
c_vo$SubjectID <- gsub(pattern = "sub-", replacement = "", x=c_vo$SubjectID)

c_md$ScanID <- gsub(pattern = "ses-", replacement = "", x=c_md$ScanID)
c_fa$ScanID <- gsub(pattern = "ses-", replacement = "", x=c_fa$ScanID)
c_vo$ScanID <- gsub(pattern = "ses-", replacement = "", x=c_vo$ScanID)

#There is an inconsistency in dataframe pt. For DETECT patients the scan ID is "ID.....", for DETECT controls (aka healthy) it is "......" (a date). I need to replicate this in the d dataframes

d_md$ScanID <- gsub(pattern = "_U-", replacement = "_", x=d_md$ScanID)
d_md        <- d_md %>% separate(ScanID, into = c("date", "id"), sep = "_")
d_md$scan   <- ifelse(d_md$Category == "patients", as.character(d_md$id), as.character(d_md$date))

d_fa$ScanID <- gsub(pattern = "_U-", replacement = "_", x=d_fa$ScanID)
d_fa        <- d_fa %>% separate(ScanID, into = c("date", "id"), sep = "_")
d_fa$scan   <- ifelse(d_fa$Category == "patients", as.character(d_fa$id), as.character(d_fa$date))

d_vo$ScanID <- gsub(pattern = "_U-", replacement = "_", x=d_vo$ScanID)
d_vo        <- d_vo %>% separate(ScanID, into = c("date", "id"), sep = "_")
d_vo$scan   <- ifelse(d_vo$Category == "patients", as.character(d_vo$id), as.character(d_vo$date))


#merging dataframes for MD, to have CENTER and DETECT MD data for paper 1 in one dataframe

#data from DETECT study
md1   <- merge(pt, d_md, 
               by.x = c("subject", "scan"), 
               by.y = c("SubjectID", "scan"), 
               all.x = FALSE, 
               all.y = FALSE)
md1   <- md1 %>% select(-Study, -Category,-date, -id )

#data from CENTER study
c_md  <- c_md %>% distinct(ScanID, .keep_all = TRUE) #only picking dti_scan1 per person
md2   <- merge(pt, c_md, 
               by.x = c("subject", "scan"), 
               by.y = c("SubjectID", "ScanID"), 
               all.x = FALSE, 
               all.y = FALSE)
md2   <- md2 %>% select(-Study, -Category, -CentreID, -DTI_set)

#Combine DETECT and CENTER data
md    <- rbind.data.frame(md1, md2) %>% as.data.frame()
dim(md)
#This reveals (ans was confirmed by checking the source data) that for 10 scans diffusion data was not obtained

#Listing patients and scans for which DTI data is missing
missing <- pt %>% filter(!scan %in% md$scan)
missing
dim(missing)
length(unique(missing$subject))

#DTI is missing for 10 scans of 9 subjects (i.e. 1 subject is missing DTI for both their scans)
#8 of the 9 subjects are patients, the other one is a control

#merging dataframes for FA, to have CENTER and DETECT FA data for paper 1 in one dataframe

#data from DETECT study
fa1     <- merge(pt, d_fa, 
                 by.x = c("subject", "scan"), 
                 by.y = c("SubjectID", "scan"), 
                 all.x = FALSE, 
                 all.y = FALSE)
fa1     <- fa1 %>% select(-Study, -Category, -date, -id)

#data from CENTER study
c_fa    <- c_fa %>% distinct(ScanID, .keep_all = TRUE) #only picking dti_scan1 per person
fa2     <- merge(pt, c_fa, 
                 by.x = c("subject", "scan"), 
                 by.y = c("SubjectID", "ScanID"), 
                 all.x = FALSE, 
                 all.y = FALSE)
fa2     <- fa2 %>% select(-Study, -Category, -CentreID, -DTI_set)

#Combine DETECT and CENTER data
fa <- rbind.data.frame(fa1, fa2) %>% as.data.frame()
dim(fa) #As expected 10 scans missing DTI data


#merging dataframes for Volume, to have CENTER and DETECT TractSeg volume data for paper 1 in one dataframe

#data from DETECT study
vo1     <- merge(pt, d_vo, 
                 by.x = c("subject", "scan"), 
                 by.y = c("SubjectID", "scan"), 
                 all.x = FALSE, 
                 all.y = FALSE)
vo1     <- vo1 %>% select(-Study, -Category, -date, -id)

#data from CENTER study
c_vo    <- c_vo %>% distinct(ScanID, .keep_all = TRUE) #only picking dti_scan1 per person
vo2     <- merge(pt, c_vo, 
             by.x = c("subject", "scan"), 
             by.y = c("SubjectID", "ScanID"), 
             all.x = FALSE, 
             all.y = FALSE)
vo2     <- vo2 %>% select(-Study, -Category, -CentreID, -DTI_set)

#Combine DETECT and CENTER data
vo      <- rbind.data.frame(vo1, vo2) %>% as.data.frame()
dim(vo)
#As expected 10 scans missing TractSeg volumes

###2. QUALITY CONTROL OF AVAILABLE DTI SCANS

##2.1 Prepare Quality Control (QC) data

c_qc <- read.csv("CENTER_QC.csv")
c_qc <- c_qc %>% select(Study:DTIFIT.num.outlier)
d_qc <- read.csv("DETECT_QC.csv")

#Making subject and scan IDs compatible for merging
c_qc$SubjectID  <- gsub(pattern = "sub-", replacement = "", x=c_qc$SubjectID)
c_qc$ScanID     <- gsub(pattern = "ses-", replacement = "", x=c_qc$ScanID)

#There is an inconsistency in dataframe pt. For DETECT patients the scan ID is "ID.....", for DETECT controls (aka healthy) it is "......" (a date). I need to replicate this in the d dataframes

d_qc$ScanID   <- gsub(pattern = "_U-", replacement = "_", x=d_qc$ScanID)
d_qc          <- d_qc %>% separate(ScanID, into = c("date", "id"), sep = "_")
d_qc$scan     <- ifelse(d_qc$Category == "patients", 
                        as.character(d_qc$id), 
                        as.character(d_qc$date))

#merging dataframes for QC, to have CENTER and DETECT QC data for paper 1 in one dataframe
#data from DETECT study
qc1     <- merge(pt, d_qc, 
                 by.x = c("subject", "scan"), 
                 by.y = c("SubjectID", "scan"), 
                 all.x = FALSE, 
                 all.y = FALSE)
qc1     <- qc1 %>% select(-Study, -Category,-date, -id )

#data from CENTER study
c_qc    <- c_qc %>% distinct(ScanID, .keep_all = TRUE) #only picking dti_scan1 per person
qc2     <- merge(pt, c_qc, 
                 by.x = c("subject", "scan"), 
                 by.y = c("SubjectID", "ScanID"), 
                 all.x = FALSE, 
                 all.y = FALSE)
qc2     <- qc2 %>% select(-Study, 
                          -Category, 
                          -CentreID, 
                          -DTI_set, 
                          -Number_DTI_Volumes, 
                          -Extra_b0_used, 
                          -ScannerManufacturer, 
                          -ScannerModel)

#combine CENTER and DETECT data
colnames(qc2) <- colnames(qc1)
qc            <- rbind.data.frame(qc1, qc2) %>% as.data.frame()
dim(qc)
colnames(qc)

##2.1 Look for outliers in motion scores

qc$group <- ifelse(qc$timepoint == "control", "control", "patient")
qc$group <- factor(qc$group, 
                   levels = c("control", "patient"), 
                   labels = c("control", "patient"))

#First metric: average total motion
ggplot(data=qc, aes(x=avg_total_motion)) +
        geom_histogram() + 
        facet_grid(study ~ group) + 
        ggtitle("avg_total_motion")

#Collect extreme values of average total motion in a dataframe
out_avg_total_motion <- qc %>% filter(avg_total_motion > 1.5) %>% select(subject:model)
out_avg_total_motion

#Second metric: maximal total motion
ggplot(data=qc, aes(x=max_total_motion)) +
  geom_histogram() + facet_grid(study ~ group) + ggtitle("max_total_motion")

#Collect extreme values of maximal total motion in a dataframe
out_max_total_motion <- qc %>% filter(max_total_motion > 2.3) %>% select(subject:model)
out_max_total_motion

#Third metric: average relative motion
ggplot(data=qc, aes(x=avg_relative_motion)) +
  geom_histogram() + facet_grid(study ~ group) + ggtitle("avg_relative_motion")

#Collect extreme values of avergae relative motion in a dataframe
out_avg_relative_motion <- qc %>% filter(avg_relative_motion > 0.6) %>% select(subject:model)
out_avg_relative_motion

#Fourth metric: maximal relative motion
ggplot(data=qc, aes(x=max_relative_motion)) +
  geom_histogram() + facet_grid(study ~ group) + ggtitle("max_relative_motion")

#Collect extreme values of maximal relative motion in a dataframe
out_max_relative_motion <- qc %>% filter(max_relative_motion > 2) %>% select(subject:model)
out_max_relative_motion

#Make a list of scans in which patients moved a lot
movers <- rbind(out_avg_relative_motion,
                out_max_relative_motion,
                out_avg_total_motion,
                out_max_total_motion) %>%
  as.data.frame() %>% 
  distinct(scan, .keep_all = TRUE) #remove duplicates
movers
#write.csv(movers, file = "movers.csv")
#After review the following scan must be removed due to movement artefact
fa <- fa %>% filter(scan != "E29776")
md <- md %>% filter(scan != "E29776")
vo <- vo %>% filter(scan != "E29776")

##2.2 Look for outliers in Tractseg volumes, FA and MD by site/scanner
# this is done separately for patients and controls

#How many patients and controls are available for each site-scanner?

fa$group <- ifelse(fa$timepoint == "control", "control", "patient")
fa$group <- factor(fa$group, levels = c("control", "patient"), labels = c("control", "patient"))

summary(fa$model)
summary(fa$site)
summary(fa$group)
summary(fa$timepoint)
fa %>% group_by(site, model, group) %>% summarise(n()) %>% as.data.frame() %>% spread("group", "n()")

#Create a new variable that incoporates both site and scanner model, since some sites have more than one scanner model
fa          <- fa %>% unite("scanner", c("site", "model"))
fa$scanner  <- factor(fa$scanner,
                      levels = unique(fa$scanner),
                      labels = unique(fa$scanner))
summary(fa$scanner) %>% as.data.frame %>% dim() #number of different scanners

md          <- md %>% unite("scanner", c("site", "model"))
md$scanner  <- factor(md$scanner,
                      levels = unique(md$scanner),
                      labels = unique(md$scanner))

vo          <- vo %>% unite("scanner", c("site", "model"))
vo$scanner  <- factor(vo$scanner,
                      levels = unique(vo$scanner),
                      labels = unique(vo$scanner))

#Identifying outliers in FA

#Exploratory plot of corpus callosum as an example ROI (ROI number X45)
ggplot(fa, aes(y=X45, x = scanner)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90)) +
  facet_grid(group~.) +
  ggtitle("FA for Corpus callosum (X45)")

fa %>% filter(X45 > 0.5) %>% select(subject:timepoint)

#Excluding the scan with gross artefact
fa <- fa %>% filter(!scan == "E37744")
md <- md %>% filter(!scan == "E37744")
vo <- vo %>% filter(!scan == "E37744")

#Making a list of control outliers for fa (having excluded one extreme outlier already)
is_outlier <- function(x) 
  {
    return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
  }

fa_out_control <-
  fa %>% 
  filter(group == "control") %>% 
  group_by(scanner) %>% 
  filter(is_outlier(X45) == TRUE) %>% 
  select(subject:scanner, X45) %>% 
  as.data.frame()
fa_out_control$metric <- "FA"
fa_out_control

#Looking at md rather than fa
md$group <- ifelse(md$timepoint == "control", "control", "patient")
md$group <- factor(md$group, 
                   levels = c("control", "patient"), 
                   labels = c("control", "patient"))

ggplot(md, aes(y=X45, x = scanner)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90))+
  facet_grid(group~.)+
  ggtitle("MD for Corpus callosum (X45)")

#Looking at patient MD outliers
md %>% filter(X45 > 0.0009) %>% select(subject:timepoint)

#Making a list of control outliers for md
md_out_control <-
  md %>% 
  filter(group == "control") %>% 
  group_by(scanner) %>% 
  filter(is_outlier(X45) == TRUE) %>% 
  select(subject:scanner, X45) %>% 
  as.data.frame()
md_out_control$metric <- "MD"
md_out_control

#Looking at vo rather than fa
vo$group <- ifelse(vo$timepoint == "control", "control", "patient")
vo$group <- factor(vo$group, 
                   levels = c("control", "patient"), 
                   labels = c("control", "patient"))


ggplot(vo, aes(y=X45, x = scanner)) +
  geom_boxplot() + theme(axis.text.x = element_text(angle = 90))+
  facet_grid(group~.) +
  ggtitle("Volume for Corpus callosum (X45)")

#Making a list of control outliers for volumes

vo_out_control <-
  vo %>% 
  filter(group == "control") %>% 
  group_by(scanner) %>% 
  filter(is_outlier(X45) == TRUE) %>% 
  select(subject:scanner, X45) %>% 
  as.data.frame()
vo_out_control$metric <- "volume"
vo_out_control

#Combing the lists of control outliers for fa, md and vo
out_control <- rbind(vo_out_control,
                     fa_out_control,
                     md_out_control) %>%
               as.data.frame() 
out_control <- out_control[order(out_control$metric), ]
out_control
#write.csv(out_control, file = "Control_diffusion_outliers_withdata.csv", row.names = FALSE)

#Removing outliers with artefact/movement that is not acceptable after visual inspection
rejected <- c("E40441", "E26156", "E11063")
md <- md %>% filter(! scan %in% rejected)
fa <- fa %>% filter(! scan %in% rejected)


#The remaining data is of high enough quality for analysis

#make variable timepoint a factor
fa$timepoint <- factor(fa$timepoint)
md$timepoint <- factor(md$timepoint)

#write csv files of raw values
write.csv(md, file = "md_allsubjects_alltracts_nooutliers.csv")
write.csv(fa, file = "fa_allsubjects_alltracts_nooutliers.csv")

###3. COMPARE MD AT MR1 WITH MR2

#Make dataframe for subjects that have paired diffusion available for MR1 and MR2
mdp <- md %>% filter(timepoint %in% c("MR1", "MR2"))
mdp <- mdp %>% group_by(subject) %>% filter(n()>1)

#Making a dataframe that shows the differences between MR1 and MR2 as
#absolute difference
#relative difference
#ln of relative difference
mdp <- mdp %>% select(subject, timepoint, X0:X71)

#Multiply all md values by 1000 to avoid small numbers
mdp <- mdp %>% 
  mutate_at(vars(X0:X71),
            .funs = list(~. * 1000))
diff <- mdp %>%  
  gather(key = tract, value = md, X0:X71)
diff <- diff %>% spread(key = timepoint, value = md)
diff <- diff %>% mutate(Absdiff = MR2-MR1)
diff <- diff %>% mutate(Reldiff = MR2/MR1)
diff <- diff %>% mutate(Logrel = log(MR2/MR1))

#spreading the table again with tracts as columns
#in this table the rows are the patients
#For each patient and tract we have one "logrel" value i.e. log(MR2/MR1)
#summarising the between-scan change in a single value per patient per tract
#allows me to use a one-sample t-test later on
logrel.df <- diff %>% select(subject, tract, Logrel)
logrel.df <- logrel.df %>% spread(key = tract, value = Logrel) %>% as.data.frame() %>% select(-subject)

#Prevent R from using exponential notation
options(scipen = 999)

#initialise a dataframe in which I will collect the summary statistics for each tract
malp_sum              <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(malp_sum)    <- c( "Tract", 
                            "MR1_median", 
                            "MR1_iqr", 
                            "MR2_median", 
                            "MR2_iqr", 
                            "N_paired", 
                            "P-value")

malp_sum <- data.frame(Tract       = character(), 
                       MR1_median  = numeric(),
                       MR1_iqr     = numeric(),
                       MR2_median  = numeric(),
                       MR2_iqr     = numeric(),
                       N_paired    = numeric(),
                      `P-value`    = numeric())


#order logrel.df columns in same order as mdp, since I will use t-test on mdp
logrel.df <- logrel.df %>% select(colnames(mdp)[3:74])


#Make a summary table of medians at MR1 and MR2 and the p-value for the t-test of the log(rel-dif)
#loop through all tracts and calculate summary statistics for each
for(i in 3:74)
  {
    name <- colnames(mdp)[[i]]
    x <- mdp[[i]]

    ttest <- t.test(logrel.df[,i-2], mu=0, alternative = "two.sided")

    mr1_median <- subset(mdp,  timepoint == "MR1", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
                  
    mr1_iqr    <- subset(mdp,  timepoint == "MR1", name) %>% 
                  unlist %>% 
                  IQR(., na.rm = TRUE) %>% 
                  round(3)
                  
    mr2_median <- subset(mdp,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
                  
    mr2_iqr    <- subset(mdp,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  IQR(., na.rm = TRUE) %>% 
                  round(3)

    n          <- mdp %>% 
                  as.data.frame %>% 
                  select(subject, name) %>% 
                  na.omit() %>% 
                  group_by(subject) %>% 
                  filter(n()>1) %>% 
                  nrow() %>% 
                  `/`(2)
                  
    #n is the number of patients in whom md can be compared between MR1 and MR2 
    #because the tract is present on both scans
 
    vec       <- c(name, mr1_median, mr1_iqr, mr2_median, mr2_iqr, n, round(ttest$p.value, 4))

    malp_sum  <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

malp_sum[,2:7] <- lapply(malp_sum[,2:7], function(x) as.numeric(as.character(x)))
colnames(malp_sum) <- c("tract","MR1_median", "MR1_iqr", "MR2_median", "MR2_iqr", "N_paired","Unadj. p-value") #to make P.value back to P-value



#To add the differences between scan I need to go back to the individual patient data 
#(note I have used individual patient data for the t-test above also)
malp_sum$MR1_iqr  <- malp_sum$MR1_iqr %>% paste0("(", ., ")")
malp_sum$MR2_iqr  <- malp_sum$MR2_iqr %>% paste0("(", ., ")")
malp_sum          <- malp_sum %>% unite(col = "MR1", MR1_median, MR1_iqr, sep = " ")
malp_sum          <- malp_sum %>% unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")

#Add medians (IQR) for absolute and relative differences between scans to the table
temp <- diff %>% 
        group_by(tract) %>% 
        summarise(Absdiff_median  = round(median(Absdiff, na.rm = TRUE),3), 
                  Absdiff_iqr     = round(IQR(Absdiff, na.rm = TRUE),3), 
                  Reldiff_median  = round(median(Reldiff, na.rm = TRUE),2), 
                  Reldiff_iqr     = round(IQR(Reldiff, na.rm = TRUE),2))

malp_sum <- merge(malp_sum, temp, by = "tract", all = TRUE)

#display medians and IQR in th same cell
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+", malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)
malp_sum$Absdiff_iqr    <- malp_sum$Absdiff_iqr %>% paste0("(", ., ")")
malp_sum$Reldiff_iqr    <- malp_sum$Reldiff_iqr %>% paste0("(", ., ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            select(Tract = tract, 
                                   MR1, 
                                   MR2, 
                                   Absdiff, 
                                   Reldiff, 
                                   N_paired, 
                                   'Unadj. p-value')

#Give tracts meaningful names
tracts    <- read.csv("TrackSeg_names_extended.csv") 
tracts    <- tracts %>% select(Number, Short_name)
malp_sum  <- merge(malp_sum, tracts, by.x = "Tract", by.y = "Number", all.x = TRUE, all.y = FALSE)
malp_sum  <- malp_sum %>% select(Short_name, everything())

#calculate the maximum p-value, among my population of p-values, 
#that is still statistically significant assuming 
#a 5% False discovery rate threshold
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant 
#according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, 
                       "significant", 
                       "not sig.")

#order tracts according to type
myorder <- c(
  "Superior_longitudinal_fascicle_I_left", "Superior_longitudinal_fascicle_I_right",
  "Superior_longitudinal_fascicle_II_left", "Superior_longitudinal_fascicle_II_right",
  "Superior_longitudinal_fascicle_III_left", "Superior_longitudinal_fascicle_III_right",
  "Uncinate_Fascicle_left", "Uncinate_Fascicle_right",
  "Arcuate_Fascicle_Left", "Arcuate_Fascicle_right",
  "Cingulum_left", "Cingulum_right",
  "Inferior_longitudinal_fascicle_left", "Inferior_longitudinal_fascicle_right",
  "Middle_longitudinal_fasicle_left", "Middle_longitudinal_fasicle_right",
  "Inferior_Occipito_frontal_fascicle_left", "Inferior_Occipito_frontal_fascicle_right",
  "Thalamo_prefrontal_left", "Thalamo_prefrontal_right",
  "Thalamo_Premotor_left", "Thalamo_Premotor_right",
  "Thalamo_Precentral_left", "Thalamo_Precentral_right",
  "Thalamo_Postcentral_left", "Thalamo_Postcentral_right",
  "Thalamo_parietal_left", "Thalamo_parietal_right",
  "Thalamo_occipital_left", "Thalamo_occpital_right",
  "Superior_Thalamic_radiation_left", "Superior_Thalamic_radiation_right",
  "Anterior_thalamic_radiation_left", "Anterior_thalamic_radiation_right",
  "Corticospinal_tract_left", "Corticospinal_tract_right",
  "CC_all", "CC_Anterior_Rostrum", "CC_Genu", "CC_Rostral_body_premotor",
  "CC_Anterior_midbody_primary_motor", "CC_isthmus", "CC_splenium",
  "Posterior_midbody_primary_somatosensory", 
  "Commissure Anterior",
  "Fornix_left", "Fornix_right",
  "Striato_fronto_orbital_left", "Striato_fronto_orbital_right",
  "Striato_prefrontal_left", "Striato_prefrontal_right",
  "Striato_premotor_left" , "Striato_premotor_right",
  "Striato_precentral_left", "Striato_precentral_right",
  "Striato_postcentral_left", "Striato_postcentral_right",
  "Striato_parietal_left", "Striato_parietal_right",
  "Striato_occipital_left", "Striato_occipital_right",
  "Parieto_occipital_pontine_left", "Parieto_occipital_pontine_right",
  "Fronto_pontine_tract_left", "Fronto_pontine_tract_right",
  "Superior_cerebellar_peduncle_left", "Superior_cerebellar_peduncle_right",
  "Inferior_cerebellar_peduncle_left", "Inferior_cerebellar_peduncle_right", 
  "Middle_cerebellar_peduncle",
  "Optic_radiation_left", "Optic_radiation_right")

malp_sum$Short_name <- factor(malp_sum$Short_name, 
                              levels = myorder)
malp_sum            <- malp_sum[order(malp_sum$Short_name), ]
malp_sum            <- malp_sum %>% select(-Tract)
names(malp_sum)[names(malp_sum) == "Short_name"] <- "White matter tract"

#Make a pretty table
malp_sum_md1 <- malp_sum

malp_sum_md1 %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  
                     "Mean diffusivity * 1000, \nMedian (IQR)" = 4, 
                     "MR1 vs MR2" = 3)) %>%
  pack_rows("Association fibres", 1, 18) %>%
  pack_rows("Thalamic radiations", 19 , 36) %>%
  pack_rows("Projection fibres", 35 , 36) %>%
  pack_rows("Commissural fibress", 37 , 49) %>%
  pack_rows("Striatal fibres", 48 , 63) %>%
  pack_rows("Brainstem", 62 , 70) %>%
  pack_rows("Other", 71 , 72) 

###4. COMPARE MD AT MR2 WITH MR3

#Make dataframe for subjects that have paired diffusion available for MR2 and MR3
mdp <- md %>% filter(timepoint %in% c("MR2", "MR3"))
mdp <- mdp %>% group_by(subject) %>% filter(n()>1)

#Making a dataframe that shows the differences between MR1 and MR2 as
#absolute difference
#relative difference
#ln of relative difference

mdp <- mdp %>% select(subject, timepoint, X0:X71)
#Multiply all md values by 1000 to avoid small numbers
mdp <- mdp %>% 
  mutate_at(vars(X0:X71),
            .funs = list(~. * 1000))

diff <- mdp %>%  
  gather(key = tract, value = md, X0:X71)
diff <- diff %>% spread(key = timepoint, value = md)
diff <- diff %>% mutate(Absdiff = MR3-MR2)
diff <- diff %>% mutate(Reldiff = MR3/MR2)
diff <- diff %>% mutate(Logrel = log(MR3/MR2))

#spreading the table again with ROIs as columns
#in this table the rows are the patients
#For each patient and tract we have one "logrel" value i.e. log(MR2/MR1)
#summarising the between-scan change in a single value per patient per tract
#allows me to use a one-sample t-test later on
logrel.df <-  diff %>% select(subject, tract, Logrel)
logrel.df <-  logrel.df %>% 
              spread(key = tract, value = Logrel) %>% 
              as.data.frame() %>% 
              select(-subject)

#Prevent R from using exponential notation
options(scipen = 999)

#Initialise a dataframe I will use to collect results in

malp_sum            <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(malp_sum)  <- c("Tract", 
                         "MR2_median", 
                         "MR2_iqr", 
                         "MR3_median", 
                         "MR3_iqr", 
                         "N_paired", 
                         "P-value")

malp_sum <- data.frame(Tract       = character(), 
                       MR1_median  = numeric(),
                       MR1_iqr     = numeric(),
                       MR2_median  = numeric(),
                       MR2_iqr     = numeric(),
                       N_paired    = numeric(),
                      `P-value`    = numeric())


#order logrel.df columns in same order as mdp, since I will use t-test on mdp
logrel.df <- logrel.df %>% select(colnames(mdp)[3:74])


#Make a summary table of medians at MR1 and MR2 for each tract
#and the p-value for the t-test of the log(rel-dif)
#by looping through all tracts

for(i in 3:74)
  {
    name  <- colnames(mdp)[[i]]
    x     <- mdp[[i]]

    ttest <- t.test(logrel.df[,i-2], mu=0, alternative = "two.sided")

    mr2_median <- subset(mdp,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
    
    mr2_iqr    <- subset(mdp,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  IQR(., na.rm = TRUE) %>% 
                  round(3)
    
    mr3_median <- subset(mdp,  timepoint == "MR3", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
    
    mr3_iqr    <- subset(mdp,  timepoint == "MR3", name) %>% 
                  unlist %>% IQR(., na.rm = TRUE) %>% 
                  round(3)

    n          <- mdp %>% 
                  as.data.frame %>% 
                  select(subject, name) %>% 
                  na.omit() %>% 
                  group_by(subject) %>% 
                  filter(n()>1) %>% 
                  nrow() %>% 
                  `/`(2)
    #n is the number of patients in whom md can be compared between MR1 and MR2 
    #because the tract is present on both scans
 
    vec <- c(name, mr2_median, mr2_iqr, mr3_median, mr3_iqr, n, round(ttest$p.value, 4))

    malp_sum <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

malp_sum[,2:7]      <- lapply(malp_sum[,2:7], function(x) as.numeric(as.character(x)))
colnames(malp_sum)  <- c("tract",
                         "MR2_median", 
                         "MR2_iqr", 
                         "MR3_median", 
                         "MR3_iqr", 
                         "N_paired",
                         "Unadj. p-value") 


#To add the differences between scan I need to go back to the individual patient data (
#(note I have used individual patient data for the t-test above also)
malp_sum$MR2_iqr  <- malp_sum$MR2_iqr %>% paste0("(", ., ")")
malp_sum$MR3_iqr  <- malp_sum$MR3_iqr %>% paste0("(", ., ")")
malp_sum          <- malp_sum %>% unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")
malp_sum          <- malp_sum %>% unite(col = "MR3", MR3_median, MR3_iqr, sep = " ")

#Add medians (IQR) for absolute and relative differences between scans to the table
temp <- diff %>% 
        group_by(tract) %>% 
        summarise(Absdiff_median = round(median(Absdiff, na.rm = TRUE),3), 
                  Absdiff_iqr    = round(IQR(Absdiff, na.rm = TRUE),3), 
                  Reldiff_median = round(median(Reldiff, na.rm = TRUE),2), 
                  Reldiff_iqr    = round(IQR(Reldiff, na.rm = TRUE),2))
malp_sum                <- merge(malp_sum, temp, by = "tract", all = TRUE)

#display medians and IQRs in the same cell
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+", malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)
malp_sum$Absdiff_iqr    <- malp_sum$Absdiff_iqr %>% paste0("(", ., ")")
malp_sum$Reldiff_iqr    <- malp_sum$Reldiff_iqr %>% paste0("(", ., ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            select(Tract = tract, 
                                   MR2, 
                                   MR3, 
                                   Absdiff, 
                                   Reldiff, 
                                   N_paired, 
                                   'Unadj. p-value')

#Give tracts meaningful names
tracts <- tracts %>% select(Number, Short_name)
malp_sum <- merge(malp_sum, tracts, by.x = "Tract", by.y = "Number", all.x = TRUE, all.y = FALSE)
malp_sum <- malp_sum %>% select(Short_name, everything())

#False discovery rate for MD
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, "significant", "not sig.")

#order tracts according to type
malp_sum$Short_name <- factor(malp_sum$Short_name, levels = myorder)
malp_sum            <-  malp_sum[order(malp_sum$Short_name), ]
malp_sum            <- malp_sum %>% select(-Tract)
names(malp_sum)[names(malp_sum) == "Short_name"] <- "White matter tract"

malp_sum_md2 <- malp_sum

#Make a pretty table
malp_sum_md2 %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  
                     "Mean diffusivity * 1000, \nMedian (IQR)" = 4, 
                     "MR2 vs MR3" = 3)) %>%
  pack_rows("Association fibres", 1, 18) %>%
  pack_rows("Thalamic radiations", 19 , 36) %>%
  pack_rows("Projection fibres", 35 , 36) %>%
  pack_rows("Commissural fibress", 37 , 49) %>%
  pack_rows("Striatal fibres", 48 , 63) %>%
  pack_rows("Brainstem", 62 , 70) %>%
  pack_rows("Other", 71 , 72) 

###5. COMPARE FA AT MR1 WITH MR2

#Make dataframe for subjects that have paired diffusion available for MR1 and MR2
fap <- fa %>% filter(timepoint %in% c("MR1", "MR2"))
fap <- fap %>% group_by(subject) %>% filter(n()>1)

#Making a dataframe that shows the differences between MR1 and MR2 as
#absolute difference
#relative difference
#ln of relative difference

fap <- fap %>% select(subject, timepoint, X0:X71)

diff <- fap %>%  
  gather(key = tract, value = fa, X0:X71)
diff <- diff %>% spread(key = timepoint, value = fa)
diff <- diff %>% mutate(Absdiff = MR2-MR1)
diff <- diff %>% mutate(Reldiff = MR2/MR1)
diff <- diff %>% mutate(Logrel = log(MR2/MR1))

#spreading the table again with ROIs as columns
logrel.df <- diff %>% select(subject, tract, Logrel)
logrel.df <- logrel.df %>% spread(key = tract, value = Logrel) %>% as.data.frame() %>% select(-subject)

#Prevent R from using exponential notation
options(scipen = 999)

#Making a summary dataframe of values per scan per tract
malp_sum            <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(malp_sum)  <- c("Tract", 
                         "MR1_median", 
                         "MR1_iqr", 
                         "MR2_median", 
                         "MR2_iqr", 
                         "N_paired", 
                         "P-value")

malp_sum <- data.frame(Tract       = character(), 
                       MR1_median  = numeric(),
                       MR1_iqr     = numeric(),
                       MR2_median  = numeric(),
                       MR2_iqr     = numeric(),
                       N_paired    = numeric(),
                      `P-value`    = numeric())


#order logrel.df columns in same order as fap, since I will use t-test on mdp
logrel.df <- logrel.df %>% select(colnames(fap)[3:74])

for(i in 3:74)
  {
    name <- colnames(fap)[[i]]
    x <- fap[[i]]

    ttest <- t.test(logrel.df[,i-2], mu=0, alternative = "two.sided")

    mr1_median <- subset(fap,  timepoint == "MR1", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
    
    mr1_iqr    <- subset(fap,  timepoint == "MR1", name) %>% 
                  unlist %>% 
                  IQR(., na.rm = TRUE) %>% 
                  round(3)
    
    mr2_median <- subset(fap,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  median(., na.rm = TRUE) %>% 
                  round(3)
    
    mr2_iqr    <- subset(fap,  timepoint == "MR2", name) %>% 
                  unlist %>% 
                  IQR(., na.rm = TRUE) %>% 
                  round(3)

    n          <- fap %>% 
                  as.data.frame %>% 
                  select(subject, name) %>% 
                  na.omit() %>% 
                  group_by(subject) %>% 
                  filter(n()>1) %>% 
                  nrow() %>% 
                  `/`(2)
#n is the number of patients in whom md can be compared between MR1 and MR2 because the tract is present on both scans
 
    vec <- c(name, mr1_median, mr1_iqr, mr2_median, mr2_iqr, n, round(ttest$p.value, 4))

    malp_sum <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

malp_sum[,2:7]      <- lapply(malp_sum[,2:7], function(x) as.numeric(as.character(x)))
colnames(malp_sum)  <- c("tract",
                         "MR1_median", 
                         "MR1_iqr", 
                         "MR2_median", 
                         "MR2_iqr", 
                         "N_paired",
                         "Unadj. p-value") 

#To add the differences between scan I need to go back to the individual patient data 
#note I have used individual patient data for the t-test above also)
malp_sum$MR1_iqr  <- malp_sum$MR1_iqr %>% paste0("(", ., ")")
malp_sum$MR2_iqr  <- malp_sum$MR2_iqr %>% paste0("(", ., ")")
malp_sum          <- malp_sum %>% unite(col = "MR1", MR1_median, MR1_iqr, sep = " ")
malp_sum          <- malp_sum %>% unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")

#calculate medians and IQRs for the absolute and relative differences
temp <- diff %>% 
        group_by(tract) %>% 
        summarise(Absdiff_median = round(median(Absdiff, na.rm = TRUE),3), 
                  Absdiff_iqr    = round(IQR(Absdiff, na.rm = TRUE),3), 
                  Reldiff_median = round(median(Reldiff, na.rm = TRUE),2), 
                  Reldiff_iqr    = round(IQR(Reldiff, na.rm = TRUE),2))

malp_sum                <- merge(malp_sum, temp, by = "tract", all = TRUE)
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+", malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)

#Displays medians and IQRs in the same cell
malp_sum$Absdiff_iqr    <- malp_sum$Absdiff_iqr %>% paste0("(", ., ")")
malp_sum$Reldiff_iqr    <- malp_sum$Reldiff_iqr %>% paste0("(", ., ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            select(Tract = tract, 
                                   MR1, 
                                   MR2, 
                                   Absdiff, 
                                   Reldiff, 
                                   N_paired, 
                                   'Unadj. p-value')

#Give tracts meaningful names
tracts    <- tracts %>% select(Number, Short_name)
malp_sum  <- merge(malp_sum, tracts, 
                   by.x = "Tract", 
                   by.y = "Number", 
                   all.x = TRUE, 
                   all.y = FALSE)
malp_sum  <- malp_sum %>% select(Short_name, everything())

#False discovery rate for MD
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, "significant", "not sig.")


#order tracts according to type
malp_sum$Short_name <- factor(malp_sum$Short_name, levels = myorder)
malp_sum            <-  malp_sum[order(malp_sum$Short_name), ]
malp_sum            <- malp_sum %>% select(-Tract)
names(malp_sum)[names(malp_sum) == "Short_name"] <- "White matter tract"

malp_sum_fa1 <- malp_sum

#Make a pretty table
malp_sum_fa1 %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  
                     "Fractional anisotropy, \nMedian (IQR)" = 4, 
                     "MR1 vs MR2" = 3)) %>%
  pack_rows("Association fibres", 1, 18) %>%
  pack_rows("Thalamic radiations", 19 , 36) %>%
  pack_rows("Projection fibres", 35 , 36) %>%
  pack_rows("Commissural fibress", 37 , 49) %>%
  pack_rows("Striatal fibres", 48 , 63) %>%
  pack_rows("Brainstem", 62 , 70) %>%
  pack_rows("Other", 71 , 72) 
###6. COMPARE FA AT MR2 WITH MR3

#Make dataframe for subjects that have paired diffusion available for MR2 and MR3
fap <- fa %>% filter(timepoint %in% c("MR2", "MR3"))
fap <- fap %>% group_by(subject) %>% filter(n()>1)

#Making a dataframe that shows the differences between MR1 and MR2 as
#absolute difference
#relative difference
#ln of relative difference

fap <- fap %>% select(subject, timepoint, X0:X71)

diff <- fap %>%  
  gather(key = tract, value = fa, X0:X71)
diff <- diff %>% spread(key = timepoint, value = fa)
diff <- diff %>% mutate(Absdiff = MR3-MR2)
diff <- diff %>% mutate(Reldiff = MR3/MR2)
diff <- diff %>% mutate(Logrel = log(MR3/MR2))

#spreading the table again with ROIs as columns
logrel.df <- diff %>% select(subject, tract, Logrel)
logrel.df <- logrel.df %>% spread(key = tract, value = Logrel) %>% as.data.frame() %>% select(-subject)

#Prevent R from using exponential notation
options(scipen = 999)

#Making a summary dataframe of values per scan per tract
malp_sum            <- data.frame(matrix(ncol = 7, nrow = 0))
colnames(malp_sum)  <- c("Tract", 
                         "MR2_median", 
                         "MR2_iqr", 
                         "MR3_median", 
                         "MR3_iqr", 
                         "N_paired", 
                         "P-value")

malp_sum <- data.frame(Tract       = character(), 
                       MR1_median  = numeric(),
                       MR1_iqr     = numeric(),
                       MR2_median  = numeric(),
                       MR2_iqr     = numeric(),
                       N_paired    = numeric(),
                      `P-value`    = numeric())


#order logrel.df columns in same order as fap, since I will use t-test on mdp
logrel.df <- logrel.df %>% select(colnames(fap)[3:74])

for(i in 3:74)
  {
  name <- colnames(fap)[[i]]
  x <- fap[[i]]

  ttest <- t.test(logrel.df[,i-2], mu=0, alternative = "two.sided")

  mr1_median <- subset(fap,  timepoint == "MR2", name) %>% 
                unlist %>% 
                median(., na.rm = TRUE) %>% 
                round(3)
  
  mr1_iqr    <- subset(fap,  timepoint == "MR2", name) %>% 
                unlist %>% 
                IQR(., na.rm = TRUE) %>% 
                round(3)
  
  mr2_median <- subset(fap,  timepoint == "MR3", name) %>% 
                unlist %>% 
                median(., na.rm = TRUE) %>% 
                round(3)
  
  mr2_iqr    <- subset(fap,  timepoint == "MR3", name) %>% 
                unlist %>% 
                IQR(., na.rm = TRUE) %>% 
                round(3)

  n          <- fap %>% 
                as.data.frame %>% 
                select(subject, name) %>% 
                na.omit() %>% 
                group_by(subject) %>% 
                filter(n()>1) %>% 
                nrow() %>% 
                `/`(2)
  #n is the number of patients in whom md can be compared between MR1 and MR2 
  #because the tract is present on both scans
 
  vec       <- c(name, mr2_median, mr2_iqr, mr3_median, mr3_iqr, n, round(ttest$p.value, 4))

  malp_sum  <- berryFunctions::addRows(malp_sum, 1, values = vec)
  }

malp_sum[,2:7]      <- lapply(malp_sum[,2:7], function(x) as.numeric(as.character(x)))
colnames(malp_sum)  <- c("tract",
                         "MR2_median", 
                         "MR2_iqr", 
                         "MR3_median", 
                         "MR3_iqr", 
                         "N_paired",
                         "Unadj. p-value") 

#To add the differences between scan I need to go back to the individual patient data 
#note I have used individual patient data for the t-test above also
malp_sum$MR2_iqr  <- malp_sum$MR2_iqr %>% paste0("(", ., ")")
malp_sum$MR3_iqr  <- malp_sum$MR3_iqr %>% paste0("(", ., ")")
malp_sum          <- malp_sum %>% unite(col = "MR2", MR2_median, MR2_iqr, sep = " ")
malp_sum          <- malp_sum %>% unite(col = "MR3", MR3_median, MR3_iqr, sep = " ")

#calculate medians and IQRs for the absolute and relative differences between scans
temp <- diff %>% 
        group_by(tract) %>% 
        summarise(Absdiff_median = round(median(Absdiff, na.rm = TRUE),3), 
                  Absdiff_iqr     = round(IQR(Absdiff, na.rm = TRUE),3), 
                  Reldiff_median = round(median(Reldiff, na.rm = TRUE),2), 
                  Reldiff_iqr     = round(IQR(Reldiff, na.rm = TRUE),2))
malp_sum                <- merge(malp_sum, temp, by = "tract", all = TRUE)
malp_sum$Absdiff_median <- ifelse(malp_sum$Absdiff_median > 0, 
                                  paste0("+", malp_sum$Absdiff_median), 
                                  malp_sum$Absdiff_median)

#display medians and IQRs in the same cell
malp_sum$Absdiff_iqr    <- malp_sum$Absdiff_iqr %>% paste0("(", ., ")")
malp_sum$Reldiff_iqr    <- malp_sum$Reldiff_iqr %>% paste0("(", ., ")")
malp_sum                <- malp_sum %>% 
                            unite(col = "Absdiff", Absdiff_median, Absdiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            unite(col = "Reldiff", Reldiff_median, Reldiff_iqr, sep = " ")
malp_sum                <- malp_sum %>% 
                            select(Tract = tract, 
                                   MR2, 
                                   MR3, 
                                   Absdiff, 
                                   Reldiff, 
                                   N_paired, 
                                   'Unadj. p-value')

#Give tracts meaningful names
tracts    <- tracts %>% select(Number, Short_name)
malp_sum  <- merge(malp_sum, tracts, 
                   by.x = "Tract", 
                   by.y = "Number", 
                   all.x = TRUE, 
                   all.y = FALSE)
malp_sum  <- malp_sum %>% select(Short_name, everything())

#False discovery rate for MD
maxp <- fdr_maxp(malp_sum$`Unadj. p-value`, 0.05)

#Add column saying whether a tract is significant according to 5% false discovery rate
malp_sum$FDR <- ifelse(malp_sum$`Unadj. p-value` <= maxp, "significant", "not sig.")


#order tracts according to type
malp_sum$Short_name <- factor(malp_sum$Short_name, levels = myorder)
malp_sum            <- malp_sum[order(malp_sum$Short_name), ]
malp_sum            <- malp_sum %>% select(-Tract)
names(malp_sum)[names(malp_sum) == "Short_name"] <- "White matter tract"

malp_sum_fa2 <- malp_sum

#Make a pretty table
malp_sum_fa2 %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) %>%
  add_header_above(c(" " = 1,  
                     "Fractional anisotropy, \nMedian (IQR)" = 4, 
                     "MR2 vs MR3" = 3)) %>%
  pack_rows("Association fibres", 1, 18) %>%
  pack_rows("Thalamic radiations", 19 , 36) %>%
  pack_rows("Projection fibres", 35 , 36) %>%
  pack_rows("Commissural fibress", 37 , 49) %>%
  pack_rows("Striatal fibres", 48 , 63) %>%
  pack_rows("Brainstem", 62 , 70) %>%
  pack_rows("Other", 71 , 72) 

###7. COMPARE PATIENTS WITH CONTROLS

##7.1 Prepare data

#Make a list of significant tracts i.e those that changed within patients between MR1 and MR2
mytracts    <- malp_sum_md1 %>% filter(FDR == "significant")

#Make a list of additional tracts based on literature
lit         <- c("CC_all")
lit_tracts  <- malp_sum_md1[which(malp_sum_md1$`White matter tract` %in% lit),]

#combine lists of tracts
mytracts    <- rbind(mytracts, lit_tracts) %>% as.data.frame()

#I will compare the above 13 tracts between patients and controls for FA and MD values, 
#at MR1 and MR2.
#I will get the column names (as X...) for the above tracts
mytracts    <- merge(mytracts, tracts, 
                     by.x = "White matter tract", 
                     by.y = "Short_name", 
                     all.x = TRUE, 
                     all.y = FALSE)
mytracts
write.csv(mytracts, file = "mytracts.csv", row.names = FALSE)
## 7.2 Testing model assumptions using diagnostic plots
#Note I am multiplying all MD values by 1000 to avoid having to operate with very tiny numbers

#Normality assumption MD
md_MR1 <- md %>% filter(timepoint %in% c("control", "MR1"))

par(mfrow=c(3,5))

for (i in 1:14) 
  { 
    select.me <-  c(mytracts$Number)
    Data      <-  md_MR1[,select.me]
    Tract     <-  colnames(Data)[[i]]
    x         <-  1000

    fit       <- lmer(Data[[i]]*x ~ group + age + sex + (1| scanner), 
                      data = md_MR1, 
                      REML = FALSE)

    qqnorm(residuals(fit))
    qqline(residuals(fit))
  }

#Normality assumption FA
fa_MR1 <- fa %>% filter(timepoint %in% c("control", "MR1"))

par(mfrow=c(3,5))

for (i in 1:14) 
  { 
    select.me <-  c(mytracts$Number)
    Data      <-  fa_MR1[,select.me]
    Tract     <-  colnames(Data)[[i]]
    x         <-  1000

    fit       <- lmer(Data[[i]]*x ~ group + age + sex + (1| scanner), 
                      data = fa_MR1, 
                      REML = FALSE)

    qqnorm(residuals(fit))
    qqline(residuals(fit))
  }
 
#Homoscedasticity assumption MD
md_MR1 <- md %>% filter(timepoint %in% c("control", "MR1"))

for (i in 1:14) 
  { 
    select.me <-  c(mytracts$Number)
    Data      <-  md_MR1[,select.me]
    Tract     <-  colnames(Data)[[i]]
    x         <-  1000

    fit       <- lmer(Data[[i]]*x ~ group + age + sex + (1| scanner), 
                      data = md_MR1, 
                      REML = FALSE)

    plot(residuals(fit))
  }
 

#Homoscedasticity assumption FA
fa_MR1 <- fa %>% filter(timepoint %in% c("control", "MR1"))

for (i in 1:14) 
  { 
    select.me <-  c(mytracts$Number)
    Data      <-  fa_MR1[,select.me]
    Tract     <-  colnames(Data)[[i]]
    x         <-  1000

    fit       <- lmer(Data[[i]]*x ~ group + age + sex + (1| scanner), 
                      data = fa_MR1, 
                      REML = FALSE)

    plot(residuals(fit))
  }
 

#7.3 MD at MR1

#Making a dataframe to collect results of pt vs ctr comparison
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract       = character(), 
                 Metric      = character(),
                 Coef_pt     = numeric(),
                 SE          = numeric(),
                 P_value     = numeric())

md_MR1 <- md %>% 
  filter(timepoint %in% c("control", "MR1")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
    Data    <- md_MR1
    Metric  <- "md_MR1"
    Tract   <- colnames(Data)[[i]]

    library(lmerTest)

    fit     <- lmer(Data[[i]]*1000 ~ group + age + sex + (1| scanner), 
                    data = Data, 
                    REML = FALSE,
                    control = lmerControl(optimizer ="Nelder_Mead"))
  
    Coef_pt <- summary(fit)$coefficients[2,1] %>% round(4)
    SE      <- summary(fit)$coefficients[2,2] %>% round(4)
    P_value <- summary(fit)$coefficients[2,5] %>% round(6)

 
    vec     <- c(Tract, Metric, Coef_pt, SE, P_value)

    df      <- berryFunctions::addRows(df, 1, values = vec)
  }

#format results table
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, 
                      paste0("+", df$Coef_pt), 
                      paste0(df$Coef_pt))

#add a column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, 
                      "significant", 
                      "not sig.")
#more formatting (incl substitute the tract number with the tract name)
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract", 
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% 
                select(-Tract) %>% 
                select(White_matter_tract = Short_name, everything())

myorder <- c(
  "Superior_longitudinal_fascicle_I_left", "Superior_longitudinal_fascicle_I_right",
  "Superior_longitudinal_fascicle_II_left", "Superior_longitudinal_fascicle_II_right",
  "Superior_longitudinal_fascicle_III_left",
  "Arcuate_Fascicle_Left",
  "Cingulum_left", 
  "Middle_longitudinal_fasicle_left", 
  "Thalamo_Precentral_left", 
  "Corticospinal_tract_right",
  "Fornix_left", 
  "Striato_precentral_left", 
  "Striato_parietal_left", 
  "CC_all"
)

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <-  df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_md_mr1 <- df
df_md_mr1%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 


#7.4 MD at MR2

#Making a dataframe to collect results of pt vs ctr comparison
df <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df) <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract      = character(), 
                 Metric      = character(),
                 Coef_pt     = numeric(),
                 SE          = numeric(),
                 P_value     = numeric())

md_MR2 <- md %>% 
  filter(timepoint %in% c("control", "MR2")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
  Data    <- md_MR2
  Metric  <- "md_MR2"
  Tract   <- colnames(Data)[[i]]

  library(lmerTest)

  fit     <- lmer(Data[[i]]*1000 ~ group + age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
  Coef_pt <- summary(fit)$coefficients[2,1] %>% round(4)
  SE      <- summary(fit)$coefficients[2,2] %>% round(4)
  P_value <- summary(fit)$coefficients[2,5] %>% round(8)

  vec     <- c(Tract, Metric, Coef_pt, SE, P_value)

  df      <- berryFunctions::addRows(df, 1, values = vec)
  }

#format dataframe
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, paste0("+", df$Coef_pt), paste0(df$Coef_pt))

#add column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

#more formatting including substituting tract numbers with tract names
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract",
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% 
                select(-Tract) %>% 
                select(White_matter_tract = Short_name, everything())

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <-  df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_md_mr2 <- df
df_md_mr2%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

#7.5 MD at MR3

#Making a dataframe to collect results of pt vs ctr comparison
df            <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df)  <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract      = character(), 
                 Metric     = character(),
                 Coef_pt    = numeric(),
                 SE         = numeric(),
                 P_value    = numeric())
md_MR3 <- md %>% 
  filter(timepoint %in% c("control", "MR3")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
  Data <- md_MR3
  Metric <- "md_MR3"
  Tract <- colnames(Data)[[i]]

  library(lmerTest)

  fit     <- lmer(Data[[i]]*1000 ~ group + age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
  Coef_pt <- summary(fit)$coefficients[2,1] %>% round(4)
  SE      <- summary(fit)$coefficients[2,2] %>% round(4)
  P_value <- summary(fit)$coefficients[2,5] %>% round(8)

  vec     <- c(Tract, Metric, Coef_pt, SE, P_value)

  df      <- berryFunctions::addRows(df, 1, values = vec)
  }

#format the table
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, 
                      paste0("+", df$Coef_pt), 
                      paste0(df$Coef_pt))

#add a column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

#more formatting including substituting tract numbers by tract names
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract", 
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% 
                select(-Tract) %>% 
                select(White_matter_tract = Short_name, everything())

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <- df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_md_mr2 <- df
df_md_mr2%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

#7.6 FA at MR1

#Making a dataframe to collect results of pt vs ctr comparison
df            <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df)  <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract       = character(), 
                 Metric      = character(),
                 Coef_pt     = numeric(),
                 SE          = numeric(),
                 P_value     = numeric())
fa_MR1 <- fa %>% 
  filter(timepoint %in% c("control", "MR1")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
  Data <- fa_MR1
  Metric <- "fa_MR1"
  Tract <- colnames(Data)[[i]]

  library(lmerTest)

  fit     <- lmer(Data[[i]] ~ group + age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
  Coef_pt <- summary(fit)$coefficients[2,1] %>% round(6)
  SE      <- summary(fit)$coefficients[2,2] %>% round(4)
  P_value <- summary(fit)$coefficients[2,5] %>% round(6)
 
  vec     <- c(Tract, Metric, Coef_pt, SE, P_value)

  df      <- berryFunctions::addRows(df, 1, values = vec)
}

#formatting the table
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, paste0("+", df$Coef_pt), paste0(df$Coef_pt))

#add column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

#more formatting including substituting tract numbers with tract names
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract", 
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% select(-Tract) %>% select(White_matter_tract = Short_name, everything())

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <-  df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_fa_mr1 <- df
df_fa_mr1%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

#7.7 FA at MR2

#Making a dataframe to collect results of pt vs ctr comparison
df            <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df)  <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract      = character(), 
                 Metric     = character(),
                 Coef_pt    = numeric(),
                 SE         = numeric(),
                 P_value    = numeric())
fa_MR2 <- fa %>% 
  filter(timepoint %in% c("control", "MR2")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
  Data      <- fa_MR2
  Metric    <- "fa_MR2"
  Tract     <- colnames(Data)[[i]]

  library(lmerTest)

  fit       <- lmer(Data[[i]] ~ group + age + sex + (1| scanner), 
                    data = Data, 
                    REML = FALSE, 
                    control = lmerControl(optimizer ="Nelder_Mead"))
  
  Coef_pt   <- summary(fit)$coefficients[2,1] %>% round(6)
  SE        <- summary(fit)$coefficients[2,2] %>% round(4)
  P_value   <- summary(fit)$coefficients[2,5] %>% round(6)

  vec       <- c(Tract, Metric, Coef_pt, SE, P_value)

  df        <- berryFunctions::addRows(df, 1, values = vec)
  }

#formatting the table
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, paste0("+", df$Coef_pt), paste0(df$Coef_pt))

#adding a column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

#more formatting including substituting tract numbers with tract names
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract", 
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% select(-Tract) %>% select(White_matter_tract = Short_name, everything())

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <-  df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_fa_mr2 <- df
df_fa_mr2%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 


#7.8 FA at MR3

#Making a dataframe to collect results of pt vs ctr comparison
df            <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(df)  <- c("Tract", "Data", "Coef_pt", "SE", "P_value")

df <- data.frame(Tract      = character(), 
                 Metric     = character(),
                 Coef_pt    = numeric(),
                 SE         = numeric(),
                 P_value    = numeric())

fa_MR3 <- fa %>% 
  filter(timepoint %in% c("control", "MR3")) %>%
  select(subject, group, age, sex, scan, scanner, timepoint, as.character(mytracts$Number))

#loop through all tracts, to fit the model onto each tract
#and collect results for each tract in a vector "vec"
#bind vectors together (as rows) of a results dataframe df
for (i in 8:21) 
  {
  Data    <- fa_MR3
  Metric  <- "fa_MR3"
  Tract   <- colnames(Data)[[i]]

  library(lmerTest)

  fit     <- lmer(Data[[i]] ~ group + age + sex + (1| scanner), 
                  data = Data, 
                  REML = FALSE, 
                  control = lmerControl(optimizer ="Nelder_Mead"))
  
  Coef_pt <- summary(fit)$coefficients[2,1] %>% round(6)
  SE      <- summary(fit)$coefficients[2,2] %>% round(4)
  P_value <- summary(fit)$coefficients[2,5] %>% round(6)

  vec     <- c(Tract, Metric, Coef_pt, SE, P_value)

  df      <- berryFunctions::addRows(df, 1, values = vec)
  }

#formatting the table
df[,3:5]    <- lapply(df[,3:5], function(x) as.numeric(as.character(x)))
df$Coef_pt  <- ifelse(df$Coef_pt > 0, paste0("+", df$Coef_pt), paste0(df$Coef_pt))

#adding a column denoting statistical significance
maxp        <- fdr_maxp(df$P_value, 0.05)
df$FDR      <- ifelse(df$P_value <= maxp, "significant", "not sig.")

#more formatting including substituting tract numbers by tract names
df          <- df %>% separate(Metric, c("Metric", "Timepoint"))
df          <- merge(df, tracts, 
                     by.x = "Tract", 
                     by.y = "Number", 
                     all.x = TRUE, 
                     all.y = FALSE)
df          <- df %>% select(-Tract) %>% select(White_matter_tract = Short_name, everything())

df$White_matter_tract <- factor(df$White_matter_tract, levels = myorder)
df                    <-  df[order(df$White_matter_tract), ]

#make a pretty version of the table
df_fa_mr2 <- df
df_fa_mr2%>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "striped"), full_width = F) 

