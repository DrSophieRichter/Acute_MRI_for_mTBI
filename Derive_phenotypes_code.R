#load libraries
library(dplyr)
library(tidyr)
library(factoextra)
library(NbClust)

##load demographic data and outcome data
  #There is one row for each patient ("subject") in this dataframe

demog <- read.csv("Paper1_demog.csv") %>% 
         select(subject,
                Age,
                Sex,
                hours, # this is the time between MR1 and MR2
                Base_RPQ,
                X2wk_RPQ,
                GOSE)


##load diffusion data
  #this data is in the long format with up to three rows per patient (MR1, MR2, MR3)
  #I only need data for MR1 and MR2
  #The columns are the white matter tracts (identified by number) 
  #plus some other columns (subject ID, scan ID, timepoint etc)

md    <- read.csv("md_allsubjects_alltracts_nooutliers.csv") %>% 
         filter(timepoint %in% c("MR1", "MR2")) 

fa    <- read.csv("fa_allsubjects_alltracts_nooutliers.csv") %>% 
         filter(timepoint %in% c("MR1", "MR2")) 


##load the list of tracts of interest
  #This dataframe contains 14 rows 
  #i.e. the 13 white matter tracts that changed significantly between MR1 and MR2
  #plus the corpus callosum
  #I will only use the 13 tracts that changed significantly
tracts <- read.csv("mytracts.csv") %>% 
          filter(FDR == "significant")

#MAKING A WHOLE_BRAIN SUMMARY MEASURE FOR FA & MD

#selecting the tracts that changed between MR1 and MR2
md <- md[,c("subject", "scanner", "timepoint", as.character(tracts$Number))]
fa <- fa[,c("subject", "scanner", "timepoint", as.character(tracts$Number))]

#calculating differences in MD between MR1 and MR2 
  #for each tract within each patient
mddif <- md %>% 
          gather(key = Tract, value = MD, X0:X50) %>%
          spread(key = timepoint, value = MD)
mddif <- mddif %>% mutate(Absdiff = MR2-MR1)
mddif <- mddif %>% mutate(Reldiff = MR2/MR1)
mddif <- mddif %>% mutate(Logrel  = log(MR2/MR1))

#Now summarising these differences per patient across 
  #all of that patient's tracts
  #to obtain a whole brain measure of change between MR1 and MR2
mdsum <- mddif %>% 
          group_by(subject) %>% 
          summarise(deltaMD     = sum(Absdiff),  
                    logratioMD  = sum(Logrel),
                    MD_sumMR1   = sum(MR1),
                    MD_sumMR2   = sum(MR2)) 

#Adding whole-brain change in md to dataframe demog
demog <- merge(demog, mdsum, 
               by = "subject", 
               all.x = TRUE, 
               all.y = FALSE)

#repeat for FA

#calculating differences between MR1 and MR2 per patient per tract
fadif <- fa %>% 
          gather(key = Tract, value = FA, X0:X50) %>%
          spread(key = timepoint, value = FA)
fadif <- fadif %>% mutate(Absdiff = MR2-MR1)
fadif <- fadif %>% mutate(Reldiff = MR2/MR1)
fadif <- fadif %>% mutate(Logrel  = log(MR2/MR1))

#Now summarising these differences per patient across all their tracts
fasum <- fadif %>% 
          group_by(subject) %>% 
          summarise(deltaFA    = sum(Absdiff),  
                    logratioFA = sum(Logrel),
                    FA_sumMR1  = sum(MR1),
                    FA_sumMR2  = sum(MR2)) 

#Adding summary measures to dataframe demog
demog <- merge(demog, fasum, 
            by = "subject", 
            all.x = TRUE, 
            all.y = FALSE)


#MAKING A MEASURE OF SYMPTOM EVOLUTION BETWEEN SCANS

#Make a variable for delta RPQ between MR1 and MR2
demog <- demog %>% mutate(deltaRPQ = X2wk_RPQ - Base_RPQ)

#number of patients with data available for both
  #symptom evolution and
  #imaging evolution
demog %>% 
  filter(is.na(demog$logratioFA) == FALSE) %>% 
  filter(is.na(.$deltaRPQ) == FALSE)  %>%
  nrow()
str(demog)

#CLUSTERING THE CHANGE IN FA AND MD BETWEEN MR1 AND MR2

#set seeed for reproducible results
set.seed(123)

#select patients with diffusion data available for both timepoints
kdata         <- demog %>% 
                 filter(is.na(demog$logratioFA) == FALSE) %>% 
                 select(subject, logratioMD, logratioFA)

##choose the right number of clusters
  #Note clinically 3 clusters would make sense
  #Let's see if the data support this
df            <- scale(kdata[,2:3])

# Elbow method
fviz_nbclust(df, kmeans, method = "wss") +
              geom_vline(xintercept = 3, linetype = 2)+
              labs(subtitle = "Elbow method")

# Silhouette method
fviz_nbclust(df, kmeans, method = "silhouette")+
              labs(subtitle = "Silhouette method")

#Derive clusters (I decided to use 3 clusters)
clusters        <- kmeans(kdata[,2:3], centers = 3, nstart = 25)
str(clusters)

#add a column that denotes the cluster each patient belongs to
kdata$phenotype <- as.factor(clusters$cluster)

#Add cluster ("imaging phenotype") to demog table 
#to allow plotting against outcome
demog           <- merge(demog, kdata[, c(1, 4)], 
                         by = "subject", 
                         all.x = TRUE,
                         all.y = TRUE)

#Inspect phenotype 
demog %>% 
  group_by(phenotype) %>% 
  summarise(FA = mean(logratioFA) %>% round(2), 
            MD = mean(logratioMD) %>% round(2))

#Name phenotypes
demog$phenotype <- factor(demog$phenotype,
                          levels = c(1, 
                                     2, 
                                     3), 
                          labels = c("No DTI change",
                                     "Pseudonormalisation",
                                     "Progressive injury"))


write.csv(demog, file = "DTI_phenotypes.csv", row.names = FALSE)

