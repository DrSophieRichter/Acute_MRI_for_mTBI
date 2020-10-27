#load libraries
library(dplyr)
library(ggplot2)
library(gridExtra)
library(grid)
library(kableExtra)
library(knitr)

##load data on symptom evolution and imaging evolution between scans
  #This data has 1 row per patient
  #It contains columns for subject ID, demographics and RPQ
  #including RPQ at MR1, at MR2 and the delta between scans
  #it also contains the imaging phenotypes for each patient
  #provided the patient had diffusion data for both scans
data      <- read.csv("DTI_phenotypes.csv") 
data      <- data %>% filter(is.na(data$phenotype) == FALSE)


#Load data on the comparison of patients vs controls
  #the table has 3 rows, one for each phenotype
  #the most important column is Num_abnormal
  #which states the number of tracts (out of 13) that
  #are significantly abnormal compared to healthy controls 
phenoMR1  <- read.csv("DTIPheno_MR1.csv")


###PREPARE DATA FOR PLOTTING

#For the barplot I want a very thin line of coloured bar
#to be displayed even when Num_abnormal = zero
phenoMR1$Num_abnormal <- ifelse(phenoMR1$Num_abnormal < 1, 
                                0.1, 
                                phenoMR1$Num_abnormal)

phenoMR1$Pheno        <- factor(phenoMR1$Pheno,
                                levels = c("prog", 
                                           "static", 
                                           "pseud"),
                                labels = c("Deteriorate", 
                                           "Minimal", 
                                           "Pseudonormalise"))

#Order levels of variable phenotype so that they are consistent in all plots
data$phenotype <- factor(data$phenotype,
                          levels =  c("Progressive injury",
                                     "No DTI change",
                                     "Pseudonormalisation"),
                          labels = c("Progressive injury",
                                     "No DTI change",
                                     "Pseudonormalisation"))


###MAKE FIGURE 1 WITH LABELLED PANELS


#Make a scatterplot of change in FA against change in MD, coloured by phenotype
p1 <- ggplot(data, aes(x=logratioFA, y=logratioMD, colour = phenotype))+
      geom_hline(yintercept = 0) +
      geom_vline(xintercept = 0) +
      geom_point(size = 1) +
      theme_bw(base_size = 6) +
      ggtitle("B. DTI evolution") + 
      theme(legend.position="none",
            plot.title = element_text(size = 9),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
      xlab("Change in FA (logratio)") + 
      ylab("Change in MD (logratio)")+  
      scale_colour_manual(values=c("red3", "steelblue", "orange"))

#Make a bar chart showing the number of abnormal tracts at MR1 for each phenotype
p2 <- ggplot(phenoMR1, aes(x = Pheno, y = Num_abnormal,group = Metric)) +
      geom_bar(colour = "black", aes(fill = Pheno), position = "dodge", stat="identity") +
      ggtitle("A. Initial DTI") +
      theme_bw(base_size = 6) +
      scale_fill_manual(values=c("red3", "steelblue", "orange"))+ 
      theme(legend.position="none", 
            plot.title = element_text(size = 9),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6)) +
      xlab("White matter changes") + 
      ylab("No. of abnormal tracts") +
      scale_y_continuous(limits = c(0, 13), breaks = seq(0, 13, 1)) +
      geom_text(aes(label=Metric),size = 2, position=position_dodge(width=0.9), vjust=-0.25)+
      scale_x_discrete(labels = c("Prog", "Min", "Pseud"))

#Make a boxplot showing the change in RPQ scores for each phenotype
p3 <- ggplot(data, aes(x=phenotype, y=deltaRPQ, fill=phenotype)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      theme_bw(base_size = 6) +
      ggtitle("C. Symptom evolution")+ 
      theme(legend.position="bottom",
            legend.title = element_text (size = 9),
            legend.text = element_text(size = 9),
            plot.title = element_text(size = 9),
            axis.title.x = element_text(size = 6),
            axis.title.y = element_text(size = 6),
            axis.text.x = element_text(size = 6),
            axis.text.y = element_text(size = 6))+
      xlab("White matter changes") + 
      ylab("delta RPQ") +
      scale_fill_manual(values=c("red3", "steelblue", "orange", "grey80"),
                        name="Imaging phenotype",
                        labels=c("Progressive injury (Prog)", 
                                 "Minimal change (Min)", 
                                 "Pseudonormalisation (Pseud)")) +
      scale_x_discrete(labels = c("Prog", "Min", "Pseud"))

get_legend  <- function(myggplot)
  {
  tmp       <- ggplot_gtable(ggplot_build(myggplot))
  leg       <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend    <- tmp$grobs[[leg]]
  return(legend)
  }

legend      <- get_legend(p3)

p3          <- p3 + theme(legend.position="none")

myplot <- grid.arrange(p2, p1, p3, 
             legend, 
             ncol = 3, 
             nrow = 2, 
             layout_matrix = rbind(c(1,2,3), 
                                   c(1,2,3),
                                   c(1,2,3),
                                   c(1,2,3),
                                   c(1,2,3),
                                   c(4,4,4)))

ggsave(plot = myplot, 
       filename = "Fig2.tiff",
       width = 18.5, height = 8, dpi = 300,
       units = "cm")



### TEST FOR ASSOCIATION BETWEEN PHENOTYPE AND DELTA RPQ

#medians and IQR for symptom evolution by phenotype
ggplot_build(p3)$data


#test for statistical significance

#simple linear model akin to Anova 
fit <- lm(deltaRPQ ~ phenotype, data = data) 

#checking model assumptions
plot(fit)

#checking if there is a difference between groups
anova(fit) 
### TEST IF BASELINE RPQs DIFFERED BETWEEN PHENOTYPES

#plot data
p4 <- ggplot(data, aes(x=phenotype, y=Base_RPQ, fill=phenotype)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      theme_bw(base_size = 10) +
      ggtitle("Baseline RPQ")+ 
      theme(legend.position="bottom")+
      xlab("White matter changes") + ylab("Baseline RPQ") +
      scale_fill_manual(values=c("red3", "steelblue", "orange", "grey80"),
                        name="Imaging phenotype",
                        labels=c("Progressive injury (Prog)", 
                                 "Minimal change (Min)", 
                                 "Pseudonormalisation (Pseud)")) +
    scale_x_discrete(labels = c("Prog", "Min", "Pseud"))
p4


#get medians (IQR) for baseline RPQ by phenotype
ggplot_build(p4)$data


#test for statistical significance

#simple linear model akin to Anova 
fit <- lm(Base_RPQ ~ phenotype, data = data) 

#checking model assumptions
plot(fit)

#checking if there is a difference between groups
anova(fit) 

### TEST IF 2-3 WEEK RPQs DIFFERED BETWEEN PHENOTYPES

#plot data
p5 <- ggplot(data, aes(x=phenotype, y=X2wk_RPQ, fill=phenotype)) +
      geom_hline(yintercept = 0) +
      geom_boxplot() +
      theme_bw(base_size = 10) +
      ggtitle("2-week RPQ")+ 
      theme(legend.position="bottom")+
      xlab("White matter changes") + ylab("Baseline RPQ") +
      scale_fill_manual(values=c("red3", "steelblue", "orange", "grey80"),
                        name="Imaging phenotype",
                        labels=c("Progressive injury (Prog)", 
                                 "Minimal change (Min)", 
                                 "Pseudonormalisation (Pseud)")) +
    scale_x_discrete(labels = c("Prog", "Min", "Pseud"))
p5

#get medians (IQR) for 2-week RPQ by phenotype
ggplot_build(p5)$data

#Test for statistical significance

#simple linear model akin to Anova 
fit <- lm(X2wk_RPQ ~ phenotype, data = data) 

#checking model assumptions
plot(fit)

#checking if there is a difference between groups
anova(fit) 

###SENSITIVITY ANALYSIS FOR MISSING RPQ DATA

#worst case delta RPQ imputation
data$deltaRPQ_worst <- ifelse(is.na(data$deltaRPQ)== TRUE, 5, data$deltaRPQ)

#best case delta RPQ imputation
data$deltaRPQ_best <- ifelse(is.na(data$deltaRPQ) == TRUE, -4.5, data$deltaRPQ)


##Calculate medians and quantiles using ggplot

#Complete case analysis
p_complete        <- ggplot(data, aes(x=phenotype, y=deltaRPQ, fill=phenotype)) +
                      geom_hline(yintercept = 0) +
                      geom_boxplot() +
                      scale_fill_manual(values=c("Prog", "Min", "Pseud", "grey8"))

d_complete        <- ggplot_build(p_complete)$data %>% 
                      as.data.frame() %>%
                      select("fill", "middle", "lower", "upper")

d_pheno           <- d_complete$fill

d_complete        <- format(round(d_complete[,-1], digits=2), nsmall = 2) 

d_complete$median <- paste0(d_complete$middle, " [", d_complete$lower, "-", d_complete$upper, "]")


#Worst case analysis
p_worst           <- ggplot(data, aes(x=phenotype, y=deltaRPQ_worst, fill=phenotype)) +
                      geom_hline(yintercept = 0) +
                      geom_boxplot() +
                      scale_fill_manual(values=c("Prog", "Min", "Pseud", "grey8"))

d_worst           <- ggplot_build(p_worst)$data %>% 
                      as.data.frame()%>%
                      select("fill", "middle", "lower", "upper")

d_worst           <- format(round(d_worst[,-1], digits=2), nsmall = 2) 

d_worst$median    <- paste0(d_worst$middle, " [", d_worst$lower, "-", d_worst$upper, "]")


#Best case analysis
p_best            <- ggplot(data, aes(x=phenotype, y=deltaRPQ_best, fill=phenotype)) +
                      geom_hline(yintercept = 0) +
                      geom_boxplot() +
                      scale_fill_manual(values=c("Prog", "Min", "Pseud", "grey8"))

d_best            <- ggplot_build(p_best)$data %>% 
                      as.data.frame() %>%
                      select("fill", "middle", "lower", "upper")

d_best            <-  format(round(d_best[,-1], digits=2), nsmall = 2) 

d_best$median     <- paste0(d_best$middle, " [", d_best$lower, "-", d_best$upper, "]")



#Combine data from the three scenarios into one dataframe
#the dataframe has 3 rows, one for each phenotype
#the dataframe has 4 columns:

sens              <- cbind(d_pheno,           #column denoting phenotype
                           d_complete$median, #contains median (Q1-Q3)
                           d_worst$median,    #contains median (Q1-Q3)
                           d_best$median) %>% #contains median (Q1-Q3)
                    as.data.frame()

colnames(sens) <- c("Phenotype", 
                    "Complete case analysis", 
                    "Worst-case scenario", 
                    "Best-case scenario")

#Make a fourth row that denotes whether the difference between phenotypes
#is significant for each of the three scenarios

fit         <- lm(deltaRPQ ~ phenotype, data = data) 
p           <- anova(fit) %>% as.data.frame()
p_complete  <- p[1,5] %>% round(3)

fit         <- lm(deltaRPQ_worst ~ phenotype, data = data) 
p           <- anova(fit) %>% as.data.frame()
p_worst     <- p[1,5] %>% round(3)

fit         <- lm(deltaRPQ_best ~ phenotype, data = data) 
p           <- anova(fit) %>% as.data.frame()
p_best      <- p[1,5] %>% round(3)

p           <- c("P-value", p_complete, p_worst, p_best)

sens        <- rbind(sens, p)
sens

#Make a pretty version of the sensitivity analysis table
sens %>% 
  kable(escape = FALSE, row.names = FALSE) %>%
  kable_styling(bootstrap_options = c("condensed", "bordered"), full_width = F) %>%
  add_header_above(c(" " = 1, "Delta RPQ score, median [Q1-Q3]" = 3))
