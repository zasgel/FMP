#Script to plot all R2 barplots:

#First we need to open the association results per p-value threshold in each tissue. This file is the output of the script to create the barplots for all thresholds in each tissue

df1 <- read.table("Tresholds_PrimeraLlengua_Whole_Blood_TRS.txt")
df1$Tissue <- "Blood"
df2 <- read.table("Tresholds_PrimeraLlengua_Brain_Amygdala_TRS.txt")
df2$Tissue <- "Amygdala"
df3 <- read.table("Tresholds_PrimeraLlengua_Brain_Anterior_cingulate_cortex_BA24_TRS.txt")
df3$Tissue <- "Anterior_cingulate_cortex_BA24"
df4 <- read.table("Tresholds_PrimeraLlengua_Brain_Caudate_basal_ganglia_TRS.txt")
df4$Tissue <- "Caudate_basal_ganglia"
df5 <- read.table("Tresholds_PrimeraLlengua_Cerebellar_Hemisphere_TRS.txt")
df5$Tissue <- "Cerebellar_Hemisphere"
df6 <- read.table("Tresholds_PrimeraLlengua_Brain_Cerebellum_TRS.txt")
df6$Tissue <- "Cerebellum"
df7 <- read.table("Tresholds_PrimeraLlengua_Brain_Cortex_TRS.txt")
df7$Tissue <- "Cortex"
df8 <- read.table("Tresholds_PrimeraLlengua_Brain_Frontal_Cortex_BA9_TRS.txt")
df8$Tissue <- "Frontal_Cortex_BA9"
df9 <- read.table("Tresholds_PrimeraLlengua_Brain_Hippocampus_TRS.txt")
df9$Tissue <- "Hippocampus"
df10 <- read.table("Tresholds_PrimeraLlengua_Brain_Hypothalamus_TRS.txt")
df10$Tissue <- "Hypothalamus"
df11 <- read.table("Tresholds_PrimeraLlengua_Brain_Nucleus_accumbens_basal_ganglia_TRS.txt")
df11$Tissue <- "Nucleus_accumbens_basal_ganglia"
df12 <- read.table("Tresholds_PrimeraLlengua_Brain_Putamen_basal_ganglia_TRS.txt")
df12$Tissue <- "Putamen_basal_ganglia"
df13 <- read.table("Tresholds_PrimeraLlengua_Spinal_cord_TRS.txt")
df13$Tissue <- "Spinal_cord_cervical_c-1"
df14 <- read.table("Tresholds_PrimeraLlengua_Brain_Substantia_nigra_TRS.txt")
df14$Tissue <- "Substantia_nigra"

#We merge them
df_all <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14)

Thresholds <- table(df_all$Threshold)

#Change the name of the Bonferroni p-value threshold 
df_all$Threshold2 <- df_all$Threshold
df_all$Threshold2[!df_all$Threshold2 %in% c("1", "0.5","0.4", "0.3", "0.2", "0.1", "0.05", "0.001")] <- "Bonf"


df_all$Threshold2 <- as.factor(df_all$Threshold2)
df_all$P_value_TRS_z <- as.numeric(df_all$P_value_TRS_z)

df_all$Threshold2 <-factor(df_all$Threshold2, levels = levels(df_all$Threshold2)[c(9,1,2,3,4,5,6,7,8)])

#We added this step to diferenciate the results from brain and blood (different colour)
df_all$category<-c()
df_all$category[df_all$Tissue == "Blood"]<-"coral"
df_all$category[df_all$Tissue != "Blood"]<-"steelblue"

#This step is necessary if you want to order the tissues in a specific way, for example blood at first position and then all the brain tissues.
positions<-c("Blood","Amygdala", "Anterior_cingulate_cortex_BA24", "Caudate_basal_ganglia","Cerebellar_Hemisphere","Cerebellum","Cortex","Frontal_Cortex_BA9","Hippocampus","Hypothalamus","Nucleus_accumbens_basal_ganglia","Putamen_basal_ganglia","Spinal_cord_cervical_c-1","Substantia_nigra")
library(plyr)  ## or dplyr (transform -> mutate)
df_all2 <- arrange(transform(df_all,Tissue=factor(Tissue,levels=positions)),Tissue)

#Plot de results
ggplot(df_all2, aes(x=Threshold2, y=R2)) + 
  geom_bar(stat="identity", fill=df_all2$category)  + 
  theme_classic() + theme(text = element_text(size = 10))+
  geom_text(aes(label=P_value_TRS_z), position=position_dodge(width=2), vjust=-0.25,hjust=-0.15, size=1.8, angle=55) +
  ggtitle("R2 TRS ADHD") +
  facet_wrap(~Tissue, nrow = 2, dir="h", scales="free_x") +
  theme(plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"), axis.text.x = element_text(angle = 45,hjust=1)) +
  ylim(0, 0.035) 

#ggplot(df_all, aes(x=Threshold2, y=R2, color=Tissue)) + geom_point() + geom_line(size=1, aes(group = Tissue)) + theme_classic()

################################################################################################################################

#Script to plot best p-value results for all tissues in R2 barplots:
#This step is the same than in the previous barplot:
df1 <- read.table("Tresholds_PrimeraLlengua_Whole_Blood_TRS.txt")
df1$Tissue <- "Blood"
df2 <- read.table("Tresholds_PrimeraLlengua_Brain_Amygdala_TRS.txt")
df2$Tissue <- "Amygdala"
df3 <- read.table("Tresholds_PrimeraLlengua_Brain_Anterior_cingulate_cortex_BA24_TRS.txt")
df3$Tissue <- "Anterior_cingulate_cortex_BA24"
df4 <- read.table("Tresholds_PrimeraLlengua_Brain_Caudate_basal_ganglia_TRS.txt")
df4$Tissue <- "Caudate_basal_ganglia"
df5 <- read.table("Tresholds_PrimeraLlengua_Cerebellar_Hemisphere_TRS.txt")
df5$Tissue <- "Cerebellar_Hemisphere"
df6 <- read.table("Tresholds_PrimeraLlengua_Brain_Cerebellum_TRS.txt")
df6$Tissue <- "Cerebellum"
df7 <- read.table("Tresholds_PrimeraLlengua_Brain_Cortex_TRS.txt")
df7$Tissue <- "Cortex"
df8 <- read.table("Tresholds_PrimeraLlengua_Brain_Frontal_Cortex_BA9_TRS.txt")
df8$Tissue <- "Frontal_Cortex_BA9"
df9 <- read.table("Tresholds_PrimeraLlengua_Brain_Hippocampus_TRS.txt")
df9$Tissue <- "Hippocampus"
df10 <- read.table("Tresholds_PrimeraLlengua_Brain_Hypothalamus_TRS.txt")
df10$Tissue <- "Hypothalamus"
df11 <- read.table("Tresholds_PrimeraLlengua_Brain_Nucleus_accumbens_basal_ganglia_TRS.txt")
df11$Tissue <- "Nucleus_accumbens_basal_ganglia"
df12 <- read.table("Tresholds_PrimeraLlengua_Brain_Putamen_basal_ganglia_TRS.txt")
df12$Tissue <- "Putamen_basal_ganglia"
df13 <- read.table("Tresholds_PrimeraLlengua_Spinal_cord_TRS.txt")
df13$Tissue <- "Spinal_cord_cervical_c-1"
df14 <- read.table("Tresholds_PrimeraLlengua_Brain_Substantia_nigra_TRS.txt")
df14$Tissue <- "Substantia_nigra"


df_all <- rbind(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14)

Thresholds <- table(df_all$Threshold)

df_all$Threshold2 <- df_all$Threshold
df_all$Threshold2[!df_all$Threshold2 %in% c("1", "0.5","0.4", "0.3", "0.2", "0.1", "0.05", "0.001")] <- "Bonf"


df_all$Threshold2 <- as.factor(df_all$Threshold2)
df_all$P_value_TRS_z <- as.numeric(df_all$P_value_TRS_z)

df_all$Threshold2 <-factor(df_all$Threshold2, levels = levels(df_all$Threshold2)[c(9,1,2,3,4,5,6,7,8)])

#Here we select the best p-value thresholds per each tissue
library(pacman)
p_load(openxlsx, ggplot2, dplyr)
df_best_pval <- df_all %>% mutate(P_value_TRS_z = as.numeric(P_value_TRS_z)) %>%  group_by(Tissue) %>% summarize(min(P_value_TRS_z))
df_all3<-df_all %>% dplyr::filter(P_value_TRS_z %in% df_best_pval$`min(P_value_TRS_z)`)

#Then we create a new column that has the name of the tissue and the p-value threshold used 
library(stringr)
df_all3$Tissue2<-str_c(df_all3$Tissue," (P=", df_all3$Threshold2,")")

#Plot the results
library(ggplot2)
ggplot(df_all3, aes(x=Tissue2, y=R2)) + 
  geom_bar(stat="identity", fill="steelblue")  + 
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label=P_value_TRS_z), position=position_dodge(width=0.9), vjust=-0.25) +
  ggtitle("R2 TRS") +
  scale_x_discrete(limits=df_all3$Tissue2)


##If you prefer to create a plot with * instead of the exact p-value of each bar (it would be very useful to identify the significant results) you need to do it manually.
#With natalia we have tried to do it in ggplot but we won't be able. So you need to add an additional column that would be: * p-value < 0.05, ** p-value < 0.01, ***p-value<0.001
setwd("~/Library/CloudStorage/OneDrive-VallHebronInstitutRecerca-G60594009/TWAS/Predixcan/Resultats/TRS i PRS TWAS brain vs microarray/TRS TWAS ditte 22 new v2.0/TRS_v6/Results_v6/2_barplots")
write.table(df_all3, file="Best_p-value_all_tissues0.txt",quote=F, col.names = T, row.names = F)
df_all3<-read.delim("Best_p-value_all_tissues.txt",h=T)

ggplot(df_all3, aes(x=Tissue2, y=R2)) + 
  geom_bar(stat="identity", fill="steelblue")  + 
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  geom_text(aes(label=sig), position=position_dodge(width=0.9), vjust=0.5, size= 6) +
  ggtitle("R2 TRS") +
  scale_x_discrete(limits=df_all3$Tissue2)

################################################################################################################
######## TRS + PRS ########
################################################################################################################

#We have prepared a file with the required information in excel: 
#R2 --> result of R2 
#TRS_Tissue --> tissue name and the best p-value threshold in brackets
#Class --> TRS, PRS or TRS+PRS

df_all<-read.delim("TRS_best-pvalue_sig_tissues_PRS.txt",h=T)

#Fem el plot

library(ggplot2)

ggplot(df_all, aes(fill=Class, x=TRS_Tissue, y=R2)) + 
  geom_bar(stat="identity", position='dodge')+  
  theme_light() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  ggtitle("R2 TRS") +
  scale_x_discrete(limits=df_all$TRS_Tissue) +
  scale_fill_manual('Results', values=c('coral2','steelblue','darkolivegreen3'))

