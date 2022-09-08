#######################################################
# TRS
# 29/03/2022
# Natalia Llonga - editat Judit Cabana
#######################################################

# Description
#----------------------------------------
# Script to perform Transcriptomic Risk Score. 

# Method to use this script in bash:

## Rscript TRS_Judit.R --twas="/home/juditc/ADHD/GWAS_TDAH/TDAH/PRS_all_demontis_GWAS1_6/TRS/TWAS_ditte/ACC_ADHD_ditte_JTI.csv" --marray="/home/juditc/ADHD/GWAS_TDAH/TDAH/GWAS_TDAH_b38/MetaXcan/software/results/JTI/ACC_ADHD_gwas1to6_b37_predict_JTI.txt"  --workingdir="/home/juditc/ADHD/GWAS_TDAH/TDAH/PRS_all_demontis_GWAS1_6/TRS/Results/ACC/p-val1/" --clumplist="/home/juditc/ADHD/GWAS_TDAH/TDAH/PRS_all_demontis_GWAS1_6/TRS/Clumpling_genes/ACC_selected_genes.txt" --covar="/home/juditc/ADHD/GWAS_TDAH/TDAH/GWAS1_6_norel_TDAH_covar_file_new.txt" --output="GWAS1_6_ACC"

#  
# Description of parameters: 

## --twas                : Specifies the file with TWAS results that we are going to use as discovery - TWAS of Demontis - they are in this path (/home/juditc/ADHD/GWAS_TDAH/TDAH/PRS_all_demontis_GWAS1_6/TRS_Demontis22_microarray_new/Results_TWAS)
## --marray              : Specifies the file with predicted expression from FE with genes in rows (genenames in ensembl id) and individuals in columns.
## --threshold           : Specifies the threshold to filter genes in TWAS
## --workingdir          : Working directory (path). The results will appear in the specified folder.
## --output              : Name of the output
## --clumplist			     : List of clumped genes that we want to keep from the marray file.
## --covar				       : File with covariates that we want to include in the logistic regression

## Thresholds:Bonferroni, 0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1

options(warn=-1) #supress warnings in r code
library(pacman)
library(rms)
library(ordinal)
library(MASS)
library(data.table)
library(rcompanion)
p_load(dplyr, GenomicRanges, tidyr, rtracklayer, tidyverse, ggplot2, data.table)

# Variables introduced with commandArgs:
args <- commandArgs(trailingOnly = TRUE)

TWAS.results      <- as.character(strsplit(grep('--twas*', args, value = TRUE), split = '=')[[1]][[2]])
marray.data       <- as.character(strsplit(grep('--marray*', args, value = TRUE), split = '=')[[1]][[2]])
workingdir        <- as.character(strsplit(grep('--workingdir*', args, value = TRUE), split = '=')[[1]][[2]])
output            <- as.character(strsplit(grep('--output*', args, value = TRUE), split = '=')[[1]][[2]])
clumped.genes	  <- as.character(strsplit(grep('--clumplist*', args, value = TRUE), split = '=')[[1]][[2]])
covar			  <- as.character(strsplit(grep('--covar*', args, value = TRUE), split = '=')[[1]][[2]])
subject <- as.character(strsplit(grep('--subject*', args, value = TRUE), split = '=')[[1]][[2]])

#   workingdir ="~/Library/CloudStorage/OneDrive-VallHebronInstitutRecerca-G60594009/TWAS/Predixcan/Resultats/Dades noves gwas ADHD/prova"
#   output ="GWAS1_6_ACC_Demontis22"
#   threshold.p=1
#   TWAS.results="~/Library/CloudStorage/OneDrive-VallHebronInstitutRecerca-G60594009/TWAS/Predixcan/Resultats/Dades noves gwas ADHD/Resultats ditte22/JTI_Brain_ACC_ADHD_ditte22_JTI.csv"
#   marray.data="~/Library/CloudStorage/OneDrive-VallHebronInstitutRecerca-G60594009/TWAS/Predixcan/Resultats/Dades noves gwas ADHD/scripts natalia/marray_f.RData"
##   clumped.genes="/home/juditc/ADHD/GWAS_TDAH/TDAH/PRS_all_demontis_GWAS1_6/TRS/Clumpling_genes/ACC_selected_genes.txt"
#	  covar="~/Library/CloudStorage/OneDrive-VallHebronInstitutRecerca-G60594009/TWAS/Predixcan/Resultats/Dades noves gwas ADHD/covariates_gwas2456_marray_mod.txt"
#

# Establish working directory
#----------------------------------------
print("Starting script.")
print(paste("Script for: ", output, sep=""))
setwd(workingdir)
print(paste("Working directory used: ", workingdir))

#Load data
#----------------------------------------


######################################################
# First part: TWAS results
######################################################
print("First part: TWAS results")
# GTEx:
TWAS.r <- read.csv(TWAS.results, header=T)
#head(TWAS.r)

#Change all strings that are numeric and appear as character 
TWAS.r = TWAS.r %>% mutate(zscore = as.numeric(as.character(zscore)),
                           effect_size = as.numeric(as.character(effect_size)),
                           pvalue = as.numeric(as.character(pvalue)),
                           var_g = as.numeric(as.character(var_g)),
                           pred_perf_r2 = as.numeric(as.character(pred_perf_r2)),
                           pred_perf_pval = as.numeric(as.character(pred_perf_pval)),
                           pred_perf_qval = as.numeric(as.character(pred_perf_qval)))

print(paste("Number of genes used in Demontis TWAS: ",TWAS.r %>% nrow, sep="")) # 7490

TWAS.r = TWAS.r %>% dplyr::filter(!is.na(pvalue))

print(paste("Number of genes with result: ",TWAS.r %>% nrow, sep="")) # 7490

print(paste("Number of significative genes: ",TWAS.r %>% dplyr::filter(pvalue < 0.05) %>% nrow, sep="")) # 764


#Plot a volcano plot to see if we have outliers

TWAS.r$delabel <- NA
pdf(paste("Volcano_plot_TWAS_results_zscore_",output, ".pdf", sep=""))
ggplot(data=TWAS.r, aes(x=zscore, y=-log10(pvalue), label=delabel)) + 
  geom_point() + 
  theme_classic() +
  geom_text()
dev.off()

list_tresh <- list()
bf_treshold <- 0.05/TWAS.r %>% nrow
tresh <- c(bf_treshold ,0.001, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1)

#Read TWAS pred file on FE 
marray.data0<-read.table(marray.data,header=T)
#dim(marray.data0)	#[1] 5146 4918


#Threshold:
for (threshold.p in tresh) {
  TWAS.r2 = TWAS.r %>% filter(pvalue < threshold.p) 
  

print(paste("Number of genes after threshold of ", threshold.p, " : ",TWAS.r2 %>% nrow, sep="")) # 7490



######################################################
# Second part: microarray data
######################################################
print("Second part: microarray data")


print(paste("Number of genes in target TWAS: ",marray.data0 %>% ncol, sep="")) # 4918

#Read the output from the clumping genes
col<-read.table(clumped.genes, header=T)
col<-col[,1]
col<-as.character(col)

#Filter the TWAS file to keep only those that are in the clumping file
marray.data1<-marray.data0[, names(marray.data0) %in% col]
#dim(marray.data1)	#[1] 5146 4633

print(paste("Number of non correlated genes in target TWAS: ",marray.data1 %>% ncol, sep="")) # 4633

#Transpose the data to have genes in rows and individuals in the columns. 
rownames(marray.data1)<-marray.data0$IID
marray.data<-as.data.frame(t(marray.data1 [,-1]))

####### Standarize TWAS target expression
######################################################
print("Header of initial microarray")
marray.data[1:5,1:5]
dim(marray.data)

standard.deviation = apply(marray.data, 1 , sd)
marray.st = (marray.data-rowMeans(marray.data))/standard.deviation
print("Header of standarized microarray data")
marray.st[1:5,1:5]
dim(marray.st)

# marray.st is the standarized microarray data

#Check
first.row = (marray.data[1,]-rowMeans(marray.data)[1])/standard.deviation[1] 
#head(first.row)
print("Check, if it's all True, is correct: ")
table(first.row == marray.st[1,]) # 556 All true

######################################################
# Third part: Calculate TRS and standarize TRS
######################################################

# En la tabla de marray.st crear una columna TRS, en la que se calcule el 
#  zscore*geneexpression por cada gen en cada individuo
#marray.st[1:5,1:5]
print("Number of genes in common between TWAS target and TWAS discovery (column TRUE): ")
table(rownames(marray.st) %in% TWAS.r2$gene) # 7476 genes in common

#Plot a VennDiagram with the genes in common
library(VennDiagram)
venn.diagram(list(marray = rownames(marray.st), GTEx8_Ditte = TWAS.r2$gene),
             paste("Venn_Genes_in_common_threshold_",threshold.p,".tiff", sep=""),
             fill=c("red", "green"), 
             alpha=c(0.5,0.5), cex = 2, 
             category.names=c("marray", "GTEx8_Ditte22"))

#Now we join the two files to obtain the same order
alldata = marray.st %>% as.data.frame
alldata$gene = rownames(alldata)
alldata2 = inner_join(alldata, TWAS.r2 %>% dplyr::select(gene, zscore, pvalue))
#Joining, by = "gene"

write.table(alldata2$gene, paste("genename_genes_common_targetTWAS_discoveryTWAS_threshold_", threshold.p,".txt", sep=""), quote=F, col.names=F, row.names=F)


#And separate them again to run TRS
marray.sorted = alldata2[,!colnames(alldata2) %in% c("gene", "zscore", "pvalue")]
#marray.sorted[1:5,1:5]
twas.sorted = alldata2[,c("gene", "zscore", "pvalue")]

TRS.zscore = colSums(marray.sorted*twas.sorted$zscore)

results = data.frame(ID = names(TRS.zscore), TRS.zscore = TRS.zscore)
results$ID<-gsub("0_","",results$ID)

# Now we standarize the TRS:

results$TRS_zscore_sd = (results$TRS.zscore-mean(results$TRS.zscore))/sd(results$TRS.zscore)

write.table(results, paste("Output_TRS_", output,"_threshold_", threshold.p,".txt", sep=""), quote=F, col.names=T, row.names=F)




#############################################################################
### THIS PART WILL NEED TO BE ADAPTED BECAUSE WE USE ORDINAL REGRESSION, SO MAYBE WE CAN STOP DE SCRIPT HERE AND RUN THE ORDINAL REGRESSION USING THE OTHER SCRIPT THAT MARIA HAS PREPARED OR COPY IT HERE.
######################################################
# Logistic regression analysis
######################################################


df <- read.table(covar,sep=",", header = T)
df$transformed_PRS_SCORE <- (df$SCORE - mean(df$SCORE))/sd(df$SCORE)
df<- merge(df,results, by.x = "subject_ID", by.y = "ID")

df$Escola <- factor(df$Escola) 

df_stage <- na.omit(df[,c(subject,"SEX_final","AGE_final","index_SES_final","Escola", "BATCH", "TRS_zscore_sd",
                          "transformed_PRS_SCORE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])

df_stage[[subject]] = factor(df_stage[[subject]], levels = c("D", "C", "B","A"), ordered = TRUE) 

vect_subj <- df_stage[[subject]]

#add Escola as random effect
model_TRS <-  clmm2(vect_subj ~ TRS_zscore_sd + SEX_final + AGE_final + index_SES_final +  BATCH + PC1 
                + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)
res_TRS <- summary(model_TRS)$coef
random.effect_TRS <- t(as.data.frame(coef(summary(model_TRS))[4,]))
colnames(random.effect_TRS) <- c("Estimate", "Std_error", "z_value", "Pvalue")
rownames(random.effect_TRS) <- output
as.data.frame(random.effect_TRS)
list_tresh[[as.character(threshold.p)]] <- c(random.effect_TRS,threshold.p)

model_r2 <-  clm(vect_subj ~ TRS_zscore_sd + SEX_final + AGE_final + index_SES_final + Escola + BATCH + PC1 
                 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)

model_null <-  clm(vect_subj ~ SEX_final + AGE_final + index_SES_final + BATCH + Escola + PC1 
                   + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)

r2<-nagelkerke(fit = model_r2,null = model_null)
r2$Pseudo.R.squared.for.model.vs.null[3]

list_tresh[[as.character(threshold.p)]] <- c(random.effect_TRS,threshold.p,r2$Pseudo.R.squared.for.model.vs.null[3])
}
thresold_df <- as.data.frame(do.call(rbind,list_tresh))
colnames(thresold_df) <- c("Estimate", "Std_error", "z_value", "Pvalue","Threshold","R2")
write.table(thresold_df,paste("Tresholds_",output,".txt",sep = ""),quote = F)
thresold_df$Threshold<-as.character(thresold_df$Threshold)

pdf(paste("Barplot_thresholds_",output, ".pdf", sep=""))
ggplot(thresold_df, aes(x=Threshold, y=R2)) + 
  geom_bar(stat="identity", fill="steelblue")  + 
  theme_classic() + 
  geom_text(aes(label=round(Pvalue, digits = 4)), position=position_dodge(width=0.5), vjust=-0.15, ) +
  ggtitle(paste("R2_",output,sep = "")) +
  scale_x_discrete(limits=thresold_df$Threshold)
dev.off()
min(thresold_df$Pvalue)
selected_thresh<-thresold_df$Threshold[thresold_df$Pvalue==min(thresold_df$Pvalue)]

results <- read.table(paste("Output_TRS_", output,"_threshold_", selected_thresh,".txt", sep=""),header = T)


df <- read.table(covar,sep=",", header = T)
df$transformed_PRS_SCORE <- (df$SCORE - mean(df$SCORE))/sd(df$SCORE)
df<- merge(df,results, by.x = "subject_ID", by.y = "ID")

df$Escola <- factor(df$Escola) 

df_stage <- na.omit(df[,c(subject,"SEX_final","AGE_final","index_SES_final","Escola", "BATCH", "TRS_zscore_sd",
                          "transformed_PRS_SCORE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])

df_stage[[subject]] = factor(df_stage[[subject]], levels = c("D", "C", "B","A"), ordered = TRUE) 

vect_subj <- df_stage[[subject]]

model_TRS <-  clmm2(vect_subj ~ TRS_zscore_sd + SEX_final + AGE_final + index_SES_final +  BATCH + PC1 
                    + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)
res_TRS <- summary(model_TRS)$coef
random.effect_TRS <- t(as.data.frame(coef(summary(model_TRS))[4,]))
colnames(random.effect_TRS) <- c("Estimate", "Std_error", "z_value", "Pvalue")
rownames(random.effect_TRS) <- output

model_all <-  clmm2(vect_subj ~ TRS_zscore_sd + transformed_PRS_SCORE + SEX_final + AGE_final + index_SES_final +  BATCH + PC1 
                + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)
res_all <- summary(model_all)$coef
random.effect_all <- as.data.frame(coef(summary(model_all))[4:5,])
rownames(random.effect_all)[1] <- output
rownames(random.effect_all)[2] <- paste(output,"_PRS_score",sep="")


model_re <-  clmm2(vect_subj ~ SEX_final + AGE_final + index_SES_final + transformed_PRS_SCORE + BATCH + PC1 
                + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)


model_PRS <-  clm(vect_subj ~ SEX_final + AGE_final + index_SES_final + transformed_PRS_SCORE + BATCH + Escola + PC1 
                   + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)

model_r2 <-  clm(vect_subj ~ TRS_zscore_sd + SEX_final + AGE_final + index_SES_final + Escola + BATCH + PC1 
                 + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)

model_null <-  clm(vect_subj ~ SEX_final + AGE_final + index_SES_final + BATCH + Escola + PC1 
                   + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)

model_PRS_re <-  clmm2(vect_subj ~ SEX_final + AGE_final + index_SES_final + transformed_PRS_SCORE + BATCH + PC1 
                + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)
res_PRS_re <- summary(model_PRS_re)$coef
random.effect_PRS_re <- t(as.data.frame(coef(summary(model_PRS_re))[7,]))
rownames(random.effect_PRS_re) <- output


r2<-nagelkerke(fit = model_r2,null = model_null)
r2$Pseudo.R.squared.for.model.vs.null[3]

r2_PRS<-nagelkerke(fit = model_PRS,null = model_null)
r2_PRS$Pseudo.R.squared.for.model.vs.null[3]

model_comp <-  clm(vect_subj ~ TRS_zscore_sd + transformed_PRS_SCORE + SEX_final + AGE_final + index_SES_final +  BATCH + PC1 
                    + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10 + Escola, data=df_stage,Hess = TRUE)


r2<-nagelkerke(fit = model_r2,null = model_null)
r2$Pseudo.R.squared.for.model.vs.null[3]

r2_comp<-nagelkerke(fit = model_comp,null = model_null)
r2_comp$Pseudo.R.squared.for.model.vs.null[3]


model_fit2<-anova(model_re,model_all)
model_fit2$`Pr(Chi)`[2]

random.effect_TRS<-as.character(random.effect_TRS)
random.effect_all_1 <- as.character(random.effect_all[1,])
random.effect_all_2 <- as.character(random.effect_all[2,])
random.effect_PRS_re<- as.character(random.effect_PRS_re)
results_TRS1 <- c(selected_thresh,random.effect_TRS,r2$Pseudo.R.squared.for.model.vs.null[3],random.effect_all_1,model_fit2$`Pr(Chi)`[2],r2_comp$Pseudo.R.squared.for.model.vs.null[3])
results_TRS2 <- c("NA",random.effect_PRS_re,r2_PRS$Pseudo.R.squared.for.model.vs.null[3],random.effect_all_2,"NA","NA")
final_table <- as.data.frame(rbind(results_TRS1,results_TRS2))
colnames(final_table) <- c("Threshold","Estimate", "Std_error", "z_value", "Pvalue","R2", "Estimate_PRS", "Std_error_PRS", "z_value_PRS", "Pvalue_PRS", "LR_Pvalue", "R2")
rownames(final_table) <- c(output,paste(output,"_PRS",sep = ""))

write.table(final_table, paste("final_results_",output,".txt",sep=""), quote = F)
library(ggplot2)
library(effects)
effect_plot <- plot(Effect(focal.predictors = c("TRS_zscore_sd"),model_TRS), main = paste("Effect plot for",output))

pdf(paste("effect_plot_TWAS_results_",output, ".pdf", sep=""))
print(effect_plot)
dev.off()

TRS_plot <- ggplot(subset(df, !is.na(df[[subject]])), aes(x=TRS_zscore_sd, color=vect_subj)) +
  geom_density() +ggtitle(paste("Plot_TWAS_result_",output))

pdf(paste("plot_TWAS_results_",output, ".pdf", sep=""))
print(TRS_plot)
dev.off()



print("Finished :) ")

sink() #Close conection to file

