
#Set options
options(warn=-1)
#Load libraries
library(ggplot2)
library(data.table)

list_res <- list()
tissues <- c("Brain_Amygdala_FE", "Brain_Cerebellum_FE", "Brain_Cortex","Brain_Cerebellar_Hemisphere_FE",
             "Brain_Caudate_basal_ganglia_FE", "Brain_Frontal_Cortex_BA9_FE","Brain_Putamen_basal_ganglia_FE",
             "Brain_Anterior_cingulate_cortex_BA24_FE","Brain_Hypothalamus_FE","Brain_Hippocampus_FE", 
             "Spinal_cord_FE", "Brain_Substantia_nigra_FE", "Brain_Nucleus_accumbens_basal_ganglia_FE", 
             "Whole_Blood_FE")
subjects <- c("Matematiques", "Angles", "PrimeraLlengua")

for (tissue in tissues) {
  pred <- paste(tissue,"_predict_JTI.txt", sep = "")
  pred <- fread(pred, header=T)
  for (subject in subjects) {
    

df <- paste(subject,"_",tissue,"_predict_JTI.txt", sep = "")
output_name <- paste(subject, "_", tissue, sep = "")

#Read tables
df <- read.table(df)
pred <- as.data.frame(pred)

df$Genes <- rownames(df)
df$log.z <- -log10(df$Pr...z..)
df4 <- df[df$log.z>10, ]
var.pred<-lapply(pred[,3:dim(pred)[2]], var)
iqr.pred <- lapply(pred[,3:dim(pred)[2]], IQR)
var <- as.data.frame(var.pred)
iqr <- as.data.frame(iqr.pred)
var <- as.data.frame(t(var))
iqr <- as.data.frame(t(iqr))


colnames(iqr) <- "IQR"
colnames(var) <- "var"
iqr$Genes <- rownames(iqr)
var$Genes <- rownames(var)
outliers <- merge(df4,var,by.x ="Genes", by.y = "Genes")
outliers <- merge(outliers,iqr,by.x ="Genes", by.y = "Genes")

pdf(paste("histogram_",output_name, ".pdf" ,sep=""))
hist(outliers$var)
hist(outliers$IQR)
dev.off()

df5 <- merge(df,var,by.x ="Genes", by.y = "Genes")
df6 <- merge(df,iqr,by.x ="Genes", by.y = "Genes")
df7 <- merge(df,var,by.x ="Genes", by.y = "Genes")
df8 <- merge(df,iqr,by.x ="Genes", by.y = "Genes")

pdf(paste("histogram_",output_name, ".pdf" ,sep=""))
hist(df5$var)
hist(df6$IQR)
hist(outliers$var)
hist(outliers$IQR)
dev.off()

df5 <- df5[df5$var>max(outliers$var),]
df6 <- df6[df6$IQR>max(outliers$IQR),]

pdf(paste("volcano_",output_name, ".pdf" ,sep=""))
ggplot(data=df, aes(x=Estimate, y=-log10(Pr...z..))) + geom_point() + ggtitle("Volcano plot for inital dataset")
ggplot(data=df5, aes(x=Estimate, y=-log10(Pr...z..))) + geom_point() + ggtitle("Volcano plot for removed outliers with var")
ggplot(data=df6, aes(x=Estimate, y=-log10(Pr...z..))) + geom_point() + ggtitle("Volcano plot for removed outliers with IQR")
dev.off()

outliers_number <- dim(outliers)[1]
outliers_number_var <- table(df7$var<max(outliers$var))[2]
outliers_number_iqr <- table(df8$IQR<max(outliers$IQR))[2]
outliers_max_var <- max(outliers$var)
outliers_max_iqr <- max(outliers$IQR)
total_genes <- dim(df)[1]

list_res[[paste(tissue,"_",subject)]] <- data.frame(total_genes,outliers_number,outliers_number_iqr,outliers_number_var,outliers_max_iqr,outliers_max_var)

  }
}
df_res <- do.call(rbind,list_res)
write.table(df_res, paste("results.txt", sep=""))


