
#Set options
options(warn=-1)
#Load libraries
library(ggplot2)
library(data.table)

tissues <- c("Brain_Amygdala_FE", "Brain_Cerebellum_FE", "Brain_Cortex","Brain_Cerebellar_Hemisphere_FE",
             "Brain_Caudate_basal_ganglia_FE", "Brain_Frontal_Cortex_BA9_FE","Brain_Putamen_basal_ganglia_FE",
             "Brain_Anterior_cingulate_cortex_BA24_FE","Brain_Hypothalamus_FE","Brain_Hippocampus_FE", 
             "Spinal_cord_FE", "Brain_Substantia_nigra_FE", "Brain_Nucleus_accumbens_basal_ganglia_FE", 
             "Whole_Blood_FE")
subjects <- c("Matematiques", "Angles", "PrimeraLlengua")

list_var <- list()
for (tissue in tissues) {
  pred <- paste(tissue,"_predict_JTI.txt", sep = "")
  pred <- fread(pred, header=T)
  for (subject in subjects) {
    
    df <- paste(subject,"_",tissue,"_predict_JTI.txt", sep = "")
    output_name <- paste(subject, "_", tissue, sep = "")
    #Read tables
    df <- read.table(df)
    pred <- as.data.frame(pred)
    
    initial_genes <- dim(df)[1]
    df$Genes <- rownames(df)
    var.pred<-lapply(pred[,3:dim(pred)[2]], var)
    var <- as.data.frame(var.pred)
    var <- as.data.frame(t(var))

    colnames(var) <- "var"
    var$Genes <- rownames(var)
    
    var_remove <- var[var$var<5*10^-5,]
    vec_genes <- var_remove$Genes
    pred_new <- pred[,!colnames(pred) %in% vec_genes]
    
    df5 <- merge(df,var,by.x ="Genes", by.y = "Genes")
    df5 <- df5[df5$var>5*10^-5,]
  
    
    filtered_genes <- dim(df5)[1]
write.table(df5,paste(subject,"_",tissue,"_predict_JTI_association_filtered.txt", sep = ""))
write.table(pred_new, paste(tissue,"_predict_JTI_filtered.txt", sep = "") )
list_var[[paste(tissue,"_",subject)]] <- data.frame(initial_genes,filtered_genes)

a <- ggplot(data=df5, aes(x=Estimate, y=-log10(Pr...z..))) + geom_point() + ggtitle(paste("Volcano plot for",subject,
                                                                                          tissue, sep = ""))
                                                                                     
pdf(paste("volcano_",output_name, "_filtered.pdf" ,sep=""))
print(a)
dev.off()

  }
}
df_var <- do.call(rbind,list_var)
write.table(df_var, paste("Gene_numbers_filtering.txt", sep=""))
