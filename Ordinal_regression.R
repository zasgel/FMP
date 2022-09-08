#https://data.library.virginia.edu/fitting-and-interpreting-a-proportional-odds-model/

options(stringsAsFactors=F)

#package to fit ordinal regression (polr)
library(MASS)

#package to fit ordinal regression with and without random effects (ordinal) with clm2 for fixed effects (equivalent to polr) and clmm2 for random effects
library(ordinal)

###########	
#read data#
###########
df <- read.table("final_covariates_gwas_4_2021_marray_PRS_corrected.csv",sep=",", header = T)
summary(df)

#need to calculate PCs and add them to the regression
df$transformed_SCORE <- (df$SCORE - mean(df$SCORE))/sd(df$SCORE)

df$Escola <- factor(df$Escola) 

list_polr <- list()
list_clmm <- list()
list_re <- list()
list_m_polr <- list()
list_m_clmm <- list()
list_m_re <- list()

#add batch to the model#################################################

for (subject in c("PrimeraLlengua","Angles","Matematiques")){
  
  df_stage <- na.omit(df[,c(subject,"SEX_final","AGE_final","index_SES_final","Escola", "BATCH",
                    "transformed_SCORE","PC1","PC2","PC3","PC4","PC5","PC6","PC7","PC8","PC9","PC10")])
  
  df_stage[[subject]] = factor(df_stage[[subject]], levels = c("D", "C", "B","A"), ordered = TRUE) 
  
  model_polr <- polr(df_stage[[subject]] ~ SEX_final + AGE_final + index_SES_final + Escola + BATCH + 
                       transformed_SCORE + PC1 + PC2 + PC3 + PC4 + PC5 + PC6+ PC7 + PC8 + PC9 + PC10, data=df_stage, Hess = TRUE)
  summary(model_polr)$coef
  polr <- coef(summary(model_polr))[50,]
  list_polr[[subject]] <- polr
  list_m_polr[[subject]] <- model_polr
  
  vect_subj <- df_stage[[subject]]
  model_clmm <-  clmm2(vect_subj ~ SEX_final + AGE_final + index_SES_final + transformed_SCORE+ BATCH + PC1 
                       + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, data=df_stage,Hess = TRUE)
  summary(model_clmm)$coef
  clmm <- coef(summary(model_clmm))[7,]
  list_clmm[[subject]] <- clmm
  list_m_clmm[[subject]] <- model_clmm
  
  #add Escola as random effect
  model <-  clmm2(vect_subj ~ SEX_final + AGE_final + index_SES_final + transformed_SCORE + BATCH + PC1 
                     + PC2 + PC3 + PC4 + PC5 + PC6 + PC7 + PC8 + PC9 + PC10, random=Escola, data=df_stage,Hess = TRUE)
  summary(model)$coef
  random.effect <- coef(summary(model))[7,]
  list_re[[subject]] <- random.effect
  list_m_re[[subject]] <- model
  
}

df_polr <- do.call(rbind,list_polr)
df_clmm <- do.call(rbind,list_clmm)
df_re <- do.call(rbind,list_re)

write.csv(df_re, "Random_effects_results.csv" )


list_m_re$Angles

library(effects)
plot(Effect(focal.predictors = c("transformed_SCORE"),list_m_re$Angles), main = "Effect plot for English")
plot(Effect(focal.predictors = c("transformed_SCORE"),list_m_re$PrimeraLlengua), main = "Effect plot for native tongue")
plot(Effect(focal.predictors = c("transformed_SCORE"),list_m_re$Matematiques), main = "Effect plot for Maths")

library(ggplot2)
ggplot(subset(df, !is.na(df$Angles)), aes(x=transformed_SCORE, color=Angles)) +
  geom_density()

ggplot(subset(df, !is.na(df$PrimeraLlengua)), aes(x=transformed_SCORE, color=PrimeraLlengua)) +
  geom_density()

ggplot(subset(df, !is.na(df$Matematiques)), aes(x=transformed_SCORE, color=Matematiques)) +
  geom_density()

