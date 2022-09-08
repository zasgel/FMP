#####################
# Script to find clumping correlated genes
#####################

# Description
#----------------------------------------
# With this script we calculated the correlation between genes that are in the same region  in the chromosome (+-kb) based on most significative pvalues order.
# We introduce TWAS results (e.g. JTI method) and a reference file in gtf format. 
# For the gtf file we use gencode.v32lift37.annotation.gtf, downloaded in http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/ in 05/10/2021.

# Method to use this script in bash:

## Rscript clump_genes.r --cov_association_JTI="/Users/nataliallongatejedor/Documents/VHIR/Analisis/TRS/Dades_Judit_TWAS/Substantia_nigra_ADHD_gwas1to6_b37_cov_association_JTI.txt"  --predict_JTI="/Users/nataliallongatejedor/Documents/VHIR/Analisis/TRS/Dades_Judit_TWAS/Substantia_nigra_ADHD_gwas1to6_b37_predict_JTI.txt" --pvals="/Users/nataliallongatejedor/Desktop/Substantia_nigra_ADHD_ditte_JTI.csv" --gtf="/Users/nataliallongatejedor/Desktop/gencode.v32lift37.annotation.gtf" --cor_threshold=0.9 --kb=500 --workingdir="/Users/nataliallongatejedor/Desktop/" --output="Substantia_nigra"

# Description of parameters: 

## --cov_association_JTI : Specifies the file that includes d'association test from TWAS, ended with cov_association_JTI (in format .txt)  e.g. Substantia_nigra_ADHD_gwas1to6_b37_cov_association_JTI.txt
## --predict_JTI         : Specifies the output file that include predicted gene expression from TWAS, ended with predict_JTI.  e.g. Substantia_nigra_ADHD_gwas1to6_b37_predict_JTI.txt
## --pvals               : Specifies the file which contains pvals column to sort and find correlated genes (in format .csv). We used the output file from TWAS used as discovery phenotype for TRS. e.g. Substantia_nigra_ADHD_ditte_JTI.csv
## --gtf                 : Specifies the gtf file which contains genomic coordinates of all genes.
## --cor_threshold       : Correlation to perform clumping. We can use 0.9 
## --kb                  : Physical distance threshold for clumping.
## --workingdir          : Working directory (path). The results will appear in the specified folder.
## --output              : Name of the output



#This script generates 2 output files:
# paste(workingdir, "df_correlations.txt", sep="")                    : Generates a txt with gene_clumping (gene with window +-kb), genes_ov (genes that overlap with gene_clumping) and correlation.
# paste(workingdir, output, "_df_correlations_collapsed.txt", sep="") : Generates a txt with gene_clumping like "df_correlations" but in collapsed table format.


# Starting script
#----------------------------------------
options(warn=-1) #supress warnings in r code
#if is not installed: install.packages("pacman")
library(pacman)
p_load(dplyr, GenomicRanges, tidyr, rtracklayer)
start.time <- Sys.time() # control time in script

#Variables:
#  results.twas <- "/Users/nataliallongatejedor/Desktop/Zeynep_clumping/Angles_Brain_Cerebellar_Hemisphere_FE_predict_JTI_association_filtered.txt"
#  pvals <- "/Users/nataliallongatejedor/Desktop/Zeynep_clumping/JTI_Brain_Cerebellar_Hemisphere_ADHD_ditte22_JTI_OR.csv"
#  gtf_file <- "/Users/nataliallongatejedor/Desktop/Zeynep_clumping/gencode.v32lift37.annotation.gtf"
#  cors.vals <- "/Users/nataliallongatejedor/Desktop/Zeynep_clumping/Brain_Cerebellar_Hemisphere_FE_predict_JTI_filtered.txt"
#  workingdir <- "/Users/nataliallongatejedor/Desktop/Zeynep_clumping/"
#  cor_threshold = 0.9
#  kb = 500
#  output = "test"

# Variables introduced with commandArgs:
args <- commandArgs(trailingOnly = TRUE)

results.twas      <- as.character(strsplit(grep('--cov_association_JTI*', args, value = TRUE), split = '=')[[1]][[2]])
cors.vals         <- as.character(strsplit(grep('--predict_JTI*', args, value = TRUE), split = '=')[[1]][[2]])
pvals             <- as.character(strsplit(grep('--pvals*', args, value = TRUE), split = '=')[[1]][[2]])
gtf_file          <- as.character(strsplit(grep('--gtf*', args, value = TRUE), split = '=')[[1]][[2]])
cor_threshold     <- as.numeric(strsplit(grep('--cor_threshold*', args, value = TRUE), split = '=')[[1]][[2]])
kb                <- as.numeric(strsplit(grep('--kb*', args, value = TRUE), split = '=')[[1]][[2]])
workingdir        <- as.character(strsplit(grep('--workingdir*', args, value = TRUE), split = '=')[[1]][[2]])
output            <- as.character(strsplit(grep('--output*', args, value = TRUE), split = '=')[[1]][[2]])

# Establish working directory
#----------------------------------------
sink(paste(workingdir, output,".log", sep="")) # Create empty txt file
print("Starting script.")
print(paste("Script for: ", output, sep=""))
setwd(workingdir)
print(paste("Working directory used: ", workingdir))

# Load functions
#----------------------------------------
#Function to calculate time spended running the script:
f <- function(start_time, end_time) {
  start_time <- as.POSIXct(start_time)
  end_time   <- as.POSIXct(end_time)
  dt <- difftime(end_time, start_time, units="secs")
  # Since you only want the H:M:S, we can ignore the date...
  # but you have to be careful about time-zone issues
  format(.POSIXct(dt,tz="GMT"), "%H:%M:%S")
}

#Load data
#----------------------------------------
genes.twas <- read.table(results.twas, header=T)[,c("Genes", "Pr...z..")]
colnames(genes.twas) <- c("gene", "pvalue")
#genes.twas %>% dim
#genes.twas %>% head
colnames(genes.twas)[colnames(genes.twas) == "pvalue"] <- "pvalue_gwas1to6"
print(paste("Number of genes in Target TWAS (gwas1to6):  ", nrow(genes.twas), sep=""))

# Open file with pvals
#----------------------------------------
pvals.data <- read.csv(pvals, header=T)[,c("gene", "gene_name", "pvalue")] #pvalue is from Ditte
list.genes <- inner_join(genes.twas, pvals.data, by="gene")
#dim(list.genes) # Tenim 4081 gens en comú
#list.genes %>% head
print(paste("Number of genes that match with Ditte TWAS: ", nrow(list.genes), sep=""))

#Obtain coordinates of the genes
#----------------------------------------
# File final: Amb aquest fan match tots
# http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/GRCh37_mapping/
# gencode.v32lift37.annotation.gtf.gz
gtf <- rtracklayer::import(gtf_file)
gtf_df=as.data.frame(gtf)[,c("seqnames", "start", "end", "width", "strand", "type", "gene_id", "gene_name")]
locs.gtf <- gtf_df %>% filter(type=="gene") %>% tidyr::separate(gene_id, into=c('ensg', 'extension'),
                                                                sep="\\.")
#Merge all information we have
all.info <- inner_join(list.genes, locs.gtf, by = c("gene" = "ensg", "gene_name"))
print("Data.frame with merged info:")
all.info %>% head
# dim(list.genes) # list initial of genes from TWAS 
# dim(locs.gtf) # number of coordinates in gtf file
# nrow(all.info) #allinfo contains info about gtf file.  if is the same as list.genes is correct.
print(paste("Number of genes including genomic location: ", nrow(all.info), sep=""))
# File of correlations
#----------------------------------------
cors <- read.table(cors.vals, header = T)
#cors[1:5,1:5]

# Prepare data for the loop
#----------------------------------------
all.info = all.info %>% arrange(pvalue) #sort genes by pvalue
all.info$gene = as.character(all.info$gene) #as.character gene
all.info2 <- all.info # we create a copy list of genes. In the loop we remove genes from all.info2

list.final.genes <- c() # list of final genes
df_correlations <- data.frame(gene_clumping=c(),
                              genes_ov = c(),
                              correlation = c())
  
#----------------------------------------
# MHC
#----------------------------------------
print("Finding the most significative gene in MHC region")
# https://www.ncbi.nlm.nih.gov/grc/human/regions/MHC?asm=GRCh37
#    chr6:28,477,797-33,448,354
# Maria : 26.000.000-33.000.000
# Pain : 28.000.000-34.000.000

#MHC = data.frame(seqnames = "chr6", start= 28477797, end=33448354, strand="*") %>% GRanges
#MHC = data.frame(seqnames = "chr6", start= 26000000, end=33000000, strand="*") %>% GRanges

# We decided to use the same parameters than: 
#  Pain, O, Glanville, KP, Hagenaars, S, Selzam, S, Fürtjes, A, Coleman, JRI, Rimfeld, K, Breen, G, Folkersen, L & Lewis, CM 2021, 
#  'Imputed gene expression risk scores: a functionally informed component of polygenic risk', 
#  Human Molecular Genetics, vol. 30, no. 8, pp. 727-738. https://doi.org/10.1093/hmg/ddab053

MHC = data.frame(seqnames = "chr6", start= 28000000, end=34000000, strand="*") %>% GRanges

locs.gr = all.info2 %>% GRanges
ovs <- findOverlaps(MHC, locs.gr,
             minoverlap=1L, # Min number of bases that overlap
             select="all", 
             ignore.strand=TRUE)
ovs.gr2 <- locs.gr[subjectHits(ovs)]
ovs.mhc <- ovs.gr2 %>% as.data.frame %>% arrange(pvalue)
ovs.mhc.sig <- ovs.mhc[1,] # 

df_correlations = rbind(df_correlations,data.frame(gene_clumping=c(ovs.mhc.sig$gene),
                              genes_ov = c(ovs.mhc$gene),
                              correlation = c(-9)))
list.final.genes <- c(list.final.genes, ovs.mhc$gene)
all.info2 = all.info2[!all.info2$gene %in% list.final.genes,]
print(paste("Number of genes identified in MHC region: ", length(list.final.genes), sep=""))
# Loop
#----------------------------------------
print("Running Loop")
for (gene in all.info$gene) {
  if (gene %in% all.info2$gene) {
    tested_gene = gene
    # Convert to GRanges object
    locs.gr = all.info2 %>% GRanges
    gene_to_test = as.data.frame(locs.gr[locs.gr$gene==tested_gene,])
    # Create 'clumping' region in the most significative gene of the list
    #----------------------------------------
    gene_to_test$start = gene_to_test$start-kb*1000
    gene_to_test$end = gene_to_test$end+kb*1000
    gene_to_test.gr = gene_to_test %>% GRanges
    # find overlaps between clumping region
    #----------------------------------------
    ov <-findOverlaps(gene_to_test.gr, locs.gr, 
                      minoverlap=1L, # Min number of bases that overlap
                      select="all", 
                      ignore.strand=TRUE) 
    ov.gr2 <- locs.gr[subjectHits(ov)] #Hits in df2 overlapping in df1
    # In case we have overlap
    #----------------------------------------
    if(nrow(as.data.frame(cors[,colnames(cors) %in% as.character(ov.gr2$gene)]))>1) {
      # Correlations
      #----------------------------------------
      correlations <- cor(as.data.frame(cors[,colnames(cors) %in% as.character(ov.gr2$gene)]))
      cors2 <- correlations[colnames(correlations) %in% gene_to_test,]^2
      #Select correlated genes:
      correlated.genes <- cors2[cors2>cor_threshold]
      list.final.genes <- c(list.final.genes, names(correlated.genes))
      #print(correlated.genes)
      df_cors <- data.frame(gene_clumping=rep(tested_gene, length(correlated.genes)),
                            genes_ov = names(correlated.genes),
                            correlation = correlated.genes)
      df_correlations = rbind(df_correlations, df_cors)
    }
    # In case we don't have overlap: 
    #----------------------------------------
    if(nrow(as.data.frame(cors[,colnames(cors) %in% as.character(ov.gr2$gene)]))==1) {
      list.final.genes <- c(list.final.genes, tested_gene)
      df_cors <- data.frame(gene_clumping=tested_gene,
                            genes_ov =tested_gene,
                            correlation = 1)
      df_correlations = rbind(df_correlations, df_cors)
    }
    # Remove those correlated genes from list and save information:
    #----------------------------------------
    all.info2 = all.info2[!all.info2$gene %in% list.final.genes,]
    #print(nrow(all.info2))
  }
  # In case we removed the gene from the list, caused by other overlaps.
  #----------------------------------------
  else {
    next
  }
}
print(paste("All genes used? ", names(table(all.info$gene %in% list.final.genes))==TRUE, sep=""))
print(paste("Loop finished. Used: ", length(list.final.genes), " genes about ", nrow(all.info)," initial genes", sep=""))

row.names(df_correlations) <- NULL

#Create collapsed table for the results
#----------------------------------------
df_correlations.red <- df_correlations %>%
  group_by(gene_clumping) %>%
  summarise(genes_ov = toString(genes_ov),
            correlation = toString(correlation)
            ) %>%
  ungroup()

#Save results
#----------------------------------------
write.table(df_correlations, paste(workingdir,"/", output, "_df_correlations.txt", sep=""),quote=F, row.names=F, col.names = T)
write.table(df_correlations.red, paste(workingdir,"/", output, "_df_correlations_collapsed.txt", sep=""),quote=F, row.names=F, col.names = T)
write.table(df_correlations.red[,1], paste(workingdir,"/", output, "_selected_genes.txt", sep=""),quote=F, row.names=F, col.names = T)

# Count time spent in script
#----------------------------------------
end.time <- Sys.time()
print(paste("Total Time running the script: ", f(start.time, end.time), sep =""))
print("Finished :) ")
sink() #Close conection to file
