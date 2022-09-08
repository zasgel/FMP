#!/bin/bash

#Set job requirements
#SBATCH --job-name=TWAS
#SBATCH -N 1
#SBATCH --tasks-per-node 1

#############
#module load#
#############
module load 2020
module load R/4.0.2-intel-2020a

###########
#variables#
###########
tissues_wd=/home/juditc/ADHD/FE/TWAS/tissues/
cd $tissues_wd 
tissues_file=$(ls | sed -n ${SLURM_ARRAY_TASK_ID}p)


#####################
#copy all to scratch#
#####################
dir_name=TWAS_${SLURM_ARRAY_TASK_ID}

mkdir "$TMPDIR"/${dir_name}
cp ${tissues_wd}${tissues_file} "$TMPDIR"/${dir_name}
cp /home/juditc/ADHD/FE/TWAS/final_covariates_gwas_4_2021_marray_PRS_corrected.csv "$TMPDIR"/${dir_name}
cp /home/juditc/ADHD/FE/TWAS/script_tissue_subject.R "$TMPDIR"/${dir_name}

cd "$TMPDIR"/${dir_name}


Rscript script_tissue_subject.R $tissues_file $subject


##############
#save results#
##############
cp  "$TMPDIR"/${dir_name}/$subject_$tissue.log /home/juditc/ADHD/FE/TWAS/
cp  "$TMPDIR"/${dir_name}/$subject_$tissue.out /home/juditc/ADHD/FE/TWAS/
cp  "$TMPDIR"/${dir_name}/$subject_$tissue.fm /home/juditc/ADHD/FE/TWAS/

