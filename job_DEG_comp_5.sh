#!/bin/sh
#PBS -lwalltime=03:00:00
#PBS -lselect=1:ncpus=10:mem=400gb
#PBS -N DEG_comp_5_new_ethnicity_new
#PBS -V
#PBS -j oe
#PBS -J 1-65
#PBS -o /rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/DEG/logs/

module load anaconda3/personal
eval "$(conda shell.bash hook)"
conda activate testing --stack

cd /rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/DEG/
input=`sed -n "${PBS_ARRAY_INDEX}p" DEG_1-5_input_new.txt`

Rscript /rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/scripts/DEG_comp_5.R $input

echo Done
