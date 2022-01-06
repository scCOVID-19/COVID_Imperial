#!/bin/bash
#BSUB -G teichlab
#BSUB -P team205
#BSUB -J dandelion
#BSUB -q normal
#BSUB -o log/%J.out
#BSUB -e log/%J.err
#BSUB -n 16
#BSUB -M 80000
#BSUB -R "select[mem>80000] rusage[mem=80000] span[hosts=1]"
### ~~~ job script below ~~~ ###

cd /lustre/scratch117/cellgen/team297/kt16/COVID_imperial_renal/processed/BCR
singularity run -B $PWD /nfs/team297/kt16/Softwares/sc-dandelion_dev.sif dandelion-preprocess --meta="bcr_meta.csv" --sep="-" --file_prefix="all" --keep_trailing_hyphen_number --clean_output --skip_format_header

