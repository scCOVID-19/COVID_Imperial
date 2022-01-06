#!/bin/bash
#BSUB -G teichlab
#BSUB -P team205
#BSUB -J glmmSeq
#BSUB -q normal
#BSUB -o %J.out
#BSUB -n 24
#BSUB -R "select[mem>80000] rusage[mem=80000] span[hosts=1]" -M80000
### ~~~ job script below ~~~ ###

Rscript --vanilla glmmSeq_jobs.R \
       -i /lustre/scratch117/cellgen/team297/kt16/COVID_imperial_renal/h5ad/df.fil3_gex_bcells_vdj_sce_B_CD11c.RDS \
       -o /lustre/scratch117/cellgen/team297/kt16/COVID_imperial_renal/h5ad/glmmSeq_B_CD11c.RData \
       -m 5 \
       -n 24

