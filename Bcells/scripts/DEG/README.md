### DEG scripts

Call `glmmSeq_jobs.R` with:

```bash
Rscript --vanilla glmmSeq_jobs.R -i [input sce .RDS file.] -m [minimum number of cells] -n [number of cores] 

# example
Rscript --vanilla glmmSeq_jobs.R \
       -i df.fil3_gex_bcells_vdj_sce_B_ASC_dividing.RDS \
       -m 5 \
       -n 24
```

The current script will run the following tests in a LMM model, accounting for age, gender and ethnicity and the individual is specified as a random effect:

#### Contrasts
1. ~case vs control (NEGATIVE > POSITIVE > RECOVERY)~
2. case vs control (POSITIVE - NEGATIVE)
3. case vs control (RECOVERY - NEGATIVE)
4. case vs control (POSITIVE - RECOVERY)
5. WHO temp severity (mild > moderate > severe > critical)
6. WHO temp severity group ((severe and critical) - (mild and moderate))
7. Date of sampling - date of first symptoms (Days from first symptoms)

# Fig1
## General figure
       UMAP
       Sampling setup
       glmm - celltype proportion
              deg ~ age + gender + ethnicity + centre/wave (all samples combined) + WHO temp severity group
       severity comparison (no negs, no recv) - all samples combined
              deg ~ age + gender + ethnicity + centre/wave (all samples combined) + WHO temp severity group
              jack's model

# Fig 2
## T and B Celltypes
       1 UMAP
       1 compositional - Milo
              *1) deg ~ age + gender + ethnicity + centre/wave (all samples combined) + WHO temp severity group + (1 | individual_id)
       wilcoxon first, then this
              *2) deg ~ age + gender + ethnicity + centre/wave + spline(time_from_first_symptoms) * WHO_severity + (1 | individual_id)
       Marker-genes - supp
       ***V(D)J - diversity vs clonality over time
       **pathway - cytoxtoxic abilities over time
       **Antigen presentation/BCRsignaling over time
       **Cytokine/chemokine dotplot - over time
       **celltype interaction - cellphonedb over time  ?query Erin Tfh - B

# Fig 3
## MNP
       1 UMAP
       1 compositional - Milo
              *1) deg ~ age + gender + ethnicity + centre/wave (all samples combined) + WHO temp severity group + (1 | individual_id)
       wilcoxon first, then this
              *2) gene ~ age + gender + ethnicity + centre/wave + spline(time_from_first_symptoms) * WHO_severity + (1 | individual_id)
       Marker-genes - supp
       **Monocyte DC, projenitor trajectories
       **Antigen presentation capacity over time
       **Cytokine/chemokine dotplot - over time
       **monoyte-platelet interaction - cellphonedb over time

# Fig 4
# wave 2 neg vs recovery - are patients back to baseline per celltype
specific samples from patients who were negative in wave 1 but positive in wave 2 
       case vs control (RECOVERY - NEGATIVE)
       case vs control (POSITIVE - NEGATIVE)
       
# Fig 5
versus Stephenson et al
renal disease specific things?
integrative analysis - cytokine/secreted molecules
