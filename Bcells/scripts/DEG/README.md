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
1. case vs control (NEGATIVE > POSITIVE > RECOVERY)
2. case vs control (POSITIVE - NEGATIVE)
3. case vs control (RECOVERY - NEGATIVE)
4. case vs control (POSITIVE - NEGATIVE)
5. WHO temp severity (mild > moderate severe > critical)
6. WHO temp severity group ((severe and critical) - (mild and moderate))
7. Date of sampling - date of first symptoms (Days from first symptoms)

These constrasts are repeated in:
1. all samples combined
2. specific samples from patients who were negative in wave 1 but positive in wave 2
3. Only samples in wave 1 cohort
4. Only samples in wave 2 cohort
