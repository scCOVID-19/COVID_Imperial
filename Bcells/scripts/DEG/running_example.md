# Example workflow
### Load the functions
```R
source('scripts/glmm_functions.R') # I store the functions here
```

### Read in sce object and perform basic setup of variables
```R
setwd('/lustre/scratch117/cellgen/team297/kt16/COVID_imperial_renal/')
sce <- readRDS('h5ad/df.fil3_gex_bcells_vdj_sce_B_mem_all.RDS')
counts(sce) <- assays(sce)[['X']] # because i'm saving from a h5ad object with anndata2ri
sce$case_control <- factor(sce$case_control, levels = c('NEGATIVE', 'POSITIVE', 'RECOVERY'))
sce$WHO_temp_severity <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'))
sce$WHO_temp_severity_group <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('NA', 'mild_moderate', 'mild_moderate', 'severe_critical', 'severe_critical'))
sce$WHO_severity_group <- factor(sce$WHO_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('NA', 'mild_moderate', 'mild_moderate', 'severe_critical', 'severe_critical')) # interpreted as peak severity
sce$grouped_temp_severity <- ifelse(sce$WHO_temp_severity %in% c("mild", "moderate"), "mild_moderate", "severe_critical")
sce$grouped_severity <- ifelse(sce$WHO_severity %in% c("mild", "moderate"), "mild_moderate", "severe_critical")
sce$age_scaled <- scale(sce$calc_age) # scale age
```

### Basic parameters to use.
```R
min_cells = 10
ncpus = 10
```

### Comparison 1: deg from ordered WHO temp severity
```R
# Remove samples with less than 10 cells
nCells <- table(sce$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
sce1 <- sce[, !sce$sample_id %in% rmSamples]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id = as.character(colData(sce1)[, c("sample_id")]))
y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))

y1 <- setupDGElist(y, "WHO_temp_severity_group", remove = "NA") # use grouped_temp_severity later
# sanity check
table(y1$samples$grouped_temp_severity, y1$samples$individual_id)
table(y1$samples$grouped_temp_severity, y1$samples$centre)
table(y1$samples$grouped_temp_severity, y1$samples$sex)
table(y1$samples$grouped_temp_severity, y1$samples$ethnicity)

res1 <- testDGElist(y1,
            formula = as.formula("~ grouped_temp_severity + sex + ethnicity + age_scaled + centre + (1|individual_id)"),
            individual_id = "individual_id",
            ncores = ncpus)
# print(colnames(res1@stats))
results1 <- degTable(res1, contrast = 'grouped_temp_severity', group = 'severe_critical')
```

### Comparison 2: deg from the interaction between WHO (peak) severity groups and time from infection (time from first symptoms)
```R
# remove samples >21 days
sce2 <- sce[,sce$time_from_infection <= 21] # prevent issues with spline
nCells <- table(sce2$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
sce2 <- sce2[, !sce2$sample_id %in% rmSamples]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce2, id = as.character(colData(sce2)[, c("sample_id")]))
y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))
y2 <- setupDGElist(y, "WHO_severity_group", remove = "NA") # use grouped_severity later
# sanity check
table(y2$samples$grouped_severity, y2$samples$individual_id)
table(y2$samples$grouped_severity, y2$samples$centre)
table(y2$samples$grouped_severity, y2$samples$sex)
table(y2$samples$grouped_severity, y2$samples$ethnicity)

res2 <- testDGElist(y2,
            formula = as.formula("~ splines::bs(time_from_infection, degree = 2) * grouped_severity + sex + ethnicity + age_scaled + centre + (1|individual_id)"),
            individual_id = "individual_id",
            modified = TRUE,
            ncores = ncpus)
# print(colnames(res2$stats))
results2 <- degTable(res2, contrast = 'splines::bs(time_from_infection, degree = 2):grouped_severity', 'severe_critical', modified = TRUE)
```

### Comparison 3: just wave 1, deg from positive vs negative
```R
sce3 <- sce[, sce$centre == 'NCL']
# Remove samples with less than MIN
nCells <- table(sce3$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
sce3 <- sce3[,!sce3$sample_id %in% rmSamples]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce3, id=as.character(colData(sce3)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y3 <- setupDGElist(y, 'case_control')
# sanity check
table(y3$samples$case_control, y3$samples$individual_id)
table(y3$samples$case_control, y3$samples$sex)
table(y3$samples$case_control, y3$samples$ethnicity)

res3 <- testDGElist(y3, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = ncpus
           )
# print(colnames(res3@stats))
results3 <- degTable(res3, contrast = 'case_control', group = 'POSITIVE')
```

### Comparison 4: just patients that were negative in wave 1 but positive in wave 2, deg from recovery vs negative
```R
sce4 <- sce[, sce$individual_id %in% c('C101', 'C108', 'C137', 'C138', 'C140', 
                                       'C141', 'C145', 'C146', 'C147', 'C168',
                                        'C169', 'C170', 'C187', 'C190', 'C33')]
# Remove samples with less than MIN
nCells <- table(sce4$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
sce4 <- sce4[,!sce4$sample_id %in% rmSamples]
# remove non-complete data (all have positive), 1 = negative, 3 = recovery
df <- table(sce4$individual_id, sce4$case_control)
keep_ids <- row.names(df)[which(df[,1] != 0 & df[,3] != 0)]
sce4 <- sce4[, sce4$individual_id %in% keep_ids]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce4, id=as.character(colData(sce4)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y4 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')
# sanity check
table(y4$samples$case_control, y4$samples$individual_id)
table(y4$samples$case_control, y4$samples$sex)
table(y4$samples$case_control, y4$samples$ethnicity)

res4 <- testDGElist(y4, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = ncpus
           )
# print(colnames(res4@stats))
results4 <- degTable(res4, contrast = 'case_control', group = 'RECOVERY')
```

### Comparison 5: just patients that were negative in wave 1 but positive in wave 2, deg from positive vs negative
```R
sce5 <- sce[, sce$individual_id %in% c('C101', 'C108', 'C137', 'C138', 'C140', 
                                       'C141', 'C145', 'C146', 'C147', 'C168',
                                        'C169', 'C170', 'C187', 'C190', 'C33')]
# Remove samples with less than MIN
nCells <- table(sce5$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
sce5 <- sce5[,!sce5$sample_id %in% rmSamples]
# remove non-complete data (all have positive), 1 = negative, 2 = positive
df <- table(sce5$individual_id, sce5$case_control)
keep_ids <- row.names(df)[which(df[,1] != 0 & df[,2] != 0)]
sce5 <- sce5[, sce5$individual_id %in% keep_ids]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce5, id=as.character(colData(sce5)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y5 <- setupDGElist(y, 'case_control', remove = 'RECOVERY')
# sanity check
table(y5$samples$case_control, y5$samples$individual_id)
table(y5$samples$case_control, y5$samples$sex)
table(y5$samples$case_control, y5$samples$ethnicity)

res5 <- testDGElist(y5, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = ncpus
           )
# print(colnames(res5@stats))
results5 <- degTable(res5, contrast = 'case_control', group = 'POSITIVE')
```