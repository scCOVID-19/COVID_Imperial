library(Seurat)
library(data.table)
setwd("/rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/DEG/")
source("/rds/general/user/emacdona/projects/covid19-transcriptome/live/sc_rnaseq/scripts/glmm_functions_4.R")

args<-commandArgs(T)

# Basic parameters to use.
min_cells<-10
cell_type<-args[1]
level<-args[2]
ncpus<-10

print(paste0(cell_type," level ",level," comparison 1"))

data<-readRDS(paste0("../prep_DEG_files/TNK_",cell_type,"_level_",level,"_gex.RDS"))
data

sce<-data

sce$individual_id<-factor(sce$individual_id)
sce$sample_id<-factor(sce$sample_id)
sce$ethnicity<-factor(sce$ethnicity)
sce$sex<-factor(sce$sex)
sce$case_control <- factor(sce$case_control, levels = c('NEGATIVE', 'POSITIVE', 'RECOVERY'))
sce$WHO_temp_severity <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'))
sce$WHO_temp_severity_group <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('NA', 'mild_moderate', 'mild_moderate', 'severe_critical', 'severe_critical'))
sce$WHO_severity_group <- factor(sce$WHO_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('NA', 'mild_moderate', 'mild_moderate', 'severe_critical', 'severe_critical')) # interpreted as peak severity
sce$grouped_temp_severity <- ifelse(sce$WHO_temp_severity %in% c("mild", "moderate"), "mild_moderate", "severe_critical")
sce$grouped_severity <- ifelse(sce$WHO_severity %in% c("mild", "moderate"), "mild_moderate", "severe_critical")
sce$age_scaled <- scale(sce$calc_age) # scale age

# Remove samples with less than 10 cells
nCells <- table(sce$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
rm_ids<-c("C116","C139")
sce1 <- sce[, !sce$sample_id %in% rmSamples & !sce$individual_id %in% rm_ids]

# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id = as.character(colData(sce1)[, c("sample_id")]))
smrzd <- smrzd[, smrzd$case_control == "POSITIVE"]

y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))

y1 <- setupDGElist(y, "WHO_temp_severity_group", remove = "NA")

# sanity check
table(y1$samples$WHO_temp_severity_group, y1$samples$individual_id)
table(y1$samples$WHO_temp_severity_group, y1$samples$centre)
table(y1$samples$WHO_temp_severity_group, y1$samples$sex)
table(y1$samples$WHO_temp_severity_group, y1$samples$corrected_ethnicity)

res1 <- testDGElist(y1,
            formula = as.formula("~ WHO_temp_severity_group + sex + PC1_nonafricanVsAfrican + PC2_asianVsEuropean + age_scaled + centre + (1|individual_id)"),
            individual_id = "individual_id",
            modified = TRUE,
            ncores = ncpus)

results1 <- degTable_modified(res1, contrast = 'WHO_temp_severity_group', group = 'severe_critical')


save(res1,results1,file=paste0("./level_",level,"/comp_1/",cell_type,"_all_deg_new_ethnicity_new.RData"))

writeLines("Done!")
