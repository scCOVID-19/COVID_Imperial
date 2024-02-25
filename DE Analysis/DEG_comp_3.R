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

print(paste0(cell_type," level ",level," comparison 3"))

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


sce3 <- sce[, sce$centre == 'NCL']
# Remove samples with less than MIN
nCells <- table(sce3$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
rm_ids<-c("C116","C139")
sce3 <- sce3[,!sce3$sample_id %in% rmSamples & !sce3$individual_id %in% rm_ids]

# Summarize Counts
smrzd <- aggregateAcrossCells(sce3, id=as.character(colData(sce3)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y3 <- setupDGElist(y, 'case_control')

# sanity check
table(y3$samples$case_control, y3$samples$individual_id)
table(y3$samples$case_control, y3$samples$sex)
table(y3$samples$case_control, y3$samples$corrected_ethnicity)


res3 <- testDGElist(y3, 
            formula = as.formula("~ case_control + sex + PC1_nonafricanVsAfrican + PC2_asianVsEuropean + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            modified = TRUE,
            ncores = ncpus
           )
results3 <- degTable_modified(res3, contrast = 'case_control', group = 'POSITIVE')

save(res3,results3,file=paste0("./level_",level,"/comp_3/",cell_type,"_all_deg_new_ethnicity_new.RData"))


writeLines("Done!")
