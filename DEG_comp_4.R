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

print(paste0(cell_type," level ",level," comparison 4"))

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
rm_ids<-c("C116","C139","C141")
sce4 <- sce4[, sce4$individual_id %in% keep_ids & !sce4$individual_id %in% rm_ids]

# Summarize Counts
smrzd <- aggregateAcrossCells(sce4, id=as.character(colData(sce4)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y4 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')

# sanity check
table(y4$samples$case_control, y4$samples$individual_id)
table(y4$samples$case_control, y4$samples$sex)
table(y4$samples$case_control, y4$samples$corrected_ethnicity)

res4 <- testDGElist(y4, 
            formula = as.formula("~ case_control + sex + PC1_nonafricanVsAfrican + PC2_asianVsEuropean + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            modified = TRUE,
            ncores = ncpus)

results4 <- degTable_modified(res4, contrast = 'case_control', group = 'RECOVERY')

save(res4,results4,file=paste0("./level_",level,"/comp_4/",cell_type,"_all_deg_new_ethnicity_new.RData"))

writeLines("Done!")
