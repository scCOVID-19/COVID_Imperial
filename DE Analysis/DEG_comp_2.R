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

print(paste0(cell_type," level ",level," comparison 2"))

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

# remove samples >21 days
sce2 <- sce[,sce$time_from_infection <= 21] # prevent issues with spline
nCells <- table(sce2$sample_id)
rmSamples <- names(nCells[nCells < min_cells])
rm_ids<-c("C116","C139")
sce2 <- sce2[, !sce2$sample_id %in% rmSamples & !sce2$individual_id %in% rm_ids]

# Summarize Counts
smrzd <- aggregateAcrossCells(sce2, id = as.character(colData(sce2)[, c("sample_id")]))
smrzd <- smrzd[, smrzd$case_control == "POSITIVE"]

y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))
y2 <- setupDGElist(y, "WHO_severity_group", remove = "NA") # use grouped_severity later

# sanity check
table(y2$samples$WHO_severity_group, y2$samples$individual_id)
table(y2$samples$WHO_severity_group, y2$samples$centre)
table(y2$samples$WHO_severity_group, y2$samples$sex)
table(y2$samples$WHO_severity_group, y2$samples$corrected_ethnicity)

res2 <- testDGElist(y2,
            formula = as.formula("~ splines::bs(time_from_infection, degree = 2) * WHO_severity_group + sex + PC1_nonafricanVsAfrican + PC2_asianVsEuropean + age_scaled + centre + (1|individual_id)"),
            individual_id = "individual_id",
            modified = TRUE,
            ncores = ncpus)


results2 <- degTable_modified2(res2, contrast = 'splines..bs.time_from_infection..degree...2.2.WHO_severity_group', contrast_p='splines..bs.time_from_infection..degree...2..WHO_severity_group',group='severe_critical')

save(res2,results2,file=paste0("./level_",level,"/comp_2/",cell_type,"_all_deg_new_ethnicity_new.RData"))

writeLines("Done!")
