#!/usr/bin/env Rscript
options(warn=-1)
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(glmmSeq))
suppressPackageStartupMessages(library(optparse))

opt_list <- list(make_option(c("-i", "--in"), dest="IN", help="input sce .RDS file."),
                 make_option(c("-m", "--min"), dest="MIN", default= 5, help="minimum number of cells."),
                 make_option(c("-n", "--ncpu"), dest="NCPU", default= 1, help="number of cores"))
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

# setup functions
setupDGElist <- function(dgelist, comparison, ordered = FALSE, remove = NULL){
    if (!is.null(remove)){
        dgelist <- dgelist[,!dgelist$samples[,comparison] %in% remove]
    }
    if (!is.factor(dgelist$samples[,comparison])){
        keep <- suppressWarnings(filterByExpr(dgelist, group=NULL, min.count=3, min.total.count=5))
    } else {
        keep <- filterByExpr(dgelist, group=dgelist$samples[,comparison], min.count=3, min.total.count=5)
    }
    dgelist <- dgelist[keep,]
    
    # ensure the factor levels are correct
    if (is.factor(dgelist$samples[,comparison])){
        dgelist$samples[,comparison] <- droplevels(dgelist$samples[,comparison])
        if (ordered){
            # ensure the factor levels are correct
            dgelist$samples[,comparison] <- droplevels(dgelist$samples[,comparison])
            dgelist$samples[,comparison] <- ordered(dgelist$samples[,comparison])
        } else {
            dgelist$samples[,comparison] <- droplevels(dgelist$samples[,comparison])
        }
    }
    return(dgelist)
}

testDGElist <- function(dgelist, comparison, ordered = FALSE, ncores = NULL) {
    if (is.null(ncores)){
        NCORES <- parallel::detectCores()-1
    } else {
        NCORES <- ncores
    }
    # Estimate Dispersion
    disp  <- suppressMessages(setNames(edgeR::estimateDisp(dgelist)$tagwise.dispersion, rownames(dgelist)))
    # Norm
    sizeFactors <- calcNormFactors(dgelist$counts)
    if (ordered){
        results <- glmmSeq(as.formula(paste0("~ ", comparison," + sex + ethnicity + calc_age + (1|individual_id)")),
                  id = "individual_id",
                  countdata = dgelist$counts,
                  metadata = dgelist$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  designMatrix = model.matrix(as.formula(paste0("~ ", comparison," + sex + ethnicity + calc_age")), data=dgelist$samples),
                  removeSingles = FALSE,
                  progress = TRUE, cores = NCORES)
    } else {
        results <- glmmSeq(as.formula(paste0("~ ", comparison," + sex + ethnicity + calc_age + (1|individual_id)")),
                  id = "individual_id",
                  countdata = dgelist$counts,
                  metadata = dgelist$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles = FALSE,
                  progress = TRUE, cores = NCORES)
    }    
    return(results)
}

sce <- readRDS(opt$IN)
counts(sce) <- assays(sce)[['X']] # because i'm saving from a h5ad object with anndata2ri
sce$case_control <- factor(sce$case_control, levels = c('NEGATIVE', 'POSITIVE', 'RECOVERY'))
sce$WHO_temp_severity <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'))
# severity groups
sce$WHO_temp_severity_group <- factor(sce$WHO_temp_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('NA', 'mild_mod', 'mild_mod', 'sev_crit', 'sev_crit'))

##### general comparison
# Remove samples with less than MIN
nCells <- table(sce$sample_id)
rmSamples <- names(nCells[nCells<opt$MIN])
sce1 <- sce[,!sce$sample_id %in% rmSamples]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id=as.character(colData(sce1)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))

# setup the comparisons
y1 <- setupDGElist(y, 'case_control', ordered = TRUE)
y2 <- setupDGElist(y, 'case_control', remove = 'RECOVERY')
y3 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')
y4 <- setupDGElist(y, 'case_control', remove = 'NEGATIVE')
y5 <- setupDGElist(y, 'WHO_temp_severity', ordered = TRUE, remove = "NA")
y6 <- setupDGElist(y, 'WHO_temp_severity_group', remove = "NA")
y7 <- setupDGElist(y, 'days_to_sampling') # days from first symptoms

# run glmmseq
res1 <- tryCatch(testDGElist(y1, 'case_control', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res2 <- tryCatch(testDGElist(y2, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res3 <- tryCatch(testDGElist(y3, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res4 <- tryCatch(testDGElist(y4, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res5 <- tryCatch(testDGElist(y5, 'WHO_temp_severity', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res6 <- tryCatch(testDGElist(y6, 'WHO_temp_severity_group', ncores = opt$NCPU), error = function(e) return(NA))
res7 <- tryCatch(testDGElist(y7, 'days_to_sampling', ncores = opt$NCPU), error = function(e) return(NA))

OUT = paste0(dirname(opt$IN), '/', gsub('.RDS', '_glmmSeq_results.RData', basename(opt$IN)))
save(res1, res1, res2, res3, res4, res5, res6, res7, file = OUT)

##### wave 1 + wave 2 specific comparison
# duplicate sce to subset to only patients that were negative in wave1 but positive in wave2
sce2 <- sce[, sce$individual_id %in% c('C101', 'C108', 'C137', 'C138', 'C140', 'C141', 'C145', 'C146', 'C147', 'C168', 'C169', 'C170', 'C187', 'C190', 'C33')]
# Remove samples with less than MIN
nCells2 <- table(sce2$sample_id)
rmSamples2 <- names(nCells2[nCells2<opt$MIN])
sce2 <- sce2[,!sce2$sample_id %in% rmSamples2]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce2, id=as.character(colData(sce2)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))

# setup the comparisons
y1 <- setupDGElist(y, 'case_control', ordered = TRUE)
y2 <- setupDGElist(y, 'case_control', remove = 'RECOVERY')
y3 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')
y4 <- setupDGElist(y, 'case_control', remove = 'NEGATIVE')
y5 <- setupDGElist(y, 'WHO_temp_severity', ordered = TRUE, remove = "NA")
y6 <- setupDGElist(y, 'WHO_temp_severity_group', remove = "NA")
y7 <- setupDGElist(y, 'days_to_sampling') # days from first symptoms

# run glmmseq
res1 <- tryCatch(testDGElist(y1, 'case_control', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res2 <- tryCatch(testDGElist(y2, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res3 <- tryCatch(testDGElist(y3, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res4 <- tryCatch(testDGElist(y4, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res5 <- tryCatch(testDGElist(y5, 'WHO_temp_severity', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res6 <- tryCatch(testDGElist(y6, 'WHO_temp_severity_group', ncores = opt$NCPU), error = function(e) return(NA))
res7 <- tryCatch(testDGElist(y7, 'days_to_sampling', ncores = opt$NCPU), error = function(e) return(NA))

OUT = paste0(dirname(opt$IN), '/', gsub('.RDS', '_glmmSeq_results_wave1and2.RData', basename(opt$IN)))
save(res1, res2, res3, res4, res5, res6, res7, file = OUT)

##### wave1 comparison
sce3 <- sce[, sce$centre == 'NCL']
# Remove samples with less than MIN
nCells3 <- table(sce3$sample_id)
rmSamples3 <- names(nCells3[nCells3<opt$MIN])
sce3 <- sce3[,!sce3$sample_id %in% rmSamples3]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce3, id=as.character(colData(sce3)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))

# setup the comparisons
y1 <- setupDGElist(y, 'case_control', ordered = TRUE)
y2 <- setupDGElist(y, 'case_control', remove = 'RECOVERY')
y3 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')
y4 <- setupDGElist(y, 'case_control', remove = 'NEGATIVE')
y5 <- setupDGElist(y, 'WHO_temp_severity', ordered = TRUE, remove = "NA")
y6 <- setupDGElist(y, 'WHO_temp_severity_group', remove = "NA")
y7 <- setupDGElist(y, 'days_to_sampling') # days from first symptoms

# run glmmseq
res1 <- tryCatch(testDGElist(y1, 'case_control', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res2 <- tryCatch(testDGElist(y2, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res3 <- tryCatch(testDGElist(y3, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res4 <- tryCatch(testDGElist(y4, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res5 <- tryCatch(testDGElist(y5, 'WHO_temp_severity', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res6 <- tryCatch(testDGElist(y6, 'WHO_temp_severity_group', ncores = opt$NCPU), error = function(e) return(NA))
res7 <- tryCatch(testDGElist(y7, 'days_to_sampling', ncores = opt$NCPU), error = function(e) return(NA))

OUT = paste0(dirname(opt$IN), '/', gsub('.RDS', '_glmmSeq_results_wave1.RData', basename(opt$IN)))
save(res1, res2, res3, res4, res5, res6, res7, file = OUT)

##### wave2 comparison
sce4 <- sce[, sce$centre == 'Cambridge']
# Remove samples with less than MIN
nCells4 <- table(sce4$sample_id)
rmSamples4 <- names(nCells4[nCells4<opt$MIN])
sce4 <- sce4[,!sce4$sample_id %in% rmSamples4]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce4, id=as.character(colData(sce4)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))

# setup the comparisons
y1 <- setupDGElist(y, 'case_control', ordered = TRUE)
y2 <- setupDGElist(y, 'case_control', remove = 'RECOVERY')
y3 <- setupDGElist(y, 'case_control', remove = 'POSITIVE')
y4 <- setupDGElist(y, 'case_control', remove = 'NEGATIVE')
y5 <- setupDGElist(y, 'WHO_temp_severity', ordered = TRUE, remove = "NA")
y6 <- setupDGElist(y, 'WHO_temp_severity_group', remove = "NA")
y7 <- setupDGElist(y, 'days_to_sampling') # days from first symptoms

# run glmmseq
res1 <- tryCatch(testDGElist(y1, 'case_control', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res2 <- tryCatch(testDGElist(y2, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res3 <- tryCatch(testDGElist(y3, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res4 <- tryCatch(testDGElist(y4, 'case_control', ncores = opt$NCPU), error = function(e) return(NA))
res5 <- tryCatch(testDGElist(y5, 'WHO_temp_severity', ordered = TRUE, ncores = opt$NCPU), error = function(e) return(NA))
res6 <- tryCatch(testDGElist(y6, 'WHO_temp_severity_group', ncores = opt$NCPU), error = function(e) return(NA))
res7 <- tryCatch(testDGElist(y7, 'days_to_sampling', ncores = opt$NCPU), error = function(e) return(NA))

OUT = paste0(dirname(opt$IN), '/', gsub('.RDS', '_glmmSeq_results_wave2.RData', basename(opt$IN)))
save(res1, res2, res3, res4, res5, res6, res7, file = OUT)
