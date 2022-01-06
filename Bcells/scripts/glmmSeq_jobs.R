#!/usr/bin/env Rscript
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
			     make_option(c("-o", "--out"), dest="OUT", help="output result .RData file."),
				 make_option(c("-m", "--min"), dest="MIN", default= 5, help="minimum number of cells."),
				 make_option(c("-n", "--ncpu"), dest="NCPU", default= 1, help="number of cores")
	)
# Parse arguments
opt <- parse_args(OptionParser(option_list=opt_list))

sce <- readRDS(opt$IN)
counts(sce) <- assays(sce)[['X']]
sce$WHO_severity <- factor(sce$WHO_severity, levels = c('NA', 'mild', 'moderate', 'severe', 'critical'), labels = c('control', 'mild', 'moderate', 'severe', 'critical'))

# Remove samples with less than MIN
nCells <- table(sce$sample_id)
rmSamples <- names(nCells[nCells<opt$MIN])
sce <- sce[,!sce$sample_id %in% rmSamples]
# Summarize Counts
smrzd <- aggregateAcrossCells(sce, id=as.character(colData(sce)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
keep <- filterByExpr(y, group=y$samples$case_control, min.count=3, min.total.count=5)
y <- y[keep,]
y$samples$case_control <- droplevels(y$samples$case_control)
y$samples$WHO_severity <- droplevels(y$samples$WHO_severity)
y$samples$WHO_severity <- ordered(y$samples$WHO_severity)
# Estimate Dispersion
disp  <- suppressMessages(setNames(edgeR::estimateDisp(y)$tagwise.dispersion, rownames(y)))
# Norm
sizeFactors <- calcNormFactors(y$counts)

results1 <- glmmSeq(~ case_control + sex + ethnicity + calc_age + (1|individual_id),
                  id = "individual_id",
                  countdata = y$counts,
                  metadata = y$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = opt$NCPU)
results2 <- glmmSeq(~ case_control + sex + ethnicity + calc_age + days_to_admission + (1|individual_id),
                  id = "individual_id",
                  countdata = y$counts,
                  metadata = y$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = opt$NCPU)
results3 <- glmmSeq(~ WHO_severity + sex + ethnicity + calc_age + (1|individual_id),
                  id = "individual_id",
                  countdata = y$counts,
                  metadata = y$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = opt$NCPU)
results4 <- glmmSeq(~ WHO_severity + sex + ethnicity + calc_age + days_to_admission + (1|individual_id),
                  id = "individual_id",
                  countdata = y$counts,
                  metadata = y$samples,
                  dispersion = disp,
                  sizeFactors = sizeFactors,
                  removeDuplicatedMeasures = FALSE,
                  removeSingles=FALSE,
                  progress=TRUE,
                  cores = opt$NCPU)
saveRDS(results1, results2, results3, results4, file = opt$OUT)
