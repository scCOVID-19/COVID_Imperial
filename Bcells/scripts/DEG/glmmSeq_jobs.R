#!/usr/bin/env Rscript
options(warn = -1)
suppressPackageStartupMessages(library(SingleCellExperiment))
suppressPackageStartupMessages(library(SummarizedExperiment))
suppressPackageStartupMessages(library(variancePartition))
suppressPackageStartupMessages(library(limma))
suppressPackageStartupMessages(library(edgeR))
suppressPackageStartupMessages(library(scran))
suppressPackageStartupMessages(library(scater))
suppressPackageStartupMessages(library(glmmSeq))
suppressPackageStartupMessages(library(lmerTest))
suppressPackageStartupMessages(library(qvalue))
suppressPackageStartupMessages(library(pbmcapply))
suppressPackageStartupMessages(library(optparse))

opt_list <- list(make_option(c("-i", "--in"), dest = "IN", help = "input sce .RDS file."),
    make_option(c("-m", "--min"), dest = "MIN", default = 5, help = "minimum number of cells."),
    make_option(c("-n", "--ncpu"), dest = "NCPU", default = 1, help = "number of cores"))
# Parse arguments
opt <- parse_args(OptionParser(option_list = opt_list))

# define functions
setupDGElist <- function(dgelist, comparison, ordered = FALSE, remove = NULL) {
    if (!is.null(remove)) {
        dgelist <- dgelist[, !dgelist$samples[, comparison] %in% remove]
    }
    if (!is.factor(dgelist$samples[, comparison])) {
        keep <- suppressWarnings(filterByExpr(dgelist, group = NULL, min.count = 3,
            min.total.count = 5))
    } else {
        keep <- filterByExpr(dgelist, group = dgelist$samples[, comparison], min.count = 3,
            min.total.count = 5)
    }
    dgelist <- dgelist[keep, ]

    # ensure the factor levels are correct
    if (is.factor(dgelist$samples[, comparison])) {
        dgelist$samples[, comparison] <- droplevels(dgelist$samples[, comparison])
        if (ordered) {
            # ensure the factor levels are correct
            dgelist$samples[, comparison] <- droplevels(dgelist$samples[, comparison])
            dgelist$samples[, comparison] <- ordered(dgelist$samples[, comparison])
        } else {
            dgelist$samples[, comparison] <- droplevels(dgelist$samples[, comparison])
        }
    }
    return(dgelist)
}

testDGElist <- function(dgelist, formula, individual_id, optimizer = 'bobyqa', designMatrix = NULL, ncores = NULL, ...) {
    if (is.null(ncores)) {
        NCORES <- parallel::detectCores() - 1
    } else {
        NCORES <- ncores
    }
    # Estimate Dispersion
    disp <- suppressMessages(setNames(edgeR::estimateDisp(dgelist)$tagwise.dispersion,
        rownames(dgelist)))
    # Norm
    sizeFactors <- calcNormFactors(dgelist$counts)
    results <- glmmSeq(formula, 
                       id = individual_id, 
                       countdata = dgelist$counts, 
                       metadata = dgelist$samples,
                       dispersion = disp,
                       sizeFactors = sizeFactors,
                       removeDuplicatedMeasures = FALSE,
                       designMatrix = designMatrix,
                       removeSingles = FALSE,
                       control = glmerControl(optimizer = optimizer, ...),
                       progress = TRUE,
                       cores = NCORES)
    return(results)
}

degTable <- function(results, contrast, group){    
    modelData <- results@modelData
    outLabels <- apply(modelData, 1, function(x) paste(x, collapse="_"))
    modelData$y <- paste0('y_', outLabels)
    cols1 = grep(group, modelData$y, value = TRUE)
    cols2 = grep(group, modelData$y, value = TRUE, invert = TRUE)
    LFC <- log2(rowMeans(results@predict[, cols2])+1) - log2(rowMeans(results@predict[, cols1])+1)
    tmp <- data.frame(results@stats[,c(paste0(contrast, group), paste0('P_',contrast), paste0('q_',contrast))])
    colnames(tmp) <- c('fixed-effects estimates', 'pval', 'qval')
    tmp$LFC <- LFC
    tmp <- tmp[order(-tmp$`fixed-effects estimates`, tmp$qval), ]
    return(tmp)
}

glmm_modified <- function(dgeList,
                             modelFormula,
                             id,
                             control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)),
                             cores = 1,
                             verbose = TRUE,
                             ...) {
    metadata <- dgeList$samples
    countdata <- dgeList$counts
    dispersion <- suppressMessages(setNames(edgeR::estimateDisp(dgeList)$tagwise.dispersion,
        rownames(dgeList)))
    sizeFactors <- calcNormFactors(dgeList$counts)

    ids <- as.character(metadata[, id])
    # Manipulate formulae
    reducedFormula <- nobars(modelFormula)
    designMatrix <- model.matrix(reducedFormula, data = dgeList$samples)
    fullFormula <- update.formula(modelFormula, count ~ ., simplify = FALSE)

    # Check numbers and alignment
    if (!all(rownames(countdata) %in% names(dispersion), nrow(countdata))) {
        stop("Dispersion length must match nrow in countdata")
    }

    if (!is.null(sizeFactors))
        offset <- log(sizeFactors) else offset <- NULL
    if (verbose)
        cat(paste0("\nn = ", length(ids), " samples, ", length(unique(ids)), " individuals\n"))

    start <- Sys.time()
    fullList <- lapply(rownames(countdata), function(i) {
        list(y = countdata[i, ], dispersion = dispersion[i])
    })

    resultList <- pbmclapply(fullList, function(geneList) {
        try(run_glmm(data = metadata, 
                     geneList = geneList, 
                     fullFormula = fullFormula,
                     control = control,
                     offset = offset, ...), silent = TRUE)
    }, mc.cores = cores)

    # Print timing if verbose
    end <- Sys.time()
    if (verbose)
        print(end - start)

    # Output
    names(resultList) <- rownames(countdata)

    # any failed?
    if (any(lapply(resultList, length) < 4)){
        failed <- names(which(lapply(resultList, length) < 4))
        errormsg <- lapply(resultList[failed], function(x) x[1])
        resultList <- resultList[-which(lapply(resultList, length) < 4)]
    }
    
    noErr <- vapply(resultList, function(x) x$tryErrors == "", FUN.VALUE = TRUE)
    if (length(which(noErr)) == 0) { 
        stop("All genes returned an error. Check sufficient data in each group")
    }
    if (sum(!noErr) != 0) {
        if (verbose) {
            if (length(failed) > 0){
                cat(paste0("Errors in ", sum(!noErr, length(failed)), " gene(s):", paste0(c(names(noErr)[! noErr], failed), collapse = ", ")))
            } else {
                cat(paste0("Errors in ", sum(!noErr), " gene(s):", paste0(names(noErr)[! noErr], collapse = ", ")))
            }
        }
        outputErrors <- vapply(resultList[!noErr], function(x) {x$tryErrors}, FUN.VALUE = c("test"))
        if (length(failed > 0)){
            outputErrors <- c(outputErrors, errormsg)
        }
    } else {
            outputErrors <- c("No errors")
    }
    optInfo <- t(vapply(resultList[noErr], function(x) {
        setNames(x$optinfo, c("Singular", "Conv"))
    }, FUN.VALUE = c(1, 1)))
  
    nCheat <- resultList[noErr][[1]]$stats
    s <- t(vapply(resultList[noErr], function(x) {x$stats},
                  FUN.VALUE = rep(1, length(nCheat))))
  
        return(list(stats = s, fit = resultList, optInfo = optInfo, errors = outputErrors))
    }

run_glmm <- function(data, geneList, fullFormula, control, offset, ...) {
    data[, "count"] <- as.numeric(geneList$y)
    fit <- lme4::glmer(fullFormula, data = data, control = control, offset = offset,
        family = MASS::negative.binomial(theta = 1/geneList$dispersion), ...)
    if (class(fit) != "try-error") {
        stats <- setNames(c(geneList$dispersion, AIC(fit),
                        as.numeric(logLik(fit))),
                      c("Dispersion", "AIC", "logLik"))
        fixedEffects <- lme4::fixef(fit)
        wald <- try(car::Anova(fit), silent = TRUE)
        if (class(wald) != "try-error") {
            waldtest <- setNames(c(wald[, "Chisq"], wald[, "Pr(>Chisq)"]),
                 c(paste0("Chisq_", rownames(wald)),
                    paste0("P_", rownames(wald))))
            singular <- as.numeric(isSingular(fit))
            conv <- length(slot(fit, "optinfo")$conv$lme4$messages)
            return(list(stats = c(stats, fixedEffects, waldtest), fit = fit, optinfo = c(singular, conv), tryErrors = ""))
        } else {
            return(list(stats = NA, fit = NA, optinfo = NA, tryErrors = wald[1]))
        }
    } else {
        return(list(stats = NA, fit = NA, optinfo = NA, tryErrors = fit[1]))
    }
}

glmm_qval <- function(result,
                      cutoff = 0.05,
                      pi0 = NULL,
                      verbose = TRUE) {
    resultStats <- data.frame(result$stats, check.names = FALSE)
    for (cn in colnames(resultStats)[grep("P_", colnames(resultStats))]) {
        q_cn <- gsub("P_", "q_", cn)
        resultStats[, q_cn] <- NA
        resultStats[!is.na(resultStats[, cn]), q_cn] <- qvalue(resultStats[!is.na(resultStats[, cn]), cn], pi0 = pi0)$qvalues
        if (verbose) {
            cat(paste0("\n", q_cn, "\n"))
            cat(paste(rep("-", nchar(q_cn)), collapse = ""))
            print(table(ifelse(resultStats[, q_cn] < cutoff, "Significant", "Not Significant")))
        }
    }
    result$stats <- as.matrix(resultStats)
    return(result)
}

### Read in object and basic formatting
sce <- readRDS(opt$IN)
counts(sce) <- assays(sce)[["X"]]  # because i'm saving from a h5ad object with anndata2ri, otherwise, use whichever is the raw counts slot
sce$case_control <- factor(sce$case_control, levels = c("NEGATIVE", "POSITIVE", "RECOVERY"))
sce$WHO_temp_severity <- factor(sce$WHO_temp_severity, levels = c("NA", "mild", "moderate",
    "severe", "critical"))
sce$WHO_severity_group <- factor(sce$WHO_temp_severity, levels = c("NA", "mild",
    "moderate", "severe", "critical"), labels = c("NA", "mild_moderate", "mild_moderate",
    "severe_critical", "severe_critical"))
sce$age_scaled <- scale(sce$calc_age)  # scale age

# Comparison 1: deg from ordered WHO temp severity Remove samples with less than MIN
nCells <- table(sce$sample_id)
rmSamples <- names(nCells[nCells < opt$MIN])
sce1 <- sce[, !sce$sample_id %in% rmSamples]
# drop unused levels
sce1$individual_id <- droplevels(sce1$individual_id)
sce1$sample_id <- droplevels(sce1$sample_id)
sce1$case_control <- droplevels(sce1$case_control)
sce1$ethnicity <- droplevels(sce1$ethnicity)
sce1$sex <- droplevels(sce1$sex)
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id = as.character(colData(sce1)[, c("sample_id")]))
y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))
y <- setupDGElist(y, "WHO_temp_severity", ordered = TRUE, remove = "NA")

res1 <- try(testDGElist(y,
                    formula = as.formula("~ WHO_temp_severity + sex + ethnicity + age_scaled + centre + (1|individual_id)"),
                    individual_id = "individual_id",
                    designMatrix = model.matrix(as.formula("~ WHO_temp_severity + sex + ethnicity + age_scaled + centre"), data = y$samples),
                    ncores = opt$NCPU
                    ), silent = TRUE)
if (class(res1) != "try-error"){
    res1 <- glmmQvals(res1, pi0 = 1)
    results1 <- res1@stats
    results1 <- cbind(results1, res1@optInfo)
} else {
    res1 <- NA
    results1 <- NA
}

# Comparison 2: deg from the interaction between WHO temp severity groups and time
y <- DGEList(counts = counts(smrzd), samples = colData(smrzd))
y <- setupDGElist(y, "WHO_severity_group", ordered = TRUE, remove = "NA")
y$samples$grouped_severity <- ifelse(y$samples$WHO_severity %in% c("mild", "moderate"),
    "mild_moderate", "severe_critical")

res2 <- try(glmm_modified(y, 
                         modelFormula = as.formula("~ splines::bs(days_to_sampling, degree = 2) * grouped_severity + sex + ethnicity + age_scaled + centre + (1|individual_id)"),
                         id = "individual_id",
                         cores = opt$NCPU), silent = TRUE)
res2 <- glmm_qval(res2, pi0 = 1)
results2 <- cbind(res2$stats, res2$optInfo)

# Comparison 3: just wave 1, deg from positive vs negative
sce2 <- sce[, sce$centre == 'NCL']
# Remove samples with less than MIN
nCells <- table(sce2$sample_id)
rmSamples <- names(nCells[nCells < opt$MIN])
sce1 <- sce2[,!sce2$sample_id %in% rmSamples]
# drop unused levels
sce1$individual_id <- droplevels(sce1$individual_id)
sce1$sample_id <- droplevels(sce1$sample_id)
sce1$case_control <- droplevels(sce1$case_control)
sce1$ethnicity <- droplevels(sce1$ethnicity)
sce1$sex <- droplevels(sce1$sex)
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id=as.character(colData(sce1)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y <- setupDGElist(y, 'case_control')
res3 <- try(testDGElist(y, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = opt$NCPU
           ), silent = TRUE)
if (class(res3) != "try-error"){
    res3 <- glmmQvals(res3, pi0 = 1)
    result3 <- degTable(res3, 'case_control', 'POSITIVE')
    result3 <- cbind(result3, res3@optInfo)
} else {
    res3 <- NA
    results3 <- NA
}

# Comparison 4: just patients that were negative in wave 1 but positive in wave 2, deg from recovery vs negative
sce2 <- sce[, sce$individual_id %in% c('C101', 'C108', 'C137', 'C138', 'C140', 
                                       'C141', 'C145', 'C146', 'C147', 'C168',
                                        'C169', 'C170', 'C187', 'C190', 'C33')]
# Remove samples with less than MIN
nCells <- table(sce2$sample_id)
rmSamples <- names(nCells[nCells < opt$MIN])
sce3 <- sce2[,!sce2$sample_id %in% rmSamples]
# remove non-complete data (all have positive), 1 = negative, 3 = recovery
df = table(sce3$individual_id, sce3$case_control)
keep_ids = row.names(df)[which(df[,1] != 0 & df[,3] != 0)]
sce1 <- sce3[, sce3$individual_id %in% keep_ids]
# drop unused levels
sce1$individual_id <- droplevels(sce1$individual_id)
sce1$sample_id <- droplevels(sce1$sample_id)
sce1$case_control <- droplevels(sce1$case_control)
sce1$ethnicity <- droplevels(sce1$ethnicity)
sce1$sex <- droplevels(sce1$sex)
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id=as.character(colData(sce1)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y <- setupDGElist(y, 'case_control', remove = 'POSITIVE')

res4 <- try(testDGElist(y, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = opt$NCPU
           ), silent = TRUE)
if (class(res4) != "try-error"){
    res4 <- glmmQvals(res4, pi0 = 1)
    result4 <- degTable(res4, 'case_control', 'RECOVERY')
    result4 <- cbind(result4, res4@optInfo)
} else {
    res4 <- NA
    results4 <- NA
}

# Comparison 5: just patients that were negative in wave 1 but positive in wave 2, deg from positive vs negative
# remove non-complete data (all have positive), 1 = negative, 2 = positive
keep_ids = row.names(df)[which(df[,1] != 0 & df[,2] != 0)]
sce1 <- sce3[, sce3$individual_id %in% keep_ids]
# drop unused levels
sce1$individual_id <- droplevels(sce1$individual_id)
sce1$sample_id <- droplevels(sce1$sample_id)
sce1$case_control <- droplevels(sce1$case_control)
sce1$ethnicity <- droplevels(sce1$ethnicity)
sce1$sex <- droplevels(sce1$sex)
# Summarize Counts
smrzd <- aggregateAcrossCells(sce1, id=as.character(colData(sce1)[,c("sample_id")]))
y <- DGEList(counts=counts(smrzd), samples=colData(smrzd))
y <- setupDGElist(y, 'case_control', remove = 'RECOVERY')

res5 <- try(testDGElist(y, 
            formula = as.formula("~ case_control + sex + ethnicity + age_scaled + (1|individual_id)"), 
            individual_id = 'individual_id',
            ncores = opt$NCPU
           ), silent = TRUE)
if (class(res5) != "try-error"){
    res5 <- glmmQvals(res5, pi0 = 1)
    result5 <- degTable(res5, 'case_control', 'POSITIVE')
    result5 <- cbind(result5, res5@optInfo)
} else {
    res5 <- NA
    results5 <- NA
}

results <- list(`ordered WHO temp severity` = list(degresults = results1, full = res1), 
                `interaction between WHO temp severity groups and time` = list(degresults = results2, full = res2),
                `wave 1 positive vs negative` = list(degresults = results3, full = res3),
                `wave 2 recovery vs negative` = list(degresults = results4, full = res4),
                `wave 2 positive vs negative` = list(degresults = results5, full = res5)
                )

OUT = paste0(dirname(opt$IN), "/", gsub(".RDS", "_glmmSeq_DEG_results.RDS", basename(opt$IN)))
saveRDS(results, OUT)