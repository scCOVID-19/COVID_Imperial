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
suppressPackageStartupMessages(library(BiocParallel))

setupDGElist <- function(dgelist, comparison, ordered = FALSE, remove = NULL, dropLevels = c("individual_id",
    "sample_id", "case_control", "ethnicity", "sex")) {
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
    if (!is.null(dropLevels)) {
        for (dL in dropLevels) {
            dgelist$samples[, dL] <- droplevels(dgelist$samples[, dL])
        }
    }
    return(dgelist)
}

testDGElist <- function(dgelist, formula, individual_id, modified = FALSE, optimizer = "bobyqa",
    designMatrix = NULL, ncores = NULL, BPPARAM = SerialParam(progress = TRUE), tpm = FALSE,
    ...) {
    if (is.null(ncores)) {
        NCORES <- parallel::detectCores() - 1
    } else {
        NCORES <- ncores
    }
    if (tpm) {
        if (modified) {
            stop("tpm=TRUE only can be run in modified=FALSE mode.")
        }
        disp <- suppressMessages(setNames(edgeR::estimateDisp(dgelist)$tagwise.dispersion,
            rownames(dgelist)))  # Estimate Dispersion
        dgelist <- calcNormFactors(dgelist)  # just to get the norm.factors
        sizeFactors <- dgelist$samples$norm.factors
    } else {
        # Norm
        dgelist <- calcNormFactors(dgelist)
        sizeFactors <- dgelist$samples$lib.size * dgelist$samples$norm.factors
        disp <- suppressMessages(setNames(edgeR::estimateDisp(dgelist)$tagwise.dispersion,
            rownames(dgelist)))  # Estimate Dispersion
    }
    if (modified) {
        results <- suppressMessages(glmm_modified(dgelist, modelFormula = formula,
            dispersion = disp, sizeFactors = sizeFactors, id = individual_id, control = glmerControl(optimizer = optimizer,
                optCtrl = list(maxfun = 2e+05)), BPPARAM = BPPARAM, ...))
        results <- glmm_qval(results, pi0 = 1)
    } else {
        results <- glmmSeq(formula, id = individual_id, countdata = dgelist$counts,
            metadata = dgelist$samples, dispersion = disp, sizeFactors = sizeFactors,
            removeDuplicatedMeasures = FALSE, removeSingles = FALSE, control = glmerControl(optimizer = optimizer,
                optCtrl = list(maxfun = 2e+05), check.conv.singular = "ignore"),
            progress = TRUE, cores = NCORES)
        results <- glmmQvals(results, pi0 = 1)
    }
    return(results)
}

degTable <- function(results, contrast, group, remove_issues = TRUE, reverse = FALSE,
    modified = FALSE) {
    if (modified) {
        contrasts <- paste0(c(gsub("):", ")1:", contrast), gsub("):", ")2:", contrast)),
            group)
        tmp <- data.frame(results$stats[, c(contrasts, paste0("P_", contrast), paste0("q_",
            contrast))], check.names = FALSE)
        colnames(tmp) <- c("beta_linear", "beta_quadratic", "pval", "qval")
        tmp <- cbind(tmp, results$optInfo)
    } else {
        tmp <- data.frame(results@stats[, c(paste0(contrast, group), paste0("P_",
            contrast), paste0("q_", contrast))], check.names = FALSE)
        colnames(tmp) <- c(beta, "pval", "qval")
        # tmp$LFC <- LFC
        tmp <- cbind(tmp, results@optInfo)
    }
    if (remove_issues) {
        tmp <- tmp[which(tmp$Conv == 0), ]
    }
    if (modified) {
        tmp <- tmp[order(-tmp$beta_quadratic, tmp$qval), ]
    } else {
        tmp <- tmp[order(-tmp$beta, tmp$qval), ]
    }
    return(tmp)
}

degTable_modified <- function(results, contrast, group, remove_issues = TRUE, reverse = FALSE) {
    tmp <- data.frame(results$stats[, c(paste0(contrast, group), paste0("P_", contrast),
        paste0("q_", contrast))], check.names = FALSE)
    colnames(tmp) <- c("beta", "pval", "qval")
    # tmp$LFC <- LFC
    tmp <- cbind(tmp, results$optInfo)
    if (remove_issues) {
        tmp <- tmp[which(tmp$Conv == 0), ]
    }
    tmp <- tmp[order(-tmp$beta, tmp$qval), ]
    return(tmp)
}

degTable_simple <- function(results, contrast, group, remove_issues = TRUE, reverse = FALSE) {
    modelData <- results@modelData
    outLabels <- apply(modelData, 1, function(x) paste(x, collapse = "_"))
    modelData$y <- paste0("y_", outLabels)
    cols1 = grep(group, modelData$y, value = TRUE)
    cols2 = grep(group, modelData$y, value = TRUE, invert = TRUE)
    if (reverse) {
        LFC <- log2(results@predict[, cols2] + 1) - log2(results@predict[, cols1] +
            1)
    } else {
        LFC <- log2(results@predict[, cols1] + 1) - log2(results@predict[, cols2] +
            1)
    }
    tmp <- data.frame(results@stats[, c(paste0(contrast, group), paste0("P_", contrast),
        paste0("q_", contrast))], check.names = FALSE)
    colnames(tmp) <- c("beta", "pval", "qval")
    tmp$LFC <- LFC
    tmp <- cbind(tmp, results@optInfo)
    if (remove_issues) {
        tmp <- tmp[which(tmp$Conv == 0), ]
    }
    tmp <- tmp[order(-tmp$beta, tmp$qval), ]
    return(tmp)
}
                           
glmm_modified <- function(dgeList, modelFormula, id, dispersion, sizeFactors = NULL,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e+05)),
    BPPARAM = SerialParam(progress = TRUE), verbose = TRUE, ...) {
    metadata <- dgeList$samples
    countdata <- dgeList$counts
    # Norm
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
    resultList <- bplapply(fullList, function(geneList) {
        try(run_glmm(data = metadata, geneList = geneList, fullFormula = fullFormula,
            control = control, offset = offset, ...), silent = TRUE)
    }, BPPARAM = BPPARAM)
    # Print timing if verbose
    end <- Sys.time()
    if (verbose)
        print(end - start)
    # Output
    names(resultList) <- rownames(countdata)
    # any failed?
    if (any(lapply(resultList, length) < 4)) {
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
            if (length(failed) > 0) {
                cat(paste0("Errors in ", sum(!noErr, length(failed)), " gene(s):",
                  paste0(c(names(noErr)[!noErr], failed), collapse = ", ")))
            } else {
                cat(paste0("Errors in ", sum(!noErr), " gene(s):", paste0(names(noErr)[!noErr],
                  collapse = ", ")))
            }
        }
        outputErrors <- vapply(resultList[!noErr], function(x) {
            x$tryErrors
        }, FUN.VALUE = c("test"))
        if (length(failed > 0)) {
            outputErrors <- c(outputErrors, errormsg)
        }
    } else {
        outputErrors <- c("No errors")
    }
    optInfo <- t(vapply(resultList[noErr], function(x) {
        setNames(x$optinfo, c("Singular", "Conv"))
    }, FUN.VALUE = c(1, 1)))
    nCheat <- resultList[noErr][[1]]$stats
    s <- t(vapply(resultList[noErr], function(x) {
        x$stats
    }, FUN.VALUE = rep(1, length(nCheat))))
    return(list(stats = s, fit = resultList, optInfo = optInfo, errors = outputErrors))
}

run_glmm <- function(data, geneList, fullFormula, control, offset, ...) {
    data[, "count"] <- as.numeric(geneList$y)
    fit <- lme4::glmer(fullFormula, data = data, control = control, offset = offset,
        family = MASS::negative.binomial(theta = 1/geneList$dispersion), ...)
    if (class(fit) != "try-error") {
        stats <- setNames(c(geneList$dispersion, AIC(fit), as.numeric(logLik(fit))),
            c("Dispersion", "AIC", "logLik"))
        fixedEffects <- lme4::fixef(fit)
        wald <- try(car::Anova(fit), silent = TRUE)
        if (class(wald) != "try-error") {
            waldtest <- setNames(c(wald[, "Chisq"], wald[, "Pr(>Chisq)"]), c(paste0("Chisq_",
                rownames(wald)), paste0("P_", rownames(wald))))
            singular <- as.numeric(isSingular(fit))
            conv <- length(slot(fit, "optinfo")$conv$lme4$messages)
            return(list(stats = c(stats, fixedEffects, waldtest), fit = fit, optinfo = c(singular,
                conv), tryErrors = ""))
        } else {
            return(list(stats = NA, fit = NA, optinfo = NA, tryErrors = wald[1]))
        }
    } else {
        return(list(stats = NA, fit = NA, optinfo = NA, tryErrors = fit[1]))
    }
}

glmm_qval <- function(result, cutoff = 0.05, pi0 = NULL, verbose = TRUE) {
    resultStats <- data.frame(result$stats, check.names = FALSE)
    for (cn in colnames(resultStats)[grep("P_", colnames(resultStats))]) {
        q_cn <- gsub("P_", "q_", cn)
        resultStats[, q_cn] <- NA
        resultStats[!is.na(resultStats[, cn]), q_cn] <- qvalue::qvalue(resultStats[!is.na(resultStats[,
            cn]), cn], pi0 = pi0)$qvalues
        if (verbose) {
            cat(paste0("\n", q_cn, "\n"))
            cat(paste(rep("-", nchar(q_cn)), collapse = ""))
            print(table(ifelse(resultStats[, q_cn] < cutoff, "Significant", "Not Significant")))
        }
    }
    result$stats <- as.matrix(resultStats)
    return(result)
}
