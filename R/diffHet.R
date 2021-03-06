#' @title Calculate Differential Heterogeneity
#' @description
#' From a user-inputted value and two subtype groups,
#' calculates the mean values for both subtypes at each
#' loci. The heterogeneity difference is calculated and the p-values
#' and adjusted p-values are calculated if the heterogeneity
#' difference is greater than a given cutoff
#' @param compare.matrix The comparison matrix generated from
#' the compMatrix() function
#' @param value The value to be used in calculations. Possible
#' values are 'read', 'pdr', 'meth', 'epipoly', and 'shannon'
#' @param group1 The first subtype group to be compared
#' @param group2 The second subtype group to be compared
#' @param subtype A dataframe containing the subtype information
#' for the samples in the comparison matrix. The row names should
#' be the names of the samples and there should be one column
#' containing the subtype information for each sample.
#' @param het.dif.cutoff A number representing the cutoff
#' for the heterogeneity difference. If the heterogeneity difference is greater
#' than the cutoff value, than the p-value and adjusted p-value are
#' calculated for the loci. If the heterogeneity difference is less than
#' the cutoff value, than the p-value and adjusted p-value are set
#' to NA. (default: 0.20)
#' @param permutations The number of permutations for the
#' EntropyExplorer function. Value must be set to 'shannon'.
#' (default: 1000)
#' @param permutationtest boolean values determining if the permutation test is 
#' applied for DEH loci identification based on customized heterogeneity metrics
#' (default: FALSE) 
#' @param p.adjust.method The method to be used as a parameter in
#' p.adjust() function. Possible methods are 'holm', 'hochberg',
#' 'hommel', 'bonferroni', 'BH', 'BY', 'fdr', and 'none'.(default: 'fdr')
#' @param cores The number of cores to be used for parallel execution.
#' Not available for 'shannon' values. (default: 5)
#' @return A dataframe containing chromosome number, loci, mean of
#' group1, mean of group2, heterogeneity difference, and the p-value and
#' adjusted p-value for the loci with a heterogeneity difference greater
#' than the cutoff
#' @export
diffHet <- function(compare.matrix, value, group1, group2,
    subtype, het.dif.cutoff = 0.2, permutations = 1000, permutationtest=FALSE,
    p.adjust.method = "fdr", cores = 5) {
    values <- c("read", "pdr", "meth", "epipoly", "shannon","myValues")#new add
    if (!(value %in% values)) {
        stop("Invalid value '", value, "': Possible values are 'read',
             'pdr', 'meth', 'epipoly', 'myValues' or 'shannon'")#new add
    }
    type.col <- which(colnames(compare.matrix) == 'type')
    subtype$tmp = seq_len(nrow(subtype))
    group1.samples <- subtype[which(subtype[, 2] == group1), 1]
    group2.samples <- subtype[which(subtype[, 2] == group2), 1]
    value.matrix <- compare.matrix[compare.matrix$type == value, -type.col]
    rownames(value.matrix) <- value.matrix$location
    val1.r <- value.matrix[, group1.samples]
    val2.r <- value.matrix[, group2.samples]
    mean1.v <- rowMeans(val1.r)
    mean2.v <- rowMeans(val2.r)
    fc.v <- mean1.v - mean2.v
    fc.v.filter <- fc.v[which(abs(fc.v) <= het.dif.cutoff)]
    loci.filter <- names(fc.v.filter)
    if (length(loci.filter)){
        nopvals <- data.frame(chromosome = splitn(loci.filter, "-", 1),
            loci = splitn(loci.filter, "-", 2),
            group1.mean = mean1.v[loci.filter],
            group2.mean = mean2.v[loci.filter],
            het.dif = fc.v.filter, p.value = NA, p.adjust = NA,
            type = value, stringsAsFactors = FALSE)
    }
    val1 <- val1.r[which(abs(fc.v) > het.dif.cutoff),]
    val2 <- val2.r[which(abs(fc.v) > het.dif.cutoff),]
    loci <- rownames(val1)
    if (value == "shannon"){
        Pertest <- TRUE
        print("using permutation test to identify DEH loci")
    } else if (value == "pdr" | value == "epipoly"){
        Pertest <- FALSE
        print("using t-test to identify DEH loci")
    } else {
        Pertest <- permutationtest
        if (Pertest){
            print("using permutation test to identify DEH loci")
        } else {
            print("using t-test to identify DEH loci")
        }
    }
    if (Pertest) {
        #seed = sample(1:1e+06, 1)
        #set.seed(seed)
        nums <- ceiling(nrow(val1)/cores)
        start <- 1
        end <- nums
        lst1 <- list()
        lst2 <- list()
        for (i in seq_len(cores)) {
            lst1[[i]] <- val1[start:end, ]
            lst2[[i]] <- val2[start:end, ]
            start <- end + 1
            end <- end + nums
            if (end > nrow(val1)) {
                end <- nrow(val1)
            }
        }
        i <- NULL
        j <- NULL
        doParallel::registerDoParallel(cores=cores)
        pval <- foreach(i = lst1, j = lst2, .combine = rbind) %dopar%
            {
                suppressWarnings(EntropyExplorer::EntropyExplorer(i,
                    j, "dse", "pr", shift = c("auto","auto"),
                    padjustmethod = p.adjust.method,
                    nperm = permutations))
            }
        if (nrow(pval) == 0) {
            print('===Note:permutation test were not excuted because of too few values!===')
            return(NA)
        } else {
            p.vals <- foreach(i = loci, .combine = rbind) %dopar%
                {
                    mean1 <- mean(as.numeric(val1[i, ]))
                    mean2 <- mean(as.numeric(val2[i, ]))
                    split <- strsplit(i, "-")
                    curChr <- split[[1]][1]
                    curLoci <- split[[1]][2]
                    fc <- mean1 - mean2
                    cur.pval <- pval[i, 1]
                    data.frame(chromosome = curChr, loci = curLoci,
                               group1.mean = mean1, group2.mean = mean2,
                               het.dif = fc, p.value = cur.pval,
                               stringsAsFactors = FALSE)
                }
        }
        
    } else {
        doParallel::registerDoParallel(cores=cores)
        p.vals <- foreach(i = loci, .combine = rbind) %dopar%
            {
                test <- t.test(val1[i, ], val2[i, ])
                mean1 <- test$estimate[1]
                mean2 <- test$estimate[2]
                split <- strsplit(i, "-")
                curChr <- split[[1]][1]
                curLoci <- split[[1]][2]
                pval <- test$p.value
                fc <- mean1 - mean2
                data.frame(chromosome = curChr, loci = curLoci,
                    group1.mean = mean1, group2.mean = mean2,
                    het.dif = fc, p.value = pval, stringsAsFactors = FALSE)
            }
    }
    print("Finish p value calculation")
    p.vals$p.adjust = p.adjust(p.vals$p.value, p.adjust.method)
    p.vals$type = value
    if (length(loci.filter)){
        p.vals <- rbind(p.vals, nopvals)
    }
    rownames(p.vals) <- seq_len(nrow(p.vals))
    column3 <- paste0(group1, ".mean")
    column4 <- paste0(group2, ".mean")
    colnames(p.vals)[3:4] <- c(column3, column4)
    p.vals
}
