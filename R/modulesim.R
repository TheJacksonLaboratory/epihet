#' @title  module comparison between two subtypes.
#' @description Compare any two modules from two subytpes based
#' on genes shared by the modules.
#' @param module.subtype1 a data frame generated from the network.construct()
#' function the module information of subtype1,the first column is module
#' nodes,the second column is module label, the third column is module color.
#' @param module.subtype2 a data frame generated from the network.construct()
#' function. The module information of subtype1, the first column is
#' module nodes, the second column is module label, the third column is
#' module color.
#' @param pdf.height An integer representing the height (in inches)
#' of the outputted boxplot pdf file (default: 10)
#' @param pdf.width An integer representing the width (in inches) of
#' the outputted boxplot pdf file (default: 10)
#' @param sve A boolean to save the plot (default: FALSE)
#' @return a matrix containing Jaccard scores
#' @examples
#' data(modulesil,package = "epihet")
#' data(moduledm,package = "epihet")
#' sim.score=epihet::modulesim(module.subtype1=modulesil,
#'                             module.subtype2=moduledm,
#'                             pdf.height = 10,pdf.width = 10,
#'                             sve = TRUE)
#' @export
modulesim = function(module.subtype1, module.subtype2,
    pdf.height = 10, pdf.width = 10, sve = FALSE) {

    jaccard.matrix = jaccard(module.subtype1, module.subtype2)
    jaccard.matrix[jaccard.matrix == 0] = NA
    rownames(jaccard.matrix) = paste("subtype1_ME",
        rownames(jaccard.matrix), sep = "")
    colnames(jaccard.matrix) = paste("subtype2_ME",
        colnames(jaccard.matrix), sep = "")
    # Heatmap
    if (sve) {
        pheatmap::pheatmap(jaccard.matrix, cluster_rows = FALSE,
            cluster_cols = FALSE, height = pdf.height,
            width = pdf.width, filename = "module_compare.pdf")
    } else {
        pheatmap::pheatmap(jaccard.matrix, cluster_rows = FALSE,
            cluster_cols = FALSE)
    }
    return(jaccard.matrix)
}

#' @title  Jaccard score calculation between modules from two subtypes.
#' @description Jaccard score calculation based on the common genes in two
#' modules from two subtypes.
#' @param module.subtype1 a data frame generated from the
#' epinetwork() function. The module information of subtype1,
#' the first column is module nodes, the second column is module label,
#' the third column is module color.
#' @param module.subtype2 a data frame generated from
#' the epinetwork() function. The module information of subtype1,
#' the first column is module nodes, the second column is module label,
#' the third column is module color.
#' @return A matrix containing Jaccard scores.
#' @examples
#' data(modulesil)
#' data(moduledm)
#' jaccard.matrix = jaccard(modulesil, moduledm)
#' @export
jaccard = function(module.subtype1, module.subtype2) {
    moduleid.1 = unique(module.subtype1[, 3])
    moduleid.2 = unique(module.subtype2[, 3])
    row.num = length(moduleid.1)
    col.num = length(moduleid.2)
    jaccard.matrix = matrix(0, nrow = row.num, ncol = col.num)
    rownames(jaccard.matrix) = moduleid.1
    colnames(jaccard.matrix) = moduleid.2
    for (i in 1:row.num) {
        for (j in 1:col.num) {
            moduleid.row = moduleid.1[i]
            moduleid.col = moduleid.2[j]
            gene.row = module.subtype1$gene[which(module.subtype1$color ==
                moduleid.row)]
            gene.col = module.subtype2$gene[which(module.subtype2$color ==
                moduleid.col)]
            gene.share = intersect(gene.row, gene.col)
            gene.union = union(gene.row, gene.col)
            jaccard.matrix[moduleid.row, moduleid.col] =
                length(gene.share)/length(gene.union)
        }
    }
    return(jaccard.matrix)
}
