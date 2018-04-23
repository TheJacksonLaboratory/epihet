#' @title Make TSNE Plot from Comparison Matrix
#' @description
#' From a user-inputted value, creates a TSNE plot from the
#' sample data and colors each point by the subtype information
#' provided.
#' @param compare.matrix The comparison matrix generated from
#' the compMatrix() function
#' @param value The value to be graphed in the PCA plot
#' @param type A dataframe containing the type information
#' for the samples in the comparison matrix. The row names should
#' be the names of the samples and there should be one column
#' containing the type information for each sample.
#' @param points.colors A vector of colors to be used as the color
#' of the individual points for each sample. One color is used per
#' subtype. (default: NULL)
#' @param theta A decimal representing the theta parameter for
#' the Rtsne() function. Represents the speed/accuracy trade-off
#' (0.0 is exact TSNE) (default: 0.5)
#' @param curTheme the theme of ggplot2 to control control the appearance of
#' all non-data components of the plot
#' @param perplexity An integer representing the perplexity
#' parameter for the Rtsne() function (default: 30)
#' @param max_iter An integer representing the max_iter parameter
#' for the Rtsne() function. Represents the number of iterations
#' (default: 1000)
#' @param pdf.height An integer representing the height (in inches)
#' of the outputted TSNE plot pdf file (default: 10)
#' @param pdf.width An integer representing the width (in inches) of
#' the outputted TSNE plot pdf file (default: 10)
#' @param sve A boolean to save the plot (default: FALSE)
#' @importFrom ggplot2 aes geom_point ggplot ggsave ggtitle
#' labs scale_color_manual xlab ylab
#' @return A T-SNE plot
#' @examples
#' comp.Matrix=data.frame(
#' p1=c(0.6,0.3,0.5,0.5,0.5,0.6,0.45,0.57,0.45,0.63,0.58,0.67,0.5,0.42,0.67),
#' p2=c(0.62,0.63,0.55,0.75,0.84,0.58,1,0.33,1,0.97,0.57,0.68,0.73,0.72,0.82),
#' p3=c(0.72,0.53,0.62,0.69,0.37,0.85,1,0.63,0.87,0.87,0.82,0.81,0.79,
#' 0.62,0.68),
#' N1=c(0.15,0.24,0.15,0.26,0.34,0.32,0.23,0.14,0.26,0.32,0.12,0.16,0.31,
#' 0.24,0.32),
#' N2=c(0.32,0.26,0.16,0.36,0.25,0.37,0.12,0.16,0.41,0.47,0.13,0.52,0.42,
#' 0.41,0.23),
#' N3=c(0.21,0.16,0.32,0.16,0.36,0.27,0.24,0.26,0.11,0.27,0.39,0.5,0.4,
#' 0.31,0.33),
#' type=rep(c("pdr","epipoly","shannon"),c(5,5,5)),
#' location=rep(c("chr22-327:350:361:364","chr22-755:761:771:773",
#' "chr22-761:771:773:781","chr22-821:837:844:849","chr22-838:845:850:858"),
#' 3),stringsAsFactors =FALSE )
#'
#' subtype = data.frame(Type= c(rep('CEBPA_sil', 3), rep('Normal', 3)),
#' row.names = colnames(comp.Matrix)[1:6],stringsAsFactors = FALSE)
#'
#' epihet::epiTSNE(compare.matrix = comp.Matrix, value = 'epipoly',
#' type = subtype, points.colors = NULL, theta = 0.5,
#' perplexity = 1, max_iter = 1000, pdf.height = 10,
#' pdf.width = 10, sve = TRUE)
#' @export
epiTSNE = function(compare.matrix, value, type, points.colors = NULL,
    theta = 0.5, curTheme = NULL, perplexity = 5, max_iter = 1000,
    pdf.height = 10, pdf.width = 10, sve = FALSE) {
    values = c("read", "pdr", "meth", "epipoly", "shannon")
    if (!(value %in% values)) {
        stop("Invalid value '", value, "': Possible values are 'read',
             'pdr', 'meth', 'epipoly', or 'shannon'")
    }
    value.matrix = compare.matrix[compare.matrix$type == value,
        -(length(compare.matrix) - 1)]
    rownames(value.matrix) = value.matrix$location
    value.matrix = value.matrix[, -length(value.matrix)]
    if (ncol(value.matrix) - 1 < 3 * perplexity) {
        perplexity = floor((ncol(value.matrix) - 1)/3)
        if (perplexity >= 1) {
            print(paste0("Perplexity too large. Setting to ", perplexity))
        } else {
            print("Not enough samples for t-SNE.Quitting.")
            return()
        }
    }
    value.matrix = t(value.matrix)
    merge.matrix = merge(type, value.matrix, by = 0, all = TRUE)
    rownames(merge.matrix) = merge.matrix[, 1]
    full.matrix = merge.matrix[, -1]
    colnames(full.matrix)[1] = "Type"
    sample.matrix = as.matrix(full.matrix[, -1])
    if (is.null(curTheme)) {
        curTheme = theme(legend.position = "right",
            legend.title = element_text(size = 10),
            legend.text = element_text(size = 10),
            axis.text.y = element_text(size = 10),
            axis.title = element_text(size = 10),
            plot.title = element_text(size = 10),
            axis.text.x = element_text(size = 10, angle = -90),
            strip.text.y = element_text(size = 10))
    }
    set.seed(42)
    tsne = Rtsne::Rtsne(sample.matrix, theta = theta,
        verbose = TRUE, perplexity = perplexity, max_iter = max_iter)
    tsne2 = tsne$Y
    title = paste0("t-SNE Plot for ", value, "(perplexity=", perplexity, ")")
    pcaRes = data.frame(Sample = rownames(sample.matrix),
        Type = factor(full.matrix[, 1]), tsne2)
    tsne.plot = ggplot(pcaRes, aes(x = X1, y = X2,
        color = Type)) + geom_point() + ggtitle(title) +
        xlab("t-SNE1") + ylab("t-SNE2")
    tsne.plot = tsne.plot + theme_linedraw() + curTheme
    if (!is.null(points.colors)) {
        tsne.plot = tsne.plot + scale_color_manual(values = points.colors)
    }
    if (sve) {
        ggsave(paste0(value, "_tsne.pdf"), plot = tsne.plot,
            height = pdf.height, width = pdf.width)
    } else {
        tsne.plot
    }
}
