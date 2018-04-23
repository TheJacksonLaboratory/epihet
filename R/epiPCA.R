#' @title Make PCA Plot from Comparison Matrix
#' @description
#' From a user-inputted value, creates a PCA plot from
#' the sample data and colors each point by the subtype
#' information provided.
#' @param compare.matrix The comparison matrix generated from
#' the compMatrix() function
#' @param value The value to be graphed in the PCA plot
#' @param type A dataframe containing the type information
#' for the samples in the comparison matrix. The row names should
#' be the names of the samples and there should be one column
#' containing the type information for each sample.
#' @param points.colors A vector to be used as the color
#' of the individual points for each sample. One color is used per
#' type. the names of vector is the types(default: NULL)
#' @param frames A boolean stating if the frames should be drawn
#' around the points for each subtype cluster. (default: F)
#' @param frames.colors A vector of colors to be used as the color
#' of the frames for each subtype cluster. (default: NULL)
#' @param probability A boolean stating if the frames should be drawn
#' as probability ellipses around the points for each subtype cluster.
#' Both 'probability' and 'frames' must be set to TRUE to have effect.
#' (default: F)
#' @param pdf.height An integer representing the height (in inches) of
#' the outputted PCA plot pdf file (default: 10)
#' @param pdf.width An integer representing the width (in inches) of
#' the outputted PCA plot pdf file (default: 10)
#' @param sve A boolean to save the plot (default: FALSE)
#' @return A PCA plot
#' @examples
#' library(ggfortify)
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
#'epihet::epiPCA(compare.matrix = comp.Matrix, value = 'epipoly',
#'               type = subtype, points.colors = NULL,
#'             frames = FALSE, frames.colors = NULL,
#'             probability = FALSE, pdf.height = 10,
#'             pdf.width = 10, sve = TRUE)
#' @export
epiPCA = function(compare.matrix, value, type, points.colors = NULL,
    frames = FALSE, frames.colors = NULL, probability = FALSE,
    pdf.height = 10, pdf.width = 10, sve = FALSE) {
    # requireNamespace(ggfortify)
    values = c("read", "pdr", "meth", "epipoly", "shannon")
    if (!(value %in% values)) {
        stop("Invalid value '", value, "': Possible values are 'read',
             'pdr', 'meth', 'epipoly', or 'shannon'")
    }
    value.matrix = compare.matrix[compare.matrix$type == value,
        -(length(compare.matrix) - 1)]
    rownames(value.matrix) = value.matrix$location
    value.matrix = value.matrix[, -length(value.matrix)]
    value.matrix = t(value.matrix)
    merge.matrix = merge(type, value.matrix, by = 0, all = TRUE)
    rownames(merge.matrix) = merge.matrix[, 1]
    full.matrix = merge.matrix[, -1]
    sample.matrix = full.matrix[, -1]
    pca.matrix = prcomp(sample.matrix)
    colnames(full.matrix)[1] = "Type"
    curTheme = theme(panel.background = element_rect(fill = "white",
        colour = "darkgrey"), plot.background = element_rect(fill = "white",
        colour = "white"), strip.text = element_text(size = 8,
        colour = "black"), strip.background = element_rect(fill = "white",
        colour = "black"))
    if (frames) {
        if (probability) {
            pca.plot = autoplot(pca.matrix, data = full.matrix,
                size = 2, shape = 19, colour = "Type",
                frame = TRUE, frame.type = "norm") +
                curTheme
        } else {
            pca.plot = autoplot(pca.matrix, data = full.matrix,
                size = 2, shape = 19, colour = "Type", frame = TRUE) +
                curTheme
        }
        if (!is.null(frames.colors)) {
            pca.plot = pca.plot + scale_fill_manual(values = frames.colors)
        }
    } else {
        pca.plot = autoplot(pca.matrix, data = full.matrix,
            size = 2, shape = 19, colour = "Type") +
            curTheme
    }
    if (!is.null(points.colors)) {
        pca.plot = pca.plot + scale_color_manual(values = points.colors)
    }
    title = paste0("PCA Plot of ", value)
    if (sve) {
        ggsave(paste0(value, "_pca.pdf"), plot = pca.plot +
            ggtitle(title), width = pdf.width, height = pdf.height)
    } else {
        pca.plot + ggtitle(title)
    }
}

