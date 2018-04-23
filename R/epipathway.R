#' @title  pathway annotation
#' @description pathway identification significantly enriched by genes
#' in one module.
#' @param gene.list a data frame generated from network.construct() function.
#' The first column is gene entrez ID, the second column is module lable,
#' the third column is module color.
#' @param cutoff Cutoff value of qvalue for pathway enrichment
#' @param prefix a prefix for PDF file name
#' @param pdf.height An integer representing the height (in inches)
#' of the outputted boxplot pdf file (default: 10)
#' @param pdf.width An integer representing the width (in inches) of
#' the outputted boxplot pdf file (default: 10)
#' @return a data frame containing pathways that are significantly enriched
#' by genes from one module
#' @examples
#' genelist=data.frame(ENTREZID=c("2902","2905","3360","286223","59338",
#' "344018","5144","55001","7410","730051","55743","6804","200634","2802",
#' "2260","651","2104","23432","10505","23194","9855","7101",
#' "389136","124857","1829","3164","3754","8614","9469","3217","9578",
#' "10516","10630"),label=rep(18,33),color=rep("lightgreen",33),
#' stringsAsFactors = FALSE)
#' pathway = epihet::epipathway(genelist,cutoff = 0.05,
#'                              prefix="CEBPA_sil",pdf.height = 10,
#'                              pdf.width = 10)
#' @export
epipathway = function(gene.list, cutoff = 0.05, prefix = NA,
    pdf.height = 10, pdf.width = 10) {
    module.id = unique(gene.list[, 3])
    pathway.result = list()
    t = 1
    for (i in module.id) {
        print(i)
        dat = gene.list[which(gene.list[, 3] == i), 1]
        x = ReactomePA::enrichPathway(gene = dat, pvalueCutoff = cutoff,
            readable = TRUE)
        if (length(x) == 0) {
            print(paste("Could not find pathways for",
                paste("ME", i, sep = "")), sep = " ")

        } else if (nrow(x@result) == 0) {
            print(paste("Could not find pathways for",
                paste("ME", i, sep = "")), sep = " ")
        } else {
            dat.na = data.frame(module = paste("ME", i, sep = ""), x@result)
            barplot(x)
            ggsave(paste(prefix, i, "pathway.pdf", sep = "_"),
                height = pdf.height, width = pdf.width)
            pathway.result[[t]] = dat.na
            t = t + 1
        }
    }
    return(pathway.result)
}
