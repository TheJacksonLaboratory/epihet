#' @title  module annotation
#' @description annotate modules using differentially expressed genes.
#' @param DEG a character vector containing up/down regulated genes
#' @param background a charactor vector containing all genes as background
#' in hypergeometric test
#' @param module.gene a data frame containing genes with genome region
#' containing DEH loci from one module, generated from epinetwork() function.
#' The first column is gene entrez ID, the second column is module lable,
#' the third column is module color.
#' @param cutoff Cutoff value of qvalue for gene enrichment
#' @param adjust.method one of 'holm', 'hochberg', 'hommel', 'bonferroni',
#' 'BH', 'BY', 'fdr', 'none'
#' @param prefix a prefix for PDF file name
#' @param pdf.height An integer representing the height (in inches)
#' of the outputted boxplot pdf file (default: 10)
#' @param pdf.width An integer representing the width (in inches) of
#' the outputted boxplot pdf file (default: 10)
#' @param sve A boolean to save the plot (default: FALSE)
#' @return a data frame showing modules that were enriched
#' by DEGs and module size, p value and q value
#' @examples
#' data(DEG,package = "epihet")
#' data(background,package = "epihet")
#' module=data.frame(gene=c("NM_000014","NM_000015","NM_000017","NM_000019",
#' "NM_052960","NR_138250","NM_000037","NM_000038","NM_000039","NM_000044",
#' "NM_000046","NM_015074","NM_183416","NM_004421","NM_001330311",
#' "NM_001145210","NM_000097","NM_000103","NM_000104",
#' "NM_000079","NM_000083","NM_000086","NM_000087","NM_000092","NM_000094",
#' "NM_000095","NM_006474"),label=rep(c(1,2),c(12,15)),
#' color=rep(c("purple","brown"),c(12,15)),
#' stringsAsFactors = FALSE)
#' module.annotation=epihet::moduleanno(DEG$refseq,background$gene,
#'                                      module.gene=module,
#'                                      cutoff=0.05,adjust.method = "fdr",
#'                                      prefix='epipoly',pdf.height = 10,
#'                                      pdf.width = 10, sve = TRUE)
#' @export
moduleanno = function(DEG, background, module.gene,
    cutoff = 0.05, adjust.method = "fdr", prefix = NA,
    pdf.height = 10, pdf.width = 10, sve = FALSE) {
    N = length(background)
    module.id = unique(module.gene[, c(2, 3)])
    gene.num = as.data.frame(table(module.gene[, 3]))
    gene.num$Var1 = as.character(gene.num$Var1)
    Kpai = length(intersect(DEG, background))
    result = foreach(i = module.id$color, .combine = rbind) %do%
        {
            n = gene.num$Freq[gene.num$Var1 == i]
            geneset = module.gene[which(module.gene[, 3] == i), 1]
            kpai = length(intersect(DEG, geneset))
            pvalue = 1 - phyper(kpai - 1, Kpai, N - Kpai, n)
            data = data.frame(moduleColors = i,
                moduleLable = module.id$label[module.id$color == i],
                background.size = N, module.size = n,
                DEG.size = Kpai, share.gene = kpai,
                pvalue = pvalue)
        }
    result$qvalue = p.adjust(result$pvalue, adjust.method)
    result.sig = result[which(result$qvalue <= cutoff), ]
    if (dim(result.sig)[1] == 0) {
        print("No module significantly enriched  by DEGs")
        result.used = result[which(result$share.gene > 0), ]
        if (dim(result.used)[1] == 0) {
            print("No common genes. Quitting.")
        } else {
            color.points = as.character(result.used$moduleColors)
            names(color.points) = color.points
            g = ggplot(data = result.used, aes(x = module.size,
                y = share.gene/module.size, color = moduleColors,
                label = share.gene)) + geom_point(shape = 19, size = 3) +
                geom_text(hjust = 0.5, vjust = -0.5) +
                scale_color_manual(values = color.points)
            g = g + labs(x = "Module size", y = "Percentage of DEGs") +
                theme(legend.title = element_blank()) +
                theme(legend.position = "none")
            g = g + theme(axis.text.x = element_text(size = 10,
                colour = "black"),
                axis.text.y = element_text(size = 10,colour = "black"),
                axis.title = element_text(size = 10,colour = "black"))+
                theme(axis.line = element_line(size = 0.5,colour = "black"),
                    panel.background = element_rect(fill = "white",
                colour = "darkgrey"),
                plot.background=element_rect(fill = "white",colour = "white"),
                strip.text = element_text(size = 12,colour = "black"),
                strip.background = element_rect(fill = "white",
                    colour = "black"))
            if (sve) {
                ggsave(paste0(prefix, "DEG_module.pdf"),plot = g,
                    height = pdf.height, width = pdf.width)
            } else {
                g
            }
        }
    } else {
        color.points = as.character(result.sig$moduleColors)
        names(color.points) = color.points
        g = ggplot(data = result.sig, aes(x = module.size,
            y = -log10(qvalue), color = moduleColors,
            label = share.gene)) + geom_point(shape = 19,
            size = 3) + geom_text(hjust = 0.5, vjust = -0.5) +
            scale_color_manual(values = color.points)
        g = g + labs(x = "Module size", y = "-log10 (adjusted p-value)") +
            theme(legend.title = element_blank()) +
            theme(legend.position = "none")
        g = g + theme(axis.text.x = element_text(size = 10,
            colour = "black"), axis.text.y = element_text(size = 10,
            colour = "black"), axis.title = element_text(size = 10,
            colour = "black")) + theme(axis.line = element_line(size = 0.5,
            colour = "black"), panel.background = element_rect(fill = "white",
            colour = "darkgrey"),
            plot.background = element_rect(fill = "white",colour = "white"),
            strip.text = element_text(size = 12,colour = "black"),
            strip.background = element_rect(fill = "white",colour = "black"))
        if (sve) {
            ggsave(paste0(prefix, "DEG_module.pdf"),
                plot = g, height = pdf.height, width = pdf.width)
        } else {
            g
        }
    }
    return(result)
}
