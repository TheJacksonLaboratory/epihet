#' @title Construct co-epihet network based on DEH loci.
#' @description Construct co-epihet network for DEH loci or for genes with
#' genome region containing DEH loci using WGCNA
#' and identify modules that are significantly associated with the measured
#' clinical traits for co-epihet DEH loci network,we identify genes with genome
#' region containing DEH loci in each module.
#' @param node.type a character suggest node type in network. Possible values
#' are 'locus','gene'.
#' @param DEH the dataframe containing the chromosome number, loci and strand
#' information of DEH loci generated from diffHet() function.
#' @param compare.matrix The comparison matrix generated from
#' the compMatrix() function.
#' @param value The value to be used to identify DEH loci
#' Possible values are 'pdr', 'epipoly',and 'shannon'.
#' @param group The subtype group to be used to construct network.
#' @param subtype A dataframe containing the subtype information
#' for the samples in the comparison matrix. The row names should
#' be the names of the samples and there should be one column
#' containing the subtype information for each sample.
#' @param datTraits a dataframe containing the clinical traits for all patients
#' from the subtype group.
#' @param annotation.obj a GRanges object containing gene annotation
#' information.
#' @param networktype network type in WGCNA.
#' @param method character string specifying the correlation to be
#' used in WGCNA.
#' @param prefix character string containing the file name base for files
#' containing the consensus
#' topological overlaps in WGCNA.
#' @param mergeCutHeight a numeric, dendrogram cut height for module merging
#' (default: 0.25).
#' @param minModuleSize a numeric, minimum module size for module detection
#' in WGCNA (default: 30).
#' @return a list, if node type is CpG site, it contains the epigenetic
#' heterogeneity matrix for patients
#' module information, gene.list which is a data frame containing genes with
#' genome region containing DEH loci from one module
#' if node type is gene,it contains the epigenetic heterogeneity matrix for
#' patients and module information.
#' @export
epinetwork = function(node.type = "locus", DEH, compare.matrix,
    value = NULL, group, subtype, datTraits = NULL,
    annotation.obj, networktype = "signed", method = "pearson",
    prefix = NULL, mergeCutHeight = 0.25, minModuleSize = 30) {
    # obtain the epigenetic heterogeneity marix on each
    # DEH loci for patients from the group subtype
    group.samples = subtype[which(subtype[, 2] == group), 1]
    compare.matrix = compare.matrix[which(compare.matrix$type == value), ]
    rownames(compare.matrix) = compare.matrix$location
    group.matrix = compare.matrix[, group.samples]
    DEH.loci = paste(DEH[, 1], DEH[, 2], sep = "-")
    value.matrix = group.matrix[DEH.loci, ]
    value.matrix = t(value.matrix)
    DEH.loci = data.frame(loci = DEH.loci, stringsAsFactors = FALSE)
    userset = userobj(DEH)
    o = GenomicRanges::findOverlaps(userset, annotation.obj)
    hit.matrix = as.matrix(o)
    if (node.type == "locus") {
        input.matrix = value.matrix
        nSamples = nrow(value.matrix)
    } else {
        geneid = unique(subjectHits(o))
        if (length(geneid) == 0) {
            print("no gene annotated by DEH loci")
        } else {
            mean.matrix = foreach(i = geneid, .combine = cbind) %do%
                {
                  gene.name = mcols(annotation.obj[i])[, 2]
                  queryhit.id = hit.matrix[which(hit.matrix[, 2] == i), 1]
                  queryhit.loci = mcols(userset[queryhit.id])
                  epivalue = value.matrix[, queryhit.loci$loci]
                  if (nrow(queryhit.loci) == 1) {
                    data.na = data.frame(epivalue)
                  } else {
                    data.na = as.data.frame(apply(epivalue, 1, mean))
                  }
                  colnames(data.na) = gene.name
                  data.na
                }
            mean.matrix = mean.matrix[, colSums(mean.matrix) != 0]
            nSamples = dim(mean.matrix)[1]
            input.matrix = mean.matrix
        }
    }
    # Choose a set of soft-thresholding powers
    powers = c(c(1:10), seq(from = 12, to = 30, by = 2))
    # Call the network topology analysis function
    sft = WGCNA::pickSoftThreshold(input.matrix, powerVector = powers,
        verbose = 5)
    # select beta:
    beta = data.frame(power = sft$fitIndices[, 1],
        R2 = -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
            stringsAsFactors = FALSE)
    if (max(beta$R2) < 0.8) {
        if (nSamples < 20) {
            if (networktype == "signed") {
                softpowers = 20
            } else {
                softpowers = 10
            }
        } else if (nSamples >= 20 & nSamples < 30) {
            if (networktype == "signed") {
                softpowers = 18
            } else {
                softpowers = 9
            }
        } else if (nSamples >= 30 & nSamples < 40) {
            if (networktype == "signed") {
                softpowers = 16
            } else {
                softpowers = 8
            }
        } else if (nSamples >= 40 & nSamples < 60) {
            if (networktype == "signed") {
                softpowers = 14
            } else {
                softpowers = 7
            }
        } else if (nSamples >= 60) {
            if (networktype == "signed") {
                softpowers = 12
            } else {
                softpowers = 6
            }
        }
    } else {
        softpowers = min(beta$power[beta$R2 > 0.8])
    }
    net = WGCNA::blockwiseModules(input.matrix, power = softpowers,
        corType = method, networkType = networktype,
        TOMType = "unsigned", minModuleSize = minModuleSize,
        reassignThreshold = 1e-06, mergeCutHeight = mergeCutHeight,
        numericLabels = TRUE, pamRespectsDendro = FALSE,
        saveTOMs = TRUE, nThreads = 1, saveTOMFileBase = prefix,
        verbose = 3, maxBlockSize = nrow(DEH))
    # To see how many modules were identified and what
    # the module sizes are
    print("the number of modules and the size of each module:")
    print(table(net$colors))
    moduleLabels = net$colors
    # Convert labels to colors for plotting
    moduleColors = WGCNA::labels2colors(net$colors)
    print("save a hierarchical clustering dendrogram for modules in PDF")
    pdf(paste("./", prefix, "dengram_module.pdf", sep = ""))
    f = WGCNA::plotDendroAndColors(net$dendrograms[[1]],
        moduleColors[net$blockGenes[[1]]], "Module colors",
        dendroLabels = FALSE, hang = 0.03, addGuide = TRUE,
        guideHang = 0.05)
    dev.off()
    if (node.type == "locus") {
        module = data.frame(DEHloci = DEH.loci$loci,
            label = moduleLabels, color = moduleColors,
            stringsAsFactors = FALSE)
    } else {
        module = data.frame(gene = colnames(input.matrix),
            label = moduleLabels, color = moduleColors,
            stringsAsFactors = FALSE)
    }
    MEs0 = WGCNA::moduleEigengenes(input.matrix, moduleColors)$eigengenes
    print(paste(prefix, "save network, module and module eigengenes",
        sep = ":"))
    save(net, module, MEs0,
        file = paste("./", prefix,"networkConstruction_auto.rda", sep = ""))
    # cilinical traits
    if (length(datTraits) >= 1) {
        print("Relating modules to external clinical traits:")
        MEs = WGCNA::orderMEs(MEs0)
        moduleTraitCor = WGCNA::cor(MEs, datTraits,use = "p")
        moduleTraitPvalue = WGCNA::corPvalueStudent(moduleTraitCor,nSamples)
        # Will display correlations and their p-values
        textMatrix = paste(signif(moduleTraitCor, 2),
            "\n(", signif(moduleTraitPvalue, 1), ")",sep = "")
        dim(textMatrix) = dim(moduleTraitCor)
        pdf(paste(prefix, "clinicaltrait_heatmap.pdf",sep = ""))
        par(mar = c(6, 8.5, 3, 3))
        # Display the correlation values within a heatmap
        # plot
        h = WGCNA::labeledHeatmap(Matrix = moduleTraitCor,
            xLabels = names(datTraits), yLabels = names(MEs),
            ySymbols = names(MEs), colorLabels = FALSE,
            colors = WGCNA::blueWhiteRed(50), textMatrix = textMatrix,
            setStdMargins = FALSE, cex.text = 0.5,
            zlim = c(-1, 1), main = paste("Module-trait relationships"))
        print(h)
        dev.off()
    } else {
        print("No clinical traits")
    }
    # annotate gene with DEH loci
    if (node.type == "locus") {
        print(" identifying genes involved in modules")
        geneset = data.frame(DEHloci = mcols(userset[hit.matrix[,1]])$loci,
            gene = mcols(annotation.obj[hit.matrix[,2]])$name,
            stringsAsFactors = FALSE)
        gene.module = merge(geneset, module)
        gene.list = unique(gene.module[, -1])
        # plot barplot showing the number of genes in each
        # module
        gene.num = as.data.frame(table(gene.list$color))
        gene.num$Var2 = paste("ME", gene.num$Var1,sep = "")
        g = ggplot(data = gene.num, aes(x = Var2, y = Freq,fill = Var2)) +
            geom_bar(stat = "identity",
            colour = "black",
            position = position_dodge(width = 0.8),width = 0.7) +
            geom_text(aes(label = Freq),position = position_dodge(width = 0.8),
            vjust = -0.5)
        g = g + labs(x = "Module", y = "Number of Gene") +
            theme(legend.title = element_blank()) +
            theme(axis.text.x = element_text(angle = 45,hjust = 1, vjust = 1))+
            theme(legend.position = "none") +
            scale_fill_manual(values = as.character(gene.num$Var1))
        g = g +
            theme(axis.text.x = element_text(size = 10,colour = "black"),
                axis.text.y = element_text(size = 10,colour = "black"),
                axis.title = element_text(size = 10,colour = "black")) +
            theme(axis.line = element_line(size = 0.5,colour = "black"),
                panel.background = element_rect(fill = "white",
                    colour = "darkgrey"),
                plot.background = element_rect(fill = "white",
                    colour = "white"),
                strip.text = element_text(size = 12,colour = "black"),
                strip.background = element_rect(fill = "white",
                    colour = "black"))
        ggsave(paste(prefix, "num_gene_module.pdf",sep = "_"),
            dpi = 300, width = 6, height = 6)
        save(gene.list, file = paste("./", prefix,"gene.module.rda", sep = ""))
    } else {
        gene.num = as.data.frame(table(module$color))
        gene.num$Var2 = paste("ME", gene.num$Var1,sep = "")
        g = ggplot(data = gene.num, aes(x = Var2, y = Freq,fill = Var2)) +
            geom_bar(stat = "identity",colour = "black",
                position = position_dodge(width = 0.8),width = 0.7) +
            geom_text(aes(label = Freq),
                position = position_dodge(width = 0.8),vjust = -0.5)
        g = g + labs(x = "Module", y = "Number of Gene") +
            theme(legend.title = element_blank()) +
            theme(axis.text.x = element_text(angle = 45,
                hjust = 1, vjust = 1)) + theme(legend.position = "none") +
            scale_fill_manual(values = as.character(gene.num$Var1))
        g = g +
            theme(axis.text.x = element_text(size = 10, colour = "black"),
                axis.text.y = element_text(size = 10,colour = "black"),
                axis.title = element_text(size = 10, colour = "black")) +
            theme(axis.line = element_line(size = 0.5, colour = "black"),
                panel.background = element_rect(fill = "white",
                    colour = "darkgrey"),
                plot.background = element_rect(fill = "white",
                    colour = "white"),
                strip.text = element_text(size = 12,colour = "black"),
                strip.background = element_rect(fill = "white",
                    colour = "black"))
        ggsave(paste(prefix, "num_gene_module.pdf",sep = "_"),
            dpi = 300, width = 6, height = 6)
    }
    if (node.type == "locus") {
        result = list(epimatrix = input.matrix, module = module,
            geneofmodule = gene.list)
    } else {
        result = list(epimatrix = input.matrix, module = module)
    }
    return(result)
}
