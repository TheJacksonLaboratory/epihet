#' @title  Modules visualization and
#' @description Visualize the modules identified by epinetwork() function.
#' @param TOM  the topological overlap matrix in WGCNA generated from
#' the network.construct() function
#' @param value.matrix A data frame generated from the epinetwork() function.
#' the row name is patients in one subtype. the column name is the DEH loci
#' the value in the matrix is epigenetic heterogeneity on one DEH loci
#' for one patient
#' @param moduleColors the module assignment generated from the epinetwork()
#' function
#' @param mymodule a character vector containing the module colors
#' you want to visulaize
#' @param cutoff adjacency threshold for including edges in the output.
#' @param prefix a character for outpuf filename
#' @param sve A boolean to save the plot (default: FALSE)
#' @return a list containing all module edge and node information for mymodule
#' @examples
#' correlation.m<-matrix(0,12,12)
#' correlation.m[1,c(2:10)]=c(0.006,0.054,0.079,0.078, 0.011,0.033,0.014,
#' 0.023,0.034)
#' correlation.m[2,c(3:10)]=c(0.026,0.014,0.045,0.037, 0.026,0.011,0.034,
#' 0.012)
#' correlation.m[3,c(4:10)]=c(0.016,0.024,0.039,0.045, 0.009,0.003,0.028)
#' correlation.m[4,c(5:10)]=c(0.039,0.002,0.053,0.066, 0.012,0.039)
#' correlation.m[5,c(6:10)]=c(0.019,0.016,0.047,0.046, 0.013)
#' correlation.m[6,c(7:10)]=c(0.017,0.057,0.029,0.056)
#' correlation.m[7,c(8:10)]=c(0.071,0.018,0.001)
#' correlation.m[8,c(9:10)]=c(0.046,0.014)
#' correlation.m[9,10]=0.054
#' correlation.m[lower.tri(correlation.m)] =
#' t(correlation.m)[lower.tri(correlation.m)]
#'
#' matrix.v<-matrix(0.5,5,12)
#' colnames(matrix.v)=c("NM_052960","NR_138250","NM_015074","NM_183416",
#' "NM_017891","NM_001330306","NM_014917","NM_001312688","NM_001330665",
#' "NM_017766","NM_001079843","NM_001040709")
#' modulecolor<-c(rep(c("yellow","cyan"),c(10,2)))
#' module.topology=epihet::modulevisual(correlation.m,
#'                                      value.matrix=matrix.v,
#'                                      moduleColors=modulecolor,
#'                                      mymodule="yellow",cutoff=0.02,
#'                                      prefix='CEBPA_sil_epipoly',sve = TRUE)
#' @export
modulevisual = function(TOM, value.matrix, moduleColors,
    mymodule, cutoff = 0.02, prefix = NULL, sve = FALSE) {
    # Select modules
    TOM = as.matrix(TOM)
    network = foreach(modules = mymodule) %do% {
        inModule = is.finite(match(moduleColors, modules))
        loci = colnames(value.matrix)
        modloci = loci[inModule]
        modTOM = TOM[inModule, inModule]
        dimnames(modTOM) = list(modloci, modloci)
        # Export the network into edge and node list files
        # Cytoscape can read
        cyt = WGCNA::exportNetworkToCytoscape(modTOM,
            edgeFile = paste(prefix, "CytoscapeInput-edges-",
                paste(modules, collapse = "-"), ".txt", sep = ""),
            nodeFile = paste(prefix, "CytoscapeInput-nodes-",
                paste(modules, collapse = "-"), ".txt", sep = ""),
            weighted = TRUE, threshold = cutoff, nodeNames = modloci,
            nodeAttr = moduleColors[inModule])
        fromNode = as.character(cyt$edgeData$fromNode)
        toNode = as.character(cyt$edgeData$toNode)
        edge = c(rbind(fromNode, toNode))
        net = igraph::graph(edge, directed = FALSE)
        l = igraph::layout_on_sphere(net)
        if (sve) {
            pdf(paste(paste(prefix, modules, "3dsphere",
                sep = "_"), ".pdf", sep = ""))
            plot(net, vertex.color = "gold", vertex.size = 10,
                vertex.frame.color = "gray", vertex.label.color = "black",
                vertex.label.cex = 0.5, vertex.label.dist = 1,
                edge.curved = 0.2, layout = l)
            dev.off()
        } else {
            plot(net, vertex.color = "gold", vertex.size = 10,
                vertex.frame.color = "gray", vertex.label.color = "black",
                vertex.label.cex = 0.5, vertex.label.dist = 1,
                edge.curved = 0.2, layout = l)
        }
        node.topology = data.frame(degree = igraph::degree(net),
            node.centrality = igraph::centr_degree(net)$res,
            node.betweenness = igraph::betweenness(net),
            node.closeness = igraph::closeness(net))

        graph.topology = data.frame(node.num = dim(node.topology)[1],
            edge.num = igraph::gsize(net),
            graph.centrality = igraph::centr_degree(net)$centralization,
            diameter = igraph::diameter(net, unconnected = FALSE))
        data.na = list(edgeData = cyt$edgeData[, c(1:4)],
            nodeData = cyt$nodeData[, c(1, 3)], node.topo = node.topology,
            graph.topo = graph.topology)
        data.na
    }
    names(network) = mymodule
    return(network)
}
