% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/epiMA.R
\name{epiMA}
\alias{epiMA}
\title{Make MA Plot}
\usage{
epiMA(pval.matrix, padjust.cutoff = 0.05, pch = ".", sve = FALSE, pointsize = 1.5)
}
\arguments{
\item{pval.matrix}{The data frame returned from the
diffHet() function that contains means, p-values,
adjusted p-values, and heterogeneity difference}

\item{padjust.cutoff}{The adjusted p-value cutoff
to confirm a significant value. (default: 0.05)}

\item{pch}{The plotting character to be used in the
MA plot (default: '.')}

\item{sve}{A boolean to save the plot (default: FALSE)}
  
\item{pointsize}{A numeric value to adjust point size (default: 1.5)}  
}
\value{
A figure
}
\description{
Creates an MA plot from the differential
heterogeneity data calculated from the diffHet()
function. For each loci, graphs the average of both
group means on the x-axis and the heterogeneity difference
on the y-axis. Graphs coordinates with significant adjusted
p-values in red.
}
\examples{
diff.het.matrix<-data.frame(chromosome=c(rep("1",10)),
loci=paste("loci",seq_len(10),sep="-"),subtype.mean=c(0.21,0.23,0.37,0.26,
0.29,0.31,0.29,0.13,0.12,0.093),Normal.mean=c(0.01,0.01,0.01,0.02,
0.02,0.01,0.01,0.79,0.73,0.79),het.dif=c(0.20,0.220,0.360,0.240,0.270,
0.300,0.280,-0.660,-0.610,-0.697),p.value=c(3.08e-03,1.43e-02,9.27e-03,
3.45e-02,2.99e-02,3.68e-02, 4.60e-02, 5.65e-10, 9.18e-10,
9.98e-11),p.adjust=c(8.84e-03,2.76e-02, 2.04e-02, 5.01e-02,
4.56e-02, 5.24e-02, 6.08e-02, 3.74e-08, 5.22e-08,
1.06e-08),type=rep("pdr",10))

epiMA(pval.matrix = diff.het.matrix, padjust.cutoff = 0.05,
pch = ".", sve = TRUE, pointsize = 1.5)
}
