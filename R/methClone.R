#' @title A wrapper function for methclone package
#'
#' @description
#' Runs methclone package on the provided two BAM files
#' file
#'
#' @param file1 File name for the first input of methclone package
#' @param file2 File name for the second input of methclone package
#' @param outfile Output file name where methclone writes its output.
#' This is a gzipped file
#' @param sampleid Name of the sample provided to methClone package
#' @param distance Maximum distance between first and forth methylation sites
#' to be considered in a single pattern
#' @param coverage Minimum coverage at a signle methylation site to
#' be considered
#' @param methdiff Threshold for difference in methylation between the two
#' different samples
#' @return Does not return any output, it writes methcone output into
#' the provided outfile name
#' @export
#'
#'


methClone = function(file1, file2, outfile, sampleid, 
    distance = 72, coverage = 60, methdiff = 0) {
    
    MethClone_two_samples(file1, file2, outfile, sampleid, 
        distance, coverage, methdiff)
}
