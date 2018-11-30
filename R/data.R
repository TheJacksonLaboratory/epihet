#' example data
#' background
#' @docType data
#' @keywords datasets
#' @name background
#' @rdname data
#' @usage data(background)
#' background
#' @description background: A data frame containing 31995 elements as background used for pathway
#' enrichment analysis
#' @format background: A data frame with 31995 rows and 1 variables:
#' \describe{
#'   \item{gene}{background gene list}
#' }
#' @return
#' A data frame
"background"

#' datTraits
#' @docType data
#' @keywords datasets
#' @name datTraits
#' @rdname data
#' @usage data(datTraits)
#' datTraits
#' @description datTraits: Clinical traits containing OS,EFS,age
#' @format datTraits: A data frame with 6 rows and 3 variables:
#' \describe{
#'   \item{OS}{overall survival time}
#'   \item{EFS}{survival time}
#'   \item{age}{ages of patients}
#' }
"datTraits"


#' DEG
#' @docType data
#' @keywords datasets
#' @name DEG
#' @rdname data
#' @usage data(DEG)
#' DEG
#' @description DEG: Differentially expressed genes compared CEBPA-sil vs.normal 
#' @format DEG: A data frame with 13685 rows and 1 variables:
#' \describe{
#'   \item{refseq}{refseq ID of DEG}
#' }
"DEG"


#' DEH
#' @docType data
#' @keywords datasets
#' @name DEH
#' @rdname data
#' @usage data(DEH)
#' DEH
#' @description DEH: DEH loci
#' @format DEH: A data frame with 479 rows and 3 variables:
#' \describe{
#'   \item{chromosome}{chromosome of DEH loci}
#'   \item{loci}{location of DEH loci}
#'   \item{strand}{strand of DEH loci}
#' }
"DEH"


#' diffhetmatrix
#' @docType data
#' @keywords datasets
#' @name diffhetmatrix
#' @rdname data
#' @usage data(diffhetmatrix)
#' diffhetmatrix
#' @description diffhetmatrix: A differentially heterogeneity matrix
#' @format diffhetmatrix: A data frame with 30 rows and 8 variables:
#' \describe{
#'   \item{chromosome}{chromosome of loci}
#'   \item{loci}{location of loci}
#'   \item{CEBPA_sil.mean}{mean value of loci in CEBPA-sil patients}
#'   \item{Normal.mean}{mean value of loci in normal patients}
#'   \item{het.dif}{CEBPA_sil.mean-Normal.mean}
#'   \item{p.value}{p value}
#'   \item{p.adjust}{adjusted p value using multiple testing correction}
#'   \item{type}{the type of epigenetic heterogeneity:pdr,epipoly,shannon}
#' }
"diffhetmatrix"


#' moduledm
#' @docType data
#' @keywords datasets
#' @name moduledm
#' @rdname data
#' @usage data(moduledm)
#' moduledm
#' @description moduledm: Module information for CEBPA-dm mutation samples
#' @format moduledm: A data frame with 140 rows and 3 variables:
#' \describe{
#'   \item{gene}{gene of module CEBPA_dm}
#'   \item{label}{label of module CEBPA_dm}
#'   \item{color}{color of module CEBPA_dm}
#' }
"moduledm"


#' modulesil
#' @docType data
#' @keywords datasets
#' @name modulesil
#' @rdname data
#' @usage data(modulesil)
#' modulesil
#' @description modulesil: Module information for CEBPA-sil mutation samples
#' @format modulesil: A data frame with 501 rows and 3 variables:
#' \describe{
#'   \item{gene}{gene of module CEBPA_sil}
#'   \item{label}{label of module CEBPA_sil}
#'   \item{color}{color of module CEBPA_sil}
#' }
"modulesil"


#' promoter
#' @docType data
#' @keywords datasets
#' @name promoter
#' @rdname data
#' @usage data(promoter)
#' promoter
#' @description promoter: The promoter region annotation file
#' @format promoter: A large GRanges with 63344 elements:
#' \describe{
#'   \item{seqnames}{seqnames}
#'   \item{ranges}{ranges}
#'   \item{strand}{strand}
#'   \item{elementMetadata}{elementMetadata}
#'   \item{seqinfo}{seqinfo}
#'   \item{metadata}{metadata}
#' }
#' @return
#' A large GRanges object
"promoter"


#' sharedmatrix
#' @docType data
#' @keywords datasets
#' @name sharedmatrix
#' @rdname data
#' @usage data(sharedmatrix)
#' sharedmatrix
#' @description sharedmatrix: Epipolymorphism values for 6 patients on DEH loci
#' @format sharedmatrix: A data frame with 479 rows and 8 variables:
#' \describe{
#'  \item{D-2238}{patients 1}
#'  \item{D-2668}{patients 2}
#'  \item{D-3314}{patients 3}
#'  \item{D-5360}{patients 4}
#'  \item{D-6947}{patients 5}
#'  \item{D-7076}{patients 6}
#'  \item{type}{type of epigenetic heterogeity}
#'  \item{location}{location of locus}
#' }
"sharedmatrix"
