#' background
#'
#' A data frame containing 31995 elements as background used for pathway
#' enrichment analysis:
#' @docType data
#' @keywords datasets
#' @name background
#' @usage data(background)
#' background
#' @format A data frame with 31995 rows and 1 variables:
#' \describe{
#'   \item{gene}{background gene list}
#' }
#' @return
#' A data frame
"background"

#' datTraits
#'
#' clinical traits containing OS,EFS,age
#'
#' @docType data
#' @keywords datasets
#' @name datTraits
#' @usage data(datTraits)
#' datTraits
#' @format A data frame with 6 rows and 3 variables:
#' \describe{
#'   \item{OS}{overall survival time}
#'   \item{EFS}{survival time}
#'   \item{age}{ages of patients}
#' }
#' @return
#' A data frame
"datTraits"

#' DEG
#'
#' differentially expressed genes compared CEBPA-sil vs.normal
#'
#' @docType data
#' @keywords datasets
#' @name DEG
#' @usage data(DEG)
#' DEG
#' @format A data frame with 13685 rows and 1 variables:
#' \describe{
#'   \item{refseq}{refseq ID of DEG}
#' }
#' @return
#' A data frame
"DEG"


#' DEH
#'
#' DEH loci
#'
#' @docType data
#' @keywords datasets
#' @name DEH
#' @usage data(DEH)
#' DEH
#' @format A data frame with 479 rows and 3 variables:
#' \describe{
#'   \item{chromosome}{chromosome of DEH loci}
#'   \item{loci}{location of DEH loci}
#'   \item{strand}{strand of DEH loci}
#' }
#' @return
#' A data frame
"DEH"


#' diffhetmatrix
#'
#' differentially heterogeneity matrix
#'
#' @docType data
#' @keywords datasets
#' @name diffhetmatrix
#' @usage data(diffhetmatrix)
#' diffhetmatrix
#' @format A data frame with 30 rows and 8 variables:
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
#' @return
#' A data frame
"diffhetmatrix"



#' moduledm
#'
#' module information for CEBPA-dm mutation samples
#'
#' @docType data
#' @keywords datasets
#' @name moduledm
#' @usage data(moduledm)
#' moduledm
#' @format A data frame with 140 rows and 3 variables:
#' \describe{
#'   \item{gene}{gene of module CEBPA_dm}
#'   \item{label}{label of module CEBPA_dm}
#'   \item{color}{color of module CEBPA_dm}
#' }
#' @return
#' A data frame
"moduledm"


#' modulesil
#'
#' module information for CEBPA-sil mutation samples
#'
#' @docType data
#' @keywords datasets
#' @name modulesil
#' @usage data(modulesil)
#' modulesil
#' @format A data frame with 501 rows and 3 variables:
#' \describe{
#'   \item{gene}{gene of module CEBPA_sil}
#'   \item{label}{label of module CEBPA_sil}
#'   \item{color}{color of module CEBPA_sil}
#' }
#' @return
#' A data frame
"modulesil"


#' promoter
#'
#' promoter region annotation file
#'
#' @docType data
#' @keywords datasets
#' @name promoter
#' @usage data(promoter)
#' promoter
#' @format A large GRanges with 63344 elements:
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
#'
#' Epipolymorphism for 6 samples on DEH loci
#'
#' @docType data
#' @keywords datasets
#' @name sharedmatrix
#' @usage data(sharedmatrix)
#' sharedmatrix
#' @format A data frame with 479 rows and 8 variables:
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
#' @return
#' A data frame
"sharedmatrix"
