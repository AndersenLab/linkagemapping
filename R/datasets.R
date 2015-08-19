#' Cross object for N2xCB4856 RIAILs
#'
#' In this object, N2 genotypes are encoded as 1 and CB4856 genotypes are
#' encoded as 2. This cross object contains no heterozyous loci. When the
#' extract genotype function is run on this cross object, N2 genotypes are
#' converted to -1 and CB4856 genotypes to 1. For mappings with this cross
#' object, negative effect sizes indicate that N2 had the greater phenotype
#' value while postive effect sizes indicate that CB4856 had a greater phenotype
#' value.
#'
#' @name N2xCB4856cross
#' @format A cross object made with the \code{qtl} package, with genotype and
#' phenotype subobjects.
NULL

#' Cross object for N2xLSJ2 RIAILs
#'
#' In this object, N2 genotypes are encoded as 1 and LSJ2 genotypes are
#' encoded as 3. This cross object had heterozyous loci removed (set to NA). If
#' heterozygous loci are introduced back into the data, they should be encoded
#' as 2. When the extract genotype function is run on this cross object, N2
#' genotypes are converted to -1 and LSJ2 genotypes to 1. If heterozygous loci
#' are introduced back into the data, they will be extracted as the value 0. For
#' mappings with this cross object, negative effect sizes indicate that N2 had
#' the greater phenotype value while postive effect sizes indicate that LSJ2 had
#' a greater phenotype value.
#'
#' @name N2xLSJ2cross
#' @format A cross object made with the \code{qtl} package, with genotype and
#' phenotype subobjects.
NULL

#' Cross object for AF16xHK104 C. brggsae RIAILs
#'
#' In this object, AF16 genotypes are encoded as 1 and HK104 genotypes are
#' encoded as 2. This cross object contains no heterozyous loci. When the
#' extract genotype function is run on this cross object, AF16 genotypes are
#' converted to -1 and HK104 genotypes to 1. For mappings with this cross
#' object, negative effect sizes indicate that Af16 had the greater phenotype
#' value while postive effect sizes indicate that HK104 had a greater phenotype
#' value.
#'
#' @name AF16xHK104cross
#' @format A cross object made with the \code{qtl} package, with genotype and
#' phenotype subobjects.
NULL

#' Marker/position lookup table for N2xCB4856 RIAILs
#'
#' This data frame contains a lookup table for the translation of marker
#' positions from genetic to physical positions. The physical positions in this
#' table come from C. elegans genome build WS244.
#'
#' @name N2xCB4856markers
#' @format A data frame with 1460 rows and 4 variables:
#' \describe{
#'   \item{marker}{the name of the marker used for genotyping}
#'   \item{chr.num}{the chromosome number in arabic numerals}
#'   \item{chr.roman}{the chromosome number in roman numerals}
#'   \item{chr.roman}{the position of the marker in genome build WS244}
#' }
NULL

#' Marker/position lookup table for N2xLSJ2 RIAILs
#'
#' This data frame contains a lookup table for the translation of marker
#' positions from genetic to physical positions. The physical positions in this
#' table come from C. elegans genome build WS195.
#'
#' @name N2xLSJ2markers
#' @format A data frame with 175 rows and 4 variables:
#' \describe{
#'   \item{marker}{the name of the marker used for genotyping}
#'   \item{chr.num}{the chromosome number in arabic numerals}
#'   \item{chr.roman}{the chromosome number in roman numerals}
#'   \item{chr.roman}{the position of the marker in genome build WS195}
#' }
NULL

#' Marker/position lookup table for AF16xLSJ2 RIAILs
#'
#' This data frame contains a lookup table for the translation of marker
#' positions from genetic to physical positions. The physical positions in this
#' table come from C. briggsae genome build cb4.
#'
#' @name AF16xLSJ2markers
#' @format A data frame with 1031 rows and 4 variables:
#' \describe{
#'   \item{marker}{the name of the marker used for genotyping}
#'   \item{chr.num}{the chromosome number in arabic numerals}
#'   \item{chr.roman}{the chromosome number in roman numerals}
#'   \item{chr.roman}{the position of the marker in genome build cb4}
#' }
NULL