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
#' table come from C. elegans genome build WS245.
#'
#' @name N2xLSJ2markers
#' @format A data frame with 175 rows and 4 variables:
#' \describe{
#'   \item{marker}{the name of the marker used for genotyping}
#'   \item{chr.num}{the chromosome number in arabic numerals}
#'   \item{chr.roman}{the chromosome number in roman numerals}
#'   \item{chr.roman}{the position of the marker in genome build WS245}
#' }
NULL

#' Marker/position lookup table for AF16xHK104 RIAILs
#'
#' This data frame contains a lookup table for the translation of marker
#' positions from genetic to physical positions. The physical positions in this
#' table come from C. briggsae genome build cb4.
#'
#' @name AF16xHK104markers
#' @format A data frame with 1031 rows and 4 variables:
#' \describe{
#'   \item{marker}{the name of the marker used for genotyping}
#'   \item{chr.num}{the chromosome number in arabic numerals}
#'   \item{chr.roman}{the chromosome number in roman numerals}
#'   \item{chr.roman}{the position of the marker in genome build cb4}
#' }
NULL

# Dataframe of all N2 fosmids in the Andersen Lab
#
# This data frame contains genomic position and laboratory stock
# information for all N2 fosmids.
#
# @name AllN2fosmids
# @format A data frame with 12468 rows and 5 variables:
# \describe{
#   \item{chr}{the chromosome (in roman numerals) on which the fosmid is found}
#   \item{clone}{the name of the fosmid}
#   \item{start}{the genetic position at which the fosmid begins}
#   \item{end}{the genetic position at which the fosmid ends}
#   \item{newpos}{the current location of the fosmid in the Andersen Lab}
# }
NULL

# Dataframe of all CB4856 fosmids in the Andersen Lab
#
# This data frame contains genomic position and laboratory stock
# information for all CB4856 fosmids.
#
# @name AllCBfosmids
# @format A data frame with 11453 rows and 5 variables:
# \describe{
#   \item{chr}{the chromosome (in roman numerals) on which the fosmid is found}
#   \item{clone}{the name of the fosmid}
#   \item{start}{the genetic position at which the fosmid begins}
#   \item{end}{the genetic position at which the fosmid ends}
#   \item{newpos}{the current location of the fosmid in the Andersen Lab}
# }
NULL

# Dataframe of all eQTL peak location information from Matt Rockman
#
# This dataframe contains all position information for eQTL peaks
#
# @name eQTLpeaks
# @format A data frame with 2447 rows and 25 variables
NULL

# Dataframe of GO terms, descriptions, and other information for
# all genes in the eQTL dataset from Matt Rockman
#
# This data frame contains descriptive information for all genes
# included in the eQTL dataset
#
# @name probe_info
# @format A data frame with 43603 rows and 12 variables
NULL

# Genotypes of RIAILs at each of the 1454 markers
#
# This data frame contains the genotype information for each strain
# of the RIAIL panel at each of the 1454 markers. "AA" denotes
# the N2 genotype and "AB" denotes the CB4856 genotype.
#
# @name RIAILgenotypes
# @format A data frame with 1454 rows and 601 variables
NULL

# Marker conversions for each of the 1454 markers used in RIAILgenotypes
#
# This data frame contains marker name and position information for
# each of the 1454 markers used in the RIAILgenotypes data frame
#
# @name RIAILmarkerconversion
# @format A data frame with 1454 rows and 2 variables
NULL

# Dataframe of indels between N2 and CB4856
#
# This data frame contains position information for all indels between
# N2 and CB4856 greater than 25 bp. Developed by the Kammenga lab.
#
# @name indels
# @format A data frame with 10257 rows and 15 variables
NULL