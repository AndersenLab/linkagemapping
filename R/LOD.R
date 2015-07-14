#' Get the LOD value for each marker based on the correlation between genotype
#' and phenotype at each marker.
#' 
#' @param npheno
#' @param pheno
#' @param geno 
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{true} on machines with an NVIDEA
#' graphics card with the gputools package installed.
#' @return The LOD scores for all markers

get_lod_by_cor <- function(npheno, pheno, gdata, doGPU=FALSE) {
    if(doGPU) {
        return((-npheno*log(1-gputools::gpuCor(pheno, gdata, use = 'pairwise.complete.obs')$coefficients ^ 2)) / (2 * log(10)))
    } else {
        return((-npheno*log(1 - cor(pheno, gdata, use = 'pairwise.complete.obs') ^ 2)) / (2 * log(10)))
    }
}

#' Map all of the traits in a given cross object
#' 
#' @param cross A complete cross object with the phenotype data merged
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{true} on machines with an NVIDEA
#' graphics card with the gputools package installed.
#' @return The LOD scores for all markers
#' @export

map <- function(cross, doGPU=FALSE) {
    # Remove interpolated SNPs
    cross <- calc.genoprob(cross, step=0)
    
    # Extract the necessary information from the cross object
    scaledpheno <- extract_scaled_phenotype(cross)
    npheno <- count_strains_per_trait(scaledpheno)
    geno <- extract_genotype(cross)
    
    # Do and return the mapping
    lods <- get_lod_by_cor(npheno, scaledpheno, geno, doGPU)
    lods <- lodmatrix2scanone(lods, cross)
    return(lods)
} 
