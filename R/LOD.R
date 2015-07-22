#' Get the LOD value for each marker based on the correlation between genotype
#' and phenotype.
#' 
#' @param npheno The number of phenotypes present in the data set
#' @param pheno The extracted phenotype data from the cross object
#' @param geno The extracted genotype data from the cross object
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
#' graphics card with the gputools package installed. Defaults to \code{FALSE}.
#' @return The LOD scores for all markers

get_lod_by_cor <- function(npheno, pheno, gdata, doGPU = FALSE) {
    if(doGPU) {
        return((-npheno * log(1 - gputools::gpuCor(pheno, gdata, use = 'pairwise.complete.obs')$coefficients ^ 2)) / (2 * log(10)))
    } else {
        return((-npheno * log(1 - cor(pheno, gdata, use = 'pairwise.complete.obs') ^ 2)) / (2 * log(10)))
    }
}

#' Map all of the traits in a given cross object
#' 
#' @param cross A complete cross object with the phenotype data merged
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
#' graphics card and the gputools package installed. Defaults to \code{FALSE}.
#' @return The LOD scores for all markers

map <- function(cross, doGPU = FALSE) {
    # Remove interpolated SNPs
    cross <- qtl::calc.genoprob(cross, step=0)
    
    # Extract the necessary information from the cross object
    scaledpheno <- extract_scaled_phenotype(cross)
    npheno <- count_strains_per_trait(scaledpheno)
    geno <- extract_genotype(cross)
    
    # Do and return the mapping
    lods <- get_lod_by_cor(npheno, scaledpheno, geno, doGPU)
    
    lods <- lodmatrix2scanone(lods, cross)
    
    return(lods)
}

#' Map all of the traits in a given cross object with forward search
#' 
#' @param cross A complete cross object with the phenotype data merged
#' @param iterations The number of iterations for the FDR calculation. Defaults
#' to \code{1000}.
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
#' graphics card with the gputools package installed. Defaults to \code{FALSE}.
#' @return The LOD scores for all markers
#' @export

fsearch <- function(cross, iterations = 1000, doGPU = FALSE) {
    iteration <- 1
    lods <- map(cross)
    fdr <- get_peak_fdr(lods, cross, iterations, doGPU)
    peaks <- get_peaks_above_thresh(lods, fdr)
    if (nrow(peaks) > 0) {
        lodslist <- list()
        while (nrow(peaks) > 0){
            lods <- data.frame(lods)
            lodswiththresh <- cbind(lods, fdr) %>%
                dplyr::mutate(marker = rownames(.), threshold = fdr, iteration = iteration)
            meltedlods <- tidyr::gather(lodswiththresh, trait, lod, -chr, -pos, -marker, -threshold, -iteration) %>%
                dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
            lodslist <- append(lodslist, list(meltedlods))
            resids <- get_pheno_resids(lods, cross, fdr)
            metapheno <- cross$pheno[1:5]
            cross$pheno <- data.frame(metapheno, resids)
            lods <- map(cross)
            fdr <- get_peak_fdr(lods, cross, iterations, doGPU)
            peaks <- get_peaks_above_thresh(lods, fdr)
            iteration <- iteration + 1
        }
        finallods <- dplyr::rbind_all(lodslist)
    } else {
        lods <- data.frame(lods)
        lodswiththresh <- cbind(lods, fdr) %>%
            dplyr::mutate(marker = rownames(.), threshold = fdr, iteration = 1)
        finallods <- tidyr::gather(lodswiththresh, trait, lod, -chr, -pos, -marker, -threshold, -iteration) %>%
            dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
    }
    cat("\nConverting marker position to physical position. This step takes a while...\n")
    finallods$pos <- vapply(finallods$marker, function(marker) {
            return(as.numeric(unlist(
                markers[markers$SNP == marker, "WS244.pos"])))
        }, numeric(1))
    
    return(finallods)
}

