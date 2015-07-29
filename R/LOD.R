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
    
    # Calculate the LODS score for each marker:phenotype combination
    if(doGPU) {
        
        # Throw an error if the gputools package is not installed
        if (!requireNamespace("gputools", quietly = TRUE)) {
            stop("The gputools package is required to complete a mapping with 
                 `doGPU` set to `TRUE`. Please install gputools and try again.")
        }
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
    
    # Do the mapping
    lods <- get_lod_by_cor(npheno, scaledpheno, geno, doGPU)
    
    # Convert the map result to a scanone object
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
#' @importFrom dplyr %>%
#' @export

fsearch <- function(cross, phenotype = NULL, iterations = 1000, doGPU = FALSE) {
    
    # Select the subset of the phenotype data frame to map
    if (!is.null(phenotype)) {
        pattern = paste(phenotype, collapse = "|")
        selectedcols <- which(grepl(pattern, colnames(cross$pheno)))
        cross$pheno <- cross$pheno %>%
            dplyr::select(id, QX, RILname, set, selectedcols)
    }
    
    # Set up the iteration count
    iteration <- 1
    
    # Complete the first mapping with FDR calculation
    lods <- map(cross)
    fdr <- get_peak_fdr(lods, cross, iterations, doGPU)
    peaks <- get_peaks_above_thresh(lods, fdr)
    
    # If there are more than 0 peaks in the map...
    if (nrow(peaks) > 0) {
        
        # Create a list to hold all of the mapping iterations
        lodslist <- list()
        
        # While the number of peaks in each iteration is greather than 0,
        # continue with the forward search
        while (nrow(peaks) > 0){
            
            # Bind the threshold and iteration information to the map data,
            # melt the resulting data frame, and add it to the lodslist
            lods <- data.frame(lods)
            lodswiththresh <- lods %>%
                dplyr::mutate(marker = rownames(.), threshold = fdr, iteration = iteration)
            meltedlods <- tidyr::gather(lodswiththresh, trait, lod, -chr, -pos, -marker, -threshold, -iteration) %>%
                dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
            lodslist <- append(lodslist, list(meltedlods))
            
            # Get the residual phenotype cvalues and put them back in the cross
            # object
            resids <- get_pheno_resids(lods, cross, fdr)
            metapheno <- cross$pheno %>% dplyr::select(id, QX, RILname, set)
            cross$pheno <- data.frame(metapheno, resids)
            
            # Repeat the mapping, FDR calculationa and peak finding
            lods <- map(cross)
            fdr <- get_peak_fdr(lods, cross, iterations, doGPU)
            peaks <- get_peaks_above_thresh(lods, fdr)
            
            # Add one to the iteration
            iteration <- iteration + 1
        }
        lodslist <- append(lodslist, list(meltedlods))
        
        # rbind the whole list of maps
        finallods <- dplyr::rbind_all(lodslist)
    } else {
        
        # If no peaks were found from the first iteration, just melt the lods
        # and return them
        lods <- data.frame(lods)
        lodswiththresh <- cbind(lods, fdr) %>%
            dplyr::mutate(marker = rownames(.), threshold = fdr, iteration = 1)
        finallods <- tidyr::gather(lodswiththresh, trait, lod, -chr, -pos, -marker, -threshold, -iteration) %>%
            dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
    }
    
    # Switch genetic position to physical position (WS244) for all of the
    # markers in the map
    cat("\nConverting marker position to physical position. This step takes a while...\n")
    finallods$pos <- vapply(finallods$marker, function(marker) {
            return(as.numeric(unlist(
                markers[markers$SNP == marker, "WS244.pos"])))
        }, numeric(1))
    
    return(finallods)
}

