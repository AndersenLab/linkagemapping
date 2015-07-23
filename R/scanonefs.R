#' Map all of the traits in a given cross object with forward search
#' 
#' @param cross A complete cross object with the phenotype data merged
#' @param phenocol The index number for the column of the phenotype of interest
#' @param nperm The number of permutations for the GWER calculation. Defaults
#' to \code{1000}.
#' @param permevery Boolean, whether to calculate genome wide error rate for
#' every iteration. Defaults to \code{TRUE}.
#' @return The LOD scores for all markers
#' @export

scanonefs <- function(cross, phenocol = 5, nperm = 1000, permevery = TRUE) {
    
    # Create the list for all of the mappings
    lodslist <- list()
    
    # Map the first iteration
    mapresult <- qtl::scanone(cross, pheno.col = phenocol)
    
    # Get the max LOD score and threshold
    maxlod <- max(mapresult)$lod
    threshold <- as.numeric(
        quantile(qtl::scanone(cross, n.perm = nperm), probs = .95))
    
    # Change the mapresult data frame into the same structure as for the Josh
    # Bloom mapping code
    mapresult <- data.frame(mapresult) %>%
        dplyr::mutate(marker = rownames(.),
                      trait = colnames(cross$pheno)[phenocol],
                      threshold = threshold,
                      iteration = 1) %>%
        dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
    
    # Append the LODs to the lodslist and create an empty covariate matrix to be
    # filled in later
    lodslist <- append(lodslist, list(mapresult))
    covars = matrix(nrow = length(cross$pheno[, phenocol]), ncol = 0)
    
    i = 2
    
    # While there is a peak LOD above the threshold...
    while (maxlod > threshold) {
        # Pull out the genotype for all individuals at the peak marker
        genoatlocus <- as.numeric(
            qtl::pull.geno(cross,
                           max(mapresult)$chr)[,rownames(max(mapresult))[1]])
        
        # Add the genotypes at peak marker to the covariate matrix
        covars <- cbind(covars, genoatlocus)
        mapresult <- qtl::scanone(cross, phenocol, addcovar = covars)
        
        # If permutations are to be completed every iteration, permute again
        if (permevery){
            threshold <- as.numeric(
                quantile(qtl::scanone(cross, n.perm = nperm, addcovar = covars),
                         probs = .95))
        }
        
        # Get the max LOD for this iteration
        maxlod <- max(mapresult)$lod
        
        # If the max LOD is above the threshold, append it to the lodslist and
        # start the next iteration
        if (maxlod > threshold) {
            mapresult <- data.frame(mapresult) %>%
                dplyr::mutate(marker = rownames(.),
                              trait = colnames(cross$pheno)[phenocol],
                              threshold = threshold,
                              iteration = i) %>%
                dplyr::select(marker, chr, pos, trait, lod, threshold, iteration)
            lodslist <- append(lodslist, list(mapresult))
            i = i+1
        }
    }
    
    # rbind all of the elements in the lodslist
    finallods <- dplyr::rbind_all(lodslist)
    
    # Convert all of the marker positions to physical position from genetic
    # position
    finallods$pos <- vapply(finallods$marker, function(marker) {
        return(as.numeric(unlist(
            markers[markers$SNP == marker, "WS244.pos"])))
    }, numeric(1))
    
    return(finallods)
}