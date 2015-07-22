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

scanonefs <- function(cross, phenocol = 5, nperm = 1000, permevery = TRUE) {
    mapresult <- scanone(cross, pheno.col = phenocol)
    threshold <- as.numeric(quantile(scanone(cross, n.perm = nperm), probs = .95))
    lods <- data.frame(mapresult)
    covars = matrix(nrow = length(cross$pheno[, phenocol]), ncol = 0)
    i = 1
    while (max(mapresult)$lod > threshold) {
        if (permevery){
            threshold <- as.numeric(quantile(scanone(cross, n.perm = nperm), probs = .95))
        }
        genoatlocus <- as.numeric(
            pull.geno(cross,
                      max(mapresult)$chr)[,rownames(max(mapresult))[1]])
        covars <- cbind(covars, genoatlocus)
        mapresult <- scanone(cross, phenocol, addcovar = covars)
        if (max(mapresult)$lod > threshold) {
            lods <- cbind(lods, mapresult$lod)
            colnames(lods)[ncol(lods)] <- paste0("lod", i)
            i = i+1
        }
    }
    markers <- linkagemapping::markers
    colnames(markers)[4] <- "phys.pos"
    colnames(lod)[1:2] <- c("chr.roman", "gen.pos")
    lods$SNP <- rownames(lods)
    lods <- left_join(lods, markers)
    return(lods)
}