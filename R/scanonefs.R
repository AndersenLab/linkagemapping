library(qtl)
library(dplyr)

scanonefs <- function(crossobject, phenocol = 5, nperm = 1000, permevery = FALSE) {
    mapresult <- scanone(crossobject, pheno.col = phenocol)
    threshold <- as.numeric(quantile(scanone(crossobject, n.perm = nperm), probs = .95))
    lods <- data.frame(mapresult)
    covars = matrix(nrow = length(crossobject$pheno[, phenocol]), ncol = 0)
    i = 1
    while (max(mapresult)$lod > threshold) {
        if (permevery){
            threshold <- as.numeric(quantile(scanone(crossobject, n.perm = nperm), probs = .95))
        }
        genoatlocus <- as.numeric(
            pull.geno(crossobject,
                      max(mapresult)$chr)[,rownames(max(mapresult))[1]])
        covars <- cbind(covars, genoatlocus)
        mapresult <- scanone(crossobject, phenocol, addcovar = covars)
        if (max(mapresult)$lod > threshold) {
            lods <- cbind(lods, mapresult$lod)
            colnames(lods)[ncol(lods)] <- paste0("lod", i)
            i = i+1
        }
    }
    load("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/markers244.Rda")
    colnames(markers)[4] <- "phys.pos"
    colnames(lod)[1:2] <- c("chr.roman", "gen.pos")
    lods$SNP <- rownames(lods)
    lods <- left_join(lods, markers)
    return(lods)
}