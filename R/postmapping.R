#' Convert LODmatrix to scanone object
#' 
#' @param lods A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @return A scanone onject with the resultant mapping data in the lod column

lodmatrix2scanone <- function(lods, cross) {
    # Throws a warning for missing phenotype data that we don't care about
    # because we're only using the fake mapping to get the scanone class object
    suppressWarnings({
        LODSm <- t(as.matrix(lods))
        LODSs <- scanone(cross, pheno.col=6, method='mr')
        LODSso <- data.frame(LODSs, LODSm)
        LODSso <- LODSso[,-3]
        class(LODSso) <- class(LODSs)
        return(LODSso)
    })
}

#' Get max LOD score and SNP index for each mapped trait
#' 
#' @param lods A data frame output by the mapping functions to be converted to a
#' @param cross An example cross object from which to extract scanone skeleton
#' @return A list consisting of two vectors, the first being the max peak LOD
#' height and the second being the index of the SNP at which the peak LOD occurs

maxpeaks <- function(lods, cross) {
    
    # Get rid of all non-LOD infomation
    lods <- data.frame(lods[,3:ncol(lods)])
    
    # Get the maximum LOD score for each trait
    maxpeaklod <- vapply(lods, max, numeric(1))
    
    # Get marker index of LOD peak per trait                         
    maxpeakindex <- vapply(lods, which.max, integer(1))
    
    # Return the LODs and the indices as a two element list
    return(list(maxpeaklod = maxpeaklod, maxpeakindex=maxpeakindex))
}

#' Get the FDR value for a particular mapping by phenotype permutation
#' 
#' @param lods A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross The original cross object used to perform the mapping
#' @param perms
#' @param doGPU Boolean, whether to use the gputools package to speed up,
#' mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
#' graphics card with the gputools package installed. Defaults to \code{FALSE}.
#' @return The value of the 5% FDR threshold
#' @importFrom foreach %do% %dopar%

get_peak_fdr <- function(lods, cross, perms=1000, doGPU=F) {
    
    # Set the appropriate divisor for printing frequency
    if (perms < 100) {
        div = 1
    } else if (perms >= 100 & perms < 1000) {
        div = 10
    } else {
        div = 100
    }
    
    # Get the maximum peak height for each trait
    peaklods <- maxpeaks(lods, cross)$maxpeaklod

    # Get the information necessary to do the permutation mapping
    pheno <- extract_scaled_phenotype(cross)
    geno <- extract_genotype(cross)
    npheno <- count_strains_per_trait(pheno)
    
    # Print a new line to the console to make the updates prettier
    cat("\n")
    
    # Permute the phenotype data and do a mapping for each round of permutation
    # Change to %dopar% for multithreaded
    permpeakLODs <- foreach::foreach(i = 1:perms) %do% {
        if (i %% div == 0) {
            cat(paste0("Permutation ", i, " of ", perms, "...\n"))
        }
        lods <- lodmatrix2scanone(get_lod_by_cor(npheno, pheno[sample(1:nrow(pheno)),], geno, doGPU), cross)
        maxpeaks(lods, cross)$maxpeaklod
    }
    
    # Get all of the permutation peak lods
    permpeakLODs <- lapply(permpeakLODs, function(x) data.frame(t(data.frame(x))))
    permpeakLODs <- dplyr::rbind_all(permpeakLODs)
    permpeakLODs <- tidyr::gather(permpeakLODs, trait, lod)
    
    # Get the obeserved number of peaks 
    obsPcnt <- sapply(seq(2, 5, .01), function(thresh) {
        sum(peaklods>thresh)
    })
    names(obsPcnt) <- seq(2,5, .01)
    
    # Expected number of QTL peaks with LOD greater than threshold
    expPcnt <- sapply(seq(2, 5, .01), function(thresh) sum(permpeakLODs$lod > thresh))
    names(expPcnt) <- seq(2, 5, .01)
    
    # Ratio of expected peaks to observed peaks
    pFDR <- expPcnt/obsPcnt
    
    # Get the threshold value such that the ratio of expected peaks to observed
    # peaks is less than .05
    belowalpha <- sapply(pFDR, function(x) x < .05)
    suppressWarnings({
            threshold <- as.numeric(names(belowalpha)[min(which(belowalpha))])
        })
    if (is.na(threshold)) {
        threshold <- as.numeric(names(belowalpha)[min(which(is.na(belowalpha)))])
    }
    return(threshold)
}

#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

get_pheno_resids = function(cross, lods, threshold, intercept=FALSE) {
    pheno <- data.frame(extract_scaled_phenotype(cross))
    lods <- data.frame(lods[,3:ncol(lods)])
    geno <- data.frame(extract_genotype(cross))
    
    traitsabovethresh <- lapply(1:ncol(lods), function(x){
        if(max(lods[x])>threshold){
            return(colnames(lods)[x])
        }
    })
    tat <- unlist(traitsabovethresh)
    
    maxsnp <- vapply(tat, function(x){
        which.max(lods[,x])
    }, numeric(1))
    
    peaks <- data.frame(trait = tat, snp = maxsnp)
    
    if (intercept) {
        presids <- data.frame(lapply(1:nrow(peaks), function(x) {
            residuals(lm(pheno[,peaks$trait[x]]~geno[,peaks$snp[x]] - 1, na.action = na.exclude))
        }))
    } else {
        presids <- data.frame(lapply(1:nrow(peaks), function(x) {
            residuals(lm(pheno[,peaks$trait[x]]~geno[,peaks$snp[x]] - 1, na.action = na.exclude))
        }))
    }
    
    colnames(presids) <- peaks$trait
    
    return(presids)
}


#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

get_peaks_above_thresh <- function(lods, threshold) {
    do.call(rbind, lapply(3:ncol(lods), function(x){
        data <- data.frame(cbind(SNP=rownames(lods), data.frame(lods[,c(1, x)])))
        data$trait <- colnames(lods)[x]
        colnames(data)[3] <- "LOD"
        peaks <- data %>%
            filter(LOD==max(LOD)) %>%
            do(data.frame(.[1,])) %>%
            filter(LOD > threshold)
        return(peaks)
    }))
}

#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

annotate_lods <- function(lods, cross) {
    peaks <- lods %>%
        dplyr::filter(lod > threshold) %>%
        dplyr::group_by(trait, iteration) %>%
        dplyr::filter(lod == max(lod))
    geno <- data.frame(extract_genotype(cross))
    pheno <- data.frame(extract_scaled_phenotype(cross))
    
    # trait chr pos LOD VE scaled_effect_size CI.L CI.R
    peaklist <- lapply(1:nrow(peaks), function(i) {
        print(i)
        peaktrait <- as.character(peaks$trait[i])
        marker <- gsub('-', '\\.', as.character(peaks$marker[i]))
        
        genotypes <- geno[, which(colnames(geno) == marker)]
        phenotypes <- pheno[, which(colnames(pheno) == peaktrait)]
        
    
        am <- lm(phenotypes ~ genotypes - 1)
        modelanova <- anova(am)
        tssq <- sum(modelanova[,2])
        variance_explained <- modelanova[1:(nrow(modelanova)-1),2]/tssq  
        effect_size <- as.vector(coefficients(am))
        
        # Calculate confidence interval bounds
        peakchr = as.character(peaks$chr[i])
        peakiteration <- peaks$iteration[i]
        traitlods <- lods %>% dplyr::filter(trait == peaktrait, iteration == peakiteration)
        confint <- cint(traitlods, peakchr)
        
        peak = peaks[i,]
        
        peak$var_exp <- variance_explained
        peak$eff_size <- effect_size
        peak$ci_l_marker <- unlist(confint[1, "marker"])
        peak$ci_l_pos <- unlist(confint[1, "pos"])
        peak$ci_r_marker <- unlist(confint[2, "marker"])
        peak$ci_r_pos <- unlist(confint[2, "pos"])
        return(peak)
    })
    ann_peaks <- do.call(rbind, peaklist)
    finallods <- dplyr::left_join(lods, ann_peaks)
    return(finallods)
}



#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

cint <- function(lods, chr, lodcolumn=5, drop=1.5){
    data <- lods[lods$chr==chr,]
    peak <- which.max(unlist(data[,lodcolumn]))
    peakLOD <- unlist(data[peak, lodcolumn])
    if(peak > 1){
        left <- peak - 1
        while(left > 1 & peakLOD - data[left, lodcolumn] < drop){
            left <- left-1
        }
    } else {
        left <- 1
    }
    if(peak < nrow(data)){
        right <- peak + 1
        while(right < nrow(data) & peakLOD - data[right, lodcolumn] < drop){
            right <- right + 1
        }
    } else {
        right <- nrow(data)
    }
    bounds <- rbind(data[left,], data[right,])
    return(bounds)
}