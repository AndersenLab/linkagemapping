#' Convert LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @return The genotype matrix, encoded as -1 or 1 for genotype

lodmatrix2scanone <- function(lods, cross) {
    suppressWarnings({
        LODSm <- t(as.matrix(lods))
        LODSs <- scanone(cross, pheno.col=10, method='mr')
        LODSso <- data.frame(LODSs, LODSm)
        LODSso <- LODSso[,-3]
        class(LODSso) <- class(LODSs)
        return(LODSso)
    })
}

#' Get max lod score and SNP index for each mapped trait
#' 
#' @param lods A data frame output by the mapping functions to be converted to a
#' @param cross An example cross object from which to extract scanone skeleton
#' @return The genotype matrix, encoded as -1 or 1 for genotype

maxpeaks <- function(lods, cross) {
    lods <- data.frame(lods[,3:ncol(lods)])
    maxpeaklod <- vapply(lods, max, numeric(1))
    
    # Get marker index of LOD peaks per chromosomes                             
    maxpeakindex <- vapply(lods, which.max, integer(1))
    return(list(maxpeaklod = maxpeaklod, maxpeakindex=maxpeakindex))
}

#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @importFrom foreach %do%
#' @export

get_peak_fdr <- function(lods, cross, perms=1000, doGPU=F) {
    
    peaklods <- maxpeaks(lods)$maxpeaklods

    pheno <- extract_scaled_phenotype(cross)
    geno <- extract_genotype(cross)
    
    npheno <- count_strains_per_trait(pheno)
    
    # Change dopar for multithreaded
    permpeakLODs <- foreach( i = 1:perms ) %do% {
        print(i)
        lods <- lodmatrix2scanone(get_lod_by_cor(npheno, pheno[sample(1:nrow(pheno)),], geno, doGPU), cross)
        maxpeaks(cross, lods)$maxpeaklod
    }
    
    permpeakLODs <- lapply(permpeakLODs, function(x) data.frame(t(data.frame(x))))
    permpeakLODs <- dplyr::rbind_all(permpeakLODs)
    permpeakLODs <- tidyr::gather(permpeakLODs, trait, lod)
    
    # Chromosome and peak-based FDR
    obsPcnt <- sapply(seq(2, 5, .01), function(thresh) {
        sum(peaklods>thresh)
    })
    names(obsPcnt) = seq(2,5, .01)
    
    # Expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(seq(2, 5, .01), function(thresh) sum(permpeakLODs$lod > thresh))
    names(expPcnt) = seq(2, 5, .01)
    pFDR = expPcnt/obsPcnt
    threshold <- as.numeric(names(pFDR)[which(is.infinite(pFDR))[1]])
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

get_peak_array = function(peaklist, threshold) {
    peaklist <- data.frame(trait = names(peaklist$maxpeaklod), lod = peaklist$maxpeaklod, peaklist$maxpeakindex)
    peaks <- peaklist %>% dplyr::filter(lod > threshold)
    return(peaks)
}

#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

###### fix QTLs and get residual phenotypes ###########################################################
get_pheno_resids = function(pdata, gdata, peakArray, intercept=FALSE) {
    presids <- pdata
    si <- c()
    for (i in 1:ncol(pdata)) {
        spA <- peakArray[peakArray$trait==i,]
        if(nrow(spA) > 0){
            si <- c(si,i)
            if(intercept) {
                rr <- residuals(lm(pdata[,i]~gdata[,spA$markerIndex]))
            }else{
                rr <- residuals(lm(pdata[,i]~gdata[,spA$markerIndex]-1))
            }
            presids[as.numeric(names(rr)),i] <- rr
        }
    }
    print(length(si))
    return(presids)
}


get_pheno_resids = function(cross, map, threshold, intercept=FALSE) {
    pheno <- extract_scaled_phenotype(cross)
    lods <- map[,3:ncol(map)]
    geno <- extract_genotype(cross)
    
    lapply(pdata, )
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
        print(x)
        data <- data.frame(cbind(SNP=rownames(lods), data.frame(lods[,c(1, x)])))
        data$trait <- colnames(lods)[x]
        colnames(data)[3] <- "LOD"
        peaks <- data %>%
            group_by(chr) %>%
            filter(LOD==max(LOD)) %>%
            do(data.frame(.[1,])) %>%
            filter(LOD > threshold)
        return(peaks)
    }))
}



annotate_lods <- function(peaks) {
    # trait chr pos LOD VE scaled_effect_size CI.L CI.R
    peakFit=list()
    for(i in 1:nrow(peaks)) {
        print(i)
        trait=as.character(peaks$trait[i])
        peak.markers=as.character(peaks$SNP[i])
        chr.vec = as.numeric(as.character(peaks$chr[i]))
        LOD.vec = LODS.01[which(rownames(LODS.01)==trait),]
        SNP.name = as.character(peaks$SNP[i]) 
        ax = paste('gdata[,', which(colnames(gdata)==peak.markers),']', sep='')
        aq = paste(ax, collapse= ' + ')
        am = lm(paste('pdata.01s[,' , which(gsub("-", "\\.", colnames(pdata.01s))==trait), ']', '~', (aq), '-1'))
        aov.a=anova(am)
        tssq = sum(aov.a[,2])
        a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
        a.eff.size= as.vector(coefficients(am))
        
        # Calculate confidence interval bounds
        lodsData <- LODS.01s[,c(1, 2, which(colnames(LODS.01s)==trait))]
        lodsData$chr <- as.numeric(as.character(lodsData$chr))
        CIs <- list()
        j <- chr.vec
        int <- cint(lodsData, chr=j, lodcolumn=3)
        CI.L.marker <- rownames(int)[1]
        CI.L.pos <- as.numeric(int[1,2])
        CI.R.marker <- rownames(int)[nrow(int)]
        CI.R.pos <- as.numeric(int[nrow(int),2])
        CI <- data.frame(CI.L.marker, CI.L.pos, CI.R.marker, CI.R.pos)
        CIs <- append(CIs, list(CI))
        CIs <- do.call(rbind, CIs)
        
        
        peakFit=append(peakFit, list(data.frame(cbind(data.frame(trait=trait, SNP=SNP.name, var.exp=a.var.exp, eff.size=a.eff.size), CIs))))
    }
    peakFit.df = do.call('rbind', peakFit)
    
    lods <- data.frame(cbind(SNP=rownames(LODS.01s), data.frame(LODS.01s)))
    lods[,c(1,2)] <- lapply(lods[,c(1,2)], as.character)
    lods[,3] <- as.numeric(lods[,3])
    meltLods <- melt(lods, id=c("SNP", "chr", "pos"), value="LOD")
    colnames(meltLods) <- c("SNP", "chr", "pos", "trait", "LOD")
    
    finalLods <- merge(meltLods, peakFit.df, by=c("trait", "SNP"), all.x=TRUE)
    return(finalLods)
}
