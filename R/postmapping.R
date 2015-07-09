#' Convert LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @return The genotype matrix, encoded as -1 or 1 for genotype

lodmatrix2scanone= function(LODS, cross) {
    LODSm = t(as.matrix(LODS))
    LODSs = scanone(cross, pheno.col=10, method='mr')
    LODSso = data.frame(LODSs, LODSm)
    LODSso= LODSso[,-3]
    class(LODSso)=class(LODSs)
    return(LODSso)
}

#' Convert LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' @param cross An example cross object from which to extract scanone skeleton
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

get_chr_peaks = function(mindex.split, chr.mindex.offset, LODS) {
    chr.peaks.lod    = sapply(mindex.split, function(markers) { apply(LODS[,markers], 1, max) })
    # get marker index of LOD peaks per chromosomes                             
    chr.peaks.index = sapply(mindex.split, function(markers)  { apply(LODS[,markers], 1, which.max) })
    # convert chromosome marker index to genome maker index                             
    chr.peaks.index = t(apply(chr.peaks.index, 1, function(x){as.numeric(x)+chr.mindex.offset}))
    return(list(chr.peaks.lod = chr.peaks.lod, chr.peaks.index=chr.peaks.index))
}

#' Convert jb LODmatrix to scanone object
#' 
#' @param LODS A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross An example cross object from which to extract scanone skeleton
#' @param LL
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' @export

get_peak_fdr = function(chromosome.peaks.lod, pdata, gdata, perms=100 ,doGPU=F) { 
    n.pheno = countStrainsPerTrait(pdata) 
    # change dopar for multithreaded
    permpeakLODs = foreach( i = 1:perms ) %do% {
        print(i)
        pLODS <- get.LOD.by.COR(n.pheno, pdata[sample(1:nrow(pdata)),], gdata, doGPU)
        sapply(mindex.split, function(markers) { apply(pLODS[,markers], 1, max) })        
    }
    permpeakLODs= abind(permpeakLODs, along=c(3))
    
    ###### CHROMOSOME and PEAK BASED FDR #################################################################
    obsPcnt = sapply(seq(2, 5, .01), function(thresh) { sum(chromosome.peaks.lod>thresh) }   )
    names(obsPcnt) = seq(2,5, .01)
    # expected number of QTL peaks with LOD greater than threshold
    expPcnt = sapply(seq(2, 5, .01),  
                     function(thresh) {
                         mean(apply(permpeakLODs, 3, function(ll) {sum(ll>thresh) }) )
                     } )
    names(expPcnt) = seq(2, 5, .01)
    pFDR = expPcnt/obsPcnt
    return(pFDR)
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
    tryCatch( {
        # trait number and marker number
        keepPeaks   = which(peaklist$chr.peaks.lod>threshold, arr.ind=T)
        #kP = data.frame(rownames(keepPeaks), peaklist$chr.peaks.index[keepPeaks])
        #names(kP)=c('trait', 'markerIndex') 
        kP = data.frame(trait=keepPeaks[,1], markerIndex=peaklist$chr.peaks.index[keepPeaks])
        kP = kP[order(kP$trait, kP$markerIndex),]
        return(kP)} ,error=function(e) {return(NULL)})
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
    presids = pdata
    si=c()
    for( i in 1:ncol(pdata) ) {
        spA = peakArray[peakArray$trait==i,]
        if(nrow(spA)>0){
            si =c(si,i)
            if(intercept) {
                rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]))
            }else{
                rr = residuals(lm(pdata[,i]~gdata[,spA$markerIndex]-1))
            }
            presids[as.numeric(names(rr)),i]=rr
        }
    }
    print(length(si))
    return(presids)
}
