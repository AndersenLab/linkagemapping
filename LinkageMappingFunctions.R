#Begin Josh Bloom mapping functions
library(qtl)
library(rrBLUP)
library(foreach)
library(doMC)
registerDoMC(cores=4)
library(abind)

# Lynch and Walsh p. 454 ##################################################################################
get.LOD.by.COR = function(n.pheno, pheno, gdata, doGPU=FALSE) {
    if(doGPU) {
        return( (-n.pheno*log(1-gpuCor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) )
    } else{
        return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10) ) )  
    }
}

###########################################################################################################


##### Extract genotype matrix from cross structure and recode as -1,1 #####################################
extractGenotype=function(impcross){ (do.call('cbind', sapply(impcross$geno, function(x) { x$argmax }))*2)-3 }
###########################################################################################################


##### Count number of strains with data for each phenotype from cross structure ###########################
countStrainsPerTrait = function(pheno) {apply(pheno, 2, function(x){sum(!is.na(x))})}
###########################################################################################################

##### Extract phenotype matrix from cross structure and mean center ... optional standardize Variance######
# NOTE!!!!!!!!!!!!! removing first five columns !!!! specific to this data set !!!!!
# badtraits = misbehaving phenos
# set = vector of set ids
# setCorrect = are you doing set correction or not?
# scaleVar = scale variance or not
extractScaledPhenotype=function(impcross, set=NULL, setCorrect=FALSE, scaleVar=TRUE){
    p = impcross$pheno[,7:ncol(impcross$pheno)]
    if(setCorrect==FALSE) { apply(p, 2, scale, scale=scaleVar) } else {
        s=apply(p, 2, function(x) { 
            xs = split(x,set)
            ms = lapply(xs, mean,na.rm=T)
            unlist(mapply(function(x,y) {x-y}, xs, ms))
        })
        apply(s, 2,scale, scale=scaleVar)
        
    }
}
###########################################################################################################

############ Convert jb LODmatrix to scanone object #######################################################
LODmatrix.2.scanone= function(LODS, cross, LL=NULL) {
    if(is.null(LL)){
        LODSm = t(as.matrix(LODS))
        LODSs = scanone(cross, pheno.col=10, method='mr')
        LODSso = data.frame(LODSs, LODSm)
        LODSso= LODSso[,-3]
        class(LODSso)=class(LODSs)
        return(LODSso)
    } else { 
        LODSso =  data.frame(LL[,c(1,2)],t(as.matrix(LODS)))
        class(LODSso)  = class(LL)
        return(LODSso)
    }
}
############################################################################################################


getChrPeaks = function(mindex.split, chr.mindex.offset, LODS) {
    chr.peaks.lod    = sapply(mindex.split, function(markers) { apply(LODS[,markers], 1, max) })
    # get marker index of LOD peaks per chromosomes                             
    chr.peaks.index = sapply(mindex.split, function(markers)  { apply(LODS[,markers], 1, which.max) })
    # convert chromosome marker index to genome maker index                             
    chr.peaks.index = t(apply(chr.peaks.index, 1, function(x){x+chr.mindex.offset}))
    return(list(chr.peaks.lod = chr.peaks.lod, chr.peaks.index=chr.peaks.index))
}


getPeakFDR = function(chromosome.peaks.lod, pdata, gdata, perms=100 ,doGPU=F) { 
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


getPeakArray = function(peaklist, threshold) {
    tryCatch( {
        # trait number and marker number
        keepPeaks   = which(peaklist$chr.peaks.lod>threshold, arr.ind=T)
        #kP = data.frame(rownames(keepPeaks), peaklist$chr.peaks.index[keepPeaks])
        #names(kP)=c('trait', 'markerIndex') 
        kP = data.frame(trait=keepPeaks[,1], markerIndex=peaklist$chr.peaks.index[keepPeaks])
        kP = kP[order(kP$trait, kP$markerIndex),]
        return(kP)} ,error=function(e) {return(NULL) })
}


###### fix QTLs and get residual phenotypes ###########################################################
getPhenoResids = function(pdata,gdata, peakArray, intercept=FALSE) {
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
    #return(presids[,si])
    return(presids)
} 
#######################################################################################################
















































map <- function(cross, phenoCol){
    maps <- list()
    residCross <- cross
    iter <- 1
    repeat{
        mapIter <- list()
        for(set in 1:3){
            if(set==1){
                newCross <- subset(residCross, ind=as.character(residCross$pheno$rep)=="A")
            } else if(set == 2){
                newCross <- subset(residCross, ind=as.character(residCross$pheno$rep) %in% c("B", "C", "D"))
            } else if(set==3){
                newCross <- subset(residCross, ind=as.numeric(residCross$pheno$set)==2)
            }
            if (!all(is.na(newCross$pheno[,phenoCol]))) {
                newCross$pheno[,phenoCol] <- as.numeric(as.character(newCross$pheno[,phenoCol]))
                mapData <- scanone(newCross, pheno.col=phenoCol, model="np")
                mapData$marker <- rownames(mapData)
                interation <- 1
                colName <- colnames(cross$pheno)[phenoCol]
                condition <- str_split(colName, "\\.", 2)[[1]][1]
                trait <- str_split(colName, "\\.", 2)[[1]][2]
                df <- as.data.frame(cbind(condition, trait, as.data.frame(mapData)))
                mapIter <- append(mapIter, list(df))
            }
        }
        data <- data.frame(do.call(rbind, mapIter)) %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod=sum(lod))
        if(max(data$lod, na.rm=TRUE)<4 | iter==15){ #Set the threshold cutoff here
            if(iter==1){
                maps <- append(maps, list(data))
            }
            break
        }
        maps <- append(maps, list(data))
        marker <- data$marker[which(data$lod==max(data$lod, na.rm=TRUE))]
        residCross$pheno[,phenoCol] <- residuals(lm(residCross$pheno[,phenoCol]~extractGenotype(residCross)[,marker], na.action=na.exclude))
        iter <- iter+1
    }
    masterMap <- do.call(rbind, maps) %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod=max(lod, na.rm=TRUE))
    return(masterMap)
}

map2 <- function(cross){
    traitMaps <- list()
    for(phenoCol in 6:30){
        print(paste0("Phenotype Column: ", phenoCol))
        maps <- list()
        residCross <- cross
        iter <- 1
        repeat{
            mapIter <- list()
            for(set in 1:3){
                if(set==1){
                    newCross <- subset(residCross, ind=as.character(residCross$pheno$rep)=="A")
                } else if(set == 2){
                    newCross <- subset(residCross, ind=as.character(residCross$pheno$rep) %in% c("B", "C", "D"))
                } else if(set==3){
                    newCross <- subset(residCross, ind=as.numeric(residCross$pheno$set)==2)
                }
                if (!all(is.na(newCross$pheno[,phenoCol]))) {
                    newCross$pheno[,phenoCol] <- as.numeric(as.character(newCross$pheno[,phenoCol]))
                    mapData <- scanone(newCross, pheno.col=phenoCol, model="np")
                    mapData$marker <- rownames(mapData)
                    interation <- 1
                    colName <- colnames(cross$pheno)[phenoCol]
                    condition <- str_split(colName, "\\.", 2)[[1]][1]
                    trait <- str_split(colName, "\\.", 2)[[1]][2]
                    df <- as.data.frame(cbind(condition, trait, as.data.frame(mapData)))
                    mapIter <- append(mapIter, list(df))
                }
            }
            data <- data.frame(do.call(rbind, mapIter)) %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod=sum(lod))
            if(max(data$lod, na.rm=TRUE)<4 | iter==15){ #Set the threshold cutoff here
                if(iter==1){
                    maps <- append(maps, list(data))
                }
                break
            }
            maps <- append(maps, list(data))
            marker <- data$marker[which(data$lod==max(data$lod, na.rm=TRUE))]
            residCross$pheno[,phenoCol] <- residuals(lm(residCross$pheno[,phenoCol]~extractGenotype(residCross)[,marker], na.action=na.exclude))
            iter <- iter+1
        }
        regressedMap <- do.call(rbind, maps) %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod=max(lod, na.rm=TRUE))
        traitMaps <- append(traitMaps, list(regressedMap))
    }
    masterMap <- data.frame(do.call(rbind, traitMaps))
    return(masterMap)
}

extractGenotype=function(impcross){(do.call(cbind, sapply(impcross$geno, function(x) { x$data }))*2)-3}

mergePheno2 <- function(cross, phenotype, set=NULL){
    if(!is.null(set)){
        phenotype = phenotype %>% filter(set == set)
    }
    cross$pheno$id <- as.numeric(cross$pheno$id)
    phenotype$id <- as.numeric(phenotype$id)
    cross$pheno <- left_join(cross$pheno, phenotype, by="id")
    cross$pheno <- cross$pheno[order(cross$pheno$id),]
    return(cross$pheno)
}

renameCols <- function(x){
    colnames(x)[which(colnames(x) == "n"):ncol(x)] <- paste0(x$drug[1], ".", colnames(x)[which(colnames(x) == "n"):ncol(x)])
    return(x)
}

structperm <- function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"), perms,
                       method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
                       addcovar=NULL, intcovar=NULL, weights=NULL,
                       use=c("all.obs", "complete.obs"), upper=FALSE,
                       ties.random=FALSE, start=NULL, maxit=4000,
                       tol=1e-4, verbose){
    
    print(paste0("Pheno.col: ", pheno.col))
    
    lastpheno <- which(colnames(cross$pheno)=="rep")
    newcol <- lastpheno +1
    res <- matrix(nrow = perms, ncol = 1)
    temp1 <- subset(cross, ind = (cross$pheno[,lastpheno] == "A"))
    temp2 <- subset(cross, ind = (cross$pheno[,lastpheno] %in% c("B", "C", "D")))
    temp3 <- subset(cross, ind = (is.na(cross$pheno[,lastpheno])))
    
    for(i in 1:perms){
        if(i %% 100 == 0){
            cat("Structured Permutation", i, "\n")
        }
        temp1$pheno[,newcol] <- sample(temp1$pheno[,pheno.col])
        temp2$pheno[,newcol] <- sample(temp2$pheno[,pheno.col])
        temp3$pheno[,newcol] <- sample(temp3$pheno[,pheno.col])
        scan1 <- scanone(temp1, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        scan2 <- scanone(temp2, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        scan3 <- scanone(temp3, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        all <- scan1+scan2+scan3
        res[i,] <- max(all[,3])
        
    }
    
    rownames(res) <- 1:perms
    colName <- colnames(cross$pheno)[pheno.col]
    condition <- str_split(colName, "\\.", 2)[[1]][1]
    trait <- str_split(colName, "\\.", 2)[[1]][2]
    return(data.frame(cbind(rownames=NULL, condition, trait, threshold=quantile(res, probs=.95)[1])))
    
}

cint <- function(lodsData, chr, lodcolumn=3, drop=1.5){
    data <- lodsData[lodsData$chr==chr,]
    peak <- which.max(data[,lodcolumn])
    peakLOD <- data[peak, 3]
    left <- peak - 1
    while(left >= 1 & peakLOD-data[left, 3] < drop){
        left <- left-1
    }
    right <- peak + 1
    while(right <= nrow(data) & peakLOD-data[right, 3] < drop){
        right <- right+1
    }
    bounds <- rbind(data[left,], data[right,])
    bounds$SNP <- rownames(bounds)
    colnames(bounds) <- c("chr", "pos", "LOD", "SNP")
    return(bounds)
}