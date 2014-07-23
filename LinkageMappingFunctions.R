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