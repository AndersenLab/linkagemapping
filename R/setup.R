#' Get the LOD value for each marker based on the correlation between genotype
#' 
#' @param x A number.
#' @param y A number.
#' @return The LOD scores for all markers

renamecols <- function(x){
    colnames(x)[which(colnames(x) == "n"):ncol(x)] <- paste0(x$condition[1], ".", colnames(x)[which(colnames(x) == "n"):ncol(x)])
    return(x)
}

#' Merge the cross object and the phenotype data frame using a dplyr left_join
#' 
#' @param cross A cross object
#' @param phenotype The phenotype data frame with the id numbers for each strain
#' @return The phenotype element of the cross object
#' @importFrom dplyr %>%

mergepheno <- function(cross, phenotype, set=NULL){
    if(!is.null(set)){
        phenotype = phenotype %>% dplyr::filter(set == set)
    }
    cross$pheno$id <- as.numeric(cross$pheno$id)
    phenotype$id <- as.numeric(phenotype$id)
    cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="id")
    cross$pheno <- cross$pheno[order(cross$pheno$id),]
    return(cross$pheno)
}


#' Extract genotype matrix from cross structure and recode as -1, 1
#' 
#' @param cross A cross object
#' @return The genotype matrix, encoded as -1 or 1 for genotype

extractGenotype=function(cross){
    (do.call('cbind', sapply(cross$geno, function(x) { x$argmax }))*2)-3
}

#' Count number of strains with data for each phenotype from cross structure
#' 
#' @param pheno The phenotype elemnt of the cross object
#' @return The genotype matrix, encoded as -1 or 1 for genotype

countStrainsPerTrait = function(pheno) {
    apply(pheno, 2, function(x){sum(!is.na(x))})
}

#' Extract phenotype matrix from cross object, mean center and optionally
#' standardize variance
#' 
#' @param cross A cross object
#' @return The genotype matrix, encoded as -1 or 1 for genotype
#' 

# NOTE!!!!!!!!!!!!! removing first five columns !!!! specific to this data set !!!!!
# badtraits = misbehaving phenos
# set = vector of set ids
# setCorrect = are you doing set correction or not?
# scaleVar = scale variance or not

extractScaledPhenotype=function(cross, set=NULL, setCorrect=FALSE, scaleVar=TRUE){
    p = cross$pheno[,7:ncol(cross$pheno)]
    if(setCorrect==FALSE) { apply(p, 2, scale, scale=scaleVar) } else {
        s=apply(p, 2, function(x) { 
            xs = split(x,set)
            ms = lapply(xs, mean,na.rm=T)
            unlist(mapply(function(x,y) {x-y}, xs, ms))
        })
        apply(s, 2, scale, scale=scaleVar)
        
    }
}







