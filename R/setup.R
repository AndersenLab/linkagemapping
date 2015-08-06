#' Get the LOD value for each marker based on the correlation between genotype
#' 
#' @param pheno A phenotype data frame output from the easysorter pipeline
#' @return A fully formatted phenotype data frame ready to be joined to the
#' cross object
#' @importFrom dplyr %>%

mapformat <- function(pheno){
    
    # Make the condensed phenotype name column (condition + trait)
    pheno$conpheno <- paste0(pheno$condition, ".", pheno$trait)
    
    pheno <- pheno %>%
        dplyr::select(strain, conpheno, phenotype)
    
    # Spread the phenotypes and condense down to one row per strain
    pheno <- pheno %>%
        tidyr::spread(conpheno, phenotype) %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise_each(funs = dplyr::funs(mean(., na.rm = TRUE)))
    
    return(pheno)
}

#' Merge the cross object and the phenotype data frame using a dplyr left_join
#' 
#' Incoming phenotype data frame must have the following columns:
#' \code{condition}
#' \code{trait}
#' \code{strain}
#' \code{phenotype}
#' 
#' 
#' 
#' @param cross A cross object
#' @param phenotype The phenotype data frame with the id numbers for each strain
#' @param set Filter the phenotype data to one specific set (Rockman=1,
#' Andersen=2) before joining to the data frame
#' @return A cross object complete with phenotype information
#' @importFrom dplyr %>%
#' @export

mergepheno <- function(cross, phenotype, set=NULL){
    phenotype <- phenotype[phenotype$strain %in% cross$pheno$strain, ]
    
    # Format the phenotype data
    phenotype <- mapformat(phenotype)
    
    # If a specific set is selected, merge only that set's information to the 
    if(!is.null(set)){
        strainset <- cross$pheno$strain[cross$pheno$set == set]
        phenotype <- phenotype[phenotype$strain %in% strainset,]
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="strain")
    } else {
        # Otherwise, merge everything
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="strain")
    }
    
    # Order the phenotype element rows by id
    cross$pheno <- cross$pheno[gtools::mixedorder(cross$pheno$strain), ]
    return(cross)
}


#' Extract genotype matrix from cross structure and recode as -1, 1
#' 
#' @param cross A cross object
#' @return The genotype matrix, encoded as -1 or 1 for genotype

extract_genotype=function(cross){
    
    # Pull out the genotypes into snp x strain matrix
    genomat <- qtl::pull.geno(cross)
    class(genomat) <- "numeric"
    
    # Handle genotype encodings with heterozygous individuals
    # (encoded as 1, 2, 3) and without (encoded 1, 2)
    if(max(genomat, na.rm = TRUE) == 3) {
        genomat <- genomat - 2
    } else {
        genomat <- (genomat * 2) - 3
    }
    return(genomat)
}

#' Count number of strains with data for each phenotype from cross structure
#' 
#' @param pheno The phenotype elemnt of the cross object
#' @return A named vector of the number of strains present for each phenotype
#' measurement

count_strains_per_trait = function(pheno) {
    apply(pheno, 2, function(x){sum(!is.na(x))})
}

#' Extract phenotype matrix from cross object, mean center and optionally
#' standardize variance
#' 
#' @param cross A cross object
#' @param set A vector of set IDs for all RIAILs
#' @param setcorrect Boolean, whether or not to correct for sets by scaling
#' phenotype values based on the set ID
#' @param scalevar Boolean, whether or not to standarize the variance
#' @return A matrix of scaled phenotype values

extract_scaled_phenotype=function(cross, set = NULL, setcorrect = FALSE,
                                  scalevar = TRUE){
    
    # Select only the phenotype columns
    p <- cross$pheno %>%
        dplyr::select(which(sapply(., class) == "numeric"), -set)
    
    # If not corrrecting for set effects...
    if(setcorrect==FALSE) {
        # Scale by the variance on all sets equally
        apply(p, 2, scale, scale=scalevar) 
    } else {
        # Otherwise, break up into individual sets and scale
        s <- apply(p, 2, function(x) { 
            xs <- split(x,set)
            ms <- lapply(xs, mean, na.rm=T)
            unlist(mapply(function(x,y) {x-y}, xs, ms))
        })
        apply(s, 2, scale, scale=scaleVar)
    }
}