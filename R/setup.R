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
    pheno <- t <- pheno %>%
        dplyr::filter(strain != "N2", strain != "CB4856") %>%
        tidyr::spread(conpheno, phenotype) %>%
        dplyr::group_by(strain) %>%
        dplyr::summarise_each(funs = dplyr::funs(mean(., na.rm = TRUE))) %>%
        
        # Change ID to QX number
        dplyr::mutate(id = stringr::str_split_fixed(.$strain, "QX", 2)[,2]) %>%
        
        # Remove non-QX strain (e.g. N2 and CB4856)
        dplyr::filter(id != "")
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
    # Format the phenotype data
    phenotype <- mapformat(phenotype)
    cross$pheno$id <- as.numeric(cross$pheno$id)
    phenotype$id <- as.numeric(phenotype$id)
    phenotype <- dplyr::left_join(phenotype, cross$pheno) %>%
        dplyr::select(-QX, -RILname, -strain)
    
    # If a specific set is selected, merge only that set's information to the 
    if(!is.null(set)){
        phenotype <- phenotype[phenotype$set == set,] %>%
            dplyr::select(-set)
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="id")
    } else {
        # Otherwise, merge everything
        phenotype <- phenotype %>%
            dplyr::select(-set)
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="id") %>%
            dplyr::select(-strain)
    }
    
    # Order the phenotype element rows by id
    cross$pheno <- cross$pheno[order(cross$pheno$id),]
    return(cross)
}


#' Extract genotype matrix from cross structure and recode as -1, 1
#' 
#' @param cross A cross object
#' @return The genotype matrix, encoded as -1 or 1 for genotype

extract_genotype=function(cross){
    (do.call('cbind', sapply(cross$geno, function(x) { x$argmax }))*2)-3
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
        dplyr::select(which(sapply(., class) == "numeric"), -id, -set)
    
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