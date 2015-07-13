#' Get the LOD value for each marker based on the correlation between genotype
#' 
#' @param pheno A phenotype data frame output from the easysorter pipeline
#' @return A fully formatted phenotype data frame ready to be joined to the
#' cross object
#' @importFrom dplyr %>%
#' @export

mapformat <- function(pheno){
    pheno$conpheno <- paste0(pheno$condition, ".", pheno$trait)
    pheno <- tidyr::spread(pheno, conpheno, phenotype)
    pheno <- pheno %>% dplyr::group_by(strain) %>%
        dplyr::summarise_each(funs = funs(mean(., na.rm = TRUE))) %>%
        dplyr::select(-date, -experiment, -round, -assay, -condition, -control,
                      -plate, -row, -col, -trait) %>%
        dplyr::mutate(id = stringr::str_split_fixed(.$strain, "QX", 2)[,2]) %>%
        dplyr::filter(id != "")
}

#' Merge the cross object and the phenotype data frame using a dplyr left_join
#' 
#' @param cross A cross object
#' @param phenotype The phenotype data frame with the id numbers for each strain
#' @param set Filter the phenotype data to one specific set (Rockman=1 or 2,
#' Andersen=3) before joining to the data frame
#' @return The phenotype element of the cross object
#' @importFrom dplyr %>%
#' @export

mergepheno <- function(cross, phenotype, set=NULL){
    cross$pheno$id <- as.numeric(cross$pheno$id)
    phenotype$id <- as.numeric(phenotype$id)
    if(!is.null(set)){
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="id") %>%
            dplyr::filter(set == set) %>% select(-contains("strain"))
    } else {
        cross$pheno <- dplyr::left_join(cross$pheno, phenotype, by="id")
    }
    cross$pheno <- cross$pheno[order(cross$pheno$id),]
    return(cross$pheno)
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
#' @return A named vector of the number of strains present for each phenotype measurement

count_strains_per_trait = function(pheno) {
    apply(pheno, 2, function(x){sum(!is.na(x))})
}

#' Extract phenotype matrix from cross object, mean center and optionally
#' standardize variance
#' 
#' @param cross A cross object
#' @param set A vector of set IDs for all RIAILs
#' @param setcorrect Boolean, whether or not to correct for sets by scaling phenotype values based on the set ID
#' @param scalevar Boolean, whether or not to standarize the variance
#' @return A matrix of scaled phenotype values

extract_scaled_phenotype=function(cross, set=NULL, setcorrect=FALSE,
                                  scalevar=TRUE){
    p <- cross$pheno %>%
        dplyr::select(which(sapply(., class) == "numeric"), -id, -set)
    if(setcorrect==FALSE) { 
        apply(p, 2, scale, scale=scalevar) 
    } else {
        s=apply(p, 2, function(x) { 
            xs = split(x,set)
            ms = lapply(xs, mean,na.rm=T)
            unlist(mapply(function(x,y) {x-y}, xs, ms))
        })
        apply(s, 2, scale, scale=scaleVar)
    }
}