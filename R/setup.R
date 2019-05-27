#' Load Cross Obj
#' 
#' Download large cross objects if need be and load or load local cross obj.
#' @param name Name of cross object ("N2xCB4856cross", "N2xLSJ2cross", "AF16xHK104cross", or "N2xCB4856cross_full")
#' @return A cross object
#' @export

load_cross_obj <- function(name) {
    dir.create(file.path("~/.linkagemapping"), showWarnings = FALSE)
    if(name == "N2xCB4856cross_full") {
        fname <- "~/.linkagemapping/N2xCB4856cross_full2.Rda"
        if(!file.exists(fname)){
            url = "https://storage.googleapis.com/linkagemapping/data/N2xCB4856cross_full2.Rda"
            res <- tryCatch(download.file(url,
                                          fname,
                                          method="auto"),
                            error=function(e) 1)
        }
        load(fname, envir = globalenv())
    } else {
        data(list = name)
    }
}

#' Get the LOD value for each marker based on the correlation between genotype
#' 
#' @param pheno A phenotype data frame output from the easysorter pipeline
#' @return A fully formatted phenotype data frame ready to be joined to the
#' cross object
#' @importFrom dplyr %>%

mapformat <- function(pheno){
    
    # Make the condensed phenotype name column (condition + trait)
    pheno$conpheno <- paste0(pheno$condition, ".", pheno$trait)
    
    pheno <- pheno[, c("strain", "conpheno", "phenotype")]
    
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
#' @export

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

#' Perform principal component analysis on traits
#' 
#' @param df A regressed dataframe of HTA traits for PCA
#' @return pca_obj, a PCA object output from princomp; 
#' loadings, a vector of loadings for each HTA trait into each PC;
#' cumsums, a matrix of cumulative sums of variance explained by each PC;
#' RIAILPCphenos, a dataframe of PC values for each strain;
#' corr_PC_trait, correlation of each trait to each PC

runPCA <- function(df){

traits_for_PCA <- df %>%
    dplyr::filter(!strain %in% c("N2","CB4856")) %>% #remove parents from the dataset
    dplyr::select(trait, strain, phenotype)%>% #select columns for PCA
    tidyr::spread(trait, phenotype) %>% #spread so each row is a strain
    na.omit()

row.names(traits_for_PCA) <- traits_for_PCA$strain

pc_traits <- traits_for_PCA %>%
    dplyr::select(-strain)

#scale the phenotypes before running princomp
scales_pc_traits <- as.data.frame(scale(pc_traits))

pca_obj <- princomp(scales_pc_traits)

#pull the loadings (amount of each trait in each PC)
loadings <- loadings(pca_obj)[]

#pull the total variance explained with each PC and call it "drug.cumsum"
cumsums <- t(as.matrix(cumsum(pca_obj$sdev^2/sum(pca_obj$sdev^2))))

#keep the PCs that explain at least 95% of the variation in phenotype
keep <- data.frame("Val" = cumsums[1,], "PC" = 1:ncol(cumsums)) %>%
    dplyr::filter(Val > .95)
keep <- as.numeric(min(keep$PC))

#make df of all loadings for the first five PCs
pcaoutput <- as.data.frame(loadings) %>%
    dplyr::mutate(trait = rownames(.)) %>%
    tidyr::gather(component, variance, -trait) %>%
    dplyr::mutate(component = as.numeric(stringr::str_split_fixed(component, "Comp.", 2)[,2])) %>%
    dplyr::filter(component <= keep) %>%
    dplyr::mutate(component = paste0("PC", component))

#calculate the phenotypic values of each strain for each of the top five PCs for mapping!
RIAIL_PCphenos <- as.data.frame(pca_obj$scores) %>%
    dplyr::mutate(strain = rownames(.), condition = "bleomycin")%>%
    tidyr::gather(trait, phenotype, -strain, -condition) %>%
    dplyr::filter(trait %in% c("Comp.1", "Comp.2","Comp.3","Comp.4","Comp.5"))

### Correlation of each PC with each of the HTA traits:
spread_PCA <- RIAIL_PCphenos %>%
    tidyr::spread(trait, phenotype) %>%
    dplyr::select(-condition)
comparephenos <- merge(spread_PCA, traits_for_PCA, by = "strain")%>%
    na.omit() %>%
    ungroup()%>%
    dplyr::select(-strain)

corr_PC_trait <- cor(comparephenos)

corr_PC_trait <- corr_PC_trait[1:keep,(keep+1):31]

return(list(pca_obj, loadings, RIAIL_PCphenos, corr_PC_trait))
}

