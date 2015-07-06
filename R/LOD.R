#' Get the LOD value for each marker based on the correlation between genotype
#' and phenotype at each marker.
#' 
#' @param x A number.
#' @param y A number.
#' @return The LOD scores for all markers

get.LOD.by.COR = function(n.pheno, pheno, gdata, doGPU=FALSE) {
    if(doGPU) {
        return( (-n.pheno*log(1-gputools::gpuCor(pheno, gdata, use='pairwise.complete.obs')$coefficients^2))/(2*log(10)) )
    } else{
        return( (-n.pheno*log(1-cor(pheno, gdata, use='pairwise.complete.obs')^2))/(2*log(10) ) )  
    }
}