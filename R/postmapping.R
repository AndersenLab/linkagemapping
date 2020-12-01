# Convert LODmatrix to scanone object
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param cross An example cross object from which to extract scanone skeleton
# @return A scanone onject with the resultant mapping data in the lod column

lodmatrix2scanone <- function(lods, cross, mod) {

    # Throws a warning for missing phenotype data that we don't care about
    # because we're only using the fake mapping to get the scanone class object
    suppressWarnings({
        LODSm <- t(as.matrix(lods))
        phenocol <- which(sapply(cross$pheno, class) == "numeric")[1]
        LODSs <- qtl::scanone(cross, pheno.col = phenocol, method='mr', model = mod)
        LODSso <- data.frame(LODSs, LODSm)
        LODSso <- LODSso[,-3]
        class(LODSso) <- class(LODSs)
        return(LODSso)
    })
}

# Get max LOD score and SNP index for each mapped trait
# 
# @param lods A data frame output by the mapping functions to be converted to a
# @param cross An example cross object from which to extract scanone skeleton
# @return A list consisting of two vectors, the first being the max peak LOD
# height and the second being the index of the SNP at which the peak LOD occurs

maxpeaks <- function(lods, cross) {
    
    # Get rid of all non-LOD infomation
    lods <- data.frame(lods[,3:ncol(lods)])
    
    # Get the maximum LOD score for each trait
    maxpeaklod <- vapply(lods, function(x) {max(x, na.rm=TRUE)}, numeric(1))
    
    # Get marker index of LOD peak per trait                         
    maxpeakindex <- vapply(lods, which.max, integer(1))
    
    # Return the LODs and the indices as a two element list
    return(list(maxpeaklod = maxpeaklod, maxpeakindex=maxpeakindex))
}

# Get the FDR value for a particular mapping by phenotype permutation
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param cross The original cross object used to perform the mapping
# @param perms
# @param doGPU Boolean, whether to use the gputools package to speed up,
# mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
# graphics card with the gputools package installed. Defaults to \code{FALSE}.
# @return The value of the 5% FDR threshold
# @importFrom foreach %do% %dopar%

get_peak_fdr <- function(lods, cross, perms=1000, doGPU=F, model) {
    library(foreach)
    # Set the appropriate divisor for printing frequency
    if (perms < 100) {
        div = 1
    } else if (perms >= 100 & perms < 1000) {
        div = 10
    } else {
        div = 100
    }
    
    # Get the maximum peak height for each trait
    peaklods <- maxpeaks(lods, cross)$maxpeaklod

    # Get the information necessary to do the permutation mapping
    pheno <- extract_scaled_phenotype(cross)
    geno <- extract_genotype(cross)
    npheno <- count_strains_per_trait(pheno)
    
    # Print a new line to the console to make the updates prettier
    cat("\n")
    
    # Permute the phenotype data and do a mapping for each round of permutation
    # Change to %dopar% for multithreaded
    permpeakLODs <- foreach::foreach(i = 1:perms) %do% {
        if (i %% div == 0) {
            cat(paste0("Permutation ", i, " of ", perms, "...\n"))
        }
        lods <- lodmatrix2scanone(
            get_lod_by_cor(npheno,
                           pheno[sample(1:nrow(pheno)),],
                           geno,
                           doGPU),
            cross,
            mod = model)
        maxpeaks(lods, cross)$maxpeaklod
    }
    
    # Get all of the permutation peak lods
    permpeakLODs <- lapply(permpeakLODs, function(x) {
            data.frame(t(data.frame(x)))
        })
    permpeakLODs <- dplyr::bind_rows(permpeakLODs)
    permpeakLODs <- tidyr::gather(permpeakLODs, trait, lod)
    
    # Get the obeserved number of peaks 
    obsPcnt <- sapply(seq(2, 10, .01), function(thresh) {
        sum(peaklods>thresh, na.rm = TRUE)
    })
    names(obsPcnt) <- seq(2, 10, .01)
    
    # Average expected number of QTL peaks with LOD greater than threshold
    expPcnt <- sapply(seq(2, 10, .01), function(thresh) {
            sum(permpeakLODs$lod > thresh, na.rm = TRUE)/perms
        })
    names(expPcnt) <- seq(2, 10, .01)
    
    # Ratio of expected peaks to observed peaks
    pFDR <- expPcnt/obsPcnt
    
    # Get the threshold value such that the ratio of expected peaks to observed
    # peaks is less than .05
    belowalpha <- sapply(pFDR, function(x) x < .05)
    suppressWarnings({
            threshold <- as.numeric(names(belowalpha)[min(which(belowalpha), na.rm = T)])
        })
    if (is.na(threshold)) {
        threshold <- as.numeric(
            names(belowalpha)[min(which(is.na(belowalpha)), na.rm = T)])
    }
    return(threshold)
}

# Get the GWER value for a particular mapping by phenotype permutation
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param cross The original cross object used to perform the mapping
# @param perms
# @param doGPU Boolean, whether to use the gputools package to speed up,
# mapping. This can only be set to \code{TRUE} on machines with an NVIDEA
# graphics card with the gputools package installed. Defaults to \code{FALSE}.
# @return The value of the 5% GWER threshold
# @importFrom foreach %do% %dopar%

get_peak_gwer <- function(lods, cross, perms=1000, doGPU=F, model) {
    library(foreach)
    # Set the appropriate divisor for printing frequency
    if (perms < 100) {
        div = 1
    } else if (perms >= 100 & perms < 1000) {
        div = 10
    } else {
        div = 100
    }
    
    # Get the maximum peak height for each trait
    peaklods <- maxpeaks(lods, cross)$maxpeaklod
    
    # Get the information necessary to do the permutation mapping
    pheno <- extract_scaled_phenotype(cross)
    geno <- extract_genotype(cross)
    npheno <- count_strains_per_trait(pheno)
    
    # Print a new line to the console to make the updates prettier
    cat("\n")
    
    # Permute the phenotype data and do a mapping for each round of permutation
    # Change to %dopar% for multithreaded
    permpeakLODs <- foreach::foreach(i = 1:perms) %do% {
        if (i %% div == 0) {
            cat(paste0("Permutation ", i, " of ", perms, "...\n"))
        }
        lods <- lodmatrix2scanone(
            get_lod_by_cor(npheno,
                           pheno[sample(1:nrow(pheno)),],
                           geno,
                           doGPU),
            cross,
            mod = model)
        maxpeaks(lods, cross)$maxpeaklod
    }
    
    # Get all of the permutation peak lods
    permpeakLODs <- lapply(permpeakLODs, function(x) {
        data.frame(t(data.frame(x)))
    })
    permpeakLODs <- dplyr::bind_rows(permpeakLODs)
    permpeakLODs <- tidyr::gather(permpeakLODs, trait, lod)
    
    threshold <- quantile(permpeakLODs$lod, probs = .95)
    return(threshold)
}

# Regress genotype from phenotype and resturn the residual phenotype values
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param cross The cross object used for the original mapping
# @param threshold The FDR threshold value used to determine significant peaks
# @param intercept Boolean stating whether or not to include intercept term in
# linear model. Defaults to \code{FALSE}.
# @return The residuals of the phenotypes

get_pheno_resids = function(lods, cross, threshold, intercept = FALSE) {
    
    # Get the scaled phenotype, the LODs, and the genotype data
    pheno <- data.frame(extract_scaled_phenotype(cross))
    clnames <- colnames(lods)
    lods <- data.frame(lods[,3:ncol(lods)])
    colnames(lods) <- clnames[3:length(clnames)]
    geno <- data.frame(extract_genotype(cross))
    
    # Get only the traits with a peak above threshold
    traitsabovethresh <- lapply(1:ncol(lods), function(x){
        if(max(lods[x], na.rm = TRUE)>threshold){
            return(colnames(lods)[x])
        }
    })
    tat <- unlist(traitsabovethresh)
    
    # Get the index of the peak marker
    maxsnp <- vapply(tat, function(x){
        which.max(lods[,x])
    }, numeric(1))
    
    # Make a data frame of traits and peak markers
    peaks <- data.frame(trait = tat, snp = maxsnp)
    
    # Do the regression with or without the intercept term
    if (intercept) {
        presids <- data.frame(lapply(1:nrow(peaks), function(x) {
            residuals(lm(pheno[,peaks$trait[x]]~geno[,peaks$snp[x]] - 1,
                         na.action = na.exclude))
        }))
    } else {
        presids <- data.frame(lapply(1:nrow(peaks), function(x) {
            residuals(lm(pheno[,peaks$trait[x]]~geno[,peaks$snp[x]] - 1,
                         na.action = na.exclude))
        }))
    }
    
    # Rename the columns to traits and return the resids as a data frame 
    colnames(presids) <- peaks$trait
    
    return(presids)
}


# Return entries from a LOD data frame that are above the threshold value
# 
# If more two or more markers from the same trait share the max LOD score, only
# the first marker will be taken to be the peak. 
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param threshold The FDR threshold value used to determine significant peaks
# @return A subset of the \code{lods} data for only markers with a LOD score
# above the threshold value

get_peaks_above_thresh <- function(lods, threshold) {
    do.call(rbind, lapply(3:ncol(lods), function(x){
        data <- data.frame(cbind(SNP=rownames(lods), data.frame(lods[,c(1, x)])))
        data$trait <- colnames(lods)[x]
        colnames(data)[3] <- "LOD"
        peaks <- data %>%
            dplyr::filter(LOD==max(LOD, na.rm = TRUE)) %>%
            
            # If there is more than one marker with the peak LOD value, only
            # take the first one 
            dplyr::do(data.frame(.[1,])) %>%
            
            dplyr::filter(LOD > threshold)
        return(peaks)
    }))
}

#' Annotate LOD peaks with variance explained, effect size, and confidence
#' interval bounds
#' 
#' @param lods A data frame output by the mapping functions to be converted to a
#' \code{scanone} object
#' @param cross The cross object used for the original mapping
#' @param annotate_all Boolean whether or not to annotate all markers with
#' variance explained and effect size. If \code{FALSE} (default), only peak lods
#' will be annotated.
#' @param bayes Boolean whether or not to calculate confidence intervals based
#' on Bayes statistics (LOD drop 1.5 used to find CI by default)
#' @param cutoff Determines strategy for setting CI (not available for Bayes CIs).
#' Either "chromosomal" (default), including leftmost and rightmost markers 
#' across peak chromosome with LOD above cutoff or "proximal" that ends the CI at 
#' the most proximal left and right marker that drop below cutoff.
#' @return The annotated lods data frame with information added for peak markers
#' of each iteration
#' @export

annotate_lods <- function(lods, cross, annotate_all = FALSE, bayes = FALSE, cutoff = "chromosomal") {
    
    if (annotate_all) {
        peaks <- lods
    } else {
        # Get the peak marker for each trait:iteration pair
        peaks <- lods %>%
            dplyr::filter(lod > threshold) %>%
            dplyr::group_by(trait, iteration) %>%
            dplyr::filter(lod == max(lod, na.rm = TRUE))%>%
            dplyr::distinct(lod, .keep_all =TRUE)
    }
    
    # Handle a situation in which no peaks are detected 
    if (nrow(peaks) == 0) {
        lods$var_exp <- NA
        lods$eff_size <- NA
        lods$ci_l_marker <- NA
        lods$ci_l_pos <- NA
        lods$ci_r_marker <- NA
        lods$ci_r_pos <- NA
        
        return(lods)
    }
    
    
    # Get the genotype and phenotype information
    geno <- data.frame(extract_genotype(cross))
    pheno <- data.frame(linkagemapping:::extract_scaled_phenotype(cross))
    
    # trait chr pos LOD VE scaled_effect_size CI.L CI.R
    peaklist <- lapply(1:nrow(peaks), function(i) {
        
        # Pretty print the progress
        if (nrow(peaks) <= 10) {
            div = 1
        } else if (nrow(peaks) > 10 & nrow(peaks) <= 100) {
            div = 10
        } else {
            div = 100
        }

        if (i %% div == 0) {
            if (!annotate_all) {
                cat(paste0("Annotating peak ", i, " of ", nrow(peaks), "...\n"))
            }
        }
        
        # Get the trait and peak marker
        peaktrait <- as.character(peaks$trait[i])
        
        marker <- gsub('-', '\\.', as.character(peaks$marker[i]))
        
        if (grepl("^[0-9]", marker)) {
            marker <- paste0("X", marker)
        }
        
        # Get the genotype and phenotype for that trait and marker
        genotypes <- geno[, which(colnames(geno) == marker)]
        phenotypes <- pheno[, which(colnames(pheno) == peaktrait)]
        
        # Calculate variance explained and effect size
        am <- lm(phenotypes ~ genotypes - 1)
        modelanova <- anova(am)
        tssq <- sum(modelanova[,2])
        variance_explained <- modelanova[1:(nrow(modelanova)-1),2]/tssq  
        effect_size <- as.vector(coefficients(am))
        
        if (!annotate_all){
            # Calculate confidence interval bounds
            peakchr = as.character(peaks$chr[i])
            peakiteration <- peaks$iteration[i]
            traitlods <- lods %>% dplyr::filter(trait == peaktrait, iteration == peakiteration)
            if(bayes == FALSE){
            confint <- linkagemapping:::cint(traitlods, peakchr, cutoff = cutoff)
            } else {
                confint <- cint_bayes(traitlods, peakchr)
            }
        }
        
        
        # Assemble everything into a row in a data frame and return it
        peak = peaks[i,]
        
        peak$var_exp <- variance_explained
        peak$eff_size <- effect_size
        
        if (!annotate_all) {
            peak$ci_l_marker <- unlist(confint[1, "marker"])
            peak$ci_l_pos <- unlist(confint[1, "pos"])
            peak$ci_r_marker <- unlist(confint[2, "marker"])
            peak$ci_r_pos <- unlist(confint[2, "pos"])
        }
        
        return(peak)
    })
    
    # rbind the list and join it to the complete lods data frame
    ann_peaks <- do.call(rbind, peaklist)
    finallods <- dplyr::left_join(lods, ann_peaks)
    return(finallods)
}

# Calculate the confidence interval bounds for each peak using a LOD drop cutoff
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param chrom The chromosome on which a peak was found
# @param lodcolumn The index of the column containing the lods scores
# @param drop The LOD drop for calculating the confidence interval
# @param cutoff The strategy for setting a threshold - either "chromosomal" (default)
# to include markers across whole chromosome above LOD drop or "proximal" to include
# only the nearest markers to the peak that are above the LOD drop
# @return The marker, chromosome, position, trait, LOD, threshold and
# iteration information for the left and right bounds of the confidence interval

cint <- function(lods, chrom, lodcolumn=5, drop=1.5, cutoff = "chromosomal"){
    
    # Get only the data for the chromosome containing the peak marker so that CI
    # doesn't overflow chromsome bounds
    data <- lods %>%
        dplyr::filter(chr == chrom) %>%
        na.omit() #added this in because there was an issue with X chromosome NAs from npr-1 markers.
    
    #Get the peak index and the peak lod score
    peak <- which.max(unlist(data[,lodcolumn]))
    peakLOD <- unlist(data[peak, lodcolumn])
    
    # If the peak is not at the end of a chromsome...
    if(peak > 1){
        if(cutoff == "chromosomal"){
        # find the leftmost peak with a LOD higher than a 1.5 drop
        left <- data %>%
            dplyr::filter(lod > (peakLOD-drop)) %>%
            dplyr::filter(pos == min(pos))
        left <- which(data[,3] == left$pos)
        } else if(cutoff == "proximal"){
            # keep moving left until you find a LOD drop below threshold
            left <- peak -1
            while(left > 1 & peakLOD - data[left, lodcolumn] < drop){
                left <- left -1
            }
        } else (stop("Cutoff strategy not recognized"))
    }
    else {
        # Otherwise, the peak LOD marker, which is at the end of the chromsome
        # is the left bound of the CI
        left <- 1
    }
    # Repeat the same process on the right side of the interval
    if(peak < nrow(data)){
        if(cutoff == "chromosomal"){
            #select rightmost marker that is below threshold as end of confidence interval
        right <- data %>%
            dplyr::filter(lod > (peakLOD-drop))%>%
            dplyr::filter(pos == max(pos))
        right <- which(data[,3]==right$pos)
        } else if(cutoff == "proximal"){
            # keep moving right until you find a LOD drop below threshold
            right <- peak +1
            while(right < nrow(data) & peakLOD - data[right, lodcolumn] < drop){
                right <- right +1
            }
        } else (stop("Cutoff strategy not recognized"))
    } else {
        right <- nrow(data)
    }
    
    # rbind the left and right bounds and return the interval bound data frame
    bounds <- rbind(data[left,], data[right,])
    return(bounds)
}

# Calculate the confidence interval bounds for each peak using a Bayes statistics
# 
# @param lods A data frame output by the mapping functions to be converted to a
# \code{scanone} object
# @param chrom The chromosome on which a peak was found
# @param prob The threshold for the confidence interval, as a decimal
# @param lodcolumn The index of the column containing the lods scores
# @return The marker, chromosome, position, trait, LOD, threshold and
# iteration information for the left and right bounds of the confidence interval

cint_bayes <- function (lods, chrom, prob = 0.95, lodcolumn = 5) {
    # Get only the data for the chromosome containing the peak marker so that CI
    # doesn't overflow chromsome bounds
    data <- lods[lods$chr==chrom,] %>%
        dplyr::arrange(pos) %>%
        na.omit
    
    loc <- data[, 3]
    #format LOD scores so area under the curve is 1
    width <- (c(loc[-1], loc[length(loc)]) - c(loc[1], loc[-length(loc)]))/2
    width[c(1, length(width))] <- width[c(1, length(width))] * 
        2
    area <- 10^data[, lodcolumn] * width
    area <- area/sum(area)
    
    #order the adjusted LOD scores by decreasing value
    o <- order(data[, lodcolumn], decreasing = TRUE)
    #state how much of the total area under the curve is described with each LOD score
    cs <- cumsum(area[o])
    #find the closest marker to the peak that explains 95% of the area under the curve
    wh <- min(which(cs >= prob))
    #output all markers that have a lod score above that of the 95% CI defining marker.
    int <- range(o[1:wh])
    
    #this was all done with midpoints between positions, expand it to the nearest marker
    markerpos <- (1:nrow(data))[-grep("^c.+\\.loc-*[0-9]+(\\.[0-9]+)*$", 
                                      rownames(data))]
    if (any(markerpos <= int[1])) 
        int[1] <- max(markerpos[markerpos <= int[1]])
    if (any(markerpos >= int[2])) 
        int[2] <- min(markerpos[markerpos >= int[2]])
    
    bounds <- rbind(data[int[1],], data[int[2],])
    return(bounds)
    
}

#' Find N2 fosmids that tile across a given interval
#' 
#' @param chrom The chromosome of interest, given as a roman numeral value in quotes
#' @param left_pos The left position of the interval
#' @param right_pos The right position of the interval
#' @return A plot of N2 fosmids that tile across the interval of interest and a data frame
#' of chromosome, clone name, fosmid start, fosmid end, and where the fosmid can
#' be found in the Andersen lab.
#' @export


findN2fosmids <- function(chrom, left_pos, right_pos) {
    region <- AllN2fosmids %>%
        dplyr::filter(chr == chrom) %>%
        dplyr::group_by(clone) %>%
        dplyr::filter(start, (between(start, left_pos, right_pos) | between(end, left_pos, right_pos) | between(left_pos, start, end) | between(right_pos, start, end))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(start) %>%
        dplyr::mutate(plot=seq(2,4*(length(start)),4))
    
    plot <- ggplot2::ggplot(region) +
        ggplot2::aes(xmin = start, xmax = end, ymin = plot, ymax = plot + 0.5) +
        ggplot2::geom_rect(fill = "orange") +
        ggplot2::geom_text(ggplot2::aes(x = start + (0.5 * (end - start)), y = plot + 1, label = clone), size = 4) +
        ggplot2::xlab("Genomic Position") +
        ggplot2::ylab("") +
        ggplot2::geom_vline(xintercept = left_pos, colour= "orange") +
        ggplot2::geom_vline(xintercept = right_pos, colour = "orange") +
        ggplot2::ggtitle(paste0("N2 Fosmids ", chrom, ":",left_pos, "-",right_pos))
    
    assign("N2fosmidsfound", region, envir = .GlobalEnv)
    return(plot)
}

#' Find CB4856 fosmids that tile across a given interval
#' 
#' @param chrom The chromosome of interest, given as a roman numeral value in quotes
#' @param left_pos The left position of the interval
#' @param right_pos The right position of the interval
#' @return A plot of CB4856 fosmids that tile across the interval of interest and a data frame
#' of chromosome, clone name, fosmid start, fosmid end, and where the fosmid can
#' be found in the Andersen lab.
#' @export

findCBfosmids <- function(chrom, left_pos, right_pos) {
    region <- AllCBfosmids %>%
        dplyr::filter(chr == chrom) %>%
        dplyr::group_by(clone) %>%
        dplyr::filter(start, (between(start, left_pos, right_pos) | between(end, left_pos, right_pos) | between(left_pos, start, end) | between(right_pos, start, end))) %>%
        dplyr::ungroup() %>%
        dplyr::arrange(start) %>%
        dplyr::mutate(plot=seq(2,4*(length(start)),4))
    
    plot <- ggplot2::ggplot(region) +
        ggplot2::aes(xmin = start, xmax = end, ymin = plot, ymax = plot + 0.5) +
        ggplot2::geom_rect(fill = "blue") +
        ggplot2::geom_text(ggplot2::aes(x = start + (0.5 * (end - start)), y = plot + 1, label = clone), size = 4) +
        ggplot2::xlab("Genomic Position") +
        ggplot2::ylab("") +
        ggplot2::geom_vline(xintercept = left_pos, colour= "blue") +
        ggplot2::geom_vline(xintercept = right_pos, colour = "blue") +
        ggplot2::ggtitle(paste0("CB4856 Fosmids ", chrom, ":",left_pos, "-",right_pos))
    
    assign("CBfosmidsfound", region, envir = .GlobalEnv)
    return(plot)
}

#' Find expression QTL from the Rockman dataset within a user-defined interval
#' 
#' @param chrom The chromosome of interest, given as an arabic numeral
#' @param left_pos The left position of the interval
#' @param right_pos The right position of the interval
#' @return Assigns two data frames to the global environment. eQTL_PositionInfo
#' contains quantitative information about detected eQTL in the interval.
#' eQTL_GeneDescriptions contains GO terms and other qualitative data.
#' @export

checkeQTLintervals <- function (chrom, left_pos, right_pos) 
{
    sigs <- eQTLpeaks %>% 
        dplyr::filter(chr == chrom,
                      ci_l_pos < right_pos,
                      ci_r_pos > left_pos) %>% 
        dplyr::arrange(trait)
    genes <- probe_info %>% 
        dplyr::filter(probe %in% sigs$trait) %>% 
        dplyr::arrange(probe)
    if (nrow(sigs) == 0) {
        print("NO eQTL")
    }
    else {
        assign("eQTL_PositionInfo", sigs, envir = .GlobalEnv)
        assign("eQTL_GeneDescriptions", genes, envir = .GlobalEnv)
    }
}

#' Find RIAILs for generating NILs to isolate a genomic region
#' 
#' @param CI.L The left marker of the interval, e.g. "UCE1-526"
#' @param CI.R The right marker of the interval, e.g. "CE1-108"
#' @param chromosome The chromosome of the interval as a quoted roman numeral
#' @return Outputs a list of four elements.
#' CompleteInterval1 - dataframe of strains with a continous genotype block across user input interval
#' strains also have a breakpoint in the surrounding region used to generate plot "complete".
#' BreakInterval - dataframe of strains with breaks in the interval and outside the interval
#' used to generate plot "complete".
#' Breaks - ggplot object of all RIAILs identified to have break points in and around the interval
#' Complete - ggplot object of all RIAILs identified to have a continous genotype stretch across the interval
#' @export

FindRIAILsforNILs <- function(CI.L, CI.R, chromosome){
    
    test <- do.call(cbind, lapply(RIAILgenotypes[,4:ncol(RIAILgenotypes)], function(x){
        temp <- data.frame(gen = x)
        temp <- mutate(temp, num = ifelse(gen =="AA",1,
                                          ifelse(gen=="AB",0,NA)))
        temp <- select(temp, -gen)
    }))
    
    test1 <- data.frame(RIAILgenotypes[,1:3], test)
    colnames(test1) <- c("id","chr","CM",seq(1,ncol(test1)-3))
    
    test1 <- left_join(test1, RIAILmarkerconversion, by = "id")
    
    interval <- c(CI.L,CI.R)
    rown <- which(test1$id%in%interval)
    
    
    startDist <- 10
    CompleteInterval1 <- data.frame()
    
    while(length(unique(CompleteInterval1$riail)) < 20){
        
        if(rown[1] > startDist & rown[2] < 1444)
        {
            
            markers <- test1 %>%
                slice(rown[1]:rown[2])%>%
                select(id)
            
            plumin <- c(rown[1]-startDist,rown[2]+startDist)
            
            leftmarks <- test1 %>%
                dplyr::slice(plumin[1]:(rown[1]-1))%>%
                dplyr::select(id)
            
            rightmarks <- test1 %>%
                dplyr::slice((rown[2]+1):plumin[2])%>%
                dplyr::select(id)
            
            intervals <- test1 %>%
                dplyr::slice(rown[1]:rown[2])%>%
                dplyr::select(pos)
            
            ints <- c(min(intervals$pos), max(intervals$pos))
            
            CompleteInterval <- test1 %>%
                dplyr::slice(plumin[1]:plumin[2])%>%
                dplyr::mutate(interval = ifelse(id%in%markers$id, "interval",
                                                ifelse(id%in%leftmarks$id, "left",
                                                       ifelse(id%in%rightmarks$id, "right", NA))))%>%
                tidyr::gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
                dplyr::group_by(riail,interval)%>%
                dplyr::mutate(block = sum(geno, na.rm=T))%>%
                dplyr::ungroup()%>%
                dplyr::group_by(riail)%>%
                dplyr::mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete",
                                            ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
                              left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL",
                                            ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
                              right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                                              ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
                # ungroup()%>%
                dplyr::mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
                              left1 = rep(unique(grep("[:alpha:]",left, value = T))),
                              right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
                dplyr::select(-inte,-left,-right)%>%
                dplyr::rename(inte = inte1, left = left1, right = right1)%>%
                dplyr::ungroup()%>%
                dplyr::mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
                dplyr::mutate(cbgen = ifelse(geno==1,0,
                                             ifelse(geno==0,1,NA)),
                              Lint = ints[1],
                              Rint = ints[2])
            
            
            
            CompleteInterval1 <- CompleteInterval %>%
                dplyr::filter(bothBreak == "breakB" & inte == "complete")
            
            startDist <- startDist + 5
            if(startDist > 40){
                break
            }
        }
        else if(rown[1] <= startDist)
        {
            print("You are at the start of Chromosome I")
            markers <- test1 %>%
                dplyr::slice(rown[1]:rown[2])%>%
                dplyr::select(id)
            
            plumin <- c(1,rown[2]+startDist)
            
            leftmarks <- test1 %>%
                dplyr::slice(plumin[1]:(rown[1]-1))%>%
                dplyr::select(id)
            
            rightmarks <- test1 %>%
                dplyr::slice((rown[2]+1):plumin[2])%>%
                dplyr::select(id)
            
            intervals <- test1 %>%
                dplyr::slice(rown[1]:rown[2])%>%
                dplyr::select(pos)
            
            ints <- c(min(intervals$pos), max(intervals$pos))
            
            CompleteInterval <- test1 %>%
                dplyr::slice(plumin[1]:plumin[2])%>%
                dplyr::mutate(interval = ifelse(id%in%markers$id, "interval",
                                                ifelse(id%in%leftmarks$id, "left",
                                                       ifelse(id%in%rightmarks$id, "right", NA))))%>%
                tidyr::gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
                dplyr::group_by(riail,interval)%>%
                dplyr::mutate(block = sum(geno, na.rm=T))%>%
                dplyr::ungroup()%>%
                dplyr::group_by(riail)%>%
                dplyr::mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete",
                                            ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
                              left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL",
                                            ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
                              right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                                              ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
                # ungroup()%>%
                dplyr::mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
                              left1 = rep(unique(grep("[:alpha:]",left, value = T))),
                              right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
                dplyr::select(-inte,-left,-right)%>%
                dplyr::rename(inte = inte1, left = left1, right = right1)%>%
                dplyr::ungroup()%>%
                dplyr::mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
                dplyr::mutate(cbgen = ifelse(geno==1,0,
                                             ifelse(geno==0,1,NA)),
                              Lint = ints[1],
                              Rint = ints[2])
            
            
            
            CompleteInterval1 <- CompleteInterval %>%
                dplyr::filter(right == "breakR" & inte == "complete")
            
            startDist <- startDist + 5
            if(startDist > 40){
                break
            }
        }
        else if(rown[2] > (1454 - startDist) )
        {
            print("You are at the end of Chromosome X")
            markers <- test1 %>%
                dplyr::slice(rown[1]:rown[2])%>%
                dplyr::select(id)
            
            plumin <- c(rown[1]-startDist,1454)
            
            leftmarks <- test1 %>%
                dplyr::slice(plumin[1]:(rown[1]-1))%>%
                dplyr::select(id)
            
            rightmarks <- test1 %>%
                dplyr::slice((rown[2]+1):plumin[2])%>%
                dplyr::select(id)
            
            intervals <- test1 %>%
                dplyr::slice(rown[1]:rown[2])%>%
                dplyr::select(pos)
            
            ints <- c(min(intervals$pos), max(intervals$pos))
            
            CompleteInterval <- test1 %>%
                dplyr::slice(plumin[1]:plumin[2])%>%
                dplyr::mutate(interval = ifelse(id%in%markers$id, "interval",
                                                ifelse(id%in%leftmarks$id, "left",
                                                       ifelse(id%in%rightmarks$id, "right", NA))))%>%
                tidyr::gather(riail, geno, -id ,-chr,-pos,  -CM,-interval)%>%
                dplyr::group_by(riail,interval)%>%
                dplyr::mutate(block = sum(geno, na.rm=T))%>%
                dplyr::ungroup()%>%
                dplyr::group_by(riail)%>%
                dplyr::mutate(inte = ifelse((interval == "interval") & (block == 0 | block == nrow(markers)), "complete",
                                            ifelse((interval == "interval") & (block != 0 | block != nrow(markers)), "break",NA)),
                              left = ifelse((interval == "left") & (block == 0 | block == nrow(leftmarks)), "completeL",
                                            ifelse((interval == "left") & (block != 0 | block != nrow(leftmarks)), "breakL", NA)),
                              right = ifelse( (interval == "right") & (block == 0 | block == nrow(rightmarks)), "completeR",
                                              ifelse((interval == "right") & (block != 0 | block != nrow(rightmarks)), "breakR",NA)))%>%
                # ungroup()%>%
                dplyr::mutate(inte1 = rep(unique(grep("[:alpha:]",inte, value = T))),
                              left1 = rep(unique(grep("[:alpha:]",left, value = T))),
                              right1 = rep(unique(grep("[:alpha:]",right, value = T))))%>%
                dplyr::select(-inte,-left,-right)%>%
                dplyr::rename(inte = inte1, left = left1, right = right1)%>%
                dplyr::ungroup()%>%
                dplyr::mutate(bothBreak = ifelse(left == "breakL" & right == "breakR", "breakB",NA))%>%
                dplyr::mutate(cbgen = ifelse(geno==1,0,
                                             ifelse(geno==0,1,NA)),
                              Lint = ints[1],
                              Rint = ints[2])
            
            
            CompleteInterval1 <- CompleteInterval %>%
                dplyr::filter(left == "breakL" & inte == "complete")
            
            startDist <- startDist + 5
            if(startDist > 40){
                break
            }
        }
    }
    
    if(length(unique(CompleteInterval1$chr)) > 1){
        CompleteInterval1 <- CompleteInterval1 %>%
            dplyr::filter(chr == chromosome)
        print("You are at the end of a chromosome")
    }
    
    BreakInterval <- CompleteInterval %>%
        dplyr::filter((left == "breakL" | right == "breakR") & inte == "break")
    
    if(length(unique(BreakInterval$chr)) > 1){
        BreakInterval <- BreakInterval %>%
            dplyr::filter(chr == chromosome)
        print("You are at the end of a chromosome")
    }
    
    breaks <- ggplot2::ggplot(BreakInterval)+
        ggplot2::aes(x = pos/1e6, y = geno)+
        ggplot2::facet_grid(riail~chr, scales = "free_y")+
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = geno, alpha=.5), fill = "orange")+
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = cbgen, alpha=.5), fill = "blue")+
        ggplot2::geom_vline(ggplot2::aes(xintercept = Lint/1e6))+
        ggplot2::geom_vline(ggplot2::aes(xintercept = Rint/1e6))+
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=0, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=20, face="bold", color="black"),
                       axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=20,face="bold", color="black"),
                       strip.text.y = ggplot2::element_text(size=20, angle =0, face="bold", color="black"),
                       plot.title = ggplot2::element_text(size=24, face="bold"),
                       legend.position = "none")+
        ggplot2::labs(x = "Genomic Position (Mb)", y = "RIAIL")
    
    complete <- ggplot2::ggplot(CompleteInterval1)+
        ggplot2::aes(x = pos/1e6, y = geno)+
        ggplot2::facet_grid(riail~chr, scales = "free_y")+
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = geno, alpha=.5), fill ="orange")+
        ggplot2::geom_ribbon(ggplot2::aes(ymin = 0, ymax = cbgen, alpha=.5), fill ="blue")+
        ggplot2::geom_vline(ggplot2::aes(xintercept = Lint/1e6))+
        ggplot2::geom_vline(ggplot2::aes(xintercept = Rint/1e6))+
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=0, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=20, face="bold", color="black"),
                       axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=20,face="bold", color="black"),
                       strip.text.y = ggplot2::element_text(size=20, angle =0, face="bold", color="black"),
                       plot.title = ggplot2::element_text(size=24, face="bold"),
                       legend.position = "none")+
        ggplot2::labs(x = "Genomic Position (Mb)", y = "RIAIL")
    
    return(list(CompleteInterval1, BreakInterval, breaks, complete))
}

#' Find indels to generate primers
#' 
#' @param left The left position of the N2 region of interest
#' @param right The right position of the N2 region of interest
#' @param chr The chromosome of interest, written as a quoted roman numeral
#' @return Outputs a dataframe with all known indel information
#' @export
#' 

findindels <- function(left, right, chr) {
    
    subset <- indels %>%
        dplyr::filter(Chromosome == chr) %>%
        dplyr::filter(Pos_N2 > left, Pos_N2 < right)
    
    return(subset)
}

#' Calculate mediation scores between phenotype and gene expression
#' 
#' @param peak peak marker of phenotype QTL in format chr:pos
#' @param probe probeID from expression data
#' @param phenodf Dataframe for phenotype in the form of strain, trait, phenotype with the strain matching RIL from set1 (where we have expression data)
#' @param cross Empty cross object, used to gather genotype information (does not need to be merged with phenotype, but can be)
#' @param scaled Boolean noting if the phenotype should be scaled for a mean of 0 and standard devitaiton of 1. Default is TRUE
#' @param lm Boolean argument. If TRUE, a set of linear models will be used to calculate mediation. If FALSE, mediation will be calculated from the R package "mediation"
#' @export

calc_mediation <- function(peak, expression_probe, phenodf, cross, scaled = T, lm = FALSE) {
    
    # load data
    data("eqtlpheno")
    data("eQTLpeaks")

    # get the genotype at the peak marker
    newpeak <- gsub(":", "_", peak)
    chrom <- stringr::str_split_fixed(newpeak, "_", 2)[,1]
    position <- as.numeric(stringr::str_split_fixed(newpeak, "_", 2)[,2])
    geno <- data.frame(cross$geno[[chrom]]$data)[,newpeak]
    strains <- cross$pheno
    geno <- cbind(strains, geno)
    
    # get the expression pheno for that probe
    probepheno <- eqtlpheno %>%
        dplyr::filter(expression_probe == probe) %>%
        dplyr::select(strain, expression)
    
    # # get the eQTL peak for that probe
    # eqtl_peak <- eQTLpeaks %>%
    #     dplyr::filter(trait == expression_probe,
    #                   chr == chrom,
    #                   ci_l_pos < position | ci_r_pos > position,
    #                   var_exp > 0.01)
    
    # if scaled = T, scale the phenotype (mean = 0, var = 1)
    if(scaled == T) {
        probepheno <- probepheno %>%
            dplyr::mutate(newpheno = (expression - mean(expression, na.rm = T)) / sd(expression, na.rm = T)) %>%
            dplyr::select(strain, expression = newpheno)
    }
    
    # merge drug phenotype and expression phenotype
    pheno <- phenodf %>%
        dplyr::left_join(probepheno, by = "strain") %>%
        dplyr::left_join(geno, by = "strain") %>%
        na.omit()
    
    if(lm) {
        ### LINEAR MODELS ###
        
        # total effect = geno estimate
        model.g <- lm(phenotype ~ geno, data = pheno)
        total <- data.frame(var = "total",
                            estimate = summary(model.g)$coef[2,1],
                            pval = summary(model.g)$coef[2,4])
        
        # direct effect = geno estimate
        model.y <- lm(phenotype ~ expression + geno, data = pheno)
        direct <- data.frame(var = "direct",
                             estimate = summary(model.y)$coef[3,1],
                             pval = summary(model.y)$coef[3,4])
        
        # mediation effect = expression estimate
        med <- data.frame(var = "med",
                          estimate = summary(model.y)$coef[2,1],
                          pval = summary(model.y)$coef[2,4])
        
        # mediation proportion = total - direct / total
        out <- rbind(total, direct, med)
    } else {
        model.m <- lm(expression ~ geno, data = pheno)
        model.y <- lm(phenotype ~ expression + geno, data = pheno)
        out <- mediation::mediate(model.m, model.y, sims = 1000, boot = T, treat = "geno", mediator = "expression")
        
        # Summarize model into data frame
        
        # causal mediation effect
        acme <- data.frame(var = "ACME", 
                           estimate = out$d0, 
                           ci_lower = out$d0.ci[[1]], 
                           ci_upper = out$d0.ci[[2]],
                           prob = out$d0.p)
        
        # direct effect
        ade <- data.frame(var = "ADE", 
                          estimate = out$z0, 
                          ci_lower = out$z0.ci[[1]], 
                          ci_upper = out$z0.ci[[2]],
                          prob = out$z0.p)
        
        # total effect
        total <- data.frame(var = "total", 
                            estimate = out$tau.coef, 
                            ci_lower = out$tau.ci[[1]], 
                            ci_upper = out$tau.ci[[2]],
                            prob = out$tau.p)
        
        # prop. mediated
        med <- data.frame(var = "MED", 
                          estimate = out$n0, 
                          ci_lower = out$n0.ci[[1]], 
                          ci_upper = out$n0.ci[[2]],
                          prob = out$n0.p)
        
        # make a dataframe
        df <- rbind(acme, ade, total, med) %>%
              dplyr::mutate(var = dplyr::case_when(var == "ADE" ~ "direct",
                                                 var == "MED" ~ "prop_med",
                                                 var == "ACME" ~ "med",
                                                 var == "total" ~ "total",
                                                 TRUE ~ "NA"),
                            trait = unique(phenodf$trait),
                            probe = expression_probe,
                            pheno_marker = peak)
                            # eqtl_marker = unique(eqtl_peak$marker))
    }
    
    return(df)
} 

