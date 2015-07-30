#' Plot all lod scores for a forward search mapping by iteration
#' 
#' @param map The mapping output for a single trait
#' @return A plot with a line for each iteration of the forward search mapping
#' @export

lodplot <- function(map){
    cis <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    maxmap <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod))
    
    maxmap <- cidefiner(cis, maxmap)
    
    plot <- ggplot(map)
    
    if(nrow(cis) != 0) {
        plot <- geom_ribbon(aes(ymin = 0, ymax = ci_lod), fill = "blue", alpha = 0.5) +
            geom_point(data = cis, aes(x=pos/1e6, y=lod + 0.75),
                       fill ="red", shape=25, size=3.2, show_guide = FALSE) +
            geom_text(data = cis,
                      aes(x=pos/1e6,
                          y=lod + 1.5,
                          label = paste0(100*round(var_exp, digits = 4),"%")),
                      colour = "black", size=5)
    }
    
    plot <- plot + 
        geom_line(aes(x = pos/1e6, y = lod, colour = as.factor(iteration)),
                  size = 1, alpha = 0.85) +
        facet_grid(.~chr, scales ="free") +
        labs(x = "Position (Mb)", y = "LOD") +
        scale_colour_discrete(name="Mapping\nIteration") +
        ggtitle(map$trait[1]) +
        
        theme_bw() +
        theme(axis.text.x = element_text(size=16, face="bold", color="black"),
              axis.text.y = element_text(size=16, face="bold", color="black"),
              axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = element_text(size=20, face="bold", color="black"),
              strip.text.x = element_text(size=20, face="bold", color="black"),
              strip.text.y = element_text(size=20, face="bold", color="black"),
              plot.title = element_text(size=24, face="bold", vjust = 1),
              panel.background = element_rect( color="black",size=1.2))
    return(plot)
}

#' Plot max lod scores for a forward search mapping
#' 
#' @param map The mapping output for a single trait
#' @return A plot with single line with the maximum LOD score at each marker
#' from a mapping
#' @export

maxlodplot <- function(map){
    map <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod))
    
    cis <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    if(nrow(cis) == 0) {
        plot <- ggplot(map, aes(x = genotype, y = pheno)) + geom_blank()
        return(plot)
    }
    
    map <- cidefiner(cis, map)
    
    plot <- ggplot(map) +
        aes(x = pos/1e6, y = lod)
    
    if(nrow(cis) != 0) {
        plot <- geom_ribbon(aes(ymin = 0, ymax = ci_lod), fill = "blue", alpha = 0.5) +
            geom_point(data = cis, aes(x=pos/1e6, y=lod + 0.75),
                       fill ="red", shape=25, size=3.2, show_guide = FALSE) +
            geom_text(data = cis,
                      aes(x=pos/1e6,
                          y=lod + 1.5,
                          label = paste0(100*round(var_exp, digits = 4),"%")),
                      colour = "black", size=5)
    }
    
    geom_line(size = 1, alpha = 0.85) +
        facet_grid(.~chr, scales ="free") +
        labs(x = "Position (Mb)", y = "LOD") +
        scale_colour_discrete(name="Mapping\nIteration") +
        ggtitle(map$trait[1]) +
        
        theme_bw() +
        theme(axis.text.x = element_text(size=16, face="bold", color="black"),
              axis.text.y = element_text(size=16, face="bold", color="black"),
              axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = element_text(size=20, face="bold", color="black"),
              strip.text.x = element_text(size=20, face="bold", color="black"),
              strip.text.y = element_text(size=20, face="bold", color="black"),
              plot.title = element_text(size=24, face="bold", vjust = 1),
              panel.background = element_rect( color="black",size=1.2))
    return(plot)
}

#' Plot the phenotype by genotype split at each peak marker in a mapping
#' 
#' @param cross The cross object used for the mapping
#' @param map The mapping output for a single trait
#' @return A boxplot of the phenotype by genotype split at each peak marker in a
#' mapping
#' @export

pxgplot <- function(cross, map) {
    peaks <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    if(nrow(peaks) == 0) {
        plot <- ggplot(map, aes(x = genotype, y = pheno)) + geom_blank()
        return(plot)
    }
    
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    
    pheno <- cross$pheno %>%
        select_(map$trait[1])
    geno <- data.frame(extract_genotype(cross)) %>%
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
        data.frame(., pheno)
    
    colnames(geno)[1:(ncol(geno)-1)] <- sapply(colnames(geno)[1:(ncol(geno)-1)],
                                               function(marker) {
                                                   paste(
                                                       unlist(
                                                           peaks[
                                                               peaks$marker == 
                                                                   gsub("\\.",
                                                                        "-",
                                                                        marker),
                                                               c("chr", "pos")]),
                                                       collapse = ":")
                                               })
    colnames(geno)[ncol(geno)] <- "pheno"
    
    split <- tidyr::gather(geno, marker, genotype, -pheno)
    
    split <- split %>%
        dplyr::mutate(genotype = ifelse(genotype == -1, "N2", "CB4856"))
    
    ggplot(split) +
        geom_boxplot(aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA) +
        scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue")) +
        geom_jitter(aes(x = genotype, y = pheno), alpha = .8) +
        facet_wrap(~ marker, ncol = 5) +
        theme_bw() +
        theme(axis.text.x = element_text(size=16, face="bold", color="black"),
              axis.text.y = element_text(size=16, face="bold", color="black"),
              axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = element_text(size=20, face="bold", color="black"),
              strip.text.x = element_text(size=20, face="bold", color="black"),
              strip.text.y = element_text(size=20, face="bold", color="black"),
              plot.title = element_text(size=24, face="bold", vjust = 1),
              legend.position="none",
              panel.background = element_rect(color="black",size=1.2)) +
        ggtitle(peaks$trait[1]) +
        labs(x = "Genotype", y = "Phenotype")
}

#' Plot effect size by marker
#' 
#' @param cross The cross object used for the mapping
#' @param map The mapping output for a single trait
#' @return A line plot of the effect size at each marker
#' @export

effectplot <- function(cross, map) {
    map <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod)) %>%
        dplyr::select(marker:iteration)
    
    annotated_map <- annotate_lods(map, cross, annotate_all = TRUE)
    
    ggplot(annotated_map) +
        geom_line(aes(x = pos/1e6, y = eff_size), size = 1, alpha = 0.85) +
        facet_grid(.~chr, scales ="free") +
        labs(x = "Position (Mb)", y = "Effect Size") +
        ggtitle(annotated_map$trait[1])+
        theme_bw() +
        theme(axis.text.x = element_text(size=16, face="bold", color="black"),
              axis.text.y = element_text(size=16, face="bold", color="black"),
              axis.title.x = element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = element_text(size=20, face="bold", color="black"),
              strip.text.x = element_text(size=20, face="bold", color="black"),
              strip.text.y = element_text(size=20, face="bold", color="black"),
              plot.title = element_text(size=24, face="bold", vjust = 1),
              panel.background = element_rect( color="black",size=1.2))
}




# A function to help in plotting the confidence intervals

cidefiner <- function(cis, map) {
    ci_lod <- vapply(1:nrow(map), function(marker) {
        pos <- map$pos[marker]
        chr <- map$chr[marker]
        lod <- map$lod[marker]
        s <- sum(chr == cis$chr & (pos >= cis$ci_l_pos & pos <= cis$ci_r_pos))
        inci <- as.logical(s)
        cilodscore <- ifelse(inci, lod, 0)
        return(cilodscore)
    }, numeric(1))
    map <- cbind(map, ci_lod)
    return(map)
}