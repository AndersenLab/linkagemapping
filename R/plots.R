#' Plot all lod scores for a forward search mapping by iteration
#' 
#' Plots the linkage mapping for each iteration of scanone
#' 
#' @param map The mapping output for a single trait
#' @return A plot with a line for each iteration of the forward search mapping
#' @export

lodplot <- function(map){
    
    maxmap <- map %>%
            dplyr::group_by(marker) %>%
            dplyr::filter(lod == max(lod))
    
    cis <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    maxmap <- cidefiner(cis, maxmap)
    
    plot <- ggplot2::ggplot(map)
    plot <- plot + 
        ggplot2::geom_line(ggplot2::aes(x = pos/1e6, y = lod, colour = as.factor(iteration)),
                           size = 1, alpha = 0.85) +
        ggplot2::facet_grid(.~chr, scales ="free") +
        ggplot2::labs(x = "Position (Mb)", y = "LOD") +
        ggplot2::scale_colour_discrete(name="Mapping\nIteration") +
        ggplot2::ggtitle(map$trait[1]) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
                       axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                       panel.background = ggplot2::element_rect( color="black",size=1.2))
    
    if(nrow(cis) != 0) {
        plot <- plot + 
            ggplot2::geom_ribbon(data = maxmap,
                        ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod, fill="gray2"),
                        alpha = 0.5, show_guide=FALSE) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=lod + 1.00, fill=as.factor(iteration)),
                       shape=25, size=3.2, show_guide = FALSE) +
            ggplot2::geom_text(data = cis,
                      ggplot2::aes(x=pos/1e6,
                          y=lod + 2.5,
                          label = paste0(100*round(var_exp, digits = 4),"%")),
                      colour = "black", size=5)
    }
    

    return(plot)
}

#' Plot max lod scores for a forward search mapping
#' 
#' @param map The mapping output for a single trait
#' @return A plot with single line with the maximum LOD score at each marker
#' from a mapping
#' @export

maxlodplot <- function(map){
    map1 <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod))
    
    cis <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::mutate(maxlod=max(lod))%>%
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    
    if(nrow(cis) == 0) {
        plot <- ggplot2::ggplot(map1, ggplot2::aes(x = genotype, y = pheno)) + ggplot2::geom_blank()
        return(plot)
    }
    
    map1 <- cidefiner(cis, map1)
    
    plot <- ggplot2::ggplot(map1) +
        ggplot2::aes(x = pos/1e6, y = lod)
    
    if(nrow(cis) != 0) {
        plot <- plot + 
            ggplot2::geom_ribbon(ggplot2::aes(x = pos/1e6, ymin = 0, ymax = ci_lod),
                        fill = "blue", alpha = 0.5) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=maxlod + 1.00),
                       fill ="red", shape=25, size=3.2, show_guide = FALSE) +
            ggplot2::geom_text(data = cis,
                               ggplot2::aes(x=pos/1e6,
                          y=lod + 2.5,
                          label = paste0(100*round(var_exp, digits = 4),"%")),
                      colour = "black", size=5)
    }
    
    plot <- plot + ggplot2::geom_line(size = 1, alpha = 0.85) +
        ggplot2::facet_grid(.~chr, scales ="free") +
        ggplot2::labs(x = "Position (Mb)", y = "LOD") +
        ggplot2::scale_colour_discrete(name="Mapping\nIteration") +
        ggplot2::ggtitle(map1$trait[1]) +
        
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.y = ggplot2::element_text(size=20, face="bold", color="black"),
              plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
              panel.background = ggplot2::element_rect( color="black",size=1.2))
    return(plot)
}

#' Plot the phenotype by genotype split at each peak marker in a mapping
#' 
#' @param cross The cross object used for the mapping
#' @param map The mapping output for a single trait
#' @return A boxplot of the phenotype by genotype split at each peak marker in a
#' mapping
#' @export

pxgplot <- function(cross, map, parent="N2/CB4856") {
    peaks <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    if(nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    
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
    
    split$genotype <- sapply(split$genotype, function(x){
        if(parent=="N2/CB4856") {
            if(x == -1) {
                "N2"
            } else {
                "CB4856"
            }
        } else {
            if(parent=="N2/LSJ2") {
                if(x == -1) {
                    "N2"
                } else {
                    "LSJ2"
                }
            } else {
                if(parent=="AF16/HK104") {
                    if(x==-1) {
                        "AF16"
                    } else {
                        "HK104"
                    }
                }
            }
        }
    })
        
    ggplot2::ggplot(split) +
        ggplot2::geom_boxplot(ggplot2::aes(x = genotype, y = pheno, fill = genotype), outlier.shape = NA) +
        ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue", "LSJ2" = "green", "AF16" = "indianred", "HK104"= "gold")) +
        ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = pheno), alpha = .8) +
        ggplot2::facet_wrap(~ marker, ncol = 5) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.y = ggplot2::element_text(size=20, face="bold", color="black"),
              plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
              legend.position="none",
              panel.background = ggplot2::element_rect(color="black",size=1.2)) +
        ggplot2::ggtitle(peaks$trait[1]) +
        ggplot2::labs(x = "Genotype", y = "Phenotype")
}

#' Plot effect size by marker
#' 
#' @param cross The cross object used for the mapping
#' @param map The mapping output for a single trait
#' @return A line plot of the effect size at each marker
#' @export

effectplot <- function(cross, map, parental = "N2/CB4856") {
    map2 <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod)) %>%
        dplyr::select(marker:iteration)
    
    cis <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::mutate(maxlod=max(lod))%>%
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        do(head(., n=1))
    
    annotated_map <- annotate_lods(map2, cross, annotate_all = TRUE)
    
    
    plot <- ggplot2::ggplot(annotated_map) +
        ggplot2::geom_line(ggplot2::aes(x = pos/1e6, y = eff_size), size = 1, alpha = 0.85) +
        ggplot2::facet_grid(.~chr, scales ="free") +
        ggplot2::labs(x = "Position (Mb)", y = "Effect Size") +
        ggplot2::ggtitle(annotated_map$trait[1])+
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
              axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
              axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
              strip.text.y = ggplot2::element_text(size=20, face="bold", color="black"),
              plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
              panel.background = ggplot2::element_rect( color="black",size=1.2))
    if(nrow(cis) != 0) {
        plot <- plot + 
            ggplot2::geom_rect(data=cis, ggplot2::aes(xmin=ci_l_pos/1e6, xmax=ci_r_pos/1e6, ymin=-Inf, ymax=Inf), fill="blue", alpha=.5, show_guide=FALSE) +
            ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=0),
                                fill ="red", shape=25, size=3.2, show_guide = FALSE)
    }
    plot
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