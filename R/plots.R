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
        dplyr::do(head(., n=1))
    
    maxmap <- cidefiner(cis, maxmap)
    
    plot <- ggplot2::ggplot(map)
    
    if(nrow(cis) !=0) {
        plot <- plot + 
            ggplot2::geom_ribbon(data=maxmap,
                                 ggplot2::aes(x=pos/1e6, ymin=0, ymax=ci_lod), fill="grey",
                                 show_guide=FALSE)
    }
    
    plot <- plot + 
        ggplot2::theme_bw() +
        ggplot2::geom_line(ggplot2::aes(x = pos/1e6, y = lod, color = as.factor(iteration)),
                           size = 1, alpha = 0.85) +
        ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
        ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
        ggplot2::ggtitle(Hmisc::capitalize(map$trait[1])) +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.text.y = ggplot2::element_text(size=16, face="bold", color="black"),
                       axis.title.x = ggplot2::element_text(size=20, face="bold", color="black", vjust=-.3),
                       axis.title.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.x = ggplot2::element_text(size=20, face="bold", color="black"),
                       strip.text.y = ggplot2::element_text(size=20, face="bold", color="black"),
                       plot.title = ggplot2::element_text(size=24, face="bold", vjust = 1),
                       panel.background = ggplot2::element_rect( color="black",size=1.2))
    
    if(nrow(cis) != 0) {
        cis$height1=NULL
        cis$height2=NULL
        for(i in 1:nrow(cis)) {
            if(cis$lod[i] < 10) {
                cis$height1[i] <- 1.3*cis$lod[i]
                cis$height2[i] <- 1.5*cis$lod[i]
            } else  if(cis$lod[i] < 20){
                cis$height1[i] <- 1.15*cis$lod[i]
                cis$height2[i] <- 1.3*cis$lod[i]
            } else {
                cis$height1[i] <- cis$lod[i]+2
                cis$height2[i] <- cis$lod[i]+4
            }
        }
            
        plot <- plot + 
            ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y= height1, fill = as.factor(iteration) ),
                       shape=25, size=3.5, alpha=1, show_guide = FALSE) +
            ggplot2::scale_fill_manual(values = c("#FF66B2", "#00CCCC", "#80FF00", "#9999FF", "#FFB266", "#0000FF","#FFFF66")) +
            ggplot2::scale_colour_manual(name="Mapping\nIteration", values = c("#FF66B2", "#00CCCC", "#80FF00", "#9999FF", "#FFB266", "#0000FF","#FFFF66")) +
            ggplot2::geom_text(data = cis,
                      ggplot2::aes(x=pos/1e6, y=height2, label = paste0(100*round(var_exp, digits = 4),"%")),
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
        dplyr::do(head(., n=1))
    
    
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
            ggplot2::geom_point(data = cis, ggplot2::aes(x=pos/1e6, y=(1.05*maxlod)),
                       fill ="red", shape=25, size=3.2, show_guide = FALSE) +
            ggplot2::geom_text(data = cis,
                               ggplot2::aes(x=pos/1e6,
                          y=(1.2*maxlod),
                          label = paste0(100*round(var_exp, digits = 4),"%")),
                      colour = "black", size=3)
    }
    
    plot <- plot + ggplot2::geom_line(size = 1, alpha = 0.85) +
        ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
        ggplot2::labs(x = "Genomic Position (Mb)", y = "LOD") +
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
#' @param parent The parental strains used to generate RIAILs. Either "N2xCB4856"
#' (default), "N2xLSJ2", or "AF16xHK104"
#' @return A boxplot of the phenotype by genotype split at each peak marker in a
#' mapping
#' @export

pxgplot <- function(cross, map, parent="N2xCB4856") {
    peaks <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n=1))
    
    if(nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    
    pheno <- cross$pheno %>%
        dplyr::select_(map$trait[1])
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
        if(is.na(x)) {
            return(NA)
        }
        if(parent=="N2xCB4856") {
            if(x == -1) {
                "N2"
            } else {
                "CB4856"
            }
        } else if(parent=="N2xLSJ2") {
            if(x == -1) {
                "N2"
            } else {
                "LSJ2"
            }
        } else if(parent=="AF16xHK104") {
            if(x==-1) {
                "AF16"
            } else {
                "HK104"
            }
        }
    })
    
    split$genotype <- factor(split$genotype, levels = c("N2","CB4856","LSJ2","AF16","HK104"))
        
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
#' @param parental The parental strains used to generate RIAILs. Either
#' "N2xCB4856" (default), "N2xLSJ2", or "AF16xHK104
#' @return A line plot of the effect size at each marker
#' @export

effectplot <- function(cross, map, parental = "N2xCB4856") {
    map2 <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::filter(lod == max(lod)) %>%
        dplyr::select(marker:iteration)
    
    cis <- map %>%
        dplyr::group_by(marker) %>%
        dplyr::mutate(maxlod=max(lod))%>%
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n=1))
    
    annotated_map <- annotate_lods(map2, cross, annotate_all = TRUE)
    
    
    plot <- ggplot2::ggplot(annotated_map) +
        ggplot2::geom_line(ggplot2::aes(x = pos/1e6, y = eff_size), size = 1, alpha = 0.85)+
        ggplot2::facet_grid(.~chr, scales ="free", space = "free") +
        ggplot2::labs(x = "Genomic Position (Mb)", y = "Effect Size") +
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

#' Plot the phenotype by genotype split at each peak marker in a mapping
#' next to parental phenotype split
#' 
#' @param cross The cross object used for the mapping
#' @param map The mapping output for a single trait
#' @param phenoframe The dataframe used to make the cross object
#' @param parent The parental cross, either "N2xCB4856" (default), "N2xLSJ2",
#' or "AF16xHK104"
#' @return A boxplot of the phenotype by genotype split at each peak marker in a
#' mapping, including parental phenotypes
#' @export
 
pxgplot_par <- function(cross, map, phenoframe, parent="N2xCB4856") {
    peaks <- map %>% 
        dplyr::group_by(iteration) %>%
        dplyr::filter(!is.na(var_exp)) %>%
        dplyr::do(head(., n=1)) %>%
        tidyr::separate(trait, into = c("condition", "trait"), extra = "merge")
    
    if(nrow(peaks) == 0) {
        stop("No QTL identified")
    }
    
    conditions <- peaks$condition
    traits <- peaks$trait
    
    uniquemarkers <- gsub("-", "\\.", unique(peaks$marker))
    
    colnames(cross$pheno) <- gsub("-", "\\.", colnames(cross$pheno))
    
    pheno <- cross$pheno %>%
        dplyr::select_(map$trait[1])
    
    geno <- data.frame(extract_genotype(cross)) %>%
        dplyr::select(which(colnames(.) %in% uniquemarkers)) %>%
        data.frame(., pheno, cross$pheno$strain)
    
    if(parent == "N2xCB4856") {
        geno <- geno %>%
            dplyr::mutate(RIAIL = as.numeric(gsub("QX", "", cross.pheno.strain)))%>%
            dplyr::filter(RIAIL > 239) %>%
            dplyr::select(-RIAIL)
    }
    
    for(i in 1:(ncol(geno)-2)) {
        for(j in 1:nrow(geno)) {
            geno[j,i] <- ifelse(parent == "N2xCB4856",
                                ifelse(geno[j,i] == -1, "N2-RIAIL", 
                                       ifelse(geno[j,i] == 1, "CB4856-RIAIL", NA)),
                                ifelse(parent == "N2xLSJ2",
                                       ifelse(geno[j,i] == -1, "N2-RIAIL", 
                                              ifelse(geno[j,i] == 1, "LSJ2-RIAIL", NA)),
                                       ifelse(parent == "AF16xHK104",
                                              ifelse(geno[j,i] == -1, "AF16-RIAIL", 
                                                     ifelse(geno[j,i] == 1, "HK104-RIAIL", NA)))))
        }
    }
    
    colnames(geno)[ncol(geno)] <- "strain"
    
    parentphe <- subset(phenoframe, strain %in% c("N2","CB4856","LSJ2","AF16","HK104")) %>%
        dplyr::filter(condition %in% conditions, trait %in% traits) %>%
        dplyr::ungroup()%>%
        dplyr::select(strain, phenotype)
    
    colnames(parentphe)[2] <- map$trait[1]
    
    new <- plyr::rbind.fill(geno, parentphe)
    
    for (i in 1:nrow(new)) {
        for (j in 1:(ncol(new)-2)) {
            new[i,j] <- ifelse(new$strain[i] %in% c("N2","CB4856","LSJ2","AF16","HK104"), new$strain[i], new[i,j])
        }
    }
    
    
    plotting <- new %>%
        tidyr::gather(marker, genotype, 1:(ncol(new)-2))
    
    ggplot2::ggplot(plotting) +
        ggplot2::geom_boxplot(ggplot2::aes(x = factor(x = genotype,
                                                      levels = c("N2" , "CB4856" , "LSJ2", "AF16", "HK104","N2-RIAIL" , "CB4856-RIAIL", "LSJ2-RIAIL", "AF16-RIAIL","HK104-RIAIL"),
                                                      labels = c("N2" , "CB4856" , "LSJ2", "AF16", "HK104","N2-RIAIL" , "CB4856-RIAIL", "LSJ2-RIAIL", "AF16-RIAIL","HK104-RIAIL"),
                                                      ordered = T),
                                           y = plotting[,1], fill = genotype),
                              outlier.shape = NA) +
        ggplot2::geom_jitter(ggplot2::aes(x = genotype, y = plotting[,1]), alpha = .8) +
        ggplot2::scale_fill_manual(values = c("N2" = "orange", "CB4856" = "blue","N2-RIAIL" = "orange", "CB4856-RIAIL" = "blue", "LSJ2" = "green", "LSJ2-RIAIL" = "green", "AF16" = "indianred", "AF16-RIAIL" = "indianred", "HK104" = "gold", "HK104-RIAIL" = "gold")) +
        ggplot2::facet_wrap(~marker, ncol = 5) +
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size = 14, face = "bold", color = "black"),
                       axis.text.y = ggplot2::element_text(size = 16, face = "bold", color = "black"),
                       axis.title.x = ggplot2::element_text(size = 20, face = "bold", color = "black", vjust = -0.3),
                       axis.title.y = ggplot2::element_text(size = 20,face = "bold", color = "black"),
                       strip.text.x = ggplot2::element_text(size = 20, face = "bold", color = "black"),
                       strip.text.y = ggplot2::element_text(size = 20,face = "bold", color = "black"),
                       plot.title = ggplot2::element_text(size = 24, face = "bold", vjust = 1),
                       legend.position = "none",
                       panel.background = ggplot2::element_rect(color = "black",size = 1.2)) +
        ggplot2::ggtitle(paste0(peaks$condition[1]," ",peaks$trait[1])) +
        ggplot2::labs(x = "Genotype", y = "Phenotype")
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
    map$ci_lod <- ci_lod
    return(map)
}

#' Test for correlation among traits in a data frame and output a plot
#' 
#' @param df A dataframe output by easysorter
#' @param gropus Number of groups to split the correlation into, usually based on number of significant PCs
#' @param phenocol When spread, the column number on which traits begin, usually 11 with HTA data.
#' @return A dendrogram and a correlation plot for the data
#' @export

corplot <- function(df, groups, phenocol = NA){
    df <- df %>%
        spread(trait, phenotype) %>%
        ungroup() 
    if(is.na(phenocol)){
        df2 <- df[,11:ncol(df)]
    } else {
        df2 <- df[,phenocolo:ncol(df)]
    }
    
    reorder_cormat <- function(cormat){
        # Use correlation between variables as distance
        dd <- as.dist((1-cormat)/2)
        hc <- hclust(dd) 
        cormat <-cormat[hc$order, hc$order]
    }
    
    traitcormat <- stats::cor(df2, use = "pairwise.complete.obs")
    reordered <- reorder_cormat(traitcormat)
    
    dd <- stats::as.dist((1-traitcormat)/2)
    hc <- stats::hclust(dd)
    
    dend <- ggdendro::ggdendrogram(hc, labels = FALSE, leaf_labels = FALSE) +
        ggplot2::theme(axis.title.y = ggplot2::element_blank(),
                       axis.text.x = ggplot2::element_blank(),
                       axis.text.y = ggplot2::element_blank(),
                       plot.title = ggplot2::element_text(size = 14, face = "bold", hjust = .5),
                       plot.margin = grid::unit(c(0,0,0,0), "cm"))
    
    #clusters1 is the Var1 trait and which cluster it is in:
    #first divide the hc tree into the number of groups you want based on how many significant PCs we have:
    clusters1 <- data.frame(level1 = stats::cutree(hc, k = groups)) %>%
        dplyr::mutate(Var1 = rownames(.))
    #then order the traits based on the way they appear in the dendrogram:
    #extract dendrogram info
    hc2 <- ggdendro::dendro_data(hc)$labels
    #factor so you can order them in the heatmap based on this order
    hc2$label <- factor(hc2$label)
    #order the traits based on their order in the dendrogram
    clusters1$Var1 <- factor(clusters1$Var1, levels = hc2$label)
    #assign cluster groups based on this order that you want the traits (they were grouped by cutree, but group number is not how we want it)
    clusters1 <- clusters1 %>%
        dplyr::arrange(Var1) %>%
        dplyr::group_by(level1)%>%
        dplyr::mutate(newlev = 1)
    for(i in 2:nrow(clusters1)){
        if(clusters1$level1[i] == clusters1$level1[i-1]){
            clusters1$newlev[i] <- clusters1$newlev[i-1]
        } else clusters1$newlev[i] <- clusters1$newlev[i-1]+1
    }
    #clusters2 is the same thing, but for clustering the y-axis traits in the heatmap
    clusters2 <- clusters1 %>%
        dplyr::rename(newlev2 = newlev, Var2 = Var1) %>%
        dplyr::ungroup()%>%
        dplyr::select(-level1)
    #add correlations and the grouping data altogether
    melted <- data.table::melt(reordered) %>%
        dplyr::left_join(clusters1)%>%
        dplyr::left_join(clusters2)
    #reorder Var2 so white line goes diagonal across whole heat map, not just within facet
    melted$Var2 <- factor(melted$Var2, levels = rev(clusters2$Var2))
    
    traitcorplot <- ggplot2::ggplot(data = melted, ggplot2::aes(x = Var1, y = Var2, fill=value^2))+
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradientn(colours=terrain.colors(6), name=expression(bold(italic("r")^2)))+
        ggplot2::scale_x_discrete(position = "top")+
        ggplot2::facet_grid(newlev2~newlev, scales = "free", space = "free", switch = "x")+
        ggplot2::theme_bw() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(size=8, color="black",angle = 90, hjust=0),
                       axis.text.y = ggplot2::element_text(size=8, color="black"),
                       panel.background = ggplot2::element_rect(fill = '#E5E8E8'),
                       strip.text = ggplot2::element_blank(),
                       strip.text.y = ggplot2::element_blank(),
                       panel.spacing = grid::unit(1, "pt"),
                       plot.margin = grid::unit(c(-1,0,0,0),"cm"))+
        ggplot2::labs(title="", x=NULL, y=NULL) 
    
    
    finaldend <- cowplot::plot_grid(NULL, dend, NULL, nrow = 1, ncol = 3, rel_widths = c(.12, 1, .13))
    finalplot <- cowplot::plot_grid(finaldend, traitcorplot, nrow = 2, ncol = 1, rel_heights = c(.5, 1))
    return(finalplot)
}  
