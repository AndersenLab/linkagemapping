########################################################################################
########################################################################################
########################################################################################
########################################################################################
#
# Figures for RIAILs0 paper
# 
# Figure 1: Schematic of high-throughput assay layout, agar plates go to 96-well plates go to worm sorter go to data of TOF by EXT
# 
# Please see examples of these plots from my recent PLoS Genetics paper
# Figure 2A: Boxplots of N2, CB4856, RIAILs for n (norm.n) trait
# Figure 2B: Non-parametric mapping of norm.n (using new code) with peaks indicated by upside-down red triangles (or something fancier), number of RIAILs scored on the plot, permutation threshold shown as dotted gray line, percent variance explained near each peak, and confidence interval denoted by dotted red line near each peak
# 
# Figure 3A: Histogram of TOF from all individual RIAIL worms, red lines denote the 10, 25, 50, 75, 90th quantiles
# Figures 3B-G: Non-parametric mappings just like Fig.2B in same size to group around the 3A figure matches with the respective quantiles
# 
# Figure 4: Clustering by correlation all traits (quantiles and means for red, green, yellow, EXT, and TOF, and n) with a heat map of correlations. Control and paraquat conditions should be separate figures. 
# 
# Figure 5A: Histogram of RIAIL TOF in paraquat vs. RIAIL TOF in control conditions
# Figure 5B: Histogram of centered RIAIL TOF paraquat and control data
# Figure 5C: Histogram of TOF after regression of assay and control
# Figure 5D: Non-parametric mapping of 5C trait with control and pq mappings underneath it
# 
# I like Raj’s figure for NILs in his Science paper
# Figure 6A: Schematic of confidence interval, and NILs with QX names
# Figure 6B: Boxplots for N2, NIL with N2 in CB background, CB, NIL with CB in N2 background for trait tested with NILs
# 
# Figure 7: Complementation tests with candidates? Candidate genes from seq experiments?
#
########################################################################################
########################################################################################
########################################################################################
########################################################################################

library(dplyr) #Version 0.2


## Figure 2A
fig2DataA <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv")
nData <- fig2DataA %>% filter(drug=="control") %>% select(strain, n) %>% mutate(group=ifelse(strain=="N2", "N2", ifelse(strain=="CB4856", "CB4856", "RIAILs")))

ggplot(nData, aes(x=factor(group), y=n)) + geom_boxplot(aes(fill=factor(group))) + xlab(NULL) + ylab("Number of offspring") + theme(legend.position="none") + ggtitle("Lifetime Fecundity") + scale_x_discrete(limits=c("N2", "CB4856", "RIAILs"), labels=c("Bristol", "Hawaii", "RIAILs")) + scale_fill_manual(values = c("N2" = "blue","CB4856" = "orange","RIAILs" = "gray"))

## Figure 2B
fig2DataBControlMap <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="n")
fig2DataBPQMap <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="paraquat", trait=="n")
fig2DataBControlRaw <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, n) %>% group_by(strain) %>% summarise(n=mean(n))
fig2DataBPQRaw <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="paraquat") %>% select(strain, n) %>% group_by(strain) %>% summarise(n=mean(n))
thresholdControl <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="n") %>% select(threshold) %>% as.numeric(.)
thresholdPQ <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="paraquat", trait=="n") %>% select(threshold) %>% as.numeric(.)

getPeaks <- function(rawData, mapData, threshold){
    strXSNPs <- read.csv("~/LinkageMapping/StrainsBySNPs.csv")
    colnames(strXSNPs)[1] <- "strain" 
    
    inPeak <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
    colnames(inPeak) <- colnames(mapData)
    peaks <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
    colnames(peaks) <- colnames(mapData)
    for(i in 1:nrow(fig2DataBControl)){
        if(mapData[i,]$lod > threshold){
            inPeak <- rbind(inPeak, mapData[i,])
        } else {
            if(nrow(inPeak)!=0){
                peaks <- rbind(peaks, inPeak[which(inPeak$lod==max(inPeak$lod)),])
            }
            inPeak <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
            colnames(inPeak) <- colnames(mapData)
        }
    }
    
    VE <- data.frame(matrix(nrow=0, ncol=(ncol(peaks)+1)))
    colnames(VE) <- c(colnames(peaks), "VE")
    for(i in 1:nrow(peaks)){
        marker <- gsub("-", "\\.", as.character(peaks[i, "marker"]))
        genos <- strXSNPs[,c("strain", marker)]
        total <- merge(genos, rawData, by="strain")
        variance <- cor(total[,2], total[,3])^2
        VE <- rbind(VE, data.frame(peaks[i,], VE=variance, leftMarker, rightMarker))
    }
    return(VE)
}


VE <- getPeaks(fig2DataBControlRaw, fig2DataBControlMap, thresholdControl)

ggplot(fig2DataBControlMap, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("Lifetime fecundity under control conditions\nn=", nrow(fig2DataBControlRaw), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig2DataBControlMap$lod)), shape=25, size = 3, fill="red", colour="red") + geom_text(data=VE, aes(x=pos, y=1.09*max(fig2DataBControlMap$lod), label=paste0(round(100*VE, 0), "%")), size = 3) + geom_segment()

ggplot(fig2DataBPQ, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdPQ, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle("Paraquat N") + xlab("Position") + ylab("LOD")


































## Figure 3A
# Read in data and summarise mean.TOF by strain
fig3DataA <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", !is.na(strain)) %>% group_by(strain) %>% summarise(mean.mean.TOF=mean(mean.TOF))

#Plot the histogram
ggplot(fig3DataA, aes(x=mean.mean.TOF)) + geom_histogram(binwidth=10) + xlab("Mean Time of Flight") + ylab("Count") + geom_vline(xintercept=quantile(fig3DataA$mean.mean.TOF, probs=c(.10, .25, .5, .75, .9), na.rm=TRUE), colour="red")

## Figure 3B
fig3DataBControlMap10 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="q10.TOF")
fig3DataBControlMap25 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="q25.TOF")
fig3DataBControlMap50 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="q50.TOF")
fig3DataBControlMap75 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="q75.TOF")
fig3DataBControlMap90 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_MasterMap.csv") %>% filter(condition=="control", trait=="q90.TOF")

fig3DataBControlRaw10 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, q10.TOF) %>% group_by(strain) %>% summarise(q10.TOF=mean(q10.TOF))
fig3DataBControlRaw25 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, q25.TOF) %>% group_by(strain) %>% summarise(q25.TOF=mean(q25.TOF))
fig3DataBControlRaw50 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, median.TOF) %>% group_by(strain) %>% summarise(median.TOF=mean(median.TOF))
fig3DataBControlRaw75 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, q75.TOF) %>% group_by(strain) %>% summarise(q75.TOF=mean(q75.TOF))
fig3DataBControlRaw90 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(strain!="N2", strain!="CB4856", drug=="control") %>% select(strain, q90.TOF) %>% group_by(strain) %>% summarise(q90.TOF=mean(q90.TOF))

thresholdControl10 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="q10.TOF") %>% select(threshold) %>% as.numeric(.)
thresholdControl25 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="q25.TOF") %>% select(threshold) %>% as.numeric(.)
thresholdControl50 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="median.TOF") %>% select(threshold) %>% as.numeric(.)
thresholdControl75 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="q75.TOF") %>% select(threshold) %>% as.numeric(.)
thresholdControl90 <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_thresholds.csv") %>% filter(condition=="control", trait=="q90.TOF") %>% select(threshold) %>% as.numeric(.)

getPeaks <- function(rawData, mapData, threshold){
    strXSNPs <- read.csv("~/LinkageMapping/StrainsBySNPs.csv")
    colnames(strXSNPs)[1] <- "strain" 
    
    inPeak <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
    colnames(inPeak) <- colnames(mapData)
    peaks <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
    colnames(peaks) <- colnames(mapData)
    for(i in 1:nrow(fig2DataBControl)){
        if(mapData[i,]$lod > threshold){
            inPeak <- rbind(inPeak, mapData[i,])
        } else {
            if(nrow(inPeak)!=0){
                peaks <- rbind(peaks, inPeak[which(inPeak$lod==max(inPeak$lod)),])
            }
            inPeak <- data.frame(matrix(ncol=ncol(mapData), nrow=0))
            colnames(inPeak) <- colnames(mapData)
        }
    }
    
    VE <- data.frame(matrix(nrow=0, ncol=(ncol(peaks)+1)))
    colnames(VE) <- c(colnames(peaks), "VE")
    for(i in 1:nrow(peaks)){
        marker <- gsub("-", "\\.", as.character(peaks[i, "marker"]))
        genos <- strXSNPs[,c("strain", marker)]
        total <- merge(genos, rawData, by="strain")
        variance <- abs(cor(total[,2], total[,3]))
        VE <- rbind(VE, data.frame(peaks[i,], VE=variance))
    }
    return(VE)
}



ggplot(fig3DataBControlMap10, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl10, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("10th quantile of time of flight under control conditions\nn=", nrow(fig3DataBControlMap10), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig3DataBControlMap10$lod)), shape=25, size = 3, fill="red", colour="red") #+ geom_text(data=getPeaks(fig3DataBControlRaw10, fig3DataBControlMap10, thresholdControl10), aes(x=pos, y=1.09*max(fig3DataBControlMap10$lod), label=paste0(round(100*VE, 0), "%")), size = 3)

ggplot(fig3DataBControlMap25, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl25, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("25th quantile of time of flight under control conditions\nn=", nrow(fig3DataBControlMap25), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig3DataBControlMap25$lod)), shape=25, size = 3, fill="red", colour="red") + geom_text(data=getPeaks(fig3DataBControlRaw25, fig3DataBControlMap25, thresholdControl25), aes(x=pos, y=1.09*max(fig3DataBControlMap25$lod), label=paste0(round(100*VE, 0), "%")), size = 3)

ggplot(fig3DataBControlMap50, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl50, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("50th quantile of time of flight under control conditions\nn=", nrow(fig3DataBControlMap50), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig3DataBControlMap50$lod)), shape=25, size = 3, fill="red", colour="red") + geom_text(data=getPeaks(fig3DataBControlRaw50, fig3DataBControlMap50, thresholdControl50), aes(x=pos, y=1.09*max(fig3DataBControlMap50$lod), label=paste0(round(100*VE, 0), "%")), size = 3)

ggplot(fig3DataBControlMap75, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl75, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("75th quantile of time of flight under control conditions\nn=", nrow(fig3DataBControlMap75), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig3DataBControlMap75$lod)), shape=25, size = 3, fill="red", colour="red") + geom_text(data=getPeaks(fig3DataBControlRaw75, fig3DataBControlMap75, thresholdControl75), aes(x=pos, y=1.09*max(fig3DataBControlMap75$lod), label=paste0(round(100*VE, 0), "%")), size = 3)

ggplot(fig3DataBControlMap90, aes(x=pos, y=lod)) + geom_line(size=1) + geom_hline(yintercept=thresholdControl90, colour = "gray", linetype="dashed") + facet_grid(.~chr) + ggtitle(paste0("90th quantile of time of flight under control conditions\nn=", nrow(fig3DataBControlMap90), " RIAILs")) + xlab("Position") + ylab("LOD") + geom_point(data=VE, aes(x=pos, y=1.05*max(fig3DataBControlMap90$lod)), shape=25, size = 3, fill="red", colour="red") + geom_text(data=getPeaks(fig3DataBControlRaw90, fig3DataBControlMap90, thresholdControl90), aes(x=pos, y=1.09*max(fig3DataBControlMap90$lod), label=paste0(round(100*VE, 0), "%")), size = 3)









## Figure 4A
# Do we want to exclude NA columns, rows (with infinites?)
library(gplots)
fig4DataAControl <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(!is.na(strain), drug=="control") %>% select(-(1:10), -contains("resid"), -contains("react"))
fig4DataAPQ <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv") %>% filter(!is.na(strain), drug=="paraquat") %>% select(-(1:10), -contains("resid"), -contains("react"))
fig4DataAControl <- as.matrix(cor(fig4DataAControl)[apply(cor(fig4DataAControl), 1, function(x){sum(is.na(x)) < (ncol(cor(fig4DataAControl))-1)}), apply(cor(fig4DataAControl), 2, function(x){sum(is.na(x)) < (nrow(cor(fig4DataAControl))-1)})])
fig4DataAPQ <- as.matrix(fig4DataAPQ[apply(fig4DataAPQ, 1, function(x){sum(is.na(x)) < (ncol(fig4DataAPQ))}), apply(fig4DataAPQ, 2, function(x){sum(is.na(x)) < (nrow(fig4DataAPQ))})])
fig4DataAPQ <- as.matrix(cor(fig4DataAPQ)[apply(cor(fig4DataAPQ), 1, function(x){sum(is.na(x)) < (ncol(cor(fig4DataAPQ))-1)}), apply(cor(fig4DataAPQ), 2, function(x){sum(is.na(x)) < (nrow(cor(fig4DataAPQ))-1)})])
heatmap.2(fig4DataAControl)
heatmap.2(fig4DataAPQ)
