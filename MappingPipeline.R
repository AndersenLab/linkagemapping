# Load required packages
library(regress)
library(dplyr)
library(qtl)
library(stringr)
library(ggplot2)
library(cluster)
library(rrBLUP)
library(foreach)
library(doMC)
library(abind)

# Source the functions
source("~/LinkageMapping/LinkageMappingFunctions.R")

# Set your phenotype data file here

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_ForMapping.csv")

# Remove wash wells (wells where the strain is NA)
reduced.pheno <- pheno[!is.na(pheno$strain),]

# Get the id number for mergePheno2 function, done by simple string split
reduced.pheno$id <- str_split_fixed(reduced.pheno$strain, "QX", 2)[,2]

# Summarize all traits by strain id, drug (really only necessary for RIAILs0, but doesn't hurt other RIAILs experiments)
trait <- reduced.pheno %>% group_by(id, drug) %>% summarise_each(funs(mean)) %>% filter(id!="")

# Rename all of the columns in the format "drug.trait"
trait2 <- trait %>% group_by(id, drug) %>% do(renameCols(.))

# Collapse down to one row per strain
trait3 <- trait2 %>% group_by(id) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>%
    filter(id != "") %>% select(id, contains("\\."), -contains("sorted"))

# Remove columns that are all NA 
trait3 = trait3[,-which(unlist(lapply(trait3, function(x){all(is.na(x))})))]
trait3$id <- as.integer(as.numeric(as.character(trait3$id)))

#-----------------------End processing, start mapping--------------------------#

# Load in the rqtl files 
load("~/HTA_Linkage/Mapping/N2xCB4856_RIAILs_Rqtlfiles.RData")

# Remove interpolated SNPs
N2xCB4856.cross <- calc.genoprob(N2xCB4856.cross, step=0)

# Make the "id" column numeric
N2xCB4856.cross$pheno$id <- as.numeric(sapply(as.character(N2xCB4856.cross$pheno$id), function(x){strsplit(x, "X")[[1]][2]}))

# Reorder by trait id
trait3 <- trait3[order(trait3$id),]

# Remove phenotype data from Matt Rockman's strains
trait3 <- trait3[trait3$id>239,]

# Set the phenotype information for the cross object
N2xCB4856.cross$pheno <- mergePheno2(N2xCB4856.cross, trait3)

# Extract the scaled and centered phenotype data
pdata.01s = extractScaledPhenotype(N2xCB4856.cross)

# Extract the genotype data and get the number of strains per trait

gdata = extractGenotype(N2xCB4856.cross)
n.pheno  = countStrainsPerTrait(pdata.01s)

# Calculate avg RIL relatedness (for VC models)
# Additive

A = A.mat(gdata, shrink=FALSE)/2

# Epistatic

AA = A*A

# Build mindex.split (list of marker indices split by chromosome)

marker.chr.cnt=sapply(N2xCB4856.cross$geno, function(x) ncol(x$data))
chr=rep(names(marker.chr.cnt), marker.chr.cnt)
mindex.split= split(seq_along(chr),chr)

# Convert from genome marker index to chromosome marker index 

chr.mindex.offset = sapply(mindex.split, min)-1

# Get the LOD scores, convert to scanone object and get peaks for each chromosome

LODS.01       = get.LOD.by.COR(n.pheno, pdata.01s, gdata, doGPU=F)
LODS.01s      = LODmatrix.2.scanone(LODS.01, N2xCB4856.cross)

load("~/LinkageMappingSummer2014/markers.Rda")
LODS.01s$pos <- sapply(rownames(LODS.01s), function(x){markers[markers$SNP==x, "WS185.pos"]})


peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01)



# Get the false discovery rate (FDR) for all traits and save immediately (this step can take several hours)
set.seed(0)
LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01s, gdata, 1000, doGPU=F)
# save(LODS.01.FDR, file="RIAILs2FDR.Rda")
# 
# load("~/RIAILs2FDR.Rda")

# Get singular LOD threshold (first infinite occurance in LODs )

threshold <- as.numeric(names(LODS.01.FDR)[which(is.infinite(LODS.01.FDR))[1]])

# Get the array of just peaks above threshold

peaks <- do.call(rbind, lapply(3:ncol(LODS.01s), function(x){
    print(x)
    data <- data.frame(cbind(SNP=rownames(LODS.01s), data.frame(LODS.01s[,c(1, x)])))
    data$trait <- colnames(data)[3]
    colnames(data)[3] <- "LOD"
    peaks <- data %>%
        group_by(chr) %>%
        filter(LOD==max(LOD)) %>%
        do(data.frame(.[1,])) %>%
        filter(LOD > threshold)
    return(peaks)
}))

# trait chr pos LOD VE scaled_effect_size CI.L CI.R
peakFit=list()
for(i in 1:nrow(peaks)) {
    print(i)
    trait=as.character(peaks$trait[i])
    peak.markers=as.character(peaks$SNP[i])
    chr.vec = as.numeric(as.character(peaks$chr[i]))
    LOD.vec = LODS.01[which(rownames(LODS.01)==trait),]
    SNP.name = as.character(peaks$SNP[i]) 
    ax = paste('gdata[,', which(colnames(gdata)==peak.markers),']', sep='')
    aq = paste(ax, collapse= ' + ')
    am = lm(paste('pdata.01s[,' , which(gsub("-", "\\.", colnames(pdata.01s))==trait), ']', '~', (aq), '-1'))
    aov.a=anova(am)
    tssq = sum(aov.a[,2])
    a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
    a.eff.size= as.vector(coefficients(am))
    
    # Calculate confidence interval bounds
    lodsData <- LODS.01s[,c(1, 2, which(colnames(LODS.01s)==trait))]
    lodsData$chr <- as.numeric(as.character(lodsData$chr))
    CIs <- list()
    j <- chr.vec
    int <- cint(lodsData, chr=j, lodcolumn=3)
    CI.L.marker <- rownames(int)[1]
    CI.L.pos <- as.numeric(int[1,2])
    CI.R.marker <- rownames(int)[nrow(int)]
    CI.R.pos <- as.numeric(int[nrow(int),2])
    CI <- data.frame(CI.L.marker, CI.L.pos, CI.R.marker, CI.R.pos)
    CIs <- append(CIs, list(CI))
    CIs <- do.call(rbind, CIs)
    
    
    peakFit=append(peakFit, list(data.frame(cbind(data.frame(trait=trait, SNP=SNP.name, var.exp=a.var.exp, eff.size=a.eff.size), CIs))))
}
peakFit.df = do.call('rbind', peakFit)

lods <- data.frame(cbind(SNP=rownames(LODS.01s), data.frame(LODS.01s)))
lods[,c(1,2)] <- lapply(lods[,c(1,2)], as.character)
lods[,3] <- as.numeric(lods[,3])
meltLods <- melt(lods, id=c("SNP", "chr", "pos"), value="LOD")
colnames(meltLods) <- c("SNP", "chr", "pos", "trait", "LOD")

finalLods <- merge(meltLods, peakFit.df, by=c("trait", "SNP"), all.x=TRUE)

write.csv(finalLods, file="~/RIAILs1Map.csv", row.names=FALSE)

h2.set=list()
for(i in 1:ncol(pdata.01s)){
    trait=i
    print(i)
    trait.num=as.numeric(i)
    trait.name=colnames(pdata.01s)[trait.num]
    rr.all=regress(pdata.01s[,trait.num]~1, ~A, pos=c(T,T) ,verbose=F) 
    
    h2.set[[trait.name]] = list(
        rr.all.sigma=rr.all$sigma, 
        rr.all.se=sqrt(diag(rr.all$sigma.cov))
    )
}

save(h2.set, file="RIAILs0_Heritability.Rda")


mean(sapply(h2.set, function(x){x$rr.all.sigma[1]}))
mean(sapply(h2.set, function(x){x$rr.all.sigma[2]}))
mean(sapply(h2.set, function(x){x$rr.all.se[1]}))
mean(sapply(h2.set, function(x){x$rr.all.se[2]}))

for(i in unique(as.character(peakFit.df$trait))){
    print(i)
    data <- finalLods[as.character(finalLods$trait)==i,]
    peaksDF <- finalLods[!is.na(finalLods$var.exp) & as.character(finalLods$trait)==i,]
    title <- i
    fileName = paste0("~/LinkagePlots/", gsub("\\.", "-", title), "_map.pdf")
    plot = ggplot(data) +
        geom_line(aes(x=pos/1e6, y=LOD), size=1) +
        facet_grid(.~chr) +
        xlab("Position (Mb)") +
        ylab("LOD") +
        ggtitle(title) +
        geom_hline(yintercept=2.97, colour="red", linetype="dashed") +
        geom_vline(data=peaksDF, aes(xintercept=CI.L.pos/1e6), colour="blue", size=1, alpha=.5) +
        geom_vline(data=peaksDF, aes(xintercept=CI.R.pos/1e6), colour="blue", size=1, alpha=.5) +
        geom_point(data=peaksDF, aes(x=pos/1e6, y = 1.15*LOD), fill="red", shape=25, size=4) +
        geom_text(data=peaksDF, aes(x=pos/1e6, y = 1.22*LOD, label=paste0(round(100*var.exp, 2), "%"))) 
    try(ggsave(plot = plot, filename = fileName, width = 11, height = 4.5))
}