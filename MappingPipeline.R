# Load required packages
library(MASS)
require(dplyr)
require(qtl)
require(stringr)
require(ggplot2)


# Source the functions
source("~/HTA_Linkage/Mapping/LinkageMappingFunctions.R")

# Set your phenotype data file here

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs2_ForMapping.csv")

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
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01)

# Get the false discovery rate (FDR) for all traits and save immediately (this step can take several hours)
set.seed(0)
LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01s, gdata, 1000, doGPU=F)
save(LODS.01.FDR, file="RIAILs2FDR.Rda")

load("~/RIAILs2FDR.Rda")

# Get singular LOD threshold (first infinite occurance in LODs )

threshold <- as.numeric(names(LODS.01.FDR)[which(is.infinite(LODS.01.FDR))[1]])

# Second argument is sig threshold, should not be hard coded (from LODS.01.FDR)

peakArray.01  = getPeakArray(peaklist.01, threshold)

#second jump run here, run it on only the traits with significant QTL that were detected in round 1
#pdata.02      = getPhenoResids(pdata.01s, gdata, peakArray.01) 
#n.pheno2 = countStrainsPerTrait(pdata.02)
#LODS.02       = get.LOD.by.COR(n.pheno2, pdata.02, gdata, doGPU=T)
#LODS.02s      = LODmatrix.2.scanone(LODS.02, N2xCB4856.cross, LODS.01s)
#peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 
#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02, gdata, 10, doGPU=T)
#peakArray.02  = getPeakArray(peaklist.02, 5)
# only paste together unique indices
#peakArray.02  = rbind(peakArray.01, peakArray.02)

peaksFDR05=data.frame(trait=as.character(colnames(pdata.01s)[peakArray.01[,1]]), marker.index=peakArray.01[,2])
peakList = split(peakArray.01$markerIndex, peakArray.01$trait)

load("~/HTA_Linkage/GenotypeProcessing/markers.Rda")

# trait chr pos LOD VE scaled_effect_size CI.L CI.R
peakFit=list()
for(i in names(peakList)) {
    print(i)
    trait=i
    trait.num=as.numeric(i)
    trait.name=colnames(pdata.01s)[trait.num]
    peak.markers= peakList[[i]]
    chr.vec = chr[peak.markers] 
    LOD.vec = LODS.01[trait.num, peak.markers]
    SNP.name = markers$SNP[peak.markers]
    SNP.pos = markers$WS185.pos[peak.markers]
    ax = paste('gdata[,', peak.markers,']', sep='')
    aq = paste(ax, collapse= ' + ')
    am = lm(paste('pdata.01s[,' , trait.num, ']', '~', (aq), '-1'))
    aov.a=anova(am)
    tssq = sum(aov.a[,2])
    a.var.exp = aov.a[1:(nrow(aov.a)-1),2]/tssq  
    a.eff.size= as.vector(coefficients(am))
    # Mixed model with additive and epistatic term assumes average value per strain
    #rr=regress(pdata.01s[,trait.num]~1, ~A+AA, pos=c(T,T,T) ,verbose=T)
    #wg.additive.trait=rr$sigma[['A']]
    # se are in sqrt(diag(rr$sigma.cov))
    #wg.epistatic.trait=rr$sigma[['AA']]
    #rr=regress(pdata.01s[,trait.num]~1, ~A, pos=c(T,T) ,verbose=T)
    #wg.additive.trait=rr$sigma[['A']]
    
    peakFit[[i]]=data.frame(trait.name, trait.num, peak.markers, chr=chr.vec, SNP.name, SNP.pos, 
                            LOD=LOD.vec, var.exp=a.var.exp, eff.size=a.eff.size)
    #, wg.additive.trait, wg.epistatic.trait)
}
peakFit.df = do.call('rbind', peakFit)


a=sapply(peakFit, function(x) x$wg.additive.trait[1])
aa=sapply(peakFit, function(x) x$wg.epistatic.trait[1])
veQTL=sapply(peakFit, function(x) sum(x$var.exp))

# save.image('~/Desktop/working_image_072514.RData')


#pick out subset of traits with more than 1 QTL
comp.traits=names(which(sapply(peakFit, function(x) nrow(x))>1))
pset= N2xCB4856.cross$pheno$set
p1 = pset ==1 
p2 = pset ==2
p12 = (pset ==1 | pset ==2 )
p3  = (pset ==3 )
A.1 = A.mat(gdata[p1,], shrink=FALSE)/2
A.2 = A.mat(gdata[p2,], shrink=FALSE)/2
A.12= A.mat(gdata[p12,], shrink=FALSE)/2
A.3 = A.mat(gdata[p3,], shrink=FALSE)/2

h2.set=list()
for(i in comp.traits){
    trait=i
    print(i)
    trait.num=as.numeric(i)
    trait.name=colnames(pdata.01s)[trait.num]
    #rr.all=regress(pdata.01s[,trait.num]~1, ~A, pos=c(T,T) ,verbose=F) 
    rr.1 = regress(pdata.01s[which(p1),trait.num]~1, ~A.1, pos=c(T,T) ,verbose=F)
    rr.2 = regress(pdata.01s[which(p2),trait.num]~1, ~A.2, pos=c(T,T) ,verbose=F)
    rr.12=regress(pdata.01s[which(p12),trait.num]~1, ~A.12, pos=c(T,T) ,verbose=F)
    rr.3=regress(pdata.01s[which(p3),trait.num]~1, ~A.3, pos=c(T,T) ,verbose=F)
    
    h2.set[[i]] = list(
        #rr.all.sigma=rr.all$sigma, 
        #rr.all.se=sqrt(diag(rr.all$sigma.cov)),
        rr.1.sigma = rr.1$sigma,
        rr.1.se=sqrt(diag(rr.1$sigma.cov)),
        rr.2.sigma = rr.2$sigma,
        rr.2.se=sqrt(diag(rr.2$sigma.cov)),
        rr.12.sigma=rr.12$sigma, 
        rr.12.se=sqrt(diag(rr.12$sigma.cov)),
        rr.3.sigma=rr.3$sigma,
        rr.3.se=sqrt(diag(rr.3$sigma.cov)) )
}

LODS.01s <- as.data.frame(LODS.01s)

for(col in 3:ncol(LODS.01s)){
    print(col)
    title <- colnames(LODS.01s)[col]
    fileName = paste0("~/LinkagePlots/", gsub("\\.", "-", title), "_map.pdf")
    plot2 = ggplot(LODS.01s, aes(x=pos, y=LODS.01s[,col])) + geom_line(size=1) + facet_grid(.~chr) + xlab("position") + ylab("LOD") + ggtitle(title) + geom_hline(yintercept=2.99, colour="red", linetype="dashed")
    try(ggsave(plot = plot2, filename = fileName, width = 11, height = 4.5))
}

lods <- LODS.01s[,-c(1:2)]

save(LODS.01s, file="RIAILs1LODs.Rda")
