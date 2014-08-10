# Load required packages
require(dplyr)
require(qtl)
require(stringr)
require(ggplot2)

# Source the functions
source("~/LinkageMapping/LinkageMappingFunctions.R")

#---------------Pick your poison-----------------#
# Set your phenotype data file here

# pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv")

# pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete.csv")

#------------------Process----------------------#

# Remove wash wells (wells where the strain is NA)
reduced.pheno <- pheno[!is.na(pheno$strain),]

# Get the id number for mergePheno2 function, done by simple string split
reduced.pheno$id <- str_split_fixed(reduced.pheno$strain, "QX", 2)[,2]

# Summarize all traits by strain id, drug and assay (really only necessary for RIAILs0, but doesn't hurt other RIAILs experiments)
trait <- reduced.pheno %>% group_by(id, drug, assay) %>% summarise_each(funs(mean)) %>% filter(id!="")

# Rename all of the columns in the format "drug.trait"
trait2 <- trait %>% group_by(id, drug, assay) %>% do(renameCols(.))

# Collapse down to one row per strain
trait3 <- trait2 %>% group_by(id, assay) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>%
    filter(id != "") %>% select(id, assay, contains("\\."), -contains("sorted"))

# Remove columns that are all NA 
trait3 = trait3[,-which(unlist(lapply(trait3, function(x){all(is.na(x))})))]
trait3$id <- as.integer(as.numeric(as.character(trait3$id)))

#-----------------------End processing, start mapping--------------------------#

# Load in the rqtl files 
load("~/LinkageMapping/N2xCB4856_RIAILs_Rqtlfiles.RData")

# Remove interpolated SNPs
N2xCB4856.cross <- calc.genoprob(N2xCB4856.cross, step=0)

N2xCB4856.cross$pheno$id <- as.integer(as.numeric(as.character(N2xCB4856.cross$pheno$id)))

# Reorder by trait id
trait3 <- trait3[order(trait3$id),] 

# Extract replicate values for Matt Rockman's set and merge in the phenotype values to the cross object
N2xCB4856.cross$pheno <- data.frame(cbind(rep = sapply(mergePheno2(N2xCB4856.cross, trait3)$RILname, function(x){x = as.character(x); tryCatch({substr(as.character(x), nchar(x), nchar(x))}, error = function(err){return(NA)})}), mergePheno2(N2xCB4856.cross, trait3)))

# Impute missing genotypes
impcross = argmax.geno(N2xCB4856.cross)

# pdata.01 = extractScaledPhenotype(impcross)

# Use this one
pdata.01s = extractScaledPhenotype(impcross, set=impcross$pheno$set, setCorrect=T)

gdata = extractGenotype(impcross)
n.pheno  = countStrainsPerTrait(pdata.01s)

# Calculate avg RIL relatedness (for VC models) 
# additive
A = A.mat(gdata, shrink=FALSE)/2
# epistatic
AA = A*A

# build mindex.split (list of marker indices split by chromosome)
marker.chr.cnt=sapply(impcross$geno, function(x) ncol(x$argmax))
chr=rep(names(marker.chr.cnt), marker.chr.cnt)
mindex.split= split(seq_along(chr),chr)
# to convert from genome marker index to chromosome marker index 
chr.mindex.offset = sapply(mindex.split, min)-1

LODS.01       = get.LOD.by.COR(n.pheno, pdata.01s, gdata, doGPU=F)
LODS.01s      = LODmatrix.2.scanone(LODS.01, impcross)
peaklist.01   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.01) 
# change number of permutations
LODS.01.FDR   = getPeakFDR(peaklist.01$chr.peaks.lod, pdata.01s, gdata, 50, doGPU=F)

# Second argument is sig threshold, should not be hard coded (from LODS.01.FDR)
peakArray.01  = getPeakArray(peaklist.01, 3.38)

#second jump run here, run it on only the traits with significant QTL that were detected in round 1
#pdata.02      = getPhenoResids(pdata.01s, gdata, peakArray.01) 
#n.pheno2 = countStrainsPerTrait(pdata.02)
#LODS.02       = get.LOD.by.COR(n.pheno2, pdata.02, gdata, doGPU=T)
#LODS.02s      = LODmatrix.2.scanone(LODS.02, impcross, LODS.01s)
#peaklist.02   = getChrPeaks(mindex.split, chr.mindex.offset, LODS.02) 
#LODS.02.FDR   = getPeakFDR(peaklist.02$chr.peaks.lod, pdata.02, gdata, 10, doGPU=T)
#peakArray.02  = getPeakArray(peaklist.02, 5)
# only paste together unique indices
#peakArray.02  = rbind(peakArray.01, peakArray.02)

peaksFDR05=data.frame(trait=as.character(colnames(pdata.01)[peakArray.01[,1]]), marker.index=peakArray.01[,2])
peakList = split(peakArray.01$markerIndex, peakArray.01$trait)

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
    am = lm(paste('pdata.01[,' , trait.num, ']', '~', (aq), '-1'))
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
pset= impcross$pheno$set
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











for(col in 294:ncol(LODS.01s)){
    print(col)
    title <- colnames(LODS.01s)[col]
    fileName = paste0("~/LinkagePlots/", gsub("\\.", "-", title), "_map.pdf")
    plot2 = ggplot(LODS.01s, aes(x=pos, y=LODS.01s[,col])) + geom_line(size=1) + facet_grid(.~chr) + xlab("position") + ylab("LOD") + ggtitle(title) + geom_hline(yintercept=3.38, colour="red", linetype="dashed")
    try(ggsave(plot = plot2, filename = fileName, width = 11, height = 4.5))
}
















###For Josh Bloom
df=N2xCB4856.cross$pheno
df$set=ifelse(as.character(df$rep)=="A", 1, ifelse(as.character(df$rep) %in% c("B", "C", "D"), 2, 3))
df$set[is.na(df$set)]=3
N2xCB4856.cross$pheno=df

save(N2xCB4856.cross, file="~/Dropbox/Coding/joshbloom/RIAILs0_CrossObject_Old.RData")









#---------Map separate three sets--------#

system.time(totalMap <- do.call(rbind, mclapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(N2xCB4856.cross, x)})))


#---------Recombine--------#

totalMap$lod <- as.numeric(as.character(totalMap$lod))
masterMap <- totalMap

#---------Calculate all of the thresholds of significance--------#

# thresholds <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){structperm(cross=N2xCB4856.cross, pheno.col=x, model="np", perms=1000)}))


#---------Save--------#
filename <- paste0(pheno$experiment[1], pheno$round[1])
# write.csv(thresholds, file.path("~/Dropbox/HTA/Results/ProcessedData", paste0(filename, "_thresholds.csv")), row.names=FALSE)
write.csv(masterMap, file.path("~/Dropbox/HTA/Results/ProcessedData", paste0(filename, "_MasterMap.csv")), row.names=FALSE)



#---------Plot--------#
traits = unique(masterMap[,1:2])

for(row in 1:nrow(traits)){
    print(row)
    colName = paste0(traits[row,1], ".", traits[row,2])
    plot1 = ggplot(trait3, aes_string(x = colName)) + geom_histogram()
    title = paste0("Condition: ", traits[row,1], ", Trait: ", traits[row,2])
    condition = traits[row,1]
    trait = traits[row,2]
    fileName1 = paste0("~/LinkagePlots/", condition, "_", gsub("\\.", "-", trait), "_hist.pdf")
    fileName2 = paste0("~/LinkagePlots/", condition, "_", gsub("\\.", "-", trait), "_map.pdf")
    traitMap = masterMap[masterMap$condition==condition & masterMap$trait==trait,]
    plot2 = ggplot(traitMap, aes(x=pos, y=lod)) + geom_line(size=1) + facet_grid(.~chr) + ggtitle(title) + geom_hline(yintercept=4.5, colour="gray", linetype="dashed")#+ geom_abline(intercept=thresholds[thresholds$condition==condition & thresholds$trait==trait,"threshold"], slope=0, colour = "red") + facet_grid(.~chr) + ggtitle(title)
    ggsave(plot = plot1, filename = fileName1, width = 8, height = 4.5)
    try(ggsave(plot = plot2, filename = fileName2, width = 11, height = 4.5))
}


