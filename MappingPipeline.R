require(dplyr)
require(qtl)
require(stringr)
require(ggplot2)

source("~/LinkageMapping/Reports/LinkageReportFunctions.R")

#---------------Pick your poison-----------------#

# pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv")

# pheno <- read.csv("~/Downloads/RIAILs0_complete.csv")

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_ForMapping.csv")

#------------------Process----------------------#

reduced.pheno <- pheno[!is.na(pheno$strain),]
reduced.pheno$id <- str_split_fixed(reduced.pheno$strain, "QX", 2)[,2]

trait <- reduced.pheno %>% group_by(id, drug) %>% summarise_each(funs(mean)) %>% select(-X) %>% filter(id!="")

trait2 <- trait %>% group_by(id, drug) %>% do(renameCols(.))

trait3 <- trait2 %>% group_by(id) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>%
    filter(id != "") %>% select(id, contains("\\."), -contains("sorted"))

trait3 = trait3[,-which(unlist(lapply(trait3, function(x){all(is.na(x))})))]

#-------------------------End processing, start mapping-------------------------------#

load("~/LinkageMapping/N2xCB4856_RIAILs_Rqtlfiles.RData")

N2xCB4856.cross <- calc.genoprob(N2xCB4856.cross, step=0)

N2xCB4856.cross$pheno <- data.frame(cbind(rep = sapply(mergePheno2(N2xCB4856.cross, trait3)$RILname, function(x){x = as.character(x); tryCatch({substr(as.character(x), nchar(x), nchar(x))}, error = function(err){return(NA)})}), mergePheno2(N2xCB4856.cross, trait3)))

#---------Map separate three sets--------#

system.time(set1Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 1)})))

system.time(set2Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 2)})))

system.time(set3Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 3)})))


#---------Recombine--------#

totalMap <- data.frame(do.call(rbind, list(set1Map, set2Map, set3Map)))
totalMap$lod <- as.numeric(as.character(totalMap$lod))
masterMap <- totalMap %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod = sum(lod))

#---------Calculate all of the thresholds of significance--------#

thresholds <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){structperm(cross=N2xCB4856.cross, pheno.col=x, model="np", perms=1000)}))


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
    plot2 = ggplot(traitMap, aes(x=pos, y=lod)) + geom_line(size=1) + geom_abline(intercept=thresholds[thresholds$condition==condition & thresholds$trait==trait,"threshold"], slope=0, colour = "red") + facet_grid(.~chr) + ggtitle(title)
    ggsave(plot = plot1, filename = fileName1, width = 8, height = 4.5)
    try(ggsave(plot = plot2, filename = fileName2, width = 11, height = 4.5))
}

#---------Save--------#
filename <- paste0(pheno$experiment[1], pheno$round[1])
write.csv(thresholds, file.path("~/Dropbox/HTA/Results/ProcessedData", paste0(filename, "_thresholds.csv")), row.names=FALSE)
write.csv(masterMap, file.path("~/Dropbox/HTA/Results/ProcessedData", paste0(filename, "_MasterMap.csv")), row.names=FALSE)
