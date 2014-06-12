require(dplyr)
require(qtl)
require(stringr)
require(ggplot2)

#--------------------Functions---------------------#

map <- function(x, set){
    data <- scanone(N2xCB4856.cross, pheno.col=x, model="np")
    data$marker <- rownames(data)
    data <- as.data.frame(data)
    colName <- colnames(N2xCB4856.cross$pheno)[x]
    condition <- str_split(colName, "\\.", 2)[[1]][1]
    trait <- str_split(colName, "\\.", 2)[[1]][2]
    df <- as.data.frame(cbind(condition, trait, data))
    return(df)
}

mergePheno2 <- function(cross, phenotype, set=NULL){
    if(!is.null(set)){
        phenotype = phenotype %>% filter(set == set)
    }
    cross$pheno$id <- as.numeric(cross$pheno$id)
    phenotype$id <- as.numeric(phenotype$id)
    cross$pheno <- left_join(cross$pheno, phenotype, by="id")
    cross$pheno <- cross$pheno[order(cross$pheno$id),]
    return(cross$pheno)
}

renameCols <- function(x){
    colnames(x)[which(colnames(x) == "n"):ncol(x)] <- paste0(x$drug[1], ".", colnames(x)[which(colnames(x) == "n"):ncol(x)])
    return(x)
}

#---------------Pick your poison-----------------#

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete.csv")

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_complete.csv")

#------------------Process----------------------#

reduced.pheno <- pheno[!is.na(pheno$strain),]
reduced.pheno$id <- str_split_fixed(reduced.pheno$strain, "QX", 2)[,2]

trait <- reduced.pheno %>% group_by(id, drug) %>% summarise_each(funs(mean))

trait2 <- trait %>% group_by(id, drug) %>% do(renameCols(.))

trait3 <- trait2 %>% group_by(id) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>%
    filter(id != "") %>% select(id, contains("\\."), -contains("sorted"))

trait3 = trait3[,-which(unlist(lapply(trait3, function(x){all(is.na(x))})))]

#-------------------------End processing, start mapping-------------------------------#

load(url("http://groups.molbiosci.northwestern.edu/andersen/Data/N2xCB4856_RIAILs_Rqtlfiles.RData"))

N2xCB4856.cross$pheno <- mergePheno2(N2xCB4856.cross, trait3)


set1Time = system.time(set1Map <- do.call(rbind, lapply(5:ncol(N2xCB4856.cross$pheno), function(x){map(x, 1)})))

set2Time = system.time(set2Map <- do.call(rbind, lapply(5:ncol(N2xCB4856.cross$pheno), function(x){map(x, 2)})))






system.time(masterMap <- do.call(rbind, lapply(5:ncol(N2xCB4856.cross$pheno), map)))

masterMap %>% filter(condition=="topotecan", trait=="resid.col") %>% ggplot(., aes(x=pos, y=lod)) + geom_point() + facet_grid(.~chr)

unique(masterMap$condition)

write.csv(masterMap, file.path("~/Dropbox/HTA/Results/ProcessedData", "RIAILs1_MasterMap.csv"), row.names=FALSE)

