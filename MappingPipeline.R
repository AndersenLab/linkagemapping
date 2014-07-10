require(dplyr)
require(qtl)
require(stringr)
require(ggplot2)

#--------------------Functions---------------------#

map <- function(x, set){
    if(set==1){
        newCross <- subset(N2xCB4856.cross, ind=as.character(N2xCB4856.cross$pheno$rep)=="A")
    } else if(set == 2){
        newCross <- subset(N2xCB4856.cross, ind=as.character(N2xCB4856.cross$pheno$rep) %in% c("B", "C", "D"))
    } else if(set==3){
        newCross <- subset(N2xCB4856.cross, ind=as.numeric(N2xCB4856.cross$pheno$set)==2)
    } else if(set==4){
        newCross <- subset(N2xCB4856.cross, ind=as.numeric(N2xCB4856.cross$pheno$set)==1)
    } else if(set==5){
        newCross <- subset(N2xCB4856.cross, ind=as.numeric(N2xCB4856.cross$pheno$set)==2)
    } else if(set==6){
        newCross <- N2xCB4856.cross
    }
    if (!all(is.na(newCross$pheno[,x]))) {
        newCross$pheno[,x] <- as.numeric(as.character(newCross$pheno[,x]))
        data <- scanone(newCross, pheno.col=x, model="np")
        data$marker <- rownames(data)
        data <- as.data.frame(data)
        colName <- colnames(N2xCB4856.cross$pheno)[x]
        condition <- str_split(colName, "\\.", 2)[[1]][1]
        trait <- str_split(colName, "\\.", 2)[[1]][2]
        df <- as.data.frame(cbind(condition, trait, data))
        return(df)
    }
}


perm <- function(x, set){
    if(set==1){
        newCross <- subset(N2xCB4856.cross, ind=as.character(N2xCB4856.cross$pheno$rep)=="A")
    } else if(set == 2){
        newCross <- subset(N2xCB4856.cross, ind=as.character(N2xCB4856.cross$pheno$rep) %in% c("B", "C", "D"))
    } else if(set==3){
        newCross <- subset(N2xCB4856.cross, ind=as.numeric(N2xCB4856.cross$pheno$set)==2)
    try({if (!all(is.na(newCross$pheno[,x]))) {
        threshold = summary(scanone(newCross, pheno.col = x, model="np", n.perm = 1000))[1]
        colName <- colnames(N2xCB4856.cross$pheno)[x]
        condition <- str_split(colName, "\\.", 2)[[1]][1]
        trait <- str_split(colName, "\\.", 2)[[1]][2]
        df <- as.data.frame(cbind(condition, trait, threshold))
        return(df)
    }})
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

# pheno <- read.csv("~/Downloads/RIAILs0_complete.csv")

# pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_complete.csv")

#------------------Process----------------------#

reduced.pheno <- pheno[!is.na(pheno$strain),]
reduced.pheno$id <- str_split_fixed(reduced.pheno$strain, "QX", 2)[,2]

trait <- reduced.pheno %>% group_by(id, drug) %>% summarise_each(funs(mean))

trait2 <- trait %>% group_by(id, drug) %>% do(renameCols(.))

trait3 <- trait2 %>% group_by(id) %>% summarise_each(funs(mean(., na.rm=TRUE))) %>%
    filter(id != "") %>% select(id, contains("\\."), -contains("sorted"))

trait3 = trait3[,-which(unlist(lapply(trait3, function(x){all(is.na(x))})))]

#-------------------------End processing, start mapping-------------------------------#

load("~/LinkageMapping/N2xCB4856_RIAILs_Rqtlfiles.RData")

N2xCB4856.cross <- calc.genoprob(N2xCB4856.cross, step=0)

N2xCB4856.cross$pheno <- data.frame(cbind(rep = sapply(mergePheno2(N2xCB4856.cross, trait3)$RILname, function(x){x = as.character(x); tryCatch({substr(as.character(x), nchar(x), nchar(x))}, error = function(err){return(NA)})}), mergePheno2(N2xCB4856.cross, trait3)))

#---------Three sets--------#

system.time(set1Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 1)})))
system.time(set1Perm <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){perm(x, 1)})))

system.time(set2Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 2)})))
system.time(set2Perm <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){perm(x, 2)})))

system.time(set3Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 3)})))
system.time(set3Perm <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){perm(x, 3)})))

totalMap <- data.frame(do.call(rbind, list(set1Map, set2Map, set3Map)))
totalMap$lod <- as.numeric(as.character(totalMap$lod))
masterMap <- totalMap %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod = sum(lod))

set1Perm <- read.csv("~/LinkageMapping/set1Perm.csv")
set2Perm <- read.csv("~/LinkageMapping/set2Perm.csv")
set3Perm <- read.csv("~/LinkageMapping/set3Perm.csv")

allPerm <- data.frame(do.call(rbind, list(set1Perm, set2Perm, set3Perm)))
allPerm$threshold <- as.numeric(as.character(allPerm$threshold))
thresholds <- allPerm %>% group_by(condition, trait) %>% summarise(threshold = sum(threshold))
    
# #---------Two sets----------#
# 
# system.time(set1Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 4)})))
# system.time(set1Perm2group <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){perm(x, 4)})))
# 
# system.time(set2Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 5)})))
# system.time(set2Perm2group <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){perm(x, 5)})))
# 
# allPerm2$threshold <- as.numeric(as.character(allPerm2$threshold))
# thresholds2 <- allPerm2 %>% group_by(condition, trait) %>% summarise(threshold = sum(threshold))
# 
# # #------One set--------#
# 
# system.time(masterMap <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 6)})))

################
################
################    2
################
################
################

# data = N2xCB4856.cross$pheno
# data$group <- ifelse(data$set==2, 3, ifelse(data$rep=="A", 1, ifelse(data$rep %in% c("B", "C", "D"), 2, NA)))
# 
# ggplot(data, aes(x = factor(group), y = control.f.ad)) + geom_violin() + geom_boxplot(width=.1)
# 
# ggplot(data, aes(x = factor(group), y = control.q10.EXT)) + geom_violin() + geom_boxplot(width=.1)

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

write.csv(masterMap, file.path("~/Dropbox/HTA/Results/ProcessedData", "RIAILs0_MasterMap.csv"), row.names=FALSE)
