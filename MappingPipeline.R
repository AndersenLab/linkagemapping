options(echo=TRUE)

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

structperm <- function(cross, chr, pheno.col=1, model=c("normal","binary","2part","np"), perms,
                       method=c("em","imp","hk","ehk","mr","mr-imp","mr-argmax"),
                       addcovar=NULL, intcovar=NULL, weights=NULL,
                       use=c("all.obs", "complete.obs"), upper=FALSE,
                       ties.random=FALSE, start=NULL, maxit=4000,
                       tol=1e-4, verbose){
    
    print(paste0("Pheno.col: ", pheno.col))
    
    lastpheno <- which(colnames(cross$pheno)=="rep")
    newcol <- lastpheno +1
    res <- matrix(nrow = perms, ncol = 1)
    temp1 <- subset(cross, ind = (cross$pheno[,lastpheno] == "A"))
    temp2 <- subset(cross, ind = (cross$pheno[,lastpheno] %in% c("B", "C", "D")))
    temp3 <- subset(cross, ind = (is.na(cross$pheno[,lastpheno])))
    
    for(i in 1:perms){
        if(i %% 100 == 0){
            cat("Structured Permutation", i, "\n")
        }
        temp1$pheno[,newcol] <- sample(temp1$pheno[,pheno.col])
        temp2$pheno[,newcol] <- sample(temp2$pheno[,pheno.col])
        temp3$pheno[,newcol] <- sample(temp3$pheno[,pheno.col])
        scan1 <- scanone(temp1, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        scan2 <- scanone(temp2, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        scan3 <- scanone(temp3, chr = chr, pheno = newcol, model = model, method = method, addcovar = addcovar, intcovar = intcovar, weights= weights, use = use, upper = upper, ties.random = ties.random, start = start, maxit = maxit, tol = tol)
        all <- scan1+scan2+scan3
        res[i,] <- max(all[,3])
        
    }
    
    rownames(res) <- 1:perms
    colName <- colnames(cross$pheno)[pheno.col]
    condition <- str_split(colName, "\\.", 2)[[1]][1]
    trait <- str_split(colName, "\\.", 2)[[1]][2]
    return(data.frame(cbind(rownames=NULL, condition, trait, threshold=quantile(res, probs=.95)[1])))
    
}

#---------------Pick your poison-----------------#

# pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs0_complete_simple.csv")

# pheno <- read.csv("~/Downloads/RIAILs0_complete.csv")

pheno <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_complete_simple.csv")

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

system.time(set2Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 2)})))

system.time(set3Map <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){map(x, 3)})))



totalMap <- data.frame(do.call(rbind, list(set1Map, set2Map, set3Map)))
totalMap$lod <- as.numeric(as.character(totalMap$lod))
masterMap <- totalMap %>% group_by(condition, trait, chr, pos, marker) %>% summarise(lod = sum(lod))

thresholds <- do.call(rbind, lapply(6:ncol(N2xCB4856.cross$pheno), function(x){structperm(cross=N2xCB4856.cross, pheno.col=x, model="np", perms=1000)}))

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


write.csv(thresholds, file.path("~/Dropbox/HTA/Results/ProcessedData", "RIAILs0_thresholds.csv"), row.names=FALSE)
write.csv(masterMap, file.path("~/Dropbox/HTA/Results/ProcessedData", "RIAILs0_MasterMap.csv"), row.names=FALSE)
