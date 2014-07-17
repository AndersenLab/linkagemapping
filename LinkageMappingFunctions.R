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