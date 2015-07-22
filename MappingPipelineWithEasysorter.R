library(tidyr)
library(dplyr)
library(stringr)
library(foreach)

# Load in the rqtl files 
load("~/Dropbox/AndersenLab/RCode/Linkage mapping/N2xCB4856_RIAILs_Rqtlfiles.RData")

pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs1_processed.rds")
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs2_processed.rds")
N2xCB4856.cross$pheno <- mergepheno(N2xCB4856.cross, pheno)
test <- fsearch(N2xCB4856.cross, 10)

debug(annotate_lods)
annlods <- annotate_lods(test, N2xCB4856.cross)
