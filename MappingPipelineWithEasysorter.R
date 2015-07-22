devtools::install_github("AndersenLab/linkagemapping")

library("linkagemapping")

data("N2xCB4856cross")

# Pick your poison
# pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs1_processed.rds")
pheno <- readRDS("~/Dropbox/AndersenLab/LabFolders/PastMembers/Tyler/ForTrip/RIAILs2_processed.rds")


N2xCB4856cross$pheno <- mergepheno(N2xCB4856cross, pheno)
map <- fsearch(N2xCB4856cross, 1000)
annlods <- annotate_lods(map, N2xCB4856cross)

test <- scanonefs(N2xCB4856cross, phenocol = 6, nperm = 10)
