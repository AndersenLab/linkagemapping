require(stringr)


data = read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_complete.csv")
load(url("http://groups.molbiosci.northwestern.edu/andersen/Data/N2xCB4856_RIAILs_Rqtlfiles.RData"))
data$id = data.frame(t(data.frame(str_split(data$strain, "X"))))[,2]
N2xCB4856.cross$pheno 

pheno2 <- mergePheno(N2xCB4856.cross, data)
