load("~/Desktop/linkagemapping/data/AF16xHK104cross.rda")
load("~/Desktop/linkagemapping/data/AF16xHK104markers.rda")
load("~/Desktop/linkagemapping/data/AllCBfosmids.rda")
load("~/Desktop/linkagemapping/data/AllN2fosmids.rda")
load("~/Desktop/linkagemapping/data/N2xCB4856cross.rda")
load("~/Desktop/linkagemapping/data/N2xCB4856markers.rda")
load("~/Desktop/linkagemapping/data/N2xLSJ2cross.Rda")
load("~/Desktop/linkagemapping/data/N2xLSJ2markers.rda")
load("~/Desktop/linkagemapping/data/RIAILgenotypes.rda")
load("~/Desktop/linkagemapping/data/RIAILmarkerconversion.rda")
load("~/Desktop/linkagemapping/data/eQTLpeaks.rda")
load("~/Desktop/linkagemapping/data/indels.rda")
load("~/Desktop/linkagemapping/data/probe_info.rda")


devtools::use_data(AF16xHK104markers, AF16xHK104cross, AllCBfosmids, AllN2fosmids, N2xCB4856cross, N2xCB4856markers, N2xLSJ2cross, N2xLSJ2markers, RIAILgenotypes, RIAILmarkerconversion, eQTLpeaks, indels, probe_info,
                   internal = TRUE, overwrite = TRUE)

