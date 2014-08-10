################################################################
################################################################
# 
# 1. RIAILs1 sets a, b, c, and d (the first two weeks) had the inappropriate concentrations of methotrexate and albendazole. We re-did RIAILs1a and 1c, so only 1b and 1d will have to have those two drugs removed.
# 
# 2. RIAILs1e and 1f used daunorubicin from 4ºC storage instead of -20ºC storage. We should look for discrepancies.
# 
# 3. RIAILs1g, h, i, j used a new dilution of daunorubicin. We should look for discrepancies.
# 
# 4. RIAILs1g used 4ºC etoposide instead of -20ºC etoposide. We should look for discrepancies.
# 
# 5. RIAILs1h, i, j used a new dilution of etoposide. We should look for discrepancies.
# 
# REDOs:
#     
# RIAILs1a (old 1) should be replaced with RIAILs1k (REDO 1)
# RIAILs1c (old 3) should be replaced with RIAILs1l (REDO 3)
# RIAILs1g (old 7) should be replaced with RIAILs1m (REDO 7)
# 
################################################################
################################################################



library(dplyr)

#Read in the data
data <- read.csv("~/Dropbox/HTA/Results/ProcessedData/RIAILs1_complete_simple.csv") %>% select(-1)

#filter out assay a, c, and g, which were redone as k, l, and m, respectively
data2 <- data %>% filter(assay!="a", assay!="c", assay!="g")

#Remove wrong concentration of methotrexate and albendazole, addresses point 1
data3 <- data2[!(grepl("methotrexate", data2$drug) & (data2$assay!="b" | data2$assay!="d")),]
data4 <- data3[!(grepl("albendazole", data3$drug) & (data3$assay!="b" | data3$assay!="d")),]

#filter out missing plates
data5 <- data4 %>% filter(drug!="missing")

write.csv(data5, "~/Dropbox/HTA/Results/ProcessedData/RIAILs1_ForMapping.csv")










