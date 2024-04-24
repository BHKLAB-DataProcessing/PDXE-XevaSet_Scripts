library(Xeva)
pdxe=readRDS("~/CXP/Xeva_dataset/data/XevaObjects/XevaSets/data/PDXE.rds")


expDF <- data.frame()
for(i in pdxe@experiment)
{
  df = i@data[, ]
  df$model.id <- i@model.id
  df$drug <- i@drug$join.name
  df <- df[, c("model.id", "drug", "time", "volume", "body.weight")]
  expDF <- rbind(expDF, df)
}

write.csv(expDF, "data/raw_data/expriment.csv", row.names = F)

mi <- unique(expDF[,c("model.id", "drug")])
rownames(mi) <- mi$model.id

minf <- pdxe@model
minf$drug <-  mi[minf$model.id, "drug"]

write.csv(minf, "data/raw_data/model_info.csv", row.names = F)

write.csv(pdxe@drug, "data/raw_data/drug_info.csv", row.names = F)

write.csv(pdxe@modToBiobaseMap[pdxe@modToBiobaseMap$mDataType!="RNASeq_FPKM",], 
          "data/raw_data/modToBiobaseMap.csv", row.names = F)

for(dt in names(pdxe@molecularProfiles))
{
  saveRDS(pdxe@molecularProfiles[[dt]], 
          sprintf("data/raw_data/molProf_%s.rds", dt))
}


