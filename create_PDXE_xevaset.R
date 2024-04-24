library(Xeva)
library(Biobase)

model <- read.csv("data/raw_data/model_info.csv")
model$patient.id <- make.names(model$patient.id)

experiment <- read.csv("data/raw_data/expriment.csv")

drugInfo <- read.csv("data/raw_data/drug_info.csv")

control.drug <- "untreated"
drugs <- unique(model$drug)
drugs <- drugs[drugs!=control.drug]
##------create expriment design list

expDesign <- list()
for(p in unique(model$patient.id))
{
  for(d in drugs)
  {
    bt <- list(batch.name = sprintf("%s.%s", p,d),
               control = model$model.id[model$patient.id==p & model$drug==control.drug],
               treatment=model$model.id[model$patient.id==p & model$drug==d])
    if(length(bt$control)>0 | length(bt$treatment)>0)
    { expDesign[[bt$batch.name]] <- bt }
  }
}

###----------------
modToBiobaseMap <- read.csv("data/raw_data/modToBiobaseMap.csv")
modToBiobaseMap$biobase.id <- make.names(modToBiobaseMap$biobase.id)

##-----read mol data ---
RNASeq <- readRDS("data/raw_data/molProf_RNASeq.rds")
sampleNames(RNASeq) <- make.names(sampleNames(RNASeq))
for(cl in c("biobase.id", "patient.id"))
{ pData(RNASeq)[, cl] <- make.names(as.character(pData(RNASeq)[, cl])) }

cnv <- readRDS("data/raw_data/molProf_cnv.rds")
sampleNames(cnv) <- make.names(sampleNames(cnv))
for(cl in c("sampleID"))
{ pData(cnv)[, cl] <- make.names(as.character(pData(cnv)[, cl])) }

mutation <- readRDS("data/raw_data/molProf_mutation.rds")
sampleNames(mutation) <- make.names(sampleNames(mutation))
for(cl in c("sampleID"))
{ pData(mutation)[, cl] <- make.names(as.character(pData(mutation)[, cl])) }


microArray <- readRDS("data/raw_data/GSE78806_Normalized_Final_geneName.Rda")
sampleNames(microArray) <- make.names(sampleNames(microArray))
for(cl in c("biobase.id", "patient.id"))
{ pData(microArray)[, cl] <- make.names(as.character(pData(microArray)[, cl])) }


##-----check if ids are matching ----

#RNASeq <- RNASeq[, sampleNames(RNASeq) %in% modToBiobaseMap$biobase.id]
#cnv <- cnv[, sampleNames(cnv)%in%modToBiobaseMap$biobase.id]
#mutation <- mutation[, sampleNames(mutation)%in%modToBiobaseMap$biobase.id]

##==== create XevaSet ====
pdxe = createXevaSet(name="PDXE xevaSet",
                     model = model,
                            drug = drugInfo,
                            experiment = experiment,
                            expDesign = expDesign,
        molecularProfiles=list(RNASeq = RNASeq,
                               mutation=mutation,
                               cnv=cnv,
                               microArray=microArray),
                            modToBiobaseMap = modToBiobaseMap)

##===== set response =====
for(res in c("mRECIST", "slope", "AUC", "angle", "abc", "TGI"))
{
  pdxe <- setResponse(pdxe, res.measure = res, verbose=TRUE)
}

saveRDS(pdxe, "data/PDXE_xevaSet/PDXE_XevaSet_All.rds")


##=========subset by tissue type ====
mi <- modelInfo(pdxe)
for(tissue in c("Breast Cancer", "Colorectal Cancer", "Cutaneous Melanoma",
                "Gastric Cancer", "Non-small Cell Lung Carcinoma",
                "Pancreatic Ductal Carcinoma"))
{
  pdxe.tissue <- subsetXeva(pdxe, ids= mi$model.id[mi$tissue.name==tissue],
                            id.name = "model.id")
  saveRDS(pdxe.tissue, sprintf("data/PDXE_xevaSet/PDXE_%s_XevaSet.rds",
                               gsub(" ", "_", tissue)))
}

