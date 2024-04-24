fl <- "data/raw_data/41591_2015_BFnm3954_MOESM10_ESM.xlsx"

subsetDrugName <- function(dr)
{
  rt <- c()
  for(di in strsplit(dr, "\\+")[[1]])
  {
    di <- trimws(di)
    tr <- strsplit(di, "")[[1]]
    trx<-sprintf("%s%s%s%s", tr[1],tr[2], tr[(length(tr)-1)], tr[length(tr)])
    rt <- c(rt, trx)
  }
  return(paste0(rt, collapse = "."))
}

pexp <- readxl::read_xlsx(fl, sheet = 4)
pexp <- data.frame(pexp)
pexp$Treatment <- gsub("\"", "",pexp$Treatment)

mt <- unique(pexp[, c("Model", "Treatment")])

expLST <- list()
for(i in 1:nrow(mt))
{
  w <- pexp[pexp$Model==mt$Model[i] & pexp$Treatment==mt$Treatment[i], ]
  if(nrow(w)>1)
  {
  md <- gsub("-", ".", mt$Model[i])
  dr <- subsetDrugName(mt$Treatment[i])
  w$model.id <- paste0(md, ".", dr)

  wl <- length(w$Days.Post.T0)
  dd <- w$Days.Post.T0[2:wl] > w$Days.Post.T0[1:(wl-1)]
  if(any(dd==FALSE))
  {
    splitIndx <- !c(TRUE,dd)
    sp <- split(1:nrow(w), cumsum(splitIndx))
    for(j in names(sp))
    {
      wr <- w[sp[[j]], ]
      wr$model.id <- paste0(wr$model.id, ".", j)
      expLST[[length(expLST)+1]] <- wr
    }
  } else { expLST[[length(expLST)+1]] <- w }
  }
}

expDF <- do.call(rbind, expLST)
colnames(expDF) <- c("Model", "Tumor.Type", "drug", "volume", "body.weight",
                     "time", "TVol.Difference", "BW.Difference", "model.id")



##--------------------------
modelin <- unique(expDF[, c("model.id", "Tumor.Type", "Model", "drug")])
modelin <- modelin[complete.cases(modelin), ]
colnames(modelin) <- c("model.id", "tissue", "patient.id", "drug")

tissueName <- c("Gastric Cancer", "Colorectal Cancer", "Pancreatic Ductal Carcinoma",
                "Breast Cancer", "Non-small Cell Lung Carcinoma", "Cutaneous Melanoma")
names(tissueName) <- c("GC", "CRC", "PDAC", "BRCA", "NSCLC", "CM")

modelin$tissue.name <- tissueName[modelin$tissue]

experiment <- expDF[, c("model.id", "drug", "time", "volume", "body.weight")]
unqMod <- intersect(unique(experiment$model.id), modelin$model.id)

experiment <- experiment[experiment$model.id %in% unqMod, ]
modelin <- modelin[modelin$model.id %in% unqMod,]

write.csv(experiment, "data/raw_data/expriment.csv", row.names = FALSE)

write.csv(modelin[,c("model.id", "tissue", "tissue.name", "patient.id", "drug")],
          "data/raw_data/model_info.csv", row.names = FALSE)

##-----------------------------

molfl <- list(RNASeq="data/raw_data/molProf_RNASeq.rds",
              cnv="data/raw_data/molProf_cnv.rds",
              mutation="data/raw_data/molProf_mutation.rds")
m2bio <- data.frame()
for(mdt in names(molfl))
{
  moldf <- readRDS(molfl[[mdt]])
  m2bio <- rbind(m2bio, data.frame(biobase.id = sampleNames(moldf),
                                   mDataType=mdt))
}


##-----------------------------------------
modelin$biobase.id <- modelin$patient.id
m2bioMap <- merge(x = modelin[, c("model.id", "biobase.id")],
                  y = m2bio, by = "biobase.id", all.x = TRUE)
m2bioMap <- m2bioMap[complete.cases(m2bioMap), ]
m2bioMap <- m2bioMap[,c("model.id", "biobase.id", "mDataType")]

###-------------------------
microarray <- readRDS("data/raw_data/GSE78806_Normalized_Final_geneName.Rda")

micrArMap <- merge(x = modelin[, c("model.id", "patient.id")],
                  y = pData(microarray), by = "patient.id", all = TRUE)
micrArMap <- micrArMap[complete.cases(micrArMap), ]
micrArMap$mDataType <- "microArray"

m2bioMap <- rbind(m2bioMap, micrArMap[,c("model.id", "biobase.id", "mDataType")])
##-------------------------------------
write.csv(m2bioMap, "data/raw_data/modToBiobaseMap.csv", row.names = FALSE)


#modToBiobaseMap <- read.csv("data/raw_data/old/modToBiobaseMap.csv")
#drugInfo <- read.csv("data/raw_data/old/drug_info.csv")
#model <- read.csv("data/raw_data/old/model_info.csv")
#experiment <- read.csv("data/raw_data/old/expriment.csv")
