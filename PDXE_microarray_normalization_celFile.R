library(affyio)
library(GEOquery)
library(affy)
library(gcrma)

BrainArrayPackageInfo<- function(celFilePath, version="19.0.0", organism="hs", annotationSource="entrezg", installPkg=FALSE)
{
  platform = cleancdfname(read.celfile.header(celFilePath,
                                              info = "full")$cdfName)
  platform = sub("cdf", "", platform)
  platform = sub("stv1", "st", platform)
  platform = sub("stv2", "st", platform)
  #pkgName = paste(platform, organism, annotationSource, "probe", sep = "")
  pkgName = paste(platform, organism, annotationSource, "cdf", sep = "")
  pkgFileName = paste(pkgName, "_", version, ".tar.gz", sep = "")
  pkgUrl = paste("http://mbni.org/customcdf/", version,
                 "/", annotationSource, ".download/", pkgFileName,
                 sep = "")

  if(installPkg==TRUE)
  {
    if( (pkgName %in% rownames(installed.packages())) == FALSE )
    { install.packages(pkgUrl, repos = NULL, type = "source") }
  }

  return(list(name=pkgName, url=pkgUrl, fileName=pkgFileName))
}

BrainArrayFixRowName <- function(datax, transpose=TRUE)
{
  datax = data.frame(datax)
  rownames(datax) = gsub("_at$",  "", rownames(datax))
  colnames(datax) = gsub(".cel$", "", colnames(datax))
  if(transpose==TRUE){datax = as.data.frame(t(as.matrix(datax)))}
  return(datax)
}


get_cel_files <- function(GSEid, downloadData=FALSE, baseDir = getwd())
{
  gseDir = path.expand(sprintf("%s/%s", baseDir, GSEid)) ## stupid R gives error in untar if ~/ is used
  if(downloadData==TRUE)
  {
    cat(sprintf("##----\nDownloading %s microarray data in\n##  %s \n##----\n", GSEid, gseDir ))
    getGEOSuppFiles(GSEid, makeDirectory = TRUE, baseDir = baseDir)
    untar(sprintf("%s/%s_RAW.tar",gseDir,GSEid), exdir=sprintf("%s/data", gseDir))
    celsGZ = list.files(sprintf("%s/data", gseDir), pattern = "CEL.gz", full.names=TRUE)
    sapply(celsGZ, gunzip,overwrite = TRUE)

  }
  cels = list.files(sprintf("%s/data",gseDir), pattern = ".cel$", full.names=TRUE)
  return(cels)
}



normalizeGSE_data <- function(GSEid, baseDir, downloadData=FALSE)
{

  celsFils = get_cel_files(GSEid, downloadData, baseDir) #[1:5]

  pkgInfo = BrainArrayPackageInfo(celsFils[1], "19.0.0", "hs", "ensg", installPkg=TRUE)

  library(pkgInfo$name, character.only = TRUE)

  raw.data=ReadAffy(verbose=TRUE, filenames=celsFils, cdfname=pkgInfo$name)

  datax = expresso(raw.data, bg.correct =FALSE, normalize.method = "quantiles",
                           normalize=TRUE, pmcorrect.method = "pmonly", summary.method = "medianpolish")

  #rmaD = exprs(data.rma.norm)
  datax = BrainArrayFixRowName(exprs(datax))

  outFilename = sprintf("%s/%s/%s_Normalized.Rda", baseDir,GSEid, GSEid)
  saveRDS(datax, file= outFilename)
  cat(sprintf("\n\n##----- Done for %s\n##----- File saved: %s\n\n", GSEid, outFilename))
}

baseDir = "~/CXP/XG/Data/Gao_2015_NatureMed/microArray"
#downloadData=TRUE  ## only once if data need to be downloaded
downloadData=FALSE
GSEid = "GSE78806"
normalizeGSE_data(GSEid, baseDir, downloadData)

if(1==1)
{
##---------------------------------------------------------------------------
##-------------------map IDs -------------------
library(BBmisc)

mdfx <- read.table("../Data/Gao_2015_NatureMed/microArray/GSE78806/E-GEOD-78806.sdrf.txt",
                  sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mdfx$Comment..ArrayExpress.FTP.file. <- NULL
mdfx$Term.Accession.Number.2 <- NULL
mdfx$Comment..Derived.ArrayExpress.FTP.file. <- NULL
mdfx$Term.Accession.Number.1 <- NULL
mdfx$Term.Accession.Number <- NULL
rownames(mdfx) <- sapply(mdfx$Source.Name, function(x)strsplit(x, " ")[[1]][1])


mdf <- mdfx[, c("Array.Data.File", "Characteristics..model.number.",
                "Characteristics..passage.", "FactorValue..primary.site.")]

colnames(mdf) <- c("array.id", "patient.id", "passage", "tumor.type")
mdf$patient.id <- paste0("X-", mdf$patient.id)

mdf$array.id <- sapply(strsplit(mdf$array.id, split = "_"), "[[", 1)
rownames(mdf) <- as.character(mdf$array.id)

mdf <- sortByCol(mdf, col = c("patient.id","passage","array.id"),
                 asc = c(TRUE, TRUE,TRUE))

datax <- readRDS("~/CXP/XG/Data/Gao_2015_NatureMed/microArray/GSE78806/GSE78806_Normalized.Rda")
colnames(datax) <- sapply(strsplit(colnames(datax), split = "_"), "[[", 1)
rownames(datax) <- sapply(strsplit(rownames(datax), split = "_"), "[[", 1)

#dfx <- list(data=datax, ids= mdf)
##saveRDS(dfx, file = "../Data/Gao_2015_NatureMed/microArray/GSE78806/GSE78806_Normalized_Final.Rda")

##------------------------------------------------------------------------------
###-----------------------------------------------------------------------------

library(Biobase)
library(biomaRt)

ensembl <- useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl",
                      #host = "www.ensembl.org",
                      #host = "http://uswest.ensembl.org",
                      GRCh=37)


if(1==2){
library(org.Hs.eg.db)
tryCatch({
       results <- select(org.Hs.eg.db, keys = gx[100:110],
                         columns=c("ENSEMBL","ENTREZID","SYMBOL","ALIAS","GENENAME"),
                         keytype="ENSEMBL")
     },error = function(err) { }
     )
}


#dx <- readRDS("~/CXP/XG/Data/Gao_2015_NatureMed/microArray/GSE78806/GSE78806_Normalized_Final.Rda")
#ddx <- dx$data
#dxExp <- Biobase::exprs(ddx)
dxExp <- t(datax)
gx <- unique(gsub("_at", "", rownames(dxExp)))
gid <- getBM(attributes=c('ensembl_gene_id', #'ensembl_transcript_id',
                          'hgnc_symbol', 'hgnc_id'),
             filters = 'ensembl_gene_id', values = gx, mart = ensembl)

gid <- unique(gid[, c("ensembl_gene_id",  "hgnc_symbol")])
gid <- gid[gid$hgnc_symbol!="", ]
gid <- gid[!duplicated(gid$ensembl_gene_id), ]

#gid$egid <- paste(gid$ensembl_gene_id, "_at", sep = "")
#rownames(gid) <- gid$egid
rownames(gid) <- as.character(gid$ensembl_gene_id)

commGene <- intersect(rownames(dxExp), rownames(gid))

dxNew <- dxExp[commGene, ]
rownames(dxNew) <- gid[commGene, "hgnc_symbol"]


featureDF <- gid[commGene, ]
rownames(featureDF) <- featureDF$hgnc_symbol
featuredata <- Biobase::AnnotatedDataFrame(data = featureDF)

colnames(mdf)[colnames(mdf)=="array.id"] <- "biobase.id"
phenodata   <- Biobase::AnnotatedDataFrame(data = mdf)


dxNew <- dxNew[rownames(featureDF), rownames(mdf)]
microAD <- Biobase::ExpressionSet(assayData=dxNew,
                                  phenoData=phenodata,
                                  featureData=featuredata)

#dw <- list(data=microAD, map=dx$ids)

saveRDS(microAD, file = "~/CXP/XG/Data/Gao_2015_NatureMed/microArray/GSE78806/GSE78806_Normalized_Final_geneName.Rda")



}
