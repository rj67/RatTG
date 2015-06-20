#' RMA normalize TG-Gate affymetrix microarray dataset
#'
#' Analyze the TG-Gate dataset under the supplied folder. Suppose the .CEL files are in "celfiles" subfolder. The experimental metadata is in "Attribute.tsv" file.
#' @param anno_db Annotation database for microarray Defaults to "rat2302"

# process expression set in the folder
RMATGEset <- function(data_dir, anno_db = "rat2302"){
  attr_file <- paste(data_dir, "Attribute.tsv", sep="/")
  celfile_dir <- paste(data_dir, "celfiles", sep="/")
  #######################################################
  # All sample pheno data 
  #######################################################
  array_design <- do.call(rbind, lapply(attr_file, function(x) read.delim(x, stringsAsFactors=F, fileEncoding="ISO-8859-1")))
  carray_design <- droplevels(subset(array_design, !is.na(ARR_DESIGN)))
  array_design$SACRI_PERIOD <- gsub(" hr", "h", array_design$SACRI_PERIOD) 
  array_design$SACRI_PERIOD <- gsub(" day", "d", array_design$SACRI_PERIOD) 
  array_design$time <- sapply(array_design$SACRI_PERIOD, function(x) ifelse(substr(x, nchar(x), nchar(x))=="d", as.numeric(substr(x, 1, nchar(x)-1))*24, as.numeric(substr(x, 1, nchar(x)-1))))
  #arrange SACRI_PERIOD factor level by time
  dummy<- arrange(array_design[c("time","SACRI_PERIOD")], time)
  array_design$SACRI_PERIOD <- factor(array_design$SACRI_PERIOD, levels=dummy[!duplicated(dummy$time),]$SACRI_PERIOD)
  array_design$DOSE_LEVEL <- factor(array_design$DOSE_LEVEL, levels=c("Control", "Low", "Middle", "High"))
  array_design$group <- apply(array_design[c("COMPOUND.Abbr.", "DOSE_LEVEL", "SACRI_PERIOD")], 1, function(x) paste0(x, collapse="."))
  array_design$ids <- apply(array_design[c("COMPOUND.Abbr.", "DOSE_LEVEL", "SACRI_PERIOD", "INDIVIDUAL_ID")], 1, function(x) paste0(x, collapse="."))
  rownames(array_design) <- array_design$ids
  if (class(array_design$BARCODE) == "numeric"){
    array_design$BARCODE <- paste("00", array_design$BARCODE, sep="")
  }
  #######################################################
  # read in the expression data, quality control
  #######################################################
  # library for affymetric
  library(affy)
  # brainarray
  #Data <- ReadAffy(cdfname="rat2302rnentrezg")
  # affy cdf
  Data <- ReadAffy(celfile.path=celfile_dir)
  ## quality control
  #library(yaqcaffy)
  #yqc <-yaqc(Data,verbose=TRUE)
  #write.table(show(yqc), file=paste(data_dir, "/yqc_summary.tsv", sep=""), quote=F, sep="\t")
  #pdf(paste(data_dir, "/yqc_summary.pdf", sep=""),width=7,height=6)
  #plot(yqc)
  #dev.off()
  #######################################################
  # RNA Normalization, MAS5 present call
  #######################################################
  eset.rma <- rma(Data)
  #eset.rma <- gcrma(Data)
  # add pheno info to sample name
  sampleNames(eset.rma) <- sapply(sampleNames(eset.rma), function(x) {barcode<-gsub(".CEL", "", x, fixed=T); array_design[array_design$BARCODE==barcode,]$ids}) 
  # phenodata 
  pData <- array_design[sampleNames(eset.rma), ][c("ids","COMPOUND_NAME", "COMPOUND.Abbr.","SIN_REP_TYPE", "SACRI_PERIOD", "DOSE",  "DOSE_LEVEL", "time", "group")]
  phenoData <- new("AnnotatedDataFrame", data=pData)
  # update eset
  eset.rma <- ExpressionSet(assayData=exprs(eset.rma), phenoData=phenoData, annotation=anno_db)
  ## MAS5 Normalization, 
  #eset.mas5 <- mas5(Data)
  ## use MAS5 to make Present/Absent call
  poa.call <- mas5calls(Data)
  sampleNames(poa.call) <- sapply(sampleNames(poa.call), function(x) {barcode=gsub(".CEL", "", x, fixed=T); array_design[array_design$BARCODE==barcode,]$ids}) 
  pData <- array_design[sampleNames(poa.call), ][c("COMPOUND_NAME","COMPOUND.Abbr.", "SIN_REP_TYPE", "SACRI_PERIOD", "DOSE",  "DOSE_LEVEL", "time", "group")]
  phenoData <- new("AnnotatedDataFrame", data=pData)
  poa.call <- ExpressionSet(assayData=exprs(poa.call), phenoData=phenoData, annotation=anno_db)
  #library(annotate)
  #Symbol <- getSYMBOL(featureNames(eset.rma), anno_db)
  #fData(eset.rma) <- data.frame(Symbol=Symbol)
  # remove probes that are not present in any
  eset.rma <- eset.rma[!apply(exprs(poa.call), 1, function(x) sum(x=="P"))==0, ]
  return(list("eset"=eset.rma, "poa"=poa.call))
}