
library(jsonlite)
library(XML)

list2df <- function(clinic){
  nms = NULL
  for(i in 1:length(clinic)){
    ci = unlist(clinic[[i]])
    nms = union(nms, names(ci))
  }
  length(nms)

  nms = sort(nms)

  clinicDF = NULL

  for(i in 1:length(clinic)){
    
    ci = unlist(clinic[[i]])
    vi = rep(NA, length(nms))
    vi[match(names(ci), nms)] = ci
    
    clinicDF = rbind(clinicDF, vi)
  }

  colnames(clinicDF) = nms
  clinicDF = as.data.frame(clinicDF, row.names=1:nrow(clinicDF),
                            stringsAsFactors=FALSE)
  clinicDF
  
}

df2df <- function(meta){
  
  nms = NULL
  for(i in 1:nrow(meta)){
    ci = unlist(meta[i,])
    nms = union(nms, names(ci))
  }
  length(nms)
  
  meta1 = NULL
  
  for(i in 1:nrow(meta)){
    ci = unlist(meta[i,])
    vi = rep(NA, length(nms))
    ni = intersect(names(ci), nms)
    vi[match(ni, nms)] = ci[match(ni, names(ci))]
    
    meta1 = rbind(meta1, vi)
  }
  
  colnames(meta1) = nms
  meta1 = as.data.frame(meta1, row.names=1:nrow(meta1))
  meta1
  
}

setwd("~/research/TCGA/COAD/GDC/original_data/gene_expression_meta")

# ----------------------------------------------------------------------------
# Read in meta information
# ----------------------------------------------------------------------------

meta = fromJSON("metadata.cart.2016-08-21T05_30_15.931767.json")
dim(meta)
meta[1:2,]

meta1 = df2df(meta)
dim(meta1)
meta1[1:2,]

write.table(meta1, file = "metadata_cart.txt", sep = "\t",
quote = FALSE, row.names = FALSE, col.names = TRUE)

# ----------------------------------------------------------------------------
# Read in clinical and biospecimen information
# ----------------------------------------------------------------------------

clinic = fromJSON("clinical.cart.2016-08-21T05_30_20.657982.json", simplifyDataFrame=FALSE)
length(clinic)
clinic[[1]]

clinic = list2df(clinic)
dim(clinic)
clinic[1:2,]

biosp = fromJSON("biospecimen.cart.2016-08-21T05_30_28.837995.json",
  simplifyDataFrame=FALSE)
length(biosp)
biosp[[1]]

biosp = list2df(biosp)
dim(biosp)
biosp[1:2,]

write.table(clinic, file = "metadata_clinic.txt", sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = TRUE)

write.table(biosp, file = "metadata_biospecimen.txt", sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = TRUE)

# ----------------------------------------------------------------------------
# prepare samples
# ----------------------------------------------------------------------------

cls = c("analysis.submitter_id", "analysis.analysis_id")
cls = c(cls, "analysis.input_files.file_name", "analysis.input_files.file_size")
cls = c(cls, "associated_entities.entity_id", "associated_entities.case_id")
cls = c(cls, "associated_entities.entity_submitter_id")

meta2 = data.frame(meta1[,cls], stringsAsFactors=FALSE)
dim(meta2)
meta2[1:2,]

barcode = as.character(meta2$associated_entities.entity_submitter_id)
barcode = strsplit(barcode, split="-")
table(sapply(barcode, length))
barcode = matrix(unlist(barcode), ncol=7, byrow=TRUE)
dim(barcode)
barcode[1:2,]

barcode = as.data.frame(barcode[,-1], stringsAsFactors=FALSE)
names(barcode) = c("tss", "participant", "sample", "portion", "plate", "center")
dim(barcode)
barcode[1:2,]

meta2 = cbind(meta2, barcode)
dim(meta2)
meta2[1:2,]

fnms = list.files("../gene_expression", pattern="htseq.counts",
recursive=TRUE)

fnmsp = fnms[grep("parcel", fnms)]
fnms  = setdiff(fnms, fnmsp)
length(fnms)
fnms[1:5]

fids = matrix(unlist(strsplit(fnms, split="/")), byrow=TRUE, ncol=2)[,2]
fids = gsub(".htseq.counts.gz", "", fids, fixed=TRUE)
fids = gsub(".htseq.counts", "", fids, fixed=TRUE)
fids[1:5]

meta2$analysis.submitter_id = gsub("_count", "", meta2$analysis.submitter_id)
table(fids %in% meta2$analysis.submitter_id)
meta2 = meta2[match(fids, meta2$analysis.submitter_id),]

dim(meta2)
meta2[1:2,]

# ----------------------------------------------------------------------------
# remove plates with 2 or less samples
# ----------------------------------------------------------------------------

tp = table(meta2$plate)
sort(tp)

plate2rm = names(tp)[which(tp <= 2)]
plate2rm

meta2 = meta2[-which(meta2$plate %in% plate2rm),]
dim(meta2)

# ----------------------------------------------------------------------------
# keep those from primary tumor samples
# ----------------------------------------------------------------------------

length(unique(meta2$participant[which(meta2$sample == "01A")]))
length(unique(meta2$participant[which(meta2$sample %in% c("01A", "01B"))]))
length(unique(meta2$participant[which(meta2$sample %in% c("01A", "01B", "01C"))]))

w2kp  = which(meta2$sample %in% c("01A", "01B", "01C"))
meta2 = meta2[w2kp,]
dim(meta2)
meta2[1:2,]

# ----------------------------------------------------------------------------
# keep one file per sample, the one with largest size of BAM file
# ----------------------------------------------------------------------------

tp = table(meta2$participant)
table(tp)
sam2check  = names(tp)[tp > 1]
ww2check   = which(meta2$participant %in% sam2check)
meta2check = meta2[ww2check,]
dim(meta2check)
meta2check[1:2,]

fun1 <- function(v){paste(sort(v), collapse="-") }
t1 = tapply(meta2check$sample, meta2check$participant, fun1)
table(t1)

sams = unique(meta2check$participant)
w2kp = NULL

for(sam1 in sams){
  ww1 = which(meta2check$participant == sam1)
  ww2 = ww1[which.max(meta2check$analysis.input_files.file_size[ww1])]
  w2kp = c(w2kp, ww2)
}

meta3 = rbind(meta2check[w2kp,], meta2[-ww2check,])
dim(meta3)
tp = table(meta3$participant)
table(tp)

dim(meta3)
meta3[1:2,]

length(unique(meta3$participant))

# ----------------------------------------------------------------------------
# collect gene expression data
# ----------------------------------------------------------------------------

setwd("~/research/TCGA/COAD/GDC/original_data/gene_expression")

eL = list()
w2rm = NULL

fnms = list.files("../gene_expression", pattern="htseq.counts",
recursive=TRUE)

fnmsp = fnms[grep("parcel", fnms)]
fnms  = setdiff(fnms, fnmsp)
length(fnms)
fnms[1:5]

for(i in 1:nrow(meta3)){
  
  gi = grep(meta3$analysis.submitter_id[i], fnms)
  
  if(length(gi) != 1){
    stop("file missing!\n")
  }
  
  fi = fnms[i]
  
  if(grepl(".gz$", fi)){
    system(sprintf("gunzip %s", fi))
    fi = gsub(".gz$", "", fi)
  }
  
  di = read.table(fi, as.is=TRUE, sep="\t")
  eL[[meta3$participant[i]]] = di
}

ng = sapply(eL, nrow)
table(ng)

eDat = matrix(NA, nrow=nrow(eL[[1]]), ncol=length(eL))

for(i in 1:length(eL)){
  if(i == 1){
    genes = eL[[i]][,1]
  }else{
    if(any(genes != eL[[i]][,1])){
      stop("gene name do not match\n")
    }
  }
  eDat[,i] = eL[[i]][,2]
}

rownames(eDat) = genes
colnames(eDat) = names(eL)

dim(eDat)
eDat[1:2,1:5]

table(colnames(eDat) == meta3$participant)

write.table(eDat, file = "../../data/expression_data_v1.txt", sep = "\t",
quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(meta3, file = "../../data/expression_sample2use_v1.txt", sep = "\t",
quote = FALSE, row.names = FALSE, col.names = TRUE)

q(save="no")
