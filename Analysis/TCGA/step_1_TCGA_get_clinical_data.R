
library(jsonlite)
library(XML)

setwd("~/research/TCGA/COAD/GDC/original_data/clinical_meta")

# ----------------------------------------------------------------------------
# Read in manifest
# ----------------------------------------------------------------------------

meta = fromJSON("metadata.cart.2016-08-21T20_16_05.132455.json")
dim(meta)
meta[1:2,]

# ----------------------------------------------------------------------------
# Read in clinical data
# ----------------------------------------------------------------------------

setwd("~/research/TCGA/COAD/GDC/original_data/clinical")

fnms = list.files("./", pattern="xml$", recursive=TRUE)
length(fnms)
fnms[1:5]

fids = matrix(unlist(strsplit(fnms, split="/")), byrow=TRUE, ncol=2)[,2]
fids = gsub("nationwidechildrens.org_clinical.", "", fids, fixed=TRUE)
fids = gsub(".xml", "", fids, fixed=TRUE)
fids[1:5]
length(fids)
length(unique(fids))

dClinic = list()

nms = NULL

for(i in 1:length(fnms)){
  fi    = fnms[i]
  xi    = xmlParse(fi)
  xitop = xmlRoot(xi)
  
  di = xmlSApply(xmlChildren(xitop)$patient, function(x) xmlSApply(x, xmlValue))
  yi = sapply(di, function(x){paste(x, collapse="; ")} )
  
  nms  = union(nms, names(yi))
  dClinic[[i]] = yi
}

length(nms)

dfClinic = NULL

for(i in 1:length(dClinic)){
  di = dClinic[[i]]
  vi = rep(NA, length(nms))
  vi[match(names(di), nms)] = di
  dfClinic = rbind(dfClinic, vi)
}

dim(dfClinic)
colnames(dfClinic) = nms
dfClinic = as.data.frame(dfClinic)
summary(dfClinic)

dfClinic[1:2,1:41]
dfClinic[1:2,42:63]
dfClinic[1:2,64:66]

write.table(dfClinic, file = "../../data/clinic_data.txt", sep = "\t",
  quote = FALSE, row.names = FALSE, col.names = TRUE)

q(save="no")
