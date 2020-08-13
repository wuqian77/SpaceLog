args=(commandArgs(TRUE))

if(length(args)==0){
    stop("No arguments supplied.")
}else{
    for(i in 1:length(args)){
        eval(parse(text=args[[i]]))
    }
}

k=as.numeric(args[[1]])

#--set environment----------------------------------
vicky.dir = file.path('~/wu_v/')
setwd(vicky.dir)

library(plyr)
library(dplyr)
library(data.table)
library(graphics)
library(PEN)
library(pcalg)
library(PenPC)
library(spacelog)

# Reference Pathway information is downloaded from https://www.pathwaycommons.org/archives/PC2/v10/
ref.pathway = fread(file.path(vicky.dir, "data","PathwayCommons10.All.hgnc.sif"), header=FALSE, stringsAsFactors = FALSE)

# GTex data can be download https://www.gtexportal.org/home/datasets
# Study accession is phs000424.v7.p2.
# expression_log_TReC_rmCohort.txt
# expression_info.txt

# After pre-processing, 
input.data  <- fread(file.path(vicky.dir,"expression_log_TReC_rmCohort.txt"),header=TRUE, na.strings = "", stringsAsFactors=FALSE) %>% as.data.frame()
subject.id=colnames(input.data)[-1]  
gene.id = input.data$id

gene.data    <- fread(file.path(out.dir,"expression_info.txt"),header=TRUE, na.strings = "", stringsAsFactors=FALSE) %>% as.data.frame()
gene.name = gene.data[,2][match(gene.id, gene.data[,1])]
which(duplicated(gene.name))

# 3. MSigDB from Broad, e.g. C6 curated oncogenic pathways 
# https://data.broadinstitute.org/gsea-msigdb/msigdb/release/5.2/
# 5.2 is an older version, which has been archived. Try v6.2 instead if needed. 
C6 = scan(file.path(base.dir,'data', 'Download_GSEA', "c6.all.v5.2.symbols.gmt.txt"), what=character(0), sep='\n')
C6 = strsplit(C6, split='\t')
C6.info = matrix(unlist(lapply(C6, function(x){x[1:2]})), ncol=2, byrow=TRUE)
C6.genes = lapply(C6, function(x){x[-c(1:2)]})
names(C6.genes)= C6.info[,1]
C6.gene.2=unique(unlist(C6.genes))
C6.select.gene=intersect(C6.gene.2,gene.name)

length(unlist(C6.genes)) # 31319 
length(C6.gene.2) # 11250
length(C6.select.gene) # 8097

# so we will use C6.genes to be the pathway
gene.name.C6 = as.character(unlist(C6.genes))
gene.num.C6  = cumsum(as.numeric(sapply(C6.genes, length)))
# check how many genes in each gene set
num.overlap=sapply(C6.genes, function(x)length(intersect(x, gene.name)))
pert.overlap=sapply(C6.genes, function(x)length(intersect(x, gene.name))/length(x))
hist(num.overlap)
hist(pert.overlap)

#select gene pathway with more overlapping genes with TCGA
cohort=which(pert.overlap > 0.98)
select.gene.list=NULL
for (it in 1:189){
cohort=it
cohort.genes=as.character(unlist(C6.genes[cohort]))
select.gene=intersect(cohort.genes, gene.name)
pathway.gene.input=match(select.gene, gene.name)
select.gene.list[[it]]=select.gene
}
save(select.gene.list, file=file.path(vicky.dir, "select.gene.name.GTEX.Rdata"))
  

name.vec=NULL
ref.list=NULL
gene.short=NULL
# short gene sets
for (ix in 1:189){
  gene.name=select.gene.list[[ix]]
  name.vec[[ix]]=gene.name
  ref.list[[ix]] = ref.pathway[which(ref.pathway$V1 %in% gene.name & ref.pathway$V3 %in% gene.name),]
  ref=ref.list[[ix]]
  gene.short[[ix]] = unique(c(ref$V1,ref$V3))
}

save(gene.short, file=file.path(out.dir, "select.gene.name.short.Rdata"))

cohort=k
cohort.genes=gene.short[[k]]
select.gene=intersect(cohort.genes, gene.name)
pathway.gene.input=match(select.gene, gene.name)
select.gene.list[[k]]=select.gene

# standardize gene expression data
Data = t(input.data[,-1])[,pathway.gene.input]
nDat = apply(Data, 2,scale)
cDat=  t(t(Data) - colMeans(Data))
n=dim(Data)[1]
p=dim(Data)[2]


nlambda=100
ntau=10
dat=nDat
meanx = apply(dat, 2, mean)
normx = sqrt(rowSums((t(dat) - meanx)^2)/n)
nDat = scale(dat, meanx, normx)
corr = abs(crossprod(nDat, nDat))
diag(corr) = 0
lamMax = max(corr)
thresh = 2 * exp(seq(log(lamMax), log(1/n), len = nlambda))
lambda = thresh/10
tau = 10^(seq(-6, 0, length.out = ntau))


lamb1.vec=lambda
lamb3.vec=tau
seed=777

iter=3
nn=n

for (imethod in 1:4){
r1 = space.joint.extBIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec=0, seed)
r2 = space.joint.BIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec=0, seed)
r3 = space.joint.extBIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)
r4 = space.joint.BIC(nDat, imethod, lamb1.vec, lam2=0, lamb3.vec, seed)
save(r1,r2,r3,r4, file=file.path(vicky.dir, paste0("space_joint_log_imethod_", imethod,"_iData_",k,"_n_",n, "_p_",p, "_GTex_short.RData")))
}

# PEN log penalty Neighborhood Selection extBIC
coef_log_extBIC  = ne.PEN(dat=Data,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=FALSE)
coef_log_BIC  = ne.PEN(dat=Data,nlambda=100,ntau=10,V=1:p,order=TRUE,verbose=FALSE, Model.selection='BIC')

# PEN log penalty Neighborhood Selection extBIC
coef_las_BIC = ne.PEN.lasso(dat=Data, V=1:p, order=TRUE, verbose=FALSE,
Model.selection="BIC")
# ExtendedBIC
coef_las_extBIC = ne.PEN.lasso(dat=Data, V=1:p, order=TRUE, verbose=FALSE,
Model.selection="ExtendedBIC")
save(coef_log_BIC, coef_log_extBIC, coef_las_BIC, coef_las_extBIC, file=file.path(vicky.dir,paste0("PEN_iData_",k,"_n_",n, "_p_",p, "_GTex_short.RData")))




