args=(commandArgs(TRUE))

if(length(args)==0){
  stop("No arguments supplied.")
}else{
  for(i in 1:length(args)){
    eval(parse(text=args[[i]]))
  }
}
k=as.numeric(args[[1]])
k

#load libraries------------------------------

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

#Loading the Input Data
input.data0  <- fread(file.path(vicky.dir,"expression_v2_log_TReC.txt"),header=TRUE, na.strings = "", stringsAsFactors=FALSE)
input.data   <- data.frame(input.data0)
gene.data    <- fread(file.path(vicky.dir,"expression_v2_info.txt"),header=TRUE, na.strings = "", stringsAsFactors=FALSE)
sample.data0 <- fread(file.path(vicky.dir,"expression_v2_sample.txt"),header=TRUE,sep='\n', sep2='\t',na.strings = "", stringsAsFactors=FALSE)
sample.data  <- data.frame(do.call(rbind, strsplit(as.character(sample.data0[[1]]), '\t', fixed=T)))
colnames(sample.data) <- strsplit(names(sample.data0), '\t')[[1]]

identical(input.data$id, gene.data$geneId)
x=cbind(colnames(input.data[,-1]), as.character(sample.data$participant))
num=which(x[,1]!=x[,2])
which(x[num,1]!=paste0('X',x[num,2]))
sample.data$newSampleID=colnames(input.data[,-1])
# the order of genes in gene.data is the same as the order in input.data
# the order of samples in sample.data is the same as the order in input.data

ref.pathway = fread(file.path(vicky.dir,'data', "PathwayCommons10.All.hgnc.sif"), header=FALSE, stringsAsFactors = FALSE)

# MSigDB from Broad, e.g. C6 curated oncogenic pathways 
C6 = scan(file.path(vicky.dir, 'data', "c6.all.v5.2.symbols.gmt.txt"), what=character(0), sep='\n')
C6 = strsplit(C6, split='\t')
C6.info = matrix(unlist(lapply(C6, function(x){x[1:2]})), ncol=2, byrow=TRUE)
C6.genes = lapply(C6, function(x){x[-c(1:2)]})
names(C6.genes)= C6.info[,1]
C6.gene.2=unique(unlist(C6.genes))

gene.id=as.character(gene.data$hgnc_symbol)
C6.select.gene=intersect(C6.gene.2,gene.id)

length(unlist(C6.genes))
length(C6.gene.2)
length(C6.select.gene)

# so we will use C6.genes to be the pathway
gene.name.C6 = as.character(unlist(C6.genes))
gene.num.C6  = cumsum(as.numeric(sapply(C6.genes, length)))
# check how many genes in each gene set
num.overlap=sapply(C6.genes, function(x)length(intersect(x, gene.id)))
pert.overlap=sapply(C6.genes, function(x)length(intersect(x, gene.id))/length(x))
hist(num.overlap)
hist(pert.overlap)

cohort=which(pert.overlap > 0.98)
select.gene.list=NULL
for (it in 1:189){
cohort=it
cohort.genes=as.character(unlist(C6.genes[cohort]))
select.gene=intersect(cohort.genes, gene.id)
pathway.gene.input=match(select.gene, gene.id)
select.gene.list[[it]]=select.gene
}
save(select.gene.list, file=file.path(vicky.dir, "select.gene.name.Rdata"))


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

sapply(gene.short,length)
sapply(select.gene.list,length)
save(gene.short, file=file.path(vicky.dir, "select.gene.name.short.Rdata"))


cohort=k
cohort.genes=gene.short[[k]]
select.gene=intersect(cohort.genes, gene.id)
pathway.gene.input=match(select.gene, gene.id)


# standardize gene expression data
Data = t(input.data[,-1])[,pathway.gene.input]
nDat = apply(Data, 2,scale)
cDat=  t(t(Data) - colMeans(Data))
n=dim(Data)[1]
p=dim(Data)[2]


# run space r1 and r2 and spacelog r3 and r4 by using extBIC and BIC, respectively
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
save(r1,r2,r3,r4, file=file.path(vicky.dir, paste0("space_joint_log_imethod_", imethod,"_iData_",k,"_n_",n, "_p_",p, "_TCGA_short.RData")))
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
save(coef_log_BIC, coef_log_extBIC, coef_las_BIC, coef_las_extBIC, file=file.path(vicky.dir,paste0("PEN_iData_",k,"_n_",n, "_p_",p, "_TCGA_short.RData")))




