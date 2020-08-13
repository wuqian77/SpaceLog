
library(knitr)
library(arsenal)
library(tidyverse)
library(readxl)
library(forcats)
library(plyr)
library(dplyr)
library(data.table)
library(graphics)
library(space)
library(ggplot2)
library(PEN)

#--set environment----------------------------------
setwd("~/dbGaP_10277_GTEx/")
out.dir=file.path("~/results")

data.file=file.path("57267/PhenoGenotypeFiles/RootStudyConsentSet_phs000424.GTEx.v7.p2.c1.GRU")
subject.path    = "ExpressionFiles/phe000020.v1.GTEx_RNAseq.sample-info.MULTI/phe000020.v1_release_manifest.txt"
express.path    = "ExpressionFiles/phe000020.v1.GTEx_RNAseq.expression-data-matrixfmt.c1"
genotype.dir=file.path("GenotypeFiles/phg000830.v1.GTEx.genotype-qc.MULTI/GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_support_files")

# subject clinical data
subject.clin <- fread("phs000424.v7.pht002742.v7.p2.c1.GTEx_Subject_Phenotypes.GRU.txt",header=FALSE, na.strings = "", stringsAsFactors=FALSE) %>% as.data.frame()
colnames(subject.clin)=subject.clin[2,]
clinical.data <- subject.clin[3:dim(subject.clin)[1],]
sample.mapping <- fread(file.path(data.geno.dir,'PhenotypeFiles',"phs000424.v7.pht002741.v7.p2.GTEx_Sample.MULTI.txt"),header=TRUE, na.strings = "", stringsAsFactors=FALSE) %>% as.data.frame()
colnames(sample.mapping)=sample.mapping[1,]
mapping <- sample.mapping[2:dim(sample.mapping)[1],]
subject.after.mapping= merge(clinical.data, mapping, by='SUBJID',all.x=TRUE)
write.table(subject.after.mapping, file.path(out.dir, "subject.after.mapping.txt"),sep = '\t', col.names=TRUE)

subject.data = fread(file.path(data.geno.dir, subject.path), sep = "\t",skip=16, header = TRUE, stringsAsFactors=TRUE) %>% as.data.frame()
head(subject.data)
subject.id = subject.data$Subject_ID 
sample.id  = subject.data$Sample_ID 

# blood type
blood.type = fread("phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt", sep = "\t",header = TRUE, stringsAsFactors=TRUE) %>% as.data.frame()
emInfo = fread('zcat GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_gene_reads.gct.gz',sep = "\t",skip=1, header = TRUE, stringsAsFactors=TRUE) %>% as.data.frame()


dim(emInfo)
emInfo[1:5,1:5]
em.sample=colnames(emInfo)
em.subject=do.call(rbind,strsplit(as.character(emInfo[,1]),'\\.'))[,1]
em.gene=emInfo$Description
length(unique(em.sample))
length(unique(em.subject))

loc=fread('GTEx_Data_20160115_v7_RNAseq_RNASeQCv1.1.8_metrics.tsv',sep='\t',header=TRUE,stringsAsFactors=TRUE) %>% as.data.frame()
loc.blood.type = merge(loc, blood.type, by.x='Sample',by.y='SAMPID')
sample.id.v2  = loc.blood.type$Sample # 12767
ind=which(loc.blood.type$Note %in% c('Whole Blood'))
table(loc.blood.type$Note, loc.blood.type$SMOMTRLTP)
sample.id.blood=sample.id.v2[ind]

all = merge(loc.blood.type, subject.after.mapping, by.x='Sample', by.y='SAMPID')
a1=table(all$SMNABTCH)
a2=table(all$SMNABTCHT)
a3=table(all$ANALYTE_TYPE)
a4=table(all$SMGEBTCH)
# SMNABTCH: Nucleic Acid Isolation Batch ID
# SMNABTCHT: Type of genotype or expression batch
# ANALYTE_TYPE: Analyte Type
# SMGEBTCH: Genotype or Expression Batch ID

# check how many unique subjects having gene expression data
subject.after.mapping[match(loc.blood.type$Sample, subject.after.mapping$SAMPID),'SUBJID'] %>% unique %>% length
clinical.data.s = all[match(sample.id.blood, all$Sample),]

# read in the 16 subjects removed due to chromosome abnormal
rm.data <- fread("GTEx_Analysis_2016-01-15_v7_WholeGenomeSeq_flagged_samples.txt", header=FALSE, stringsAsFactors = FALSE) %>% as.data.frame()
rm.id = subject.after.mapping$SUBJID[match(rm.data$V1,subject.after.mapping$SAMPID)]
clinical.data.GE = clinical.data.s[-which(clinical.data.s$SUBJID %in% rm.id),]

# ----------------------------------------------------------------------------
# Read in gene expression data
# ----------------------------------------------------------------------------
emInfo.s = emInfo[,which(colnames(emInfo)%in% clinical.data.GE$Sample)]
dim(emInfo.s)
rownames(emInfo.s)= emInfo[,1]
infoE = emInfo[,1:2]
# ----------------------------------------------------------------------------
# Read in gene location information
# ----------------------------------------------------------------------------
# ----------------------------------------------------------------------------
# Find a cutoff to filter out low expressed genes
# ----------------------------------------------------------------------------

datEA = data.matrix(emInfo.s)

rMin = apply(datEA, 1, min)
rMax = apply(datEA, 1, max)
rMed = apply(datEA, 1, median)
r75  = apply(datEA, 1, quantile, probs = 0.75)
r90  = apply(datEA, 1, quantile, probs = 0.90)

summary(rMin)
summary(rMax)
summary(rMed)
summary(r75)
summary(r90)

cor(rMin, rMed)
cor(r75,  rMed)
cor(r90,  rMed)

png(file.path(out.dir, "expression_cts_summary.png"), width = 6, height = 6,
    units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(5, 4, 1, 1), bty = "n")
hist(log10(1 + rMin), xlab = "log10(min + 1)", main = "")
hist(log10(1 + rMed), xlab = "log10(median + 1)", main = "")
hist(log10(1 + r75),  xlab = "log10(75 percentile + 1)", main = "")
hist(log10(1 + r90),  xlab = "log10(90 percentile + 1)", main = "")
dev.off()

summary(rMin[rMed >= 10])
summary(rMed[r75 >= 20])

table(rMed >= 10)
table(r75 >= 20)

r80  = apply(datEA, 1, quantile, probs = 0.80)

table(r80 >= 6)

w2kp = which(r75 >= 20)

dim(datEA)
datEA = datEA[w2kp, ]
dim(datEA)

dim(infoE)
infoE = infoE[w2kp, ]
dim(infoE)

if(! all(rownames(datEA) == infoE$gene) ){
  stop("gene name mismatch\n")
}

# ----------------------------------------------------------------------------
# Normalize gene expression by read-depth
# ----------------------------------------------------------------------------

tot = colSums(datEA)
s75 = apply(datEA, 2, quantile, prob = 0.75)

png(file.path(out.dir, "hist_75_percentile.png"), width = 4, height = 4,
    units = "in", res = 400)
par(mar = c(5, 4, 1, 1), bty = "n")
hist(log10(tot),breaks = 50)
dev.off()

cor(tot, s75)

png(file.path(out.dir, "expression_total_vs_75_percentile.png"), width = 4, height = 4,
    units = "in", res = 400)
par(mar = c(5, 4, 1, 1), bty = "n")
plot(tot / 1e6, s75 / 1000, xlab = "total reads (million)",
     ylab = "75 percentile (thousand)", cex = 0.5)
dev.off()

nDat = t(log10(t((datEA + 1)) / (s75+1)))
dim(nDat)
png(file.path(out.dir, "histogram.png"), width = 4, height = 4,
    units = "in", res = 400)
par(mfrow = c(2, 2), mar = c(5, 4, 1, 1), bty = "n")
par(mar = c(5, 4, 1, 1), bty = "n")
hist(nDat[,1],cex = 0.5)
hist(nDat[,100],cex = 0.5)
hist(nDat[1,],cex = 0.5)
hist(nDat[100,],cex = 0.5)
dev.off()

# ----------------------------------------------------------------------------
# Run PCA using gene expression data, check possible outlier
# these PCs do include many batch effect information
# ----------------------------------------------------------------------------

datR14Pr = nDat - rowMeans(nDat, na.rm = TRUE)

datR14Pr[is.na(datR14Pr)] = 0
covdatR1 = t(datR14Pr) %*% datR14Pr / nrow(datR14Pr)
dim(covdatR1)
prdatR1 = eigen(covdatR1)

prdatR1$values[1:20]

PC1 = prdatR1$vectors[, 1]
PC2 = prdatR1$vectors[, 2]
PC3 = prdatR1$vectors[, 3]

sample.id.s = colnames(covdatR1)
do.call(rbind, strsplit(sample.id.s, '-'))[,1] %>% table
do.call(rbind, strsplit(sample.id.s, '-'))[,2] %>% table
do.call(rbind, strsplit(sample.id.s, '-'))[,3] %>% table
do.call(rbind, strsplit(sample.id.s, '-'))[,4] %>% table

pca.check = do.call(rbind, strsplit(sample.id.s, '-'))[,3] 

require(gridExtra)
require(ggplot2)

pc.dat=data.frame(subject.id=rownames(covdatR1),PC1, PC2, PC3, pca.check)
pc.dat$pca.check=as.factor(pc.dat$pca.check)
pc.dat.v2 = merge(pc.dat, clinical.data.GE[,c('Sample','COHORT','SEX','RACE','ETHNCTY',"AGE","WGHT","BMI", "SMNABTCH", "SMNABTCHT", "SMGEBTCH")], by.x='subject.id',by.y='Sample')
pc.dat.v2$SEX = as.factor(ifelse(pc.dat.v2$SEX==1,'Male','Female'))
pc.dat.v2$COHORT = as.factor(pc.dat.v2$COHORT)
pc.dat.v2$RACE = as.factor(ifelse(pc.dat.v2$RACE==1,'Asian',
                                  ifelse(pc.dat.v2$RACE==2,'Black',
                                         ifelse(pc.dat.v2$RACE==3,'White',
                                                ifelse(pc.dat.v2$RACE==4, "American Indian or Alaska","Unknown")))))
pc.dat.v2$ETHNCTY = as.factor(ifelse(pc.dat.v2$ETHNCTY==0,'Not Latino',
                                     ifelse(pc.dat.v2$ETHNCTY==1,'Latino',
                                            ifelse(pc.dat.v2$ETHNCTY==98,'Not Reported',
                                                   ifelse(pc.dat.v2$ETHNCTY==99, "Unknown","CEPH")))))


a1=table(pc.dat.v2$SMNABTCH)
a2=table(pc.dat.v2$SMNABTCHT)
a3=table(pc.dat.v2$SMGEBTCH)
table(a1)
table(a2)
table(a3)

ind.a3=names(which(a3>6))

pc.dat.v2$SMGEBTCH.fac = as.factor(ifelse(pc.dat.v2$SMGEBTCH %in% ind.a3,pc.dat.v2$SMGEBTCH,'<= 4'))


#------ by subject name ---------
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=pca.check)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=pca.check)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=pca.check)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'by.sample.name.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#------ by race ---------
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=RACE)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=RACE)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=RACE)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.RACE.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#------ by ethnicity ---------
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=ETHNCTY)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=ETHNCTY)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=ETHNCTY)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.ethnicity.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#------ by gender
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=SEX)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=SEX)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=SEX)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.SEX.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#------ by cohort
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=COHORT)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=COHORT)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=COHORT)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.COHORT.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()


# more variables might be batch
#PhenotypeFiles/phs000424.v7.pht002743.v7.p2.c1.GTEx_Sample_Attributes.GRU.txt
#column 20 SMNABTCH: Nucleic Acid Isolation Batch ID
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=SMNABTCH)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=SMNABTCH)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=SMNABTCH)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.SMNABTCH.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()
#column 21 SMNABTCHT: Type of genotype or expression batch
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=SMNABTCHT)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=SMNABTCHT)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=SMNABTCHT)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.SMNABTCHT.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#column 26 SMGEBTCH: Genotype or Expression Batch ID
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=SMGEBTCH)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=SMGEBTCH)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=SMGEBTCH)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.SMGEBTCH.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()

#column 26 SMGEBTCH.fac: Genotype or Expression Batch ID
p1 <- ggplot(pc.dat.v2, aes(PC1,PC2,color=SMGEBTCH.fac)) + geom_point() + theme(legend.position='none')
p2 <- ggplot(pc.dat.v2, aes(PC1,PC3,color=SMGEBTCH.fac)) + geom_point() + theme(legend.position='none')
p3 <- ggplot(pc.dat.v2, aes(PC2,PC3,color=SMGEBTCH.fac)) + geom_point() + theme(legend.position=c(0.8, 0.8),legend.key.size = unit(0.4,'cm'),legend.text=element_text(size=7))

png(filename=file.path(out.dir, 'PCA.by.SMGEBTCH.fac.png'), width=12, height=4, unit="in", res=300)
grid.arrange(p1, p2, p3, ncol=3)
dev.off()



a1=table(pc.dat.v2$SMNABTCH)
a2=table(pc.dat.v2$SMNABTCHT)
a3=table(pc.dat.v2$SMGEBTCH)

pvls = matrix(NA, nrow=5, ncol=5)
for(i in 1:5){
  PCi = prdatR1$vectors[, i]
  
  a1  = anova(lm(PCi ~ pc.dat.v2$SMNABTCH))
  pvls[i,1] = a1$`Pr(>F)`[1]
  
  a1  = anova(lm(PCi ~ pc.dat.v2$SMNABTCHT))
  pvls[i,2] = a1$`Pr(>F)`[1]
  
  a1  = anova(lm(PCi ~ pc.dat.v2$SMGEBTCH))
  pvls[i,3] = a1$`Pr(>F)`[1]
  
  a1  = anova(lm(PCi ~ pc.dat.v2$COHORT))
  pvls[i,4] = a1$`Pr(>F)`[1]
  
}

colnames(pvls) = c("COHORT", "SEX", "RACE", "ETHNCTY", 
                   "pca.check")
rownames(pvls) = paste0("PC", 1:5)

signif(pvls, 2)

write.csv(pvls, file.path(out.dir, "pvalue.summary.byPCs.csv"))

# ----------------------------------------------------------------------------
# Write out data and information
# ----------------------------------------------------------------------------

dim(datEA)
datEA[1:2, 1:5]

write.table(datEA, file.path(data.dir, "expression_counts.txt"), append = FALSE,
            quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

table(colnames(nDat) == emInfo$patient_id)
pDat = data.frame(id = rownames(nDat), nDat)

dim(pDat)
pDat[1:2, 1:5]

write.table(pDat, file.path(data.dir, "expression_log_TReC.txt"), append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(infoE)
infoE[1:5, ]

write.table(infoE, file.path(data.dir,  "expression_info.txt"), append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

pcDat = data.frame(pc.dat.v2[,1:4],prdatR1$vectors[,4:5],pc.dat.v2[,5:dim(pc.dat.v2)[2]])
names(pcDat) = c(colnames(pc.dat.v2)[1:4], paste0("PC", 4:5), colnames(pc.dat.v2)[5:dim(pc.dat.v2)[2]])
dim(pcDat)
pcDat[1:5, ]

write.table(pcDat, file.path(data.dir, "expression_PCs.txt"), 
            append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, 
            col.names = TRUE)

# remove the batch effect COHORT from cohort and 
which(colnames(nDat)!=pc.dat.v2$subject.id) 
nDat.new=NULL
for (j in 1:dim(nDat)[1]){
  x=pc.dat.v2$COHORT
  y=nDat[j,]
  m1=lm(y~x)
  y.new=summary(m1)$residuals
  nDat.new=rbind(nDat.new, y.new)
}

pDat = data.frame(id = rownames(nDat), nDat.new)

dim(pDat)
pDat[1:2, 1:5]

write.table(pDat, file.path(out.dir, "expression_log_TReC_rmCohort.txt"), append = FALSE,
            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)


