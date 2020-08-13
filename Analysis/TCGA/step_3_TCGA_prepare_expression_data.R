
setwd("~/research/TCGA/COAD/GDC/data/")

# ----------------------------------------------------------------------------
# Read in sample information
# ----------------------------------------------------------------------------

info = read.table("../original_data/gene_expression_meta/metadata_clinic.txt",
  sep="\t", header=TRUE)
dim(info)
info[1:2,]

sam = read.table("expression_sample2use_v1.txt", sep="\t", header=TRUE)
dim(sam)
sam[1:2,]

info$bcr_patient_barcode = substr(info$demographic.submitter_id, 1, 12)
barcode = substr(sam$associated_entities.entity_submitter_id, 1, 12)
table(barcode %in% info$bcr_patient_barcode)

barcode2kp = intersect(barcode, info$bcr_patient_barcode)

length(unique(barcode))
info = info[match(barcode2kp, info$bcr_patient_barcode),]
dim(info)
info[1:2,]

sam = sam[match(barcode2kp, barcode),]

sam = cbind(sam, info)
dim(sam)
sam[1:2,]

table(sam$tss)
table(sam$demographic.race)
table(sam$demographic.gender)
table(sam$demographic.ethnicity)
table(sam$diagnoses.tumor_grade)
table(sam$diagnoses.tumor_stage)

# ----------------------------------------------------------------------------
# Read in expression data
# ----------------------------------------------------------------------------

dat = read.table("expression_data_v1.txt", sep="\t", header=TRUE)
dim(dat)
dat[1:2,1:5]

samIDs = gsub("^X", "", colnames(dat))
table(sam$participant %in% samIDs)

dat = dat[,match(sam$participant, samIDs)]
table(sam$participant == gsub("^X", "", colnames(dat)))

# ----------------------------------------------------------------------------
# Read in gene location information
# ----------------------------------------------------------------------------

ff2 = "~/research/data/human/gencode/gencode.v22.genes.txt"

infoE = read.table(ff2, sep = "\t", header = TRUE, as.is = TRUE, quote="")
dim(infoE)
infoE[1:2, ]
length(unique(infoE$geneId))
length(unique(infoE$ensembl_gene_id))

table(rownames(dat) %in% infoE$geneId)
rownames(dat)[which(! rownames(dat) %in% infoE$geneId)]

features = intersect(rownames(dat), infoE$geneId)

datE  = dat[match(features, rownames(dat)), ]
infoE = infoE[match(features, infoE$geneId), ]
dim(datE)
dim(infoE)

infoE[1:5, ]

table(rownames(datE) == infoE$geneId)
table(infoE$chr, useNA="ifany")
table(infoE$strand, useNA="ifany")
table(is.na(infoE$ensembl_gene_id))

# ----------------------------------------------------------------------------
# Find a cutoff to filter out low expressed genes
# ----------------------------------------------------------------------------

setwd("~/research/TCGA/COAD/GDC/")

datEA = data.matrix(datE)

rMin = apply(datEA, 1, min)
rMed = apply(datEA, 1, median)
r75  = apply(datEA, 1, quantile, probs = 0.75)
r90  = apply(datEA, 1, quantile, probs = 0.90)

summary(rMin)
summary(rMed)
summary(r75)
summary(r90)

cor(rMin, rMed)
cor(r75,  rMed)
cor(r90,  rMed)

cbind(infoE, rMin, rMed, r75)[which(rMin > 6000),]

png("figures/expression_cts_summary.png", width = 6, height = 6,
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
table(r75 >= 10)
table(r75 >= 20)

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

cor(tot, s75)

png("figures/expression_total_vs_75_percentile.png", width = 4, height = 4,
units = "in", res = 400)
par(mar = c(5, 4, 1, 1), bty = "n")
plot(tot / 1e6, s75 / 1000, xlab = "total reads (million)",
ylab = "75 percentile (thousand)", cex = 0.5)
dev.off()

nDat = t(log10(t((datEA + 1)) / s75))
dim(nDat)

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

r1 = rgb(0.8,0.2,0.2,0.6)
b1 = rgb(0.2,0.2,0.8,0.6)

png("figures/expression_PCs_log_TReC.png", width = 6, height = 6,
units = "in", res = 400)
par(mar = c(5, 4, 1, 1), mfrow = c(2, 2), bty = "n")
barplot(prdatR1$values[1:20], main = "", xlab = "Index", ylab = "Eigen-value")

plot(PC1, PC2, cex = 0.8)
plot(PC1, PC3, cex = 0.8)
plot(PC2, PC3, cex = 0.8)

dev.off()


# ----------------------------------------------------------------------------
# Write out data and information
# ----------------------------------------------------------------------------

dim(datEA)
datEA[1:2, 1:5]

write.table(datEA, file = "data/expression_v2_counts.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = TRUE, col.names = TRUE)

pDat = data.frame(id = rownames(nDat), nDat)

dim(pDat)
pDat[1:2, 1:5]

write.table(pDat, file = "data/expression_v2_log_TReC.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(infoE)
infoE[1:5, ]

write.table(infoE, file = "data/expression_v2_info.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

dim(sam)
sam[1:2,]

write.table(sam, file = "data/expression_v2_sample.txt", append = FALSE,
quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

q(save="no")
