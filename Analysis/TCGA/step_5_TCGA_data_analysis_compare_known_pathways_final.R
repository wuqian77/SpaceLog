vicky.dir = file.path('~/wu_v/')
setwd(vicky.dir)

library(igraph)
library(ggplot2)
library(reshape2)
library(data.table)

methods = c("log_extBIC_or","log_BIC_or","lasso_extBIC_or","lasso_BIC_or", 
            "log_extBIC_and","log_BIC_and","lasso_extBIC_and","lasso_BIC_and", 
                "NS_L1_extBIC","space_extBIC","space_res_extBIC", "space_df_extBIC", "NS_L1_BIC","space_BIC","space_res_BIC", "space_df_BIC", 
                  "NS_L1_extBIC","space_log_extBIC","space_log_res_extBIC", "space_log_df_extBIC", "NS_L1_BIC","space_log_BIC","space_log_res_BIC", "space_log_df_BIC")

ref.pathway = fread(file.path(vicky.dir,"data", "PathwayCommons10.All.hgnc.sif"), header=FALSE, stringsAsFactors = FALSE)

result.list.space = list.files(file.path(vicky.dir), pattern="_n_450_")
result.list.space
fnms=result.list.space

n=450
nD = TP = TN = FP = FN = nD0 = matrix(nrow=189, ncol=length(methods))
lbs = NULL
truegraph.list=NULL
estimate.list=NULL
ref.list=NULL
name.vec=NULL
p.vec=NULL
# ix =37 
for (ix in 1:189)
{
  (load(file=file.path(vicky.dir, "select.gene.name.short.Rdata")))
  gene.name=gene.short[[ix]]
  name.vec[[ix]]=gene.name
  
  
  ref.list[[ix]] = ref.pathway[which(ref.pathway$V1 %in% gene.name & ref.pathway$V3 %in% gene.name),]
  ref=ref.list[[ix]]
  Ak=matrix(0, nrow=length(gene.name),ncol=length(gene.name))
  colnames(Ak)=gene.name
  rownames(Ak)=gene.name
  for (i in 1:dim(ref)[1]){
    x=as.character(ref[i,1])
    y=as.character(ref[i,3])
    ind.x=which(colnames(Ak) ==x)
    ind.y=which(rownames(Ak) ==y)
    Ak[ind.x, ind.y]=1
    Ak[ind.y, ind.x]=1
  }
  
  print(table(rowSums(Ak!=0))) # number of parents per node
  print(table(colSums(Ak!=0))) # number of parents/children per node
  
  trueGraph  = Ak
  table(trueGraph)
  
  table(rowSums( trueGraph !=0))
  table(colSums( trueGraph !=0))
  
  truegraph.list[[ix]] = trueGraph
  p=dim(trueGraph)[1]
  p.vec[ix]=p
  
  cidx = 0
 
 try(load(file=file.path(vicky.dir,paste0("PEN_iData_",ix,"_n_",n, "_p_",p, "_TCGA_short.RData"))),silent=TRUE)

  for (indc in c('coef_log_extBIC', 'coef_log_BIC', 'coef_las_extBIC', 'coef_las_BIC')){
    cidx=cidx+1
    coefk = eval(parse(text=indc))
     if (dim(coefk)[1]==dim(trueGraph)[1]){
         estimateGraph = matrix(0,p,p)
         estimateGraph[coefk!=0 | t(coefk)!=0] =1
         
    diff = c(estimateGraph) - c(trueGraph)
    estimate.list[[cidx]] = estimateGraph
    
    TP[ix, cidx] = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
    TN[ix, cidx] = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
    FP[ix, cidx] = sum(diff == 1)/2
    FN[ix, cidx] = sum(diff == -1)/2
    nD[ix, cidx] = sum(c(estimateGraph))/2
    nD0[ix, cidx] = sum(c(trueGraph))/2
     }
  }
    
 for (indc in c('coef_log_extBIC', 'coef_log_BIC', 'coef_las_extBIC', 'coef_las_BIC')){
        cidx=cidx+1
        coefk = eval(parse(text=indc))
        if (dim(coefk)[1]==dim(trueGraph)[1]){
        estimateGraph = matrix(0,p,p)
        estimateGraph[coefk!=0 & t(coefk)!=0] =1
        diff = c(estimateGraph) - c(trueGraph)
        estimate.list[[cidx]] = estimateGraph
        
        TP[ix, cidx] = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
        TN[ix, cidx] = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
        FP[ix, cidx] = sum(diff == 1)/2
        FN[ix, cidx] = sum(diff == -1)/2
        nD[ix, cidx] = sum(c(estimateGraph))/2
        nD0[ix, cidx] = sum(c(trueGraph))/2
        }
    }
  
  
  for(i in 1:4){
       cidx=cidx+1
    name=list.files(file.path(vicky.dir), pattern=paste0("space_joint_log_imethod_",i,"_iData_",ix,"_n_450_p_"))
    ind=grep('_short.RData',name)
    try(load(file=file.path(vicky.dir, name[ind])),silent=TRUE)
    kk=0
    for (indc in c('r2','r1','r4','r3')){
      kk=kk+1
      r = eval(parse(text=indc))
      sigma = as.numeric(r$result$sig.fit)
      sigma.mat = matrix(NA, nrow=length(sigma), ncol=length(sigma))

      for (j in 1:length(sigma)){
      sigma.mat[,j] = sqrt(sigma/sigma[j])
      }
    
      coef.r=r$result$ParCor*sigma.mat
      estimateGraph = (r$result$ParCor != 0)
      diag(estimateGraph) = 0
      estimate.list[[(kk-1)*4+cidx]] = estimateGraph
      
      if (dim(estimateGraph)[1]==dim(trueGraph)[1]){

        diff = c(estimateGraph) - c(trueGraph)
        
        TP[ix, (kk-1)*4+cidx] = sum(estimateGraph ==1 & trueGraph==1 & diff == 0)/2
        TN[ix, (kk-1)*4+cidx] = sum(estimateGraph ==0 & trueGraph==0 & diff == 0)/2
        FP[ix, (kk-1)*4+cidx] = sum(diff == 1)/2
        FN[ix, (kk-1)*4+cidx] = sum(diff == -1)/2
        nD[ix, (kk-1)*4+cidx] = sum(c(estimateGraph))/2
        nD0[ix,(kk-1)*4+cidx] = sum(c(trueGraph))/2
      }
    }
  }
}

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
FDR=FP/nD
Error=FP+FN
precision = TP/(TP+FP)
recall = TPR
F1 = 2*precision*recall/(precision+recall)
save(TP, TN, FP, FN, nD, nD0, FDR, F1, TPR, methods,p.vec,truegraph.list, estimate.list, file=file.path(vicky.dir,paste0("TCGA_log_extBIC.RData")))

load(file=file.path(vicky.dir,paste0("TCGA_log_extBIC.RData")))

TPR=TP/(TP+FN)
FPR=FP/(FP+TN)
FDR=FP/nD
Error=FP+FN
precision = TP/(TP+FP)
recall = TPR
F1 = 2*precision*recall/(precision+recall)

nD.mat = sapply(truegraph.list,function(x)apply(x,2,function(x)sum(x!=0)))
a=sapply(nD.mat,max)/sapply(nD.mat,median)
summary(a)
b=sapply(nD.mat,sd)/sapply(nD.mat,median)
old_b=sapply(nD.mat, sd)
p.vec=sapply(truegraph.list,dim)[1,]


method.vec=methods
lab.BIC = rep("BIC", ncol(FP))
lab.BIC[grep("extBIC", method.vec)] = "extBIC"

lab.intersect = rep('OR',ncol(FP))
lab.intersect[grep('intersect',method.vec)] ='intersect'

lab.method = rep("log", ncol(FP))
lab.method[grep("lasso", method.vec)] = "NS_LASSO"
lab.method[grep("log", method.vec)] = "NS_log"
lab.method[grep("NS_L1", method.vec)] = "NS_L1"
lab.method[grep("space_", method.vec)] = "space"
lab.method[grep("space_res", method.vec)] = "space"
lab.method[grep("space_df", method.vec)] = "space"
lab.method[grep("space_log", method.vec)] = "sp_log"

methods = c("NS_log","NS_log","NS_lasso","NS_lasso", 
            "NS_log","NS_log","NS_lasso","NS_lasso", 
                "NS_L1_extBIC","space_no","space_res", "space_df", "NS_L1_BIC","space_no","space_res", "space_df", 
                  "NS_L1_extBIC","space_log_no","space_log_res", "space_log_df", "NS_L1_BIC","space_log_no","space_log_res", "space_log_df")


colnames(TP) = colnames(TN) = colnames(nD) = colnames(nD0) =colnames(FP) = colnames(FN) = colnames(F1) = colnames(F1) = colnames(FDR)= colnames(TPR) = methods


cols = rep("blue", ncol(FP))
cols[which(lab.method=="NS_LASSO")] = "orange"
cols[which(lab.method=="NS_log")] = "darkgrey"
cols[which(lab.method=="space")] = "palegreen"



cols.bd = rep("blue", ncol(FP))
cols.bd[which(lab.method=="NS_LASSO")] = "orange"
cols.bd[which(lab.method=="NS_log")] = "grey"
cols.bd[which(lab.method=="space")] = "darkgreen"



FDR.median = apply(FDR, 2, function(x)median(x, na.rm=TRUE))
recall.median = apply(recall,2, function(x)median(x, na.rm=TRUE))
F1.median = apply(F1,2, function(x)median(x, na.rm=TRUE))

FDR.Q1 = apply(FDR, 2, function(x)quantile(x,0.25, na.rm=TRUE))
recall.Q1 = apply(recall,2, function(x)quantile(x,0.25, na.rm=TRUE))
F1.Q1 = apply(F1,2, function(x)quantile(x,0.25, na.rm=TRUE))

FDR.Q3 = apply(FDR, 2, function(x)quantile(x,0.75, na.rm=TRUE))
recall.Q3 = apply(recall,2, function(x)quantile(x,0.75, na.rm=TRUE))
F1.Q3 = apply(F1,2, function(x)quantile(x,0.75, na.rm=TRUE))


plot.mat = data.frame(names=method.vec, BIC=lab.BIC,
                      method=lab.method,
                      FDR.Q1, recall.Q1, F1.Q1,
                      FDR.median, recall.median, F1.median,
                      FDR.Q3, recall.Q3, F1.Q3)

write.csv(plot.mat, file=file.path(vicky.dir,"TCGA_recall_precision.csv"))

pdf(file.path(vicky.dir,paste0("TCGA_ALL_extBIC_BA_v2.pdf")), width=12, height=10)
par(mfrow=c(2,3), mar=c(10,5,3,1), bty="n", cex=0.85, las=1)
w2kp = which(lab.BIC %in% c('extBIC')& lab.method !='NS_L1' & lab.intersect=='OR')
keep = which(b > quantile(b,0.25) & p.vec<100)

boxplot(FP[keep,w2kp], main="False Positives", col=cols[w2kp],border=cols.bd[w2kp],las=2,srt=45,adj=1, xpd=TRUE)
legend("topright", c("space_log","space","NS_lasso","NS_log"), fill=c("blue","palegreen","orange","darkgrey"), bty="n")

boxplot(FN[keep,w2kp], main="False Negatives", col=cols[w2kp],border=cols.bd[w2kp],las=2,srt=45, xpd=TRUE)
boxplot(FP[keep,w2kp] + FN[keep,w2kp], main="FP + FN", col=cols[w2kp],
        border=cols.bd[w2kp], las=2)

boxplot(F1[keep,w2kp], main="F1",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)
boxplot(nD[keep,w2kp], main="# of discovery",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)
boxplot(TPR[keep,w2kp], main="Power / TPR",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)

dev.off()


pdf(file.path(vicky.dir,paste0("TCGA_ALL_extBIC_ALL.pdf")), width=12, height=10)
par(mfrow=c(2,3), mar=c(10,5,3,1), bty="n", cex=0.85, las=1)
w2kp = which(lab.BIC %in% c('extBIC')& lab.method !='NS_L1' & lab.intersect=='OR')
keep = which(b > 0 )

boxplot(FP[keep,w2kp], main="False Positives", col=cols[w2kp],border=cols.bd[w2kp],las=2,srt=45,adj=1, xpd=TRUE)
legend("topright", c("space_log","space","NS_lasso","NS_log"), fill=c("blue","palegreen","orange","darkgrey"), bty="n")

boxplot(FN[keep,w2kp], main="False Negatives", col=cols[w2kp],border=cols.bd[w2kp],las=2,srt=45, xpd=TRUE)
boxplot(FP[keep,w2kp] + FN[keep,w2kp], main="FP + FN", col=cols[w2kp],
        border=cols.bd[w2kp], las=2)

boxplot(F1[keep,w2kp], main="F1",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)
boxplot(nD[keep,w2kp], main="# of discovery",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)
boxplot(TPR[keep,w2kp], main="Power / TPR",col=cols[w2kp],
        border=cols.bd[w2kp],  las=2)

dev.off()
