#Chang Chen   BioE582 Project 2

#Load packages
library(ALL)
library(hgu95av2.db)
library(simpleaffy)
library(genefilter)
library(multtest)
library(GOstats)
library(scatterplot3d)

#Subsetting and non-specific filtering
data(ALL)
ALL<-get.array.subset(ALL,"mol.biol",c("BCR/ABL","NEG"))
filist<-filterfun(pOverA(0.25,100),cv(0.7,10))
filter<-genefilter(2^exprs(ALL),filist)
filtALL<-ALL[filter,]

#Standard t-tests 
BCR<-rep(0,length(pData(ALL)$mol.bio))
BCR[grep("BCR",as.character(pData(ALL)$mol.bio))]<-1
BCR<-factor(BCR)
f <- function(x, factor) {
    d<-data.frame(x,factor)
    t.test(x~factor,factor)$p.value
}
pvs <- esApply(filtALL,1,f,BCR)
length(pvs[pvs<0.05])

#marginal adjustment
adjp <- mt.rawp2adjp(pvs,c("Bonferroni","Holm","Hochberg"))
bonf <- adjp$adjp[,2]
length(bonf[bonf<0.05])
holm <- adjp$adjp[,3]
length(holm[holm<0.05])
hoch <- adjp$adjp[,4]
length(hoch[hoch<0.05])

adjp2 <- mt.rawp2adjp(pvs,c("BH","BY"))
bh<- adjp2$adjp[,2]
length(bh[bh<0.05])
by <- adjp2$adjp[,3]
length(by[by<0.05])

#Biological annotation
index<-adjp$index[1:length(bh[bh<0.05])]
probeset.ID<-names(pvs)[index]
Gene.Symbol <- mget(probeset.ID,envir=hgu95av2SYMBOL)
Gene.Name <-  mget(probeset.ID,envir=hgu95av2GENENAME)
raw.p.value<- pvs[probeset.sig]
adjusted.p.value <- bh[bh<0.05]
tb <- cbind(probeset.ID,Gene.Symbol,Gene.Name,raw.p.value,adjusted.p.value)
write.table(tb,file = "de.gene.report.txt",sep=",",row.names = FALSE)

#GO term enrichment analysis
sigids <- mget(probeset.ID,envir=hgu95av2ENTREZID)
sigids <- unique(sigids)
sigids <- as.vector(unlist(sigids))
sigids <- sigids[!is.na(sigids)]

params <- new("GOHyperGParams",geneIds=sigids,annotation="hgu95av2",ontology="MF",
pvalueCutoff=0.05,conditional=FALSE,testDirection="over")
resultMF <- hyperGTest(params)
params <- new("GOHyperGParams",geneIds=sigids,annotation="hgu95av2",ontology="BP",
pvalueCutoff=0.05,conditional=FALSE,testDirection="over")
resultBP <- hyperGTest(params)
params <- new("GOHyperGParams",geneIds=sigids,annotation="hgu95av2",ontology="CC",
pvalueCutoff=0.05,conditional=FALSE,testDirection="over")
resultCC <- hyperGTest(params)

htmlReport(resultMF,file="GO-hypergeo.html",append=TRUE)
htmlReport(resultBP,file="GO-hypergeo.html",append=TRUE)
htmlReport(resultCC,file="GO-hypergeo.html",append=TRUE)

#Priciple component analysis
filtX <- t(exprs(filtALL))
pp <- prcomp(filtX)
summary(pp)
pdf("pcaplot.pdf")
plot(pp$x[, 1:2],col = ifelse(BCR == "0", "green","red"), pch = 16, xlab = "PC.1", ylab = "PC.2")
scatterplot3d(pp$x[, 1:3],color = ifelse(BCR == "0", "green","red"),xlab = "PC.1", ylab = "PC.2",zlab = "PC.3")
dev.off()

#end