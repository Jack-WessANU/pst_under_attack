if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("edgeR")
BiocManager::install("limma")
BiocManager::install("tximport")
#importing the Salmon files
library(edgeR)
library(limma)
library(tximport)


penx_attack1 <-tximport("./pst_salmon_inpla1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
penx_attack2 <-tximport("./pst_salmon_inpla2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
penx_attack3 <-tximport("./pst_salmon_inpla3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day10_1 <-tximport("./pst_salmon_10days1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day10_2 <-tximport("./pst_salmon_10days2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day10_3 <-tximport("./pst_salmon_10days3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day12_1 <-tximport("./pst_salmon_12days1/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day12_2 <-tximport("./pst_salmon_12days2/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)
pst_alone_day12_3 <-tximport("./pst_salmon_12days3/quant.sf", type = c("salmon"), txIn = TRUE, txOut = TRUE, countsFromAbundance = c("lengthScaledTPM"),dropInfReps=TRUE)

#combining the Salmon count information into a single matrix
countmatrix<-data.frame(penx_attack1$counts, penx_attack2$counts, penx_attack3$counts, pst_alone_day12_1$counts, pst_alone_day12_2$counts, pst_alone_day12_3$counts)
row.names(countmatrix)<-row.names(penx_attack1$counts)
plotMDS(countmatrix)

#import design matrix (generated seperately in excel)
design<-read.csv("./design.csv")

#filtering genes by expression
dge <- DGEList(counts=countmatrix)
group <- as.factor(c("penxattack","penxattack","penxattack","day12","day12","day12"))
keep <- filterByExpr(dge, group=group)
dge <- dge[keep,,keep.lib.sizes=FALSE]

#checking the results of filtering
dim(dge)
dim(design)

#Use voom to determine weights for each gene to be passed to limma, voom normalises the data
v<-voom(dge, design, plot=TRUE, normalize="quantile")

#fit linear model for each gene
fit <- lmFit(v,design)

#define the comparisons you want to examine
contrast.matrix <- makeContrasts(penxattack - day12, levels = design)

#estimate contrasts for each gene
fit2 <- contrasts.fit(fit, contrast.matrix)

#Smoothing of standard error
fit2 <- eBayes(fit2)

#Generate results tables
pst_underattackbypenx_vs_day12 <- topTable(fit2,coef=1, number=Inf)

#Export transcriptome profiles
write.csv(pst_underattackbypenx_vs_day12, file = "pst_underattackbypenx_vs_day12.csv")

pst_underattackbypenx_vs_day12
