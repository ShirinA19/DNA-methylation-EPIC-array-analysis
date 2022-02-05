# load packages required for analysis
library("limma")
library("minfi")
library("RColorBrewer")
library("missMethyl") # Can take a short time...
library("minfiData")
library("Gviz")
library("DMRcate")
library("DMRcatedata")
library("stringr")
library("mCSEA")

#BiocManager::install("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")

library("IlluminaHumanMethylationEPICanno.ilm10b2.hg19")


annepic <- getAnnotation(IlluminaHumanMethylationEPICanno.ilm10b2.hg19)
head(annepic)

#ann450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
# Use the head command to get a quick overview of the data and see what types of annotations are available
#head(ann450k)

dataDirectory <- "/Users/sher0001/Documents/Epic/Malaria_DNA_methylation/"
list.files(dataDirectory, recursive = TRUE)


# read in the sample sheet for the experiment
targets <- read.metharray.sheet(dataDirectory, pattern="BEA21P163_AM_sample_sheet.csv")
targets


# read in the raw data from the IDAT files; warnings can be ignored.
rgSet <- read.metharray.exp(targets=targets)


# Get an overview of the data
rgSet
pData(rgSet)
getManifest(rgSet)


# give the samples descriptive names
targets$ID <- paste(targets$Sample_Group,targets$Sample_Name,sep=".")
sampleNames(rgSet) <- targets$ID


# Check the names have been updated by looking at the rownames of the phenoData
pData(rgSet)

MSet <- preprocessRaw(rgSet)
MSet

# Compare to previous object
rgSet
head(getMeth(MSet)[,1:3])
head(getUnmeth(MSet)[,1:3])


ratioSet <- ratioConvert(MSet, what = "both", keepCN = TRUE)
# Observe the change of the assays
ratioSet

gset <- mapToGenome(ratioSet)
gset

## Get beta, m and copy number values
beta <- getBeta(gset)
head(beta)
m <- getM(gset)
head(m)
cn <- getCN(gset)
head(cn)

##Quality controls

qc <- getQC(MSet)
plotQC(qc)

# calculate the detection p-values
detP <- detectionP(rgSet)
head(detP)

# examine mean detection p-values across all samples to identify any failed samples
barplot(colMeans(detP), las=2, cex.names=0.8, ylab="Mean detection p-values")
abline(h=0.05,col="red")

phenoData <- pData(MSet)
densityPlot(MSet, sampGroups = phenoData$Sample_Group)


# normalize the data; this results in a GenomicRatioSet object
mSetSq <- preprocessQuantile(rgSet)


## Normalize by noob
mSetSq.noob <- preprocessFunnorm(rgSet)

#Compare with the unnormalized data to visualize the effect of the normalization. First a comparison of the Beta distributions for the different probe designs. This will give an indication of the effectiveness of the within-array normalization.

par(mfrow=c(1,2))
# Plot distributions prior to normalization for sample 1
plotBetasByType(MSet[,1],main="Raw")
# The normalized object is a GenomicRatioSet which does not contain
# the necessary probe info, we need to extract this from the MethylSet first.
typeI <- getProbeInfo(MSet, type = "I")[, c("Name","nCpG")]
typeII <- getProbeInfo(MSet, type = "II")[, c("Name","nCpG")]
probeTypes <- rbind(typeI, typeII)
probeTypes$Type <- rep(x = c("I", "II"), times = c(nrow(typeI), nrow(typeII)))
# Now plot the distributions of the normalized data for sample 1
plotBetasByType(getBeta(mSetSq)[,1], probeTypes = probeTypes, main="Normalized",)


# visualise what the data looks like before and after normalisation
par(mfrow=c(1,2))
densityPlot(rgSet, sampGroups=targets$Sample_Group,main="Raw", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))

densityPlot(getBeta(mSetSq), sampGroups=targets$Sample_Group,
            main="Normalized", legend=FALSE)
legend("top", legend = levels(factor(targets$Sample_Group)),
       text.col=brewer.pal(8,"Dark2"))


# MDS plots to look at largest sources of variation
# Create color panel
pal <- brewer.pal(8,"Dark2")
# Plot figures
par(mfrow=c(1,2))
plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.5,dim=c(1,2))
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       bg="white", cex=0.5)

plotMDS(getM(mSetSq), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)], cex=0.5,dim=c(1,2))
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       bg="white", cex=0.5)



# ensure probes are in the same order in the mSetSq and detP objects
detP <- detectionP(rgSet)
detP <- detP[match(featureNames(mSetSq),rownames(detP)),]

##For the noob

detP.noob <- detectionP(rgSet)
detP.noob <- detP.noob[match(featureNames(mSetSq.noob),rownames(detP.noob)),]


# remove any probes that have failed in one or more samples; this next line
# checks for each row of detP whether the number of values < 0.01 is equal
# to the number of samples (TRUE) or not (FALSE)
keep <- rowSums(detP < 0.01) == ncol(mSetSq)
table(keep)
##For noob
keep.noob <- rowSums(detP.noob < 0.01) == ncol(mSetSq.noob)

# Subset the GenomicRatioSet
mSetSqFlt <- mSetSq[keep,]
mSetSqFlt

##For noob
mSetSqFlt.noob <- mSetSq.noob[keep.noob,]

#There is a function in minfi that provides a simple interface for the removal of probes where common SNPs may affect the CpG. You can either remove all probes affected by SNPs (default), or only those with minor allele frequencies greater than a specified value.


mSetSqFlt <- dropLociWithSnps(mSetSqFlt)
mSetSqFlt




#For noob

mSetSqFlt.noob <- dropLociWithSnps(mSetSqFlt.noob)
mSetSqFlt.noob

##Check quality check after filtering and normalization

par(mfrow=c(1,2))
plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Group)], cex=0.5)
legend("top", legend=levels(factor(targets$Sample_Group)), text.col=pal,
       cex=0.5, bg="white")

plotMDS(getM(mSetSqFlt), top=1000, gene.selection="common",
        col=pal[factor(targets$Sample_Source)],cex=0.5)
legend("top", legend=levels(factor(targets$Sample_Source)), text.col=pal,
       cex=0.5, bg="white")
# Close double plotting window
dev.off()

# calculate M-values for statistical analysis: as previously mentioned, M-values have nicer statistical properties and are thus better for use in statistical analysis of methylation data
mVals <- getM(mSetSqFlt)

write.csv(mVals,file = "mVals.csv")

##For noob
mVals.noob <- getM(mSetSqFlt.noob)

# Set up the design matrix for the Differential Methylation analysis
# Define the factor of interest
Type <- factor(targets$Sample_Group)
# Define is the individual effect that we need to account for
individual <- factor(targets$Sample_Source)
# use the above to create a design matrix
#design1 <- model.matrix(~0+Type+individual, data=targets)
#colnames(design1) <- c(levels(Type),levels(individual)[-1])


#design <- model.matrix(~0+Type, data=targets)
#colnames(design) <- c(levels(Type))


design <- model.matrix(~0+Type+individual, data=targets)
colnames(design) <- c(levels(Type),levels(individual)[-1])


##To check the group in the design
make.names(colnames(design))


##For Mattias
#design2<-design
#IndMatrix<-matrix(0,nrow=dim(design2)[1],ncol=8)
#colnames(IndMatrix)<-paste("Ind",1:8,sep="")
#design2<-cbind(design2,IndMatrix)
#colnames(design2)<-c(colnames(design2),paste("Ind",1:8,sep=""))
#for(Ind in 1:8)
#{
# Col<-4+Ind
# Rows<-1:4+(Ind-1)*4
#design2[Rows,Col]<-1
#}
#design3<-design2
#design3<-design3[,-5]
#design3<-cbind(design2,sample(c(0,1),32,replace=TRUE))
#design3<-cbind(design3,design2[,5])
#colnames(design3)[13]<-"Ind1A"


# fit the actual linear model to the data
fit <- lmFit(mVals, design)


# create a contrast matrix for specific comparisons
contMatrix <- makeContrasts(control-before_treatment,
                            control-after_treatment,
                            before_treatment-after_treatment,
                            control-recovery,
                            levels=design)
contMatrix

contMatrix1 <- makeContrasts(before_treatment-control,
                             after_treatment-control,
                             after_treatment-before_treatment,
                             recovery-control,
                             before_treatment-recovery,
                             after_treatment-recovery,
                             levels=design)
contMatrix1





# fit the contrasts
fit2 <- contrasts.fit(fit, contMatrix)
# Rank genes
fit2 <- eBayes(fit2)


# fit the contrasts
fit3 <- contrasts.fit(fit, contMatrix1)
# Rank genes
fit3 <- eBayes(fit3)


#get the table of results from contrasts

DMPs1.contmat1 <- topTable(fit3, num=Inf, coef=1)

DMPs1.contmat1.sig<-subset(DMPs1.contmat1,adj.P.Val<0.05)

write.csv(DMPs1.contmat1.sig,file = "DMPs1.contmat1.sig.csv")

DMPs2.contmat1 <- topTable(fit3, num=Inf, coef=2)

DMPs2.contmat1.sig<-subset(DMPs2.contmat1,adj.P.Val<0.05)
write.csv(DMPs2.contmat1.sig,file = "DMPs2.contmat1.sig.csv")


DMPs3.contmat1 <- topTable(fit3, num=Inf, coef=3)

DMPs3.contmat1.sig<-subset(DMPs3.contmat1,adj.P.Val<0.05)

DMPs4.contmat1 <- topTable(fit3, num=Inf, coef=4)

DMPs4.contmat1.sig<-subset(DMPs4.contmat1,adj.P.Val<0.05)

DMPs5.contmat1 <- topTable(fit3, num=Inf, coef=5)
DMPs5.contmat1.sig<-subset(DMPs5.contmat1,adj.P.Val<0.05)

write.csv(DMPs5.contmat1.sig,file = "DMPs5.contmat1.sig.csv")

DMPs6.contmat1 <- topTable(fit3, num=Inf, coef=6)
DMPs6.contmat1.sig<-subset(DMPs6.contmat1,adj.P.Val<0.05)


# get the table of results for the first contrast, control-before_treatment
DMPs1 <- topTable(fit2, num=Inf, coef=1)
head(DMPs1)
DMPs1.sig<-subset(DMPs1,adj.P.Val<0.05)

write.csv(DMPs1.sig,file = "DMPs1.sig.csv")




write.csv(DMPs1,file = "DMPs1.csv")

DMPs2 <- topTable(fit2, num=Inf, coef=2)

DMPs2.sig<-subset(DMPs2,adj.P.Val<0.05)

write.csv(DMPs2.sig,file = "DMPs2.sig.csv")

write.csv(DMPs2,file = "DMPs2.csv")


DMPs3 <- topTable(fit2, num=Inf, coef=3)
DMPs3.sig<-subset(DMPs3,adj.P.Val<0.05)
write.csv(DMPs3,file = "DMPs3.csv")

DMPs4 <- topTable(fit2, num=Inf, coef=4)
DMPs4.sig<-subset(DMPs4,adj.P.Val<0.05)


write.csv(DMPs4,file = "DMPs4.csv")

# Retrieve data from the array annotation package; this is array-specific
annepicSub <- annepic[match(rownames(mVals),annepic$Name),
                      c(1:4,12:19,24:ncol(annepic))]
DMPs1.contmat1.ann <- topTable(fit3, num=Inf, coef=1, genelist=annepicSub)
head(DMPs1.contmat1.ann)

write.table(DMPs1.contmat1.ann, file="DMPs1.contmat1.ann.10.11.21.csv", sep=",", row.names=FALSE)






DMPs2.contmat1.ann <- topTable(fit3, num=Inf, coef=2, genelist=annepicSub)

DMPs5.contmat1.ann <- topTable(fit3, num=Inf, coef=5, genelist=annepicSub)

DMPs6.contmat1.ann <- topTable(fit3, num=Inf, coef=6, genelist=annepicSub)


write.table(DMPs1.ann, file="DMPs1.ann.new.csv", sep=",", row.names=FALSE)


# The resulting data.frame can easily be written to a CSV file, which can be opened in Excel.
# write.table(DMPs, file="DMPs.csv", sep=",", row.names=FALSE)

# eXtract Beta-values
bVals <- getBeta(mSetSqFlt)

write.csv(bVals,file = "bVals.csv")

##For noob
bVals.noob <- getBeta(mSetSqFlt.noob)


BiocManager::install("ENmix")

library(ENmix)
##Methyl age
mage=methyAge(bVals)

mage.noob=methyAge(bVals.noob)


write.csv(mage.noob,file = "methylage.csv")



# Plot most significant differentially methylated CpG
#Plus logFC

plotCpg(bVals, cpg="cg07361629", pheno=targets$Sample_Group, ylab = "Beta values")

##minus logFC
plotCpg(bVals, cpg="cg10773372", pheno=targets$Sample_Group, ylab = "Beta values")

#Plus logFC

plotCpg(bVals, cpg="cg05304729", pheno=targets$Sample_Group, ylab = "Beta values")

##minus logFC
plotCpg(bVals, cpg="cg05332943", pheno=targets$Sample_Group, ylab = "Beta values")



##
myAnnotation1 <- cpg.annotate(object = mVals,
                              datatype = "array",
                              what = "M",
                              analysis.type = "differential",
                              design = design,
                              contrasts = TRUE,
                              cont.matrix = contMatrix1,
                              coef = "before_treatment - control" ,
                              arraytype = "EPIC")
myAnnotation1

DMRs.ann1 <- dmrcate(myAnnotation1, lambda=1000, C=2)
DMRs.ann1




# Create GRanges object; create directory when prompted
results.ranges1 <- extractRanges(DMRs.ann1)
results.ranges1

write.table(results.ranges1, file="results.ranges1.new.csv", sep=",", row.names=FALSE)



myAnnotation2 <- cpg.annotate(object = mVals,
                              datatype = "array",
                              what = "M",
                              analysis.type = "differential",
                              design = design,
                              contrasts = TRUE,
                              cont.matrix = contMatrix1,
                              fdr = 0.05,
                              coef = "after_treatment - control" ,
                              arraytype = "EPIC")



DMRs.ann2 <- dmrcate(myAnnotation2, lambda=1000, C=2)
#DMRs.ann2.test <- dmrcate(myAnnotation2, lambda=500, C=2)
DMRs



# Create GRanges object; create directory when prompted
results.ranges2.test <- extractRanges(DMRs.ann2.test)
results.ranges2

results.ranges2 <- extractRanges(DMRs.ann2)

write.table(results.ranges2, file="results.ranges2.new.csv", sep=",", row.names=FALSE)
write.table(results.ranges2.test, file="results.ranges2.test.fdr.0.1.csv", sep=",", row.names=FALSE)



myAnnotation5 <- cpg.annotate(object = mVals,
                              datatype = "array",
                              what = "M",
                              analysis.type = "differential",
                              design = design,
                              contrasts = TRUE,
                              cont.matrix = contMatrix1,
                              coef = "before_treatment - recovery" ,
                              arraytype = "EPIC")

DMRs.ann5 <- dmrcate(myAnnotation5, lambda=1000, C=2)


results.ranges5 <- extractRanges(DMRs.ann5)

write.table(results.ranges5, file="results.ranges5.new.csv", sep=",", row.names=FALSE)



myAnnotation6 <- cpg.annotate(object = mVals,
                              datatype = "array",
                              what = "M",
                              analysis.type = "differential",
                              design = design,
                              contrasts = TRUE,
                              cont.matrix = contMatrix1,
                              coef = "after_treatment - recovery" ,
                              arraytype = "EPIC")

DMRs.ann6 <- dmrcate(myAnnotation6, lambda=1000, C=2)


results.ranges6 <- extractRanges(DMRs.ann6)

write.table(results.ranges6, file="results.ranges6.new.csv", sep=",", row.names=FALSE)








#Just as for the single CpG analysis, it is a good idea to visually inspect the results to make sure they make sense. For this, use the DMR.plot function. By default, this plot draws the location of the DMR in the genome, the position of nearby genes, the positions of the CpG probes, the Beta value levels of each sample as a heatmap and the mean methylation levels for the various sample groups in the experiment.

# set up the grouping variables and colours
pal <- brewer.pal(8,"Dark2")
groups <- pal[1:length(unique(targets$Sample_Group))]
names(groups) <- levels(factor(targets$Sample_Group))
cols <- groups[as.character(factor(targets$Sample_Group))]
# draw the plot for the second DMR - first gives error for some reason...
DMR.plot(ranges = results.ranges1,
         dmr = 2,
         CpGs = mSetSqFlt,
         phen.col = cols,
         genome = "hg19")


DMR.plot(ranges = results.ranges2,
         dmr = 1,
         CpGs = mSetSqFlt,
         phen.col = cols,
         genome = "hg19")


DMR.plot(ranges = results.ranges5,
         dmr = 2,
         CpGs = mSetSqFlt,
         phen.col = cols,
         genome = "hg19")

DMR.plot(ranges = results.ranges6,
         dmr = 2,
         CpGs = mSetSqFlt,
         phen.col = cols,
         genome = "hg19")



# Create a named vector containing the rank metric (here: logFC)
myRank <- DMPs1.contmat1.ann$logFC
names(myRank) <- rownames(DMPs1)

# Reshape the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
myResults <- mCSEATest(myRank,
                       bVals,
                       pheno,
                       regionsTypes = "promoters",
                       platform = "EPIC")
head(myResults$promoters)



write.table(myResults$promoters, file="myResults.mCSEA.DMPs1.csv", sep=",", row.names=TRUE)

write.table(myResults, file="myResults.mCSEA.DMPs1.csv")

mCSEAPlot(myResults,
          regionType = "promoters",
          dmrName = "AURKB",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)


# Run the mCSEA for genes
myResults.gene <- mCSEATest(myRank,
                            bVals,
                            pheno,
                            regionsTypes = "genes",
                            platform = "EPIC")
head(myResults.gene$genes)


write.table(myResults.gene$genes, file="myResults.gene.mCSEA.DMPs1.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.gene,
          regionType = "genes",
          dmrName = "WARS",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)


# Run the mCSEA for CpG
myResults.CGI <- mCSEATest(myRank,
                           bVals,
                           pheno,
                           regionsTypes = "CGI",
                           platform = "EPIC")
head(myResults.CGI$CGI)

write.table(myResults.CGI$CGI, file="myResults.CGI.mCSEA.DMPs1.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.CGI,
          regionType = "CGI",
          dmrName = "chrX:48754500-48755329",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)

mCSEAPlot(myResults.CGI, regionType, dmrName, extend = 1000,
          chromosome = TRUE, leadingEdge = TRUE, CGI = FALSE, genes = TRUE,
          transcriptAnnotation = "transcript", makePDF = TRUE,
          col = c("blue", "magenta", "green", "red", "black"))


##new test with DMPs1.ann
# Create a named vector containing the rank metric (here: logFC)
myRank1 <- DMPs1.contmat1.ann$logFC
names(myRank1) <- rownames(DMPs1.contmat1.ann)

# Reshape the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
myResults.promoters.ann <- mCSEATest(myRank1,
                                     bVals,
                                     pheno,
                                     regionsTypes = "promoters",
                                     platform = "EPIC")
head(myResults.promoters.ann$promoters)


write.table(myResults.promoters.ann$promoters, file="myResults.promoters.ann.mCSEA.DMPs1.contmat1.ann.new.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.promoters.ann,
          regionType = "promoters",
          dmrName = "PCNA",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)


myResults.genes.ann <- mCSEATest(myRank1,
                                 bVals,
                                 pheno,
                                 regionsTypes = "genes",
                                 platform = "EPIC")

head(myResults.genes.ann$genes)

write.table(myResults.promoters.ann$genes, file="myResults.genes.ann.mCSEA.DMPs1.contmat1.ann.new.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.genes.ann,
          regionType = "genes",
          dmrName = "PCNA",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)



myResults.CGI.ann <- mCSEATest(myRank1,
                               bVals,
                               pheno,
                               regionsTypes = "CGI",
                               platform = "EPIC")

head(myResults.CGI.ann$CGI)


write.table(myResults.CGI.ann$CGI, file="myResults.CGI.ann.mCSEA.DMPs1.contmat1.ann.new.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.CGI.ann,
          regionType = "CGI",
          dmrName = "chrX:48754500-48755329",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)


##2nd contrast


##new test with DMPs1.ann
# Create a named vector containing the rank metric (here: logFC)
myRank2 <- DMPs2.contmat1.ann$logFC
names(myRank2) <- rownames(DMPs2.contmat1.ann)

# Reshape the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
myResults.myRank2.promoters.ann <- mCSEATest(myRank2,
                                             bVals,
                                             pheno,
                                             regionsTypes = "promoters",
                                             platform = "EPIC")
head(myResults.myRank2.promoters.ann$promoters)



write.table(myResults.myRank2.promoters.ann$promoters, file="myResults.myRank2.promoters.ann.mCSEA.DMPs2.contmat1.ann.new.csv", sep=",", row.names=TRUE)

mCSEAPlot(myResults.myRank2.promoters.ann,
          regionType = "promoters",
          dmrName = "PCNA",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)


myResults.myRank2.genes.ann <- mCSEATest(myRank2,
                                         bVals,
                                         pheno,
                                         regionsTypes = "genes",
                                         platform = "EPIC")

head(myResults.myRank2.genes.ann$genes)

write.table(myResults.myRank2.genes.ann$genes, file="myResults.myRank2.genes.ann.mCSEA.DMPs2.contmat1.ann.new.csv", sep=",", row.names=TRUE)


mCSEAPlot(myResults.myRank2.genes.ann,
          regionType = "genes",
          dmrName = "PCNA",
          transcriptAnnotation = "symbol",
          makePDF = FALSE)



myResults.myRank2.CGI.ann <- mCSEATest(myRank2,
                                       bVals,
                                       pheno,
                                       regionsTypes = "CGI",
                                       platform = "EPIC")

head(myResults.myRank2.CGI.ann$CGI)

write.table(myResults.myRank2.CGI.ann$CGI, file="myResults.myRank2.CGI.ann.ann.mCSEA.DMPs2.contmat1.ann.new.csv", sep=",", row.names=TRUE)


##5th contrast


##new test with DMPs5.ann
# Create a named vector containing the rank metric (here: logFC)
myRank5 <- DMPs5.contmat1.ann$logFC
names(myRank5) <- rownames(DMPs5.contmat1.ann)

# Reshape the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
myResults.myRank5.promoters.ann <- mCSEATest(myRank5,
                                             bVals,
                                             pheno,
                                             regionsTypes = "promoters",
                                             platform = "EPIC")
head(myResults.myRank5.promoters.ann$promoters)



write.table(myResults.myRank5.promoters.ann$promoters, file="myResults.myRank5.promoters.ann.mCSEA.DMPs2.contmat1.ann.new.csv", sep=",", row.names=TRUE)


##new test with DMPs6.ann
# Create a named vector containing the rank metric (here: logFC)
myRank6 <- DMPs6.contmat1.ann$logFC
names(myRank6) <- rownames(DMPs6.contmat1.ann)

# Reshape the phenotype data to a format suitable for mCSEA
pheno <- as.data.frame(pData(mSetSqFlt))
pheno <- pheno[,"Sample_Group", drop=FALSE]

# Run the mCSEA
myResults.myRank6.promoters.ann <- mCSEATest(myRank6,
                                             bVals,
                                             pheno,
                                             regionsTypes = "promoters",
                                             platform = "EPIC")
head(myResults.myRank6.promoters.ann$promoters)



write.table(myResults.myRank6.promoters.ann$promoters, file="myResults.myRank6.promoters.ann.mCSEA.DMPs2.contmat1.ann.new.csv", sep=",", row.names=TRUE)




#Gene Ontology Testing
# Get the significant CpG sites at less than 5% FDR
sigCpGs <- DMPs1.contmat1.ann$Name[DMPs1.contmat1.ann$adj.P.Val<0.05]

write.table(sigCpGs, file="DMPs1.contmat1.ann.csv", sep=",", row.names=TRUE)


# First 10 significant CpGs
sigCpGs[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs)
# Get all the CpG sites used in the analysis to form the background
all <- DMPs1.contmat1.ann$Name
# Total number of CpG sites tested
length(all)


# Run enrichment - Can take a bit of time...
gst <- gometh(sig.cpg=sigCpGs, all.cpg=all)

write.table(gst, file="gst.DMPs1.contmat1.ann.csv", sep=",", row.names=TRUE)

# Top 10 GO categories
topGSA(gst, number=10)

#Contrast2
sigCpGs2 <- DMPs2.contmat1.ann$Name[DMPs2.contmat1.ann$adj.P.Val<0.05]

write.table(sigCpGs2, file="DMPs2.contmat1.ann.csv", sep=",", row.names=TRUE)


# First 10 significant CpGs
sigCpGs2[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs2)
# Get all the CpG sites used in the analysis to form the background
all2 <- DMPs2.contmat1.ann$Name
# Total number of CpG sites tested
length(all2)


# Run enrichment - Can take a bit of time...
gst2 <- gometh(sig.cpg=sigCpGs2, all.cpg=all2)

write.table(gst2, file="gst2.DMPs2.contmat1.ann.csv", sep=",", row.names=TRUE)

# Top 10 GO categories
topGSA(gst2, number=10)

#Contrast5
sigCpGs5 <- DMPs5.contmat1.ann$Name[DMPs5.contmat1.ann$adj.P.Val<0.05]

write.table(sigCpGs5, file="DMPs5.contmat1.ann.csv", sep=",", row.names=TRUE)


# First 10 significant CpGs
sigCpGs5[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs5)
# Get all the CpG sites used in the analysis to form the background
all5 <- DMPs5.contmat1.ann$Name
# Total number of CpG sites tested
length(all5)


# Run enrichment - Can take a bit of time...
gst5 <- gometh(sig.cpg=sigCpGs5, all.cpg=all5)

write.table(gst5, file="gst5.DMPs5.contmat1.ann.csv", sep=",", row.names=TRUE)

# Top 10 GO categories
topGSA(gst5, number=10)


#Contrast6
sigCpGs6 <- DMPs6.contmat1.ann$Name[DMPs6.contmat1.ann$adj.P.Val<0.05]

write.table(sigCpGs6, file="DMPs6.contmat1.ann.csv", sep=",", row.names=TRUE)


# First 10 significant CpGs
sigCpGs6[1:10]
# Total number of significant CpGs at 5% FDR
length(sigCpGs6)
# Get all the CpG sites used in the analysis to form the background
all6 <- DMPs6.contmat1.ann$Name
# Total number of CpG sites tested
length(all6)


# Run enrichment - Can take a bit of time...
gst6 <- gometh(sig.cpg=sigCpGs6, all.cpg=all6)

write.table(gst6, file="gst6.DMPs5.contmat1.ann.csv", sep=",", row.names=TRUE)

# Top 10 GO categories
topGSA(gst6, number=10)
