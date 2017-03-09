##Differential Expression Assignment
##Ecological Genomics
##Lauren Ashlock
##3_06_2017


library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)
ncol(countData)
dim(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)
dim(colData)


#comparing healthy vs sick individuals within the intertidal


colInt<-subset(colData,location=="int")
colSub <- subset(colData, location=="sub")

nrow(colInt)
#output
#48

countDataInt<-countData[, which(colnames(countData) %in% row.names(colInt))]
countDataSub<-countData[, -which(colnames(countData) %in% row.names(colInt))]
dim(countDataInt)
dim(countDataSub)

ddsInt <- DESeqDataSetFromMatrix(countData = countDataInt, colData = colInt ,design = ~ health)

ddsInt <- ddsInt[ rowSums(counts(ddsInt)) > 100, ]

colData(ddsInt)$health <- factor(colData(ddsInt)$health, levels=c("H","S"))

ddsInt <- DESeq(ddsInt) 

resInt <- results(ddsInt)

resInt <- resInt[order(resInt$padj),]

head(resInt)

summary(resInt)

#output
# out of 12399 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 205, 1.7% 
# LFC < 0 (down)   : 37, 0.3% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 9381, 76% 
# (mean count < 41)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results


#Comparing healthy vs sick individuals in subtidal

ddsSub <- DESeqDataSetFromMatrix(countData = countDataSub, colData = colSub ,design = ~ health)

ddsSub <- ddsSub[ rowSums(counts(ddsSub)) > 100, ]

colData(ddsSub)$health <- factor(colData(ddsSub)$health, levels=c("H","S"))

ddsSub <- DESeq(ddsSub) 

resSub <- results(ddsSub)

resSub <- resSub[order(resSub$padj),]

head(resSub)

summary(resSub)


#output
# out of 12392 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 20, 0.16% 
# LFC < 0 (down)   : 113, 0.91% 
# outliers [1]     : 647, 5.2% 
# low counts [2]   : 4289, 35% 
# (mean count < 13)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#Comparing healthy vs sick controlling for location 

ddsLoc <- DESeqDataSetFromMatrix(countData = countData, colData = colData ,design = ~ location + health)

ddsLoc <- ddsLoc[ rowSums(counts(ddsLoc)) > 100, ]

colData(ddsLoc)$health <- factor(colData(ddsLoc)$health, levels=c("H","S"))

ddsLoc <- DESeq(ddsLoc) 

resLoc <- results(ddsLoc)

resLoc <- resLoc[order(resLoc$padj),]

head(resLoc)

summary(resLoc)

#output
# out of 12947 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 209, 1.6% 
# LFC < 0 (down)   : 65, 0.5% 
# outliers [1]     : 400, 3.1% 
# low counts [2]   : 7679, 59% 
# (mean count < 23)
# [1] see 'cooksCutoff' argument of ?results
# [2] see 'independentFiltering' argument of ?results

#Table of up/down regulated genes in three models
vec<- c(205, 37, 20, 113, 209, 65)

summaryMatrix<- matrix(data=vec, nrow=2, ncol=3, byrow=FALSE)
summaryMatrix
row.names(summaryMatrix)<- c("Significantly Upregulated", "Significantly Downregulated")
summaryMatrix
colnames(summaryMatrix) <- c("Intertidal", "Subtidal", "Full Model")
summaryMatrix

print(summaryMatrix, type="html")

#PCA for all three models

#Intertidal

vsdInt <- varianceStabilizingTransformation(ddsInt, blind=FALSE)

plotPCA(vsdInt, intgroup=c("health"))

#Subtidal

vsdSub <- varianceStabilizingTransformation(ddsSub, blind=FALSE)

plotPCA(vsdSub, intgroup=c("health"))

#Location

vsdLoc <- varianceStabilizingTransformation(ddsLoc, blind=FALSE)

plotPCA(vsdLoc, intgroup=c("health"))

##Plot counts

#Intertidal

## Check out one of the genes to see interaction between score, health and expression....
dInt <- plotCounts(ddsInt, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
dInt
pInt <- ggplot(dInt, aes(x= score, y=count, color = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pInt <- pInt + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pInt

pInt <- ggplot(dInt, aes(x=score, y=count, color=health, group=health)) 
pInt <- pInt +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pInt

#Subtidal

dSub <- plotCounts(ddsSub, gene="TRINITY_DN42073_c0_g1_TRINITY_DN42073_c0_g1_i1_g.12173_m.12173", intgroup=(c("health","score","location")), returnData=TRUE)
dSub
pSub <- ggplot(dSub, aes(x= score, y=count, color = health)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pSub <- pSub + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pSub

pSub <- ggplot(dSub, aes(x=score, y=count, color=health, group=health)) 
pSub <- pSub +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pSub

#Full model

dLoc <- plotCounts(ddsLoc, gene="TRINITY_DN43080_c1_g1_TRINITY_DN43080_c1_g1_i3_g.14110_m.14110", intgroup=(c("health","score","location")), returnData=TRUE)
dLoc
pLoc <- ggplot(dLoc, aes(x= score, y=count, shape = health, color = location)) + theme_minimal() + theme(text = element_text(size=20), panel.grid.major = element_line(colour = "grey"))
pLoc <- pLoc + geom_point(position=position_jitter(w=0.3,h=0), size = 3) 
pLoc

pLoc <- ggplot(dLoc, aes(x=score, y=count, color=health, group=health)) 
pLoc <- pLoc +  geom_point() + stat_smooth(se=FALSE,method="loess") +  scale_y_log10()
pLoc

citation("DESeq2")
