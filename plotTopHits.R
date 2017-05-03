countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)


conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
colData

colDay3 <- subset(colData, day=="day03")

colDay3
# colInt<-subset(colDay3,location=="int")
# colSub <- subset(colDay3, location=="sub")
# 
# colInt
# colSub
# 
# countDataInt<-countData[, which(colnames(countData) %in% row.names(colInt))]
# countDataSub<-countData[, which(colnames(countData) %in% row.names(colSub))]
# dim(countDataInt)
# dim(countDataSub)

countDataDay3<-countData[, which(colnames(countData) %in% row.names(colDay3))]
head(countDataDay3)


ddsLoc<- DESeqDataSetFromMatrix(countData = countDataDay3, colData = colDay3 ,design = ~ location)

ddsLoc

ddsLoc <- ddsLoc[ rowSums(counts(ddsLoc)) > 100, ]

colData(ddsLoc)$location <- factor(colData(ddsLoc)$location, levels=c("int","sub"))

colData(ddsLoc)

ddsLoc <- DESeq(ddsLoc) 

names(ddsLoc)
str(ddsLoc)

resLoc <- results(ddsLoc)
resLoc



#sort results by padj
resLoc <- resLoc[order(resLoc$padj),]

head(resLoc)

summary(resLoc)

## Merge with normalized count data
resSumm <- merge(as.data.frame(resLoc), as.data.frame(counts(ddsLoc, normalized=TRUE)), by="row.names", sort=FALSE)
names(resSumm)[1] <- "Gene"
resSumm <- head(resSumm)
resSumm

row.names(resSumm) <- c("DN41041_c3_g2", "DN47102_c1_g1", "DN19042_c0_g1", "DN39048_c5_g1", "DN39233_c0_g1", "DN39328_c5_g1")

resSumm <- resSumm[,2:7]

resSumm

## Write results
write.csv(resSumm, file="diffexpr-results.csv")

plotCounts(ddsLoc, gene="TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543", intgroup="location")

vsdata <- vst(ddsLoc, blind=FALSE)

plotPCA(vsdata, intgroup="location")

#plot counts for top 6 genes

normcounts<-data.frame(genes=rownames(counts(ddsLoc, normalized=TRUE)),counts(ddsLoc, normalized=TRUE))
head(normcounts)

#names.to.keep <- as.factor(c("TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613", "TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581", "TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710", "TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316", "TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447 ", "TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543"))
#length(names.to.keep)
#rows.to.keep<-which(rownames(normcounts) %in% names.to.keep) 
normcounts2 <- normcounts[normcounts$genes %in% as.factor(rownames(head(resLoc))),]

normcounts2
dim(normcounts2)



install.packages("tidyr")

library(tidyr)

countsDF <- gather(normcounts,individuals,normcounts,I03_5.08_S_2:I38_6.12_H_0)
dim(countsDF)
head(countsDF)

summary(countsDF)

head(conds)
conds$individuals <- as.factor(rownames(conds))
str(conds)
library(dplyr)

topHits <- inner_join(countsDF,conds,by="individuals")

ggplot(topHits, aes(x=location,y=log(normcounts+1)))+geom_boxplot()+facet_grid(.~genes)

