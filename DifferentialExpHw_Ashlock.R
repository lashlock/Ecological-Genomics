##Differential Expression Assignment
##Ecological Genomics
##Lauren Ashlock
##3_06_2017


library("DESeq2")

library("ggplot2")

countsTable <- read.delim('countsdata_trim2.txt', header=TRUE, stringsAsFactors=TRUE, row.names=1)
countData <- as.matrix(countsTable)
head(countData)

conds <- read.delim("cols_data_trim.txt", header=TRUE, stringsAsFactors=TRUE, row.names=1)
head(conds)
colData <- as.data.frame(conds)
head(colData)


#comparing healthy vs sick individuals within the intertidal

colInt<-subset(colData,location=="int")

nrow(colInt)
#output
#47

countInt <- countData[,1:48]

head(countInt)
nrow(countInt)


intdds <- DESeqDataSetFromMatrix(countData = countInt, colData = colInt ,design = ~ health)





