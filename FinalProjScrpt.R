
#FinalProj Script

#load in libraries

library("phyloseq")
library("DESeq2")
library(XML)
library("ggplot2")
theme_set(theme_bw())

#Import the OTU table
otutable <- import_biom(BIOMfilename = 'otu_table_mc2_w_tax_no_pynast_failures_no_chimeras_frequency_filtered.biom', 
                        treefilename = 'rep_set_no_chimeras.tre', 
                        parseFunction = parse_taxonomy_greengenes)

#The warnings are ok. There is 1 warning for every OTU that doesn't have a taxonomy assignment

#Import the mapping file
mapping <- import_qiime_sample_data(mapfilename = 'R_map.txt')

#Merge the mapping file to the OTU table

phylo <- merge_phyloseq(otutable, mapping)


#Check to make sure the imports worked
phylo

head(otu_table(phylo))

tail(otu_table(phylo))

phylo_subset = subset_samples(phylo, Day == 3 & Phenotype !="Dead")
phylo_int = subset_samples(phylo, Day==3 & Tide=="intertidal")
phylo_int
phylo_sub = subset_samples(phylo, Day==3 & Tide=="subtidal")
phylo_sub
phylo_subset
sample_data(phylo_subset)
sample_data(phylo_int)
sample_data(phylo_sub)

##We want to numbers of each individual to be a factor. Right now it's an integer so we have to change that.
class(sample_data(phylo_subset)$individual)
sample_data(phylo_subset)$individual<-factor(sample_data(phylo_subset)$individual)

##Phyloseq's wrapper to get OTU data into DESeq

loc <- phyloseq_to_deseq2(phylo_subset, ~  Tide)
names(loc)

##Run DESeq. This command takes a bit of time

loc_deseq_test <- DESeq(loc, test="Wald")

str(loc_deseq_test)
##Get results from DESeq

loc_results <- results(loc_deseq_test)

##Have a look at the summary
summary(loc_results)
head(loc_results)


#Subtidal is the reference
# out of 1063 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)     : 73, 6.9% 
# LFC < 0 (down)   : 52, 4.9% 
# outliers [1]     : 0, 0% 
# low counts [2]   : 248, 23% 


##Make table with OTUs padj<0.05
alpha <- 0.05

loc_sigtab <- loc_results[which(loc_results$padj < alpha), ]

loc_sigtab

##Add taxa info to that table

loc_sigtab <- cbind(as(loc_sigtab, "data.frame"), as(tax_table(phylo)[rownames(loc_sigtab), ], "matrix"))
tail(loc_sigtab)
loc_sigtab

##Lets look at the things that look different

sphingomonadales <- subset(loc_sigtab, Order=="Sphingomonadales")
sphingomonadales

oceanospirillales <- subset(loc_sigtab, Order=="Oceanospirillales")
oceanospirillales

vibrionales <- subset(loc_sigtab, Order=="Vibrionales")
vibrionales


pseudomonadales <- subset(loc_sigtab, Order=="Pseudomonadales")
pseudomonadales



#sort by padj
loc_sigtab <- loc_sigtab[order(resLoc$padj),1:12]
head(loc_sigtab)

loc_sigtab <- loc_sigtab[1:6, 1:12]
loc_sigtab

##Save the table to your desktop

write.table(loc_sigtab, "DE_OTU_intertidal_vs_subtidal.txt", sep="\t")

###Diversity analyses

richness.df <- estimate_richness(phylo_subset, split=TRUE, measures=NULL)

str(richness.df)
richness.df
#make two separate data frames to run richness analyses separately, add on a column denoting tide and then merge back together

int_rich.df <- estimate_richness(phylo_int, split=TRUE, measures=NULL)
int_rich.df$Tide <-as.factor("Intertidal")
str(int_rich.df)
int_rich.df

sub_rich.df <- estimate_richness(phylo_sub, split=TRUE, measures=NULL)
sub_rich.df$Tide <-as.factor("Subtidal")
str(sub_rich.df)
sub_rich.df

richByTide.df <- rbind(int_rich.df, sub_rich.df)

richByTide.df

#removing standard error columns

richByTide.df$se.chao1 <- NULL

richByTide.df$se.ACE <- NULL

richByTide.df

# t1<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$p.value})
# t2<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$estimate})
# t3<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$conf.int[1:2]})
# #t1<-t.test(Fisher~Tide,data=richByTide.df)$p.value
# #t2<-t.test(Fisher~Tide,data=richByTide.df)$estimate
# #t3<-t.test(Fisher~Tide,data=richByTide.df)$conf.int[1:2]
# 
# data.frame(cbind(Pvalue=t1,Estimate=t2,ConfidenceInterval=t3))
# data.frame(cbind(t1,t2,t3))


t1<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$p.value})
t2<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$estimate})
t3<-apply(richByTide.df[,-8],2, function(x){t.test(x~richByTide.df$Tide)$conf.int[1:2]})

str(t1)
t2
t3
t1 <- as.data.frame(t1)
t1

t1 <- t(as.data.frame(t1))
t2 <- t(as.data.frame(t2))
t3 <- t(as.data.frame(t3))

t1
t2
t3

SummTable <- rbind(t1,t2,t3)

SummTable

row.names(SummTable) <- c("pvalue", "Mean in Intertidal", "Mean in Subtidal", "LCI", "UCI")

SummTable

write.table(SummTable, "RichnessSummary.txt", sep="\t")

########################################################

##Plot abundances of taxa in different pheno numbers
##First we must rarefy
set.seed(28132)
phyloR = rarefy_even_depth(phylo_subset, sample.size = min(sample_sums(phylo_subset)))

#Then we check to make sure we rarefied corectly.
title = "Sum of reads for each sample, phyloR"
plot(sort(sample_sums(phyloR), TRUE), type = "h", main = title, ylab = "reads", 
     ylim = c(0, 20000))


##We merge on basis of pheno number so we can make relative abundance graph
#phyloRm = merge_samples(phyloR, "Pheno_num")

##Merge on basis of phenotype
phyloRmD = merge_samples(phyloR, "Tide")

phyloRmD

phyloRt = transform_sample_counts(phyloRmD, function(x) 100 * x/sum(x))


p = plot_bar(phyloRt,x="Tide",facet_grid=NULL, fill="Order")
p + geom_bar(aes(color=Order, fill=Order), stat="identity", position="stack")

########Potential Ordination Plot#########




###################DGE Analysis###################

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
dim(countDataDay3)

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

normcounts2$genes <- c("DN41041_c3_g2", "DN47102_c1_g1", "DN19042_c0_g1", "DN39048_c5_g1", "DN39233_c0_g1", "DN39328_c5_g1")
dim(normcounts2)
str(normcounts2)
normcounts2

install.packages("tidyr")

library(tidyr)

countsDF <- gather(normcounts2,individuals,normcounts,I03_5.08_S_2:I38_6.12_H_0)
dim(countsDF)
head(countsDF)

summary(countsDF)

head(conds)
conds$individuals <- as.factor(rownames(conds))
str(conds)
library(dplyr)

topHits <- inner_join(countsDF,conds,by="individuals")
topHits
dim(topHits)
str(topHits)

ggplot(topHits, aes(x=location,y=log(normcounts+1), color=location))+
  geom_boxplot()+facet_grid(.~genes)+theme(legend.position="none")

#####################Annotation of topHits#############################
#Blastx results
#>TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613 type:5prime_partial len:413 TRINITY_DN41041_c3_g2_i1:2720-1482(-)
#CCTAATGAGGTGCTTGTGCGCATGGACCCAAGTGTTAAGCAGCCCAAGTGTGTCAAGAAAGAGTTTCCTCCATGCGTTTACAGAGATCGGGCCTGTCCCAATTATCGAGTTGTGGACACATTGACGGATAACATTGAGCATCGTCAGTTAGCTGCTAGCCATTACGTCATGAAGAGGACTCGAACCTGTAACATGAGCCAGGCATACGATGATGCATATATGCCACTTTACAATTACCTCAAAGGAGCAAACAGCGAGTCTGTGACTATGAGCCCTACAGCTCCGGACATGATCCTGGTCCACAAATCAACAGACAACCCACCCGAGGGGTGTGACTTTGCTTACAGTCTGTTGTTTTACCTTCCCGTTGACAGACAGCTCAGTGCCCCAGCACCGACTGAGGATGGTGTTGAACTCATCAAGCAAGAAGCCAAGGAAGCGTTTGTCATTACCTTCACAGATTCAACAACGGACCGAGTACTTGCAGCCAAGGAGACCGAGCTACACAATGAACTGGACAAACGTGAACTCTGTTATAGCACCGATGAATTCTTCTTTGCAAGCTATGATGCATCGTGGAGACCAGAGCCCCATCGCAAAGAGATTTGGATCCCAGAAGCAATGTGCACTACACCAAGTGTCCAAGGTCCTGGTTACACCATCCTGGACAATGGATGTGGGGAAGAGGTCAACTGTCCATCTTATCAGAAGGGTCAATCCCACGGGACATTCGAGGAGAGGAACTACGAGGAGAGCATGTGGATTTGTAAGACGACCGTGAGCTGTAATGTTGAGCAGGCCTTCAGTCAGTCCATCCTGCCGCTGTACCAGTATTTCACCAGCAGCTATGACATCCAACCTGTAGCCAGACCTATCATTTCCTATATGAGAATTGCCGACTTACTGAGCACGGAGTGCAACAAGGAGATTAAGACCTGTGCCTACCTACCAGAAACTGCTTTCAGACTGCCCCTGTTGACCCCACCTAGTCAATTGACATTGCTGAGCGTTCCCAGTGATGACCATCCATGGCAGGCAGCTTATGTAACCGCCCTCCAAGGTCCACCGACAGTCACCAACATCCAGTCAAGTATTACAACACTGCTGCAACAGATCATGGATGCAGACTTGTTCCACAGTGGCCGTGTGTTTGTGGCTCGTTACACCCTTCCGATGGAGGATGGACAGCAATACATTGAAGTCGGAGCTCTGACTGACACTATATCTTATGACGTCTAA
#Can promote mitochondrial permeability transition and facilitate necrotic cell death under different types of stress conditions
#Heme binding protein 2
#human protein

#"Induction of necrotic cell death and mitochondrial permeabilization by heme binding protein 2/SOUL."
#Szigeti A., Bellyei S., Gasz B., Boronkai A., Hocsak E., Minik O., Bognar Z., Varbiro G., Sumegi B., Gallyas F. Jr.
#FEBS Lett. 580:6447-6454(2006) [PubMed] [Europe PMC] [Abstract]

#>TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581 type:5prime_partial len:127 TRINITY_DN47102_c1_g1_i3:420-40(-)
#GATCGATGTGCACCATCTTCGTACCGGAAATGGCCTGCAGGCCGGGACACTCAGAATGCCACTCAAAAGCATCTTCAGCGGGATTTTTGTACACAGAGCTGCAAAAAAGAAAAAAGAAGAAGAAACGAAATAAAATACAAGGAGGGTCAAATCACACCCATAACAACATCGGATTTGCACAAGTTTTCCACCACTGGAAAGCCACGATTCAGGGCTCTCACCTCGTTGTGGATCCGTGTCTCTGACCATTTTCTAGGCGGGGCGAGCGGCCCCGGAATAGGGCTACAACCACCGGCAAATTTGCCCACGCGCCGGGGCCGAGCAACAACAAAGGCAAGACATGTCCATCCGAAGAATGATTTGAAACAACACCCGGGCTAA
#No significant similarity found

#>TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710 type:internal len:168 TRINITY_DN19042_c0_g1_i1:2-502(+)
#TACACTATGGTTGATGCGGTAGCTAAACACTCCTTCAAAGCAACTAGCGAAGATGAGCTGTCCTTCGAGAAAGGAGATACAATCGTTGTTAACGTTTTAGAGAACTCGGAGCAGCACTGGTACAAAGCCTCTTTGAACGGCAAAAAAGGATTAGTCCCAGCAAACTACATTGAAATAAAACCTTACCCTTTCTTACATGGAAACATCTCAAGAGACGGCGCTGTTTGGAAGCTTAACGACCAGCCAGACGGAGCTTTCCTTATCCGTAGAAGTGAGAGTGACATCAGCGAAATCAGTTACTCACTGTCTGTCAAGTACGGTGATGGAGTACAGCACTTCAAGATCTTAACAGATGGCAGCTACAAATACTTCCTGTGGGTAGTGAAGTTTAGATCGCTAAATGAACTGGTGGAATACTACAAGACGTCCTCTGTTAGCCGGACACAAACAATTTTATTGAAGGATATGGAAGTGGAGCACGAAGTAGCGGTGGCAGAATTT
#Adapter protein that provides a critical link between cell surface growth factor receptors and the Ras signaling pathway.
#Isoform 2 does not bind to phosphorylated epidermal growth factor receptor (EGFR) but inhibits EGF-induced transactivation of a RAS-responsive element. Isoform 2 acts 
#as a dominant negative protein over GRB2 and by suppressing proliferative signals, may trigger active programmed cell death.
#Associate with host-virus interaction
#human protein


#>TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316 type:5prime_partial len:107 TRINITY_DN39048_c5_g1_i1:361-41(-)
#TTTCGAGGGCCGGCGGGAGCTCGCCGGACGCCGCCGGAGACGCGGCGCTCTACGGGGCCGCAGCCCCTAGCTCCGGACAAGCCGATTCCAGGGACCCCGCCCCTTACCAAGAAAAGAAAACTCTTCCCGGGGCCCCGGCCGGCGTCTGCGGGCTCGCTCGCTTTACAGCTCGGACCGCCACGAAGACGGAGGGCTCCGCCTCCGGGTTCGGGAATATTGACCCGATTCCCTTTCGCCCGGGGAGCGAGGCCCCCGCTCACATGGAGAGGGGCGCCCGCCGCCGCTCGGAGCGGAACTTCCCTGGGGCTTAGGATCTACTGA
#No Blastx results
#Blast Asterias forbesii sea star 5.8S rRNA gene

#>TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447 type:internal len:150 TRINITY_DN39233_c0_g1_i1:2-448(+)
#GGACGCGACGGATCAACTGGTTCTCCTGGTAATAAGGGATTCTCTGGTGATGTTGGACGAATTGGAATCCCAGGACTTCCTGGAGCCCCAGGTGACCGTGGTGGTGTGGGACAATCAGGTATTTCTGGAGCACCTGGTGCCCAAGGTTCATCAGGCAAGCGAGGAATTTCTGGAAAAGCCGGATCTAGAGGATCTCCAGGTCCACAAGGACCAGCTGGGTCTACAGGTCTGAGCGGTGAGCCGGGTGATATTGGAAGATCAGGAAGGACTGGTCCCCCTGGACCTAAAGGAAACCCGGGAGAGAGAGGAGGCTCTGGTCAACGAGGAGCTCAGGGTATTCCAGGAGGACAGGGACTAGCAGGAGGTGCCGGAGCCCGTGGGGAGAAGGGACCAAAAGGCTCCATCGGTGAGGCTGGAGTCCAGGGTGGCTCTGGCCCAACTGGACGA
#Human
#Collagen alpha-1(XI) chain
#May play an important role in fibrillogenesis by controlling lateral growth of collagen II fibrils.



#>TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543 type:internal len:103 TRINITY_DN39328_c5_g1_i1:3-308(+)
#TCCTCGAGGACCGAGTGGCGGAGAAGGGTTCCATGTGAACAGCAGTTGTACATGGGTCAGTCGATCCTAAGCCCCAGGGAAGTTCCGCTCCGAGCGGAGGCGGGCGCCCCTCTCCATGTGAGCGGGGGCCTCGCTCCCCGGGCGAAAGGGAATCGGGTCAATATTCCCGAACCCGGAGGCGGAGCCCTCCGTCTTCGTGGCGGTCCGAGCTGTAAAGCGAGCGAGCCCGCAGACGCCGGCCGGGGCCCCGGGAAGAGTTTTCTTCTCTTGGTAAGGGGCGGGATCCCTGGAATCGGCGTGTCCGGA
#BlastX results: uncharacterized protein
#blast results 5.8s rRNA


##################Plotcounts for DGE####################


par(mfrow=c(2,3))

plotCounts(ddsLoc, gene="TRINITY_DN41041_c3_g2_TRINITY_DN41041_c3_g2_i1_g.10613_m.10613", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN47102_c1_g1_TRINITY_DN47102_c1_g1_i3_g.24581_m.24581", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN19042_c0_g1_TRINITY_DN19042_c0_g1_i1_g.1710_m.1710", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39048_c5_g1_TRINITY_DN39048_c5_g1_i1_g.8316_m.8316", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39233_c0_g1_TRINITY_DN39233_c0_g1_i1_g.8447_m.8447", intgroup="location")
plotCounts(ddsLoc, gene="TRINITY_DN39328_c5_g1_TRINITY_DN39328_c5_g1_i1_g.8543_m.8543", intgroup="location")




