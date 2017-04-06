#OutFLANK
#ecological genomics
#4 3 17

#Install packages
install.packages("devtools")
library(devtools)
source("http://bioconductor.org/biocLite.R")
biocLite("qvalue")
install_github("whitlock/OutFLANK")
install.packages("tcltk")
install.packages("stringi")
library("stringi")
a

#load these packages
library(OutFLANK)
library(vcfR)
library(adegenet)
library(tcltk)
library(tibble)
library(qvalue)


#Read in your geno file

ssw.geno_in <- read.fwf("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf.geno", 
                        width=rep(1,24))
ssw.geno <- t(ssw.geno_in)


#read in the metadata

ssw_meta <- read.table("ssw_healthloc.txt", header=TRUE)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),]
ssw_meta$Trajectory[which(ssw_meta$Trajectory=='MM')] = NA

OF_SNPs <- MakeDiploidFSTMat(ssw.geno,locusNames=seq(1,5317,1) , 
                             popNames=ssw_meta$Trajectory)

dim(OF_SNPs)

head(OF_SNPs)

OF_out <- OutFLANK(FstDataFrame = OF_SNPs, LeftTrimFraction = 0.05, 
                   RightTrimFraction = 0.05, Hmin = 0.1, NumberOfSamples = 3, 
                   qthreshold = 0.1)

OutFLANKResultsPlotter(OF_out, withOutliers = T, NoCorr = T, Hmin = 0.1, 
                       binwidth = 0.005, titletext = "Scan for local selection")
#find your outliers

outliers <- which(OF_out$results$OutlierFlag=="TRUE")
outliers

#we can extract info about the outliers by reading in the vcf file and looking at the annotations
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")
vcfann <- as.data.frame(getFIX(vcf1))
vcfann[outliers,]
