##Ecological Genomics Assignment 3
##Lauren Ashlock

install.packages("vcfR")
install.packages("adegenet")

library(adegenet)
library(vcfR)

#read in vcf files

vcf1.0 <- read.vcfR("filteredSNPS1.0.recode.vcf")

vcf2.0 <- read.vcfR("filteredSNPS2.0.recode.vcf")

gl1.0 <- vcfR2genlight(vcf1.0)

gl2.0 <- vcfR2genlight(vcf2.0)

ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort it by Individual ID

#data check
gl1.0$ind.names
ssw_meta$Individual
ssw_meta$Trajectory
str(gl1.0$other)


gl1.0$pop <- ssw_meta$Location # assign locality info
gl1.0$other <- as.list(ssw_meta$Trajectory) # assign disease status

gl2.0$pop <- ssw_meta$Location # assign locality info
gl2.0$other <- as.list(ssw_meta$Trajectory) # assign disease status

# Now, let's compute the PCA on the SNP genotypes and plot it:
pca1.0 <- glPca(gl1.0, nf=4, parallel=F) # nf = number of PC axes to retain (here, 4)

pca2.0 <- glPca(gl2.0, nf=4, parallel=F) # nf = number of PC axes to retain (here, 4)


# Perhaps we want to show disease status instead of locality:
plot(pca1.0$scores[,1], pca1.0$scores[,2], 
     cex=2, pch=20, col=as.factor(ssw_meta$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on Dataset 1.0 (Freq missing=15%; 1945 SNPs)")
legend("topleft", 
       legend=unique(ssw_meta$Trajectory), 
       pch=20, 
       col=as.factor(unique(ssw_meta$Trajectory)))

plot(pca2.0$scores[,1], pca2.0$scores[,2], 
     cex=2, pch=20, col=as.factor(ssw_meta$Trajectory), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on Dataset 2.0 (Freq missing=15%; 1861 SNPs)")
legend("topleft", 
       legend=unique(ssw_meta$Trajectory), 
       pch=20, 
       col=as.factor(unique(ssw_meta$Trajectory)))



# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1.0$loadings[,1]),
            threshold=quantile(abs(pca1.0$loadings), 0.999), main="Loading Plot PCA 1.0")

# Get their locus names
gl1.0$loc.names[which(abs(pca1.0$loadings)>quantile(abs(pca1.0$loadings), 0.999))]

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca2.0$loadings[,1]),
            threshold=quantile(abs(pca2.0$loadings), 0.999), main="Loading Plot PCA 2.0")

# Get their locus names
gl2.0$loc.names[which(abs(pca2.0$loadings)>quantile(abs(pca2.0$loadings), 0.999))]

