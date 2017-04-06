##Lauren Ashlock
##Ecological Genomics
##3 27 17
##Tutorial number 4

install.packages("vcfR")

install.packages("adegenet")

library(vcfR)
library(adegenet)


#Read the vcf SNP data into R
vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names
gl1$loc.names[1:10]
gl1$chromosome[1:3]

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
str(ssw_meta)
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort by Individual ID, just like the VCF file

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names
ssw_meta$Individual

gl1$pop <- ssw_meta$Location # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory) # assign disease status

# WE can explore the structure of our SNP data using the glPlot function, which gives us a sample x SNP view of the VCF file
glPlot(gl1, posi="bottomleft")

#0 count reference homozygote
#1 count heterozygote
#2 count alternative homozygote

# Now, let's compute the PCA on the SNP genotypes and plot it:
pca1 <- glPca(gl1, nf=4, parallel = FALSE) # nf = number of PC axes to retain (here, 4)
pca1 # prints summary

# Plot the individuals in SNP-PCA space, with locality labels:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=gl1$pop, 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(gl1$pop), 
       pch=20, 
       col=c("black", "red"))
# Not any clear structure here

# Perhaps we want to show disease status instead of locality:
plot(pca1$scores[,1], pca1$scores[,2], 
     cex=2, pch=20, col=as.factor(unlist(gl1$other)), 
     xlab="Principal Component 1", 
     ylab="Principal Component 2", 
     main="PCA on SSW data (Freq missing=20%; 5317 SNPs)")
legend("topleft", 
       legend=unique(as.factor(unlist(gl1$other))), 
       pch=20, 
       col=as.factor(unique(gl1$other$Trajectory)))

# Which SNPs load most strongly on the 1st PC axis?
loadingplot(abs(pca1$loadings[,1]),
            threshold=quantile(abs(pca1$loadings), 0.999))

# Get their locus names
threshold = quantile(abs(pca1$loadings), 0.999)
gl1$loc.names[which(abs(pca1$loadings)>threshold)]


#May need to debug and update some of the code above (PCA and DAPC)







