##Ecological Genomics
##Lauren Ashlock
## 3 22 17
## PopGen Tutorial 4

list.files() # Do you see your downloaded files there? If not, double check to make sure you've set your working directory to the right spot

# We'll need to install 2 packages to work with the SNP data:
install.packages("vcfR") # reads in vcf files and proides tools for file conversion 
install.packages("adegenet") # pop-genetics package with some handy routines, including PCA and other multivariate methods (DAPC)

# ...and load the libraries
library(adegenet)
library(vcfR)

#Read the vcf SNP data into R
download.file("https://raw.githubusercontent.com/stephenrkeller/PBIO381_srkeller_labnotebook/master/data/SNP_data/SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf",dest="test.vcf")

vcf1<-read.vcfR("test.vcf")
#vcf1 <- read.vcfR("SSW_all_biallelic.MAF0.02.Miss0.8.recode.vcf")

# The adegenet package uses a highly efficient way of storing large SNP datasets in R called a "genlight" object. The following function creates a genlight object from your vcf:
gl1 <- vcfR2genlight(vcf1)
print(gl1) # Looks good! Right # of SNPs and individuals!

# For info, try:
gl1$ind.names
gl1$loc.names[1:10]

# Notice there's nothing in the field that says "pop"? Let's fix that...
ssw_meta <- read.table("ssw_healthloc.txt", header=T) # read in the metadata
ssw_meta <- ssw_meta[order(ssw_meta$Individual),] # sort it by Individual ID

# Confirm the ID's are ordered the same in gl1 and ssw_meta:
gl1$ind.names
ssw_meta$Individual

gl1$pop <- ssw_meta$Locality # assign locality info
gl1$other <- as.list(ssw_meta$Trajectory) # assign disease status

# Now, let's compute the PCA on the SNP genotypes and plot it:
PCA1 <- glPca(gl1, nf=4) 
PCA1 # prints summary
scatter(PCA1, label=gl1$pop) # plots PCA scores and labels by locale
loadingplot(pca1) # which SNPs load most strongly on the 1st PC axis?