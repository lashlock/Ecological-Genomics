##Looking at allele frequencies 
##Ecological Genomics
##Lauren Ashlock
##3 20 17

#### read in your text files (remember to edit header first)

H_freq <- read.table("H_AlleleFreqs.frq", header=T)
S_freq <- read.table("S_AlleleFreqs.frq", header=T)

str(H_freq)
str(S_freq)

# Since these files have identical numbers of SNPs in the exact same order, 
# we can concatenate them together into one large dataframe:

All_freq <- merge(H_freq, S_freq, by=c("CHROM", "POS"))

# Looks good, now let's calculate the difference in minor allele frequency
# at each SNP and plot as a histogram

All_freq$diff <- (All_freq$H_ALT - All_freq$S_ALT)

hist(All_freq$diff, breaks=50, col="red", main="Allele frequency difference (H-S)")


# Looks like most loci show little difference (i.e., likely drift), 
# but perhaps a few show very large differences between healthy and sick 
# (drift or selection?)

# How do these highly divergent frequenices compare to Fst at the same SNPs?
fst <- read.table("HvS_Fst.weir.fst", header=T)

All_freq.fst <- merge(All_freq, fst, by=c("CHROM", "POS"))

plot(All_freq.fst$diff, All_freq.fst$WEIR_AND_COCKERHAM_FST, xlab="Allele frequency difference (H-S)", ylab="Fst", main="Healthy vs. Sick SNP divergence")

#As I guess we would expect.... with greater allele frequency difference
#there appears to be a greater Fst
#Each point on this plot is a SNP

# Which are the genes that are showing the highest divergence between Healthy and Sick?
All_freq.fst[which(All_freq.fst$WEIR_AND_COCKERHAM_FST>0.2),]

points(0.2500000,0.2415350,col="red", fill="red", cex=2)

#This modifies a point on the plot you currently have open




