# Read in the Romiguier data:
Rom <- read.csv("Romiguier_nature13685-s3.csv", header=T)

# Import OK?
str(Rom) 
head(Rom)

# Looks good!
# Now let's look at how the strength of purifying selection (piN/piS) 
# compares to the size of Ne (piS). We'll plot these on a log scale to 
# linearize the relationship.
plot(log(Rom$piS), log(Rom$piNpiS), pch=21, bg="blue", xlab="log Synonymous Nucleotide Diversity (piS)", ylab="log Ratio of Nonysn to Syn Diversity (piN/piS)", main="Purifying Selection vs. Effective Population Size")

# Now let's add our SSW points to the existing plot and give them a different symbol
points(log(0.00585312), log(0.264041), pch=24, cex=1.5, bg="red") 

# We can also add a regression line to the plot to see how far off the SSW estimates are from expectation
reg <- lm(log(Rom$piNpiS) ~ log(Rom$piS)) # Fits a linear regression
abline(reg) # adds the regression line to the plot

# It would be useful to highlight the other echinoderms in the dataset...do our seastars behave similarly?
echino <- Rom[which(Rom$Phylum=="Echinodermata"),] # subsets the data
points(log(echino$piS), log(echino$piNpiS), pch=21, bg="red") # adds the points

# Lastly, let's add a legend:
legend("bottomleft", cex=1, legend=c("Metazoans", "Echinoderms", "P. ochraceus"), pch=c(21,21,24), col=c("blue", "red", "red"))

# Pisaster seems to be in a group with other echinoderms that have relaxed 
#purifying selection (high piN/piS), given their Ne...Interesting!
#Can we hypothesize why this might be?