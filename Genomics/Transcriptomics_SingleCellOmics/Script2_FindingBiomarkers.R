##############
# PRACTICE 2 #
##############

# Loading data
data <- read.csv("TableS5.csv", row.names = 1, header = TRUE)

##############
# Exercise 1 #
##############
boxplot(data, las = 2)

##############
# Exercise 2 #
###############
# Calculating PCA
PCA <- princomp(data, cor = FALSE, scores = TRUE)

# Calculating and plotting variance
variance<-((PCA$sdev^2/sum(PCA$sdev^2))*100)
plot(variance, type="s", ylab="% of variance", xlab="PCA component")

# Displaying the variance for the first 4 components
variance[1]
variance[2]
variance[3]
variance[4]


comp1 <- PCA$loadings[,1]
comp2 <- PCA$loadings[,2]
comp3 <- PCA$loadings[,3]

sample_names <- colnames(data)

# Comp1 vs Comp2
plot(comp1, comp2, xlab = "Comp. 1", ylab = "Comp. 2")
text(comp1, comp2, labels = sample_names, cex = 0.75, pos = c(1,3,3,4,1,3,1,2,3,3,4,4,1,3,3))

# Comp1 vs Comp3
plot(comp1, comp3, xlab = "Comp. 1", ylab = "Comp. 3")
text(comp1, comp3, labels = sample_names, cex = 0.75, pos = c(3,4,3,4,3,3,1,2,3,4,4,4,1,3,3))


##############
# Exercise 3 #
##############

markers <- subset(data, WOX5>1 & WOL<1 & S4<1 & SCR<1 & WER<1 & PET111<1 & S17<1 & S18<1 & S32<1 & COR<1 & E30<1 & GL2<1 & APL<1 & CO2<1 & COBL9<1)

markers <- data.matrix(markers)
heatmap(markers)

##############
# Exercise 4 #
##############

markersPCA <- princomp(markers, cor = FALSE, scores = TRUE)

# Variance calculation and plotting
markersVar <- ((markersPCA$sdev^2/sum(markersPCA$sdev^2))*100)
plot(markersVar, type="s", ylab="% of variance", xlab="PCA component")
markersVar[1]
markersVar[2]


# Plotting loadings
plot(markersPCA$loadings)
text(markersPCA$loadings, labels = row.names(markersPCA$loadings), pos = c(2,2,2,4,3,2,2,2,1,4,2,2,4,2,4))

# Plotting scores
plot(markersPCA$scores)
text(markersPCA$scores, labels = row.names(markersPCA$scores), pos = 4)

##############
# Exercise 5 #
##############
comp1 <- PCA$scores[,1]
comp2 <- PCA$scores[,2]
comp3 <- PCA$scores[,3]
comp4 <- PCA$scores[,4]

par(mfrow=c(1,2))
plot(comp1, comp2, xlab = "Comp. 1", ylab = "Comp. 2")
text(comp1, comp2, labels = row.names(PCA$scores), pos = 3)
plot(comp3, comp4, xlab = "Comp. 3", ylab = "Comp. 4")
text(comp3, comp4, labels = row.names(PCA$scores), pos = 3)




