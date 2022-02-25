##############
# PRACTICE 1 #
##############

# import data
all_data <- read.csv("table.csv", header = TRUE, row.names = 1)

# Selection of columns which contain the transcriptome data of the wild type,
# the mutant and the complemented lines (columns 1 to 8)
data <- all_data[,1:8]

###############
# Exercise 1 #
##############

# PCA Calculation
PCA1 <- princomp(data, cor = FALSE, scores = TRUE)

# Calculating and plotting the variance for each component
variance1<-((PCA1$sdev^2/sum(PCA1$sdev^2))*100)
plot(variance1, type="s", ylab="% of variance", xlab="PCA component")

# Variance of the first component
variance1_component1 <- ((PCA1$sdev[1]^2/sum(PCA1$sdev^2))*100)
variance1_component1

# Plotting loadings 
plot(PCA1$loadings)
text(PCA1$loadings, labels = row.names(PCA1$loadings), pos = c(2,2,1,4,1,4,1,4))

###############
# Exercise 2 #
##############

# Creation of intermediate transcriptomes 
int25 <- (0.25*data$J0571) + (0.75*data$shrJ0571)
int50 <- (0.50*data$J0571) + (0.50*data$shrJ0571)
int75 <- (0.75*data$J0571) + (0.25*data$shrJ0571)

# Creation of a new data frame containing the intermediate transcriptomes
# and the other transcriptome profiles
genes_names <- row.names(data)
int <- data.frame(genes_names, int25,int50,int75, data, row.names = 1)

# Calculating PCA for intermediate transcriptomes, wt and mutant
PCA2 <- princomp(int, cor = FALSE, scores = TRUE)

# Calculating and plotting PCA variance
int_var <-((PCA2$sdev^2/sum(PCA2$sdev^2))*100)
plot(int_var, type="s", ylab="% of variance", xlab="PCA component")

# Variance explained by the first component
int_var[1]

# Plotting PCA loadings
plot(PCA2$loadings)
text(PCA2$loadings, labels = c("25%", "50%", "75%", "shr", "wt", "BLJ", "JKD", "MGP", "NUC", "IME", "SCR"), pos = 4)


##############
# Exercise 3 #
##############

# Select transcriptome profile of SCR and including it in a new dataframe
# with the intermediate and all samples transcriptomes
SCR_domain <- all_data$SCRdomain
int_scr <- data.frame(int, SCR_domain)

# Calculating PCA
PCA3 <- princomp(int_scr, cor = FALSE, scores = TRUE)

# Plotting loadings
plot(PCA3$loadings)
text(PCA3$loadings, labels = c("25%", "50%", "75%", "shr", "wt", "BLJ", "JKD", "MGP", "NUC", "IME", "SCR", "SCR domain"), pos = 4)

var3 <-((PCA3$sdev^2/sum(PCA3$sdev^2))*100)
plot(var3, type="s", ylab="% of variance", xlab="PCA component")
var3[1]

##############
# Exercise 4 #
##############

# Select the first component and sorting according it 
comp1_scores <- PCA1$scores[,1]
comp1_sorted <- sort(comp1_scores)

# Select the highest and the lowest score
highest <- tail(comp1_sorted, n = 10L)
lowest <- head(comp1_sorted, n = 10L)

# Select the genes and it is name
significant_genes <- c(lowest, highest)
significant_genes_names <- names(significant_genes)

# Define the matrix containing the transcriptome profile for the most 
# significant genes
A <- data.matrix(data[significant_genes_names,])

# Plotting heatmap
heatmap(A)
