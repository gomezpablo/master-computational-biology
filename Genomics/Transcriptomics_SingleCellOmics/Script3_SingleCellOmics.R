##############
# Practice 3 #
##############

# Package requirements
library(dplyr)
library(Seurat)
library(Matrix)

# Question 1

# Load data
planaria <- read.table("SCOplanaria.txt", header = FALSE, row.names = 1, sep = "\t")
# Convert it to matrix
planaria <- data.matrix(planaria)
# Convert it to sparse matrix
planaria <- Matrix(data = planaria, sparse = TRUE)

# - A - Create Seurat Object
planaria <- CreateSeuratObject(counts = planaria, project = "planaria", min.cells = 0, min.features = 0)

# - B - Quality control, normalization, variable features selection and scaling of the data.

# Plotting features vs counts to detect outliers
FeatureScatter(planaria, feature1 = "nFeature_RNA", feature2 = "nCount_RNA") 

planaria <- subset(planaria, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

pre_norm <- VlnPlot(object = planaria, features= c("nCount_RNA", "Normalized.Reads"), ncol = 2)
# Normalization
planaria <- NormalizeData(planaria, normalization.method = "LogNormalize", scale.factor = 10000)

# Visualize effects of normalizalition in total cell expression
                            # Normalized values are stored in planaria[["RNA"]]@data
Normalized.Reads <- colSums(planaria[["RNA"]]@data)
planaria <- AddMetaData (planaria, Normalized.Reads, col.name="Normalized.Reads")
after_norm <- VlnPlot(object = planaria, features= c("nCount_RNA", "Normalized.Reads"), ncol = 2)

# Before and after normalization
pre_norm
after_norm

# Feature selection
planaria <- FindVariableFeatures(planaria, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes and plot their label
top10 <- head(VariableFeatures(planaria), 10)
LabelPoints(plot = VariableFeaturePlot(planaria), points = top10, repel = TRUE)

# Scaling data
genes <- rownames(planaria) # storing genes n
planaria <- ScaleData(planaria, features = genes)

# -C- Perform a PCA analysis and obtain a t-SNE clustering plot.
planaria <- RunPCA(planaria, features = VariableFeatures(object = planaria))

# Finding Neighbors and Cluster cells
planaria <- FindNeighbors(planaria, dims = 1:10)
planaria <- FindClusters(planaria, resolution = 0.5)

# Visualize cells through non-linear dimensional reduction (tSNE)
planaria <- RunTSNE(planaria, dims = 1:10)
DimPlot(planaria, reduction = "tsne", label = TRUE)

# Question 2
new_planaria <- read.table("SCOplanaria.txt", header = FALSE, row.names = 1, sep = "\t")
# Convert it to matrix
new_planaria <- data.matrix(new_planaria)
# Convert it to sparse matrix
new_planaria <- Matrix(data = new_planaria, sparse = TRUE)
# Seurat object: new min cells and min features parameter applied
new_planaria <- CreateSeuratObject(counts = new_planaria, project = "new_planaria", min.cells = 3, min.features = 200)

# Subset for further analysis 
new_planaria <- subset(new_planaria, subset = nFeature_RNA > 200 & nFeature_RNA < 2500)

#  Normalization
new_planaria <- NormalizeData(new_planaria, normalization.method = "LogNormalize", scale.factor = 10000)

# Feature selection (3000 genes)
new_planaria <- FindVariableFeatures(new_planaria, selection.method = "vst", nfeatures = 3000)

# Scaling 
new_genes_names <- rownames(new_planaria) #Storing gene names.
new_planaria <- ScaleData(new_planaria, vars.to.regress = "nFeature_RNA")

# PCA analysis.
new_planaria <- RunPCA(new_planaria, features = VariableFeatures(object = new_planaria))

# Clustering.
new_planaria <- FindNeighbors(new_planaria, dims = 1:5)
new_planaria <- FindClusters(new_planaria, resolution = 0.6)

# Running and plotting t-SNE.
new_planaria <- RunTSNE(new_planaria, dims = 1:5)
DimPlot(new_planaria, reduction = "tsne", label = TRUE)

# Question 3
# Finding top5 biomarkers
planaria.markers <- FindAllMarkers(new_planaria, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top5 <- planaria.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)

# Heatmap with top5 biomarkers for each
DoHeatmap(new_planaria, features = top5$gene,hjust = TRUE) 

# Question 4
# Violine plot
VlnPlot(new_planaria, features = c("dd-Smed-v6-61-0", "dd-Smed-v6-2178-0", 
                                   "dd-Smed-v6-298-0", "dd-Smed-v6-1410-0",
                                   "dd-Smed-v6-702-0", "dd-Smed-v6-2548-0", 
                                   "dd-Smed-v6-9977-0", "dd-Smed-v6-48-0", 
                                   "dd-Smed-v6-175-0", "dd-Smed-v6-1161-1"))
# Features plot
FeaturePlot(new_planaria, features = c("dd-Smed-v6-61-0", "dd-Smed-v6-2178-0", 
                               "dd-Smed-v6-298-0", "dd-Smed-v6-1410-0", 
                               "dd-Smed-v6-702-0", "dd-Smed-v6-2548-0", 
                               "dd-Smed-v6-9977-0", "dd-Smed-v6-48-0",
                               "dd-Smed-v6-175-0", "dd-Smed-v6-1161-1"))
# Rename clusters
cluster.names <- c("Neural progenitors - 0", "Neural progenitors - 1", 
                   "Early epidermal progenitors", "Epidermis", 
                   "GABA neurons", "Parenchymal cells", "Muscle body", 
                   "Late epidermal progenitors", "Muscle progenitors", 
                   "Phagocytes", "Pigment")
names(cluster.names) <- levels(new_planaria)
new_planaria <- RenameIdents(new_planaria, cluster.names)

# Question 5

# Question 6
VlnPlot(new_planaria, features = c("dd-Smed-v6-1999-0"))
FeaturePlot(new_planaria, features = c("dd-Smed-v6-1999-0"))
