source("~yangri13/SHC/R/AllGenerics.R")
source("~yangri13/SHC/R/null_eigval.R")
source("~yangri13/SHC/R/package-sigclust2.R")
source("~yangri13/SHC/R/shc-constructor.R")
source("~yangri13/SHC/R/shc-fwer_control.R")
source("~yangri13/SHC/R/shc-plot.R")
source("~yangri13/SHC/R/shc-print.R")
source("~yangri13/SHC/R/shc-summary.R")
source("~yangri13/SHC/R/shcutree.R")
source("~yangri13/SHC/R/sigclust-constructor.R")

library(WGCNA)
library(stats)
library(ggplot2)
library(ggdendro)
library(dplyr)
library(dynamicTreeCut)
library(RANN)

# library(sigclust2)
library(Seurat)
library(Rcpp)
library(HGC)
set.seed(1)
seurat_object <- readRDS("~yangri13/projects/def-wainberg/yangri13/pca_green_50k.rds") #~yangri13/projects/def-wainberg/yangri13/Green_500k.rds
pca_comp_50 <- Embeddings(seurat_object, reduction = "pca", dims = 1:50)
# # seurat_object@reductions$PCs@cell.embeddings
# 
time_taken <- system.time({
  shc_result <- shc(pca_comp_50, metric="euclidean", linkage="ward", rcpp=TRUE, alpha=0.05)
})
print(time_taken)
# shc_result <- readRDS("SHC_norm_mini.rds")

# png("shc_p500_mini.png", 800, 600)
# plot(shc_result, hang=.1)
# dev.off()

hc_dat <- shc_result$hc_dat
p_norm <- shc_result$p_norm

num_clusters <- length(shc_result$clusters)
clusters <- shc_result$clusters
# Set cluster names and add them to the Seurat object
names(clusters) <- Cells(seurat_object)
clusters <- as.factor(clusters)
seurat_object$shc_clusters <- clusters

# Run UMAP on the Seurat object
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50)

# Generate a discrete palette with the correct number of colors
num_colors <- num_clusters

# Function to generate random hexadecimal color
generate_random_color <- function() {
  sprintf("#%02X%02X%02X", sample(0:255, 1), sample(0:255, 1), sample(0:255, 1))
}

# Generate 62 random colors
random_colors <- replicate(num_colors, generate_random_color())
# random_colors <- DiscretePalette(n=5)
# Plot UMAP with clusters
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "shc_clusters", cols = random_colors) +
  labs(title = "UMAP Plot with SHC Clusters") +
  theme_minimal()

# Add cluster labels to the UMAP plot
umap_plot <- LabelClusters(plot = umap_plot, id = "shc_clusters", repel = TRUE)

# Save the plot with labels
png("UMAP_green_50k_inv.png", 1200, 800)
print(umap_plot)
dev.off()

# Load necessary libraries
library(ComplexHeatmap)
library(circlize)
library(clue)
library(ComplexHeatmap)
library(circlize)
library(clue)

true_labels <- seurat_object@meta.data[["state"]] # Subclass, state
predicted_labels <- seurat_object$shc_clusters

# Create confusion matrix
confusion_matrix <- table(True = true_labels, Predicted = predicted_labels)

# Calculate row sums
row_sums <- rowSums(confusion_matrix)

# Avoid division by zero by replacing 0s in row_sums with 1
row_sums[row_sums == 0] <- 1

# Calculate percentages
confusion_matrix_percent <- confusion_matrix / row_sums * 100

# Transpose the confusion matrix to swap the axes
confusion_matrix_percent_transposed <- t(confusion_matrix_percent)

# Check for and replace NA or NaN values with 0
confusion_matrix_percent_transposed[is.na(confusion_matrix_percent_transposed)] <- 0

# Make the matrix square by adding dummy columns
n_rows <- nrow(confusion_matrix_percent_transposed)
n_cols <- ncol(confusion_matrix_percent_transposed)

if (n_rows > n_cols) {
  # Add dummy columns with zeros
  confusion_matrix_percent_transposed <- cbind(
    confusion_matrix_percent_transposed, 
    matrix(0, nrow = n_rows, ncol = n_rows - n_cols)
  )
} else if (n_cols > n_rows) {
  # Add dummy rows with zeros (just in case, though unlikely)
  confusion_matrix_percent_transposed <- rbind(
    confusion_matrix_percent_transposed, 
    matrix(0, nrow = n_cols - n_rows, ncol = n_cols)
  )
}

# Use the Hungarian algorithm to find the optimal assignment
assignment <- solve_LSAP(confusion_matrix_percent_transposed, maximum = TRUE)
row_order <- order(assignment)

# Reorder the matrix (only include original columns)
confusion_matrix_percent_transposed_reordered <- confusion_matrix_percent_transposed[row_order, 1:n_cols]

# Define the color palette
col_fun <- colorRamp2(c(0, 50, 100), c("white", "lightblue", "blue"))

# Create the heatmap
heatmap <- Heatmap(
  confusion_matrix_percent_transposed_reordered,
  name = "Percentage",
  col = col_fun,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_dend = FALSE,
  show_column_dend = FALSE,
  row_names_side = "left",
  column_names_side = "top",
  heatmap_legend_param = list(title = "Percentage", at = c(0, 50, 100)),
  border = TRUE,  # Enable borders
  rect_gp = gpar(col = "grey", lwd = 1)  
)

# Save the heatmap to a PNG file
png("cm_green_50k_inv_or.png", width = 1800, height = 800)
draw(heatmap)
dev.off()



# # TRANSPOSED CONFUSION MATRIX
# library(gplots)
# library(clue)  # For the Hungarian algorithm
# 
# true_labels <- seurat_object@meta.data[["state"]] # Subclass
# predicted_labels <- seurat_object$shc_clusters
# 
# # Create confusion matrix
# # Create confusion matrix
# confusion_matrix <- table(True = true_labels, Predicted = predicted_labels)
# 
# # Calculate row sums
# row_sums <- rowSums(confusion_matrix)
# 
# # Avoid division by zero by replacing 0s in row_sums with 1
# row_sums[row_sums == 0] <- 1
# 
# # Calculate percentages
# confusion_matrix_percent <- confusion_matrix / row_sums * 100
# 
# # Transpose the confusion matrix to swap the axes
# confusion_matrix_percent_transposed <- t(confusion_matrix_percent)
# 
# # Check for and replace NA or NaN values with 0
# confusion_matrix_percent_transposed[is.na(confusion_matrix_percent_transposed)] <- 0
# 
# # Use the Hungarian algorithm to find the optimal assignment
# assignment <- solve_LSAP(confusion_matrix_percent_transposed, maximum = TRUE)
# row_order <- order(assignment)
# 
# # Reorder the matrix
# confusion_matrix_percent_transposed_reordered <- confusion_matrix_percent_transposed[row_order, ]
# 
# # Plot the normalized and transposed confusion matrix with grid lines
# png("cm_p500_mini.png", width = 1800, height = 800)
# heatmap.2(confusion_matrix_percent_transposed_reordered, 
#           col = colorRampPalette(c("white", "blue"))(100),
#           Rowv = FALSE,
#           Colv = FALSE,
#           scale = "none",
#           main = "Confusion Matrix: Green Mini Sample",
#           key.title = "Percentage",
#           trace = "none",
#           dendrogram = "none",
#           sepwidth = c(0.01, 0.01),
#           sepcolor = "grey",
#           colsep = 1:ncol(confusion_matrix_percent_transposed_reordered),
#           rowsep = 1:nrow(confusion_matrix_percent_transposed_reordered))
# dev.off()


# V2, with SingleCell
library(gplots)

dendrogram <- shc_result$hc_dat
hc <- as.hclust(dendrogram)
clusters <- cutreeDynamic(dendro = hc, distM = as.matrix(dist(shc_result$in_mat)), verbose = 0)
num_clusters <- length(unique(clusters))

names(clusters) <- Cells(seurat_object)
clusters <- as.factor(clusters)

true_labels <- seurat_object@meta.data$Subclass
predicted_labels <- clusters

# Create confusion matrix
confusion_matrix <- table(True = true_labels, Predicted = predicted_labels)

# Normalize the confusion matrix by columns
col_sums <- colSums(confusion_matrix)
confusion_matrix_percent <- t(t(100 * confusion_matrix) / col_sums)

# Plot the normalized confusion matrix
png("cm_shc_sc_mini.png", width = 800, height = 600)
heatmap.2(confusion_matrix_percent, 
          col = colorRampPalette(c("white", "blue"))(100),
          scale = "none",
          main = "Normalized Confusion Matrix for SHC Clusters",
          key.title = "Percentage",
          trace = "none",
          dendrogram = "none")
dev.off()

# 
# cluster_10_seurat_object <- subset(seurat_object, subset = scSHC == 10)
# avg_nFeature_RNA_cluster_10 <- mean(cluster_10_seurat_object$nFeature_RNA)
# avg_nCount_RNA_cluster_10 <- mean(cluster_10_seurat_object$nCount_RNA)

# Calculate average nFeature_RNA and nCount_RNA for all clusters
avg_nFeature_RNA_all <- mean(seurat_object$nFeature_RNA)
avg_nCount_RNA_all <- mean(seurat_object$nCount_RNA)
# ---------------------------------------------------------------
cluster_ids <- unique(seurat_object$shc_clusters)
avg_nFeature_RNA <- numeric(length(cluster_ids))
avg_nCount_RNA <- numeric(length(cluster_ids))

# Loop through each cluster and calculate the averages
for (i in seq_along(cluster_ids)) {
  cluster_id <- cluster_ids[i]
  cluster_data <- subset(seurat_object, subset = shc_clusters == cluster_id)

  avg_nFeature_RNA[i] <- mean(cluster_data$nFeature_RNA)
  avg_nCount_RNA[i] <- mean(cluster_data$nCount_RNA)
}

# Combine results into a data frame
results <- data.frame(
  cluster_id = cluster_ids,
  avg_nFeature_RNA = avg_nFeature_RNA,
  avg_nCount_RNA = avg_nCount_RNA
)

# Print the results
print(results)
seurat_object <- readRDS("pca_seurat_mini.rds")
set.seed(1)

cluster_colors <- ifelse(clusters == 16, "red", "blue")
x <- Embeddings(seurat_object, reduction = "pca")[, 1]

# y <- seurat_object$nCount_RNA
y <- seurat_object$nFeature_RNA
correlation <- cor(x, y)
png("logCP10k_nfeature_pc1.png", 800, 600)
plot(x, y, pch = 19, col = cluster_colors, xlab = "PC1", ylab = "nFeature_RNA", main = "nFeature vs PC1")
text(x = min(x), y = max(y), labels = paste("Correlation: ", round(correlation, 2)), pos = 4)
dev.off()

y <- seurat_object$nCount_RNA
correlation <- cor(x, y)
png("logCP10k_ncount_pc1.png", 800, 600)
plot(x, y, pch = 19, col = "lightblue", xlab="PC1", ylab="nCount_RNA", main="nCount vs PC1")
text(x = min(x), y = max(y), labels = paste("Correlation: ", round(correlation, 2)), pos = 4)
dev.off()

x <- Embeddings(seurat_object, reduction = "pca")[, 2]
y <- seurat_object$nFeature_RNA
correlation <- cor(x, y)
png("logCP10k_nfeature_pc2.png", 800, 600)
plot(x, y, pch = 19, col = "lightblue", xlab="PC2", ylab="nFeature_RNA", main="nFeature vs PC2")
text(x = min(x), y = max(y), labels = paste("Correlation: ", round(correlation, 2)), pos = 4)
dev.off()

y <- seurat_object$nCount_RNA
correlation <- cor(x, y)
png("logCP10k_ncount_pc2.png", 800, 600)
plot(x, y, pch = 19, col = "lightblue", xlab="PC2", ylab="nCount_RNA", main="nCount vs PC2")
text(x = min(x), y = max(y), labels = paste("Correlation: ", round(correlation, 2)), pos = 4)
dev.off()

