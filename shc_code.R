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
options(error=recover)
set.seed(1)

seurat_object <- readRDS("~yangri13/projects/def-wainberg/yangri13/pca_green_50k.rds") #~yangri13/projects/def-wainberg/yangri13/Green_500k.rds
pca_comp_50 <- Embeddings(seurat_object, reduction = "pca", dims = 1:50)
# # seurat_object@reductions$PCs@cell.embeddings
# 
source("~yangri13/SHC/R/shc-constructor.R")
options(error=recover)
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

child_clusters <- shc_result$child_clusters
child_borders <- shc_result$child_borders

num_clusters <- length(shc_result$clusters)
clusters <- shc_result$clusters
# Set cluster names and add them to the Seurat object
names(clusters) <- Cells(seurat_object)
clusters <- as.factor(clusters)
seurat_object@meta.data$shc_clusters <- clusters
# seurat_object@meta.data$child_borders <- child_borders
# seurat_object@meta.data$latest_borders <- latest_borders

# Run UMAP on the Seurat object
seurat_object <- RunUMAP(seurat_object, reduction = "pca", dims = 1:50)

# Generate a discrete palette with the correct number of colors
num_colors <- num_clusters

# Function to generate random hexadecimal color
generate_random_color <- function() {
  sprintf("#%02X%02X%02X", sample(0:255, 1), sample(0:255, 1), sample(0:255, 1))
}

random_colors <- replicate(num_colors, generate_random_color())
# random_colors <- DiscretePalette(n=5)
# Plot UMAP with clusters
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "shc_clusters", cols = random_colors) +
  labs(title = "UMAP Plot with SHC Clusters") +
  theme_minimal()

# Add cluster labels to the UMAP plot
umap_plot <- LabelClusters(plot = umap_plot, id = "shc_clusters", repel = TRUE)

# Save the plot with labels
png("UMAP_green_50k_simple.png", 1200, 800)
print(umap_plot)
dev.off()

# CHILD DISPLAY HERE ------------------------------------------------------------------------------------------
# Define the variable to choose between full dataset or subset
seurat_object@meta.data$shc_clusters <- clusters
seurat_object@meta.data$child_borders <- child_borders
seurat_object@meta.data$child_clusters <- child_clusters

display_option <- "subset"  # Change to "full" to display the entire dataset

# Define the child clusters you want to subset (relevant only if display_option is "subset")
clusters_to_display <- c("0a", "0b")  # Replace with your desired cluster IDs

if (display_option == "subset") {
  # Subset the Seurat object based on the selected child clusters
  seurat_plot_data <- subset(seurat_object, subset = child_clusters %in% clusters_to_display)
  # Assuming child_clusters is a factor or character vector in seurat_object@meta.data
  unique_clusters <- unique(seurat_object@meta.data$child_clusters)
  
  # Generate random colors
  random_child_colors <- setNames(replicate(length(unique_clusters), generate_random_color(), simplify = FALSE), unique_clusters)
  
  # Generate the UMAP plot for the subsetted Seurat object
  umap_plot <- DimPlot(
    seurat_plot_data, 
    reduction = "umap", 
    group.by = "child_clusters", 
    cols = random_child_colors  # Reuse or generate new colors if necessary
  ) +
    labs(title = "UMAP Plot with Selected Child Clusters") +  # Add title
    theme_minimal()
  
} else if (display_option == "full") {
  # Use the full Seurat object
  seurat_plot_data <- seurat_object
  # Assuming child_clusters is a factor or character vector in seurat_object@meta.data
  unique_clusters <- unique(seurat_object@meta.data$child_clusters)
  
  # Generate random colors
  random_child_colors <- setNames(replicate(length(unique_clusters), generate_random_color(), simplify = FALSE), unique_clusters)
  
  # Generate the UMAP plot for the full Seurat object
  umap_plot <- DimPlot(
    seurat_plot_data, 
    reduction = "umap", 
    group.by = "child_clusters", 
    cols = random_child_colors  # Reuse or generate new colors if necessary
  ) +
    labs(title = "UMAP Plot with All Child Clusters") +  # Add title
    theme_minimal()
} else {
  stop("Invalid display_option. Choose 'full' or 'subset'.")
}

# Optionally add labels to the plot
umap_plot <- LabelClusters(
  plot = umap_plot, 
  id = "child_clusters", 
  repel = TRUE
)

# Save the plot to a PNG file
png(paste0("UMAP_", display_option, "_child_clusters.png"), width = 1200, height = 800)
print(umap_plot)
dev.off()

# CHILD DISPLAY END ---------------------------------------------------------------------------------------------
library(colorspace)  # For the `darken` function

# Assuming child_clusters and child_borders are already assigned
child_clusters <- shc_result$child_clusters
child_borders <- shc_result$child_borders

# Assign cell names to child_clusters
names(child_clusters) <- Cells(seurat_object)  
child_clusters <- as.factor(child_clusters)  # Convert to a factor if not already
seurat_object@meta.data$child_clusters <- child_clusters  # Add child_clusters to the Seurat object
seurat_object@meta.data$child_borders <- child_borders

# Extract UMAP coordinates
umap_coords <- Embeddings(seurat_object, "umap")

# Create a dataframe for cells and whether they are in child_borders (based on index matching)
cell_info <- data.frame(
  cell_name = Cells(seurat_object),
  cluster = child_clusters,
  umap_1 = umap_coords[, 1],
  umap_2 = umap_coords[, 2],
  is_border = !is.na(child_borders)  # Border if child_borders at the same index is not NA
)

# Specify clusters to display
selected_clusters <- c("Cluster1", "Cluster2")  # Modify with your desired clusters

# Filter cell_info for only the selected clusters
cell_info <- subset(cell_info, cluster %in% selected_clusters)

# Generate a discrete color palette for the number of unique clusters in the filtered data
num_child_clusters <- length(unique(cell_info$cluster))
random_child_colors <- replicate(num_child_clusters, generate_random_color())
names(random_child_colors) <- unique(cell_info$cluster)  # Map colors to filtered cluster labels

# Darken the colors for the border cells
darkened_colors <- random_child_colors
border_cells <- subset(cell_info, is_border == TRUE)
darkened_colors[border_cells$cluster] <- sapply(darkened_colors[border_cells$cluster], darken, amount = 0.3)

# Plot UMAP with child clusters using the normal colors for non-border cells
umap_plot_child_clusters <- DimPlot(seurat_object, reduction = "umap", group.by = "child_clusters", cols = random_child_colors) +
  labs(title = "UMAP Plot with Selected Child Clusters") +
  theme_minimal()

# Add darker colors for the border cells
umap_plot_child_clusters <- umap_plot_child_clusters + 
  geom_point(data = border_cells, aes(x = umap_1, y = umap_2, color = cluster), 
             shape = 16, size = 1.5) +
  scale_color_manual(values = darkened_colors)

# Add cluster labels to the UMAP plot for child_clusters
umap_plot_child_clusters <- LabelClusters(plot = umap_plot_child_clusters, id = "child_clusters", repel = TRUE)

# Save the UMAP plot with child cluster labels
png("UMAP_selected_child_clusters_with_darkened_borders.png", 1200, 800)
print(umap_plot_child_clusters)
dev.off()

# END OF CHILD WITH BORDER HIGHLIGHT -------------------------------------------------------------------------------
# BEGINNING OF DISPLAYING LATEST BORDERS ---------------------------------------------------------------------------
# Load necessary libraries
library(Seurat)
library(ggplot2)
library(colorspace)  # For the `darken` function

# Function to generate random hexadecimal color
generate_random_color <- function() {
  sprintf("#%02X%02X%02X", sample(0:255, 1), sample(0:255, 1), sample(0:255, 1))
}

# Assuming you already have shc_result with latest_borders and the seurat_object
latest_borders <- shc_result$latest_borders  # List or vector where border cells are indicated
shc_clusters <- seurat_object$shc_clusters  # Keep original clusters

# Assign cell names to shc_clusters
names(shc_clusters) <- Cells(seurat_object)  
shc_clusters <- as.factor(shc_clusters)  # Convert to a factor if not already
seurat_object$shc_clusters <- shc_clusters  # Add shc_clusters to the Seurat object

# Extract UMAP coordinates
umap_coords <- Embeddings(seurat_object, "umap")

# Create a dataframe for cells and check if they are in the latest_borders
cell_info <- data.frame(
  cell_name = Cells(seurat_object),
  cluster = shc_clusters,
  umap_1 = umap_coords[, 1],
  umap_2 = umap_coords[, 2],
  is_border = sapply(latest_borders, function(x) !is.na(x))  # Check if the cell is a border (not NULL)
)

# Generate a discrete color palette for the number of unique clusters
num_clusters <- length(unique(cell_info$cluster))
random_colors <- replicate(num_clusters, generate_random_color())
names(random_colors) <- unique(cell_info$cluster)  # Map colors to the cluster labels

# Darken the colors for the border cells
darkened_colors <- random_colors
border_cells <- subset(cell_info, is_border == TRUE)
darkened_colors[border_cells$cluster] <- sapply(darkened_colors[border_cells$cluster], darken, amount = 0.3)

# Plot UMAP with the original clusters (shc_clusters)
umap_plot <- DimPlot(seurat_object, reduction = "umap", group.by = "shc_clusters", cols = random_colors) +
  labs(title = "UMAP Plot with SHC Clusters and Darkened Borders") +
  theme_minimal()

# Add darker colors for the border cells
umap_plot <- umap_plot + 
  geom_point(data = border_cells, aes(x = umap_1, y = umap_2, color = cluster), 
             shape = 16, size = 1.5) +
  scale_color_manual(values = darkened_colors)

# Add cluster labels to the UMAP plot for shc_clusters
umap_plot <- LabelClusters(plot = umap_plot, id = "shc_clusters", repel = TRUE)

# Save the UMAP plot with cluster labels and darkened border cells
png("UMAP_shc_clusters_latest_borders.png", 1200, 800)
print(umap_plot)
dev.off()
# ------------------------------------------------------------------------------------------


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

