---
title: "quality_control"
author: "Anushree Hosahalli Nagarajarao"
date: "2024-05-05"
output: html_document
---

```{r}
install.packages("SingleCellExperiment")
install.packages("Seurat")
install.packages("tidyverse")
install.packages("Matrix")
install.packages("scales")
install.packages("cowplot")
install.packages("RCurl")

```
```{r}
# Load libraries
library(SingleCellExperiment)
library(Seurat)
library(tidyverse)
library(Matrix)
library(scales)
library(cowplot)
library(RCurl)
```
```{r}

# Load single cell seq data and create Seurat objects
for (file in c("ctrl_raw_feature_bc_matrix", "stim_raw_feature_bc_matrix")) {
  seurat_data <- Read10X(data.dir = paste0("data/", file))
  seurat_obj <- CreateSeuratObject(counts = seurat_data,
                                   min.features = 100,
                                   project = file)
  assign(file, seurat_obj)
}

# Check the metadata in the new Seurat objects
head(get("ctrl_raw_feature_bc_matrix")@meta.data)
head(get("stim_raw_feature_bc_matrix")@meta.data)

```
```{r}
# Create a merged Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))
```
```{r}
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
```
```{r}
# Get the index of the cell with barcode AAACATACATTTCC-1
ctrl_cell_index <- which(colnames(ctrl_raw_feature_bc_matrix) == "AAACATACATTTCC-1")

# Get the metadata for the cell with barcode AAACATACATTTCC-1
ctrl_cell_metadata <- ctrl_raw_feature_bc_matrix@meta.data[ctrl_cell_index, ]

# Print nCount_RNA and nFeature_RNA
print(paste("nCount_RNA for AAACATACATTTCC-1:", ctrl_cell_metadata$nCount_RNA))
print(paste("nFeature_RNA for AAACATACATTTCC-1:", ctrl_cell_metadata$nFeature_RNA))

```
```{r}
# Create a merged Seurat object
merged_seurat <- merge(x = ctrl_raw_feature_bc_matrix, 
                       y = stim_raw_feature_bc_matrix, 
                       add.cell.id = c("ctrl", "stim"))
```
```{r}
# Check that the merged object has the appropriate sample-specific prefixes
head(merged_seurat@meta.data)
tail(merged_seurat@meta.data)
```
```{r}
# Explore merged metadata
View(merged_seurat@meta.data)
```
```{r}
# Add number of genes per UMI for each cell to metadata
merged_seurat$log10GenesPerUMI <- log10(merged_seurat$nFeature_RNA) / log10(merged_seurat$nCount_RNA)
```
```{r}
# Compute percent mito ratio
merged_seurat$mitoRatio <- PercentageFeatureSet(object = merged_seurat, pattern = "^MT-")
merged_seurat$mitoRatio <- merged_seurat@meta.data$mitoRatio / 100

```
```{r}
# Create metadata dataframe
metadata <- merged_seurat@meta.data
```
```{r}
# Add cell IDs to metadata
metadata$cells <- rownames(metadata)
```
```{r}
# Create sample column
metadata$sample <- NA
metadata$sample[which(str_detect(metadata$cells, "^ctrl_"))] <- "ctrl"
metadata$sample[which(str_detect(metadata$cells, "^stim_"))] <- "stim"
```
```{r}
# Rename columns
metadata <- metadata %>%
        dplyr::rename(seq_folder = orig.ident,
                      nUMI = nCount_RNA,
                      nGene = nFeature_RNA)

```
```{r}
# Add metadata back to Seurat object
merged_seurat@meta.data <- metadata

# Create .RData object to load at any time
save(merged_seurat, file = "data/merged_filtered_seurat.RData")

```
```{r}
# Visualize the number of cell counts per sample
metadata %>% 
  	ggplot(aes(x=sample, fill=sample)) + 
  	geom_bar() +
  	theme_classic() +
  	theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  	theme(plot.title = element_text(hjust=0.5, face="bold")) +
  	ggtitle("NCells")
```
```{r}
# Visualize the number UMIs/transcripts per cell
metadata %>% 
  	ggplot(aes(color=sample, x=nUMI, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	ylab("Cell density") +
  	geom_vline(xintercept = 500)
```
```{r}
# Visualize the distribution of genes detected per cell via histogram
metadata %>% 
  	ggplot(aes(color=sample, x=nGene, fill= sample)) + 
  	geom_density(alpha = 0.2) + 
  	theme_classic() +
  	scale_x_log10() + 
  	geom_vline(xintercept = 300)

```
```{r}
# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI (novelty score)
metadata %>%
  	ggplot(aes(x=log10GenesPerUMI, color = sample, fill=sample)) +
  	geom_density(alpha = 0.2) +
  	theme_classic() +
  	geom_vline(xintercept = 0.8)
```
```{r}
# Visualize the distribution of mitochondrial gene expression detected per cell
metadata %>% 
  	ggplot(aes(color=sample, x=mitoRatio, fill=sample)) + 
  	geom_density(alpha = 0.2) + 
  	scale_x_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 0.2)
```
```{r}
# Visualize the correlation between genes detected and number of UMIs and determine whether strong presence of cells with low numbers of genes/UMIs
metadata %>% 
  	ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  	geom_point() + 
	scale_colour_gradient(low = "gray90", high = "black") +
  	stat_smooth(method=lm) +
  	scale_x_log10() + 
  	scale_y_log10() + 
  	theme_classic() +
  	geom_vline(xintercept = 500) +
  	geom_hline(yintercept = 250) +
  	facet_wrap(~sample)
```
```{r}
# Filter out low quality cells using selected thresholds - these will change with experiment
filtered_seurat <- subset(x = merged_seurat, 
                         subset= (nUMI >= 500) & 
                           (nGene >= 250) & 
                           (log10GenesPerUMI > 0.80) & 
                           (mitoRatio < 0.20))
```
```{r}
# Extract count matrix from the RNA assay using counts
counts <- GetAssayData(object = filtered_seurat, assay = "RNA", slot = "counts.ctrl_raw_feature_bc_matrix")

# Extract count matrix from the RNA assay using scounts
scounts <- GetAssayData(object = filtered_seurat, assay = "RNA", slot = "counts.ctrl_raw_feature_bc_matrix")

# Output a logical matrix specifying for each gene whether there are more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

```
```{r}
# Reassign to filtered Seurat object
filtered_seurat <- CreateSeuratObject(filtered_counts, meta.data = filtered_seurat@meta.data)

```
```{r}
# Save filtered subset to new metadata
metadata_clean <- filtered_seurat@meta.data

```
```{r}
dim(filtered_counts)
dim(metadata_clean)

```
```{r}
# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")
```

```{r}
filtered_seurat <- NormalizeData(filtered_seurat, normalization.method = "LogNormalize", scale.factor = 10000)


```
```{r}
filtered_seurat <- FindVariableFeatures(filtered_seurat, selection.method = "vst", nfeatures = 2000)

```
```{r}
filtered_seurat <- ScaleData(filtered_seurat, features = rownames(filtered_seurat), vars.to.regress = "mitoRatio")


```
```{r}
filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(object = filtered_seurat))

```
```{r}

filtered_seurat <- FindNeighbors(filtered_seurat, dims = 1:10)
filtered_seurat <- FindClusters(filtered_seurat, resolution = 0.5)



```
```{r}
filtered_seurat <- RunUMAP(filtered_seurat, dims = 1:10)
DimPlot(filtered_seurat, reduction = "umap")

```
```{r}
# Perform PCA if not already done
if (!"pca" %in% names(filtered_seurat@reductions)) {
  filtered_seurat <- RunPCA(filtered_seurat, features = VariableFeatures(object = filtered_seurat))
}

```
```{r}
# Visualize variance explained by each principal component
ElbowPlot(filtered_seurat)

```
```{r}
# Plot PCA results, coloring by sample to see if different conditions cluster separately
PCAPlot <- DimPlot(filtered_seurat, reduction = "pca", group.by = "sample")
PCAPlot

```
```{r}
# Check if any outliers are identified
if (length(outliers) > 0) {
  print(paste("Number of outliers identified:", length(outliers)))
  
  # Plot PCA without outliers
  PCA_plot <- DimPlot(filtered_seurat, group.by = "outliers", cells = outliers, label = TRUE)
  
  # Highlight outliers in the plot
  PCA_plot + geom_point(data = outliers, aes(color = "Outliers")) +
    ggtitle("PCA plot with identified outliers") +
    theme(legend.position = "bottom")
} else {
  print("No outliers were identified based on the given thresholds.")
}


```


```{r}
# Check for statistically significant principal components using JackStraw
if (!"jackStraw" %in% names(filtered_seurat@misc)) {
    filtered_seurat <- JackStraw(filtered_seurat, num.replicate = 100)
}

```


```{r}
# Make sure to score the JackStraw results to test the significance of each PC
filtered_seurat <- ScoreJackStraw(filtered_seurat, dims = 1:20)

# Now you can plot the JackStraw results
JackStrawPlot(filtered_seurat, dims = 1:20)



```
```{r}
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)

```

```{r}
# Load cell cycle markers
load("data/cycle.rda")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)
# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data) 

```
```{r}
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)

```
```{r}
# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)

```
```{r}
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))

```
```{r}
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")

```
```{r}
# Check which assays are stored in objects
split_seurat$ctrl@assays

```
```{r}
# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")

```






```{r}
# Load required libraries
library(ggplot2)

# Perform QC plots using the filtered data
# For example, let's plot Number of Genes vs Number of UMIs
plot1 <- ggplot(metadata_clean, aes(x = nGene, y = nUMI)) +
  geom_point() +
  labs(x = "Number of Genes", y = "Number of UMIs") +
  ggtitle("QC Plot: Number of Genes vs Number of UMIs")

# Display the plot
print(plot1)

```





```{r}
# Save filtered subset to new metadata
 metadata_clean <- filtered_seurat@meta.data


```
```{r}
# Example QC plot using filtered data
plot_filtered <- ggplot(metadata_clean, aes(x = nGene, y = nUMI)) +
  geom_point() +
  labs(x = "Number of Genes", y = "Number of UMIs") +
  ggtitle("Filtered Data: QC Plot - Number of Genes vs Number of UMIs")


```
```{r}
# Calculate the number of cells left for each sample
num_cells_left <- table(metadata_clean$sample)

# Calculate the total number of cells left
total_cells_left <- sum(num_cells_left)

# Compare with the total number of cells loaded for the experiment
total_cells_loaded <- 12000  # Total number of cells loaded for the experiment

# Print the number of cells left for each sample and the comparison
print(num_cells_left)
if (total_cells_left > total_cells_loaded) {
  cat("The number of cells removed is relatively low compared to the total number of cells loaded.\n")
} else {
  cat("The number of cells removed is relatively high compared to the total number of cells loaded.\n")
}

```
```{r}
# Create .RData object to load at any time
save(filtered_seurat, file="data/seurat_filtered.RData")

```
```{r}
# Load necessary libraries
library(Seurat)

# Load your merged and filtered Seurat object
load("data/merged_filtered_seurat.RData")

# Assuming you have already calculated additional metrics and updated metadata

merged_seurat <- FindVariableFeatures(merged_seurat)
merged_seurat <- NormalizeData(merged_seurat)
merged_seurat <- ScaleData(merged_seurat, features = VariableFeatures(object = merged_seurat))
memory.limit(size = xxx)  # Set the memory limit to a higher value

# Perform PCA
merged_seurat <- RunPCA(merged_seurat, features = VariableFeatures(object = merged_seurat))

# Find neighbors
merged_seurat <- FindNeighbors(merged_seurat, dims = 1:30)

# Cluster cells
merged_seurat <- FindClusters(merged_seurat, resolution = 0.6)

# Compute UMAP
merged_seurat <- RunUMAP(merged_seurat, dims = 1:30)

# Visualize clusters using UMAP
DimPlot(merged_seurat, reduction = "umap", label = TRUE)


```
```{r}
# Assuming you have loaded your Seurat object
# Scale data
# Scale data
seurat_data <- ScaleData(seurat_data)

# Perform PCA
pca_result <- RunPCA(seurat_object)


```
```{r}
# Normalize the counts
seurat_phase <- NormalizeData(filtered_seurat)
```
```{r}
# Load cell cycle markers
load("data/cycle.rda")

# Verify and update cell cycle marker gene lists if necessary
# Example of updating gene lists:
# g2m_genes <- c("gene1", "gene2", "gene3")
# s_genes <- c("gene4", "gene5", "gene6")

# Score cells for cell cycle
seurat_phase <- CellCycleScoring(seurat_phase, 
                                 g2m.features = g2m_genes, 
                                 s.features = s_genes)

# View cell cycle scores and phases assigned to cells                                 
View(seurat_phase@meta.data)  
 
```
```{r}
# Identify the most variable genes
seurat_phase <- FindVariableFeatures(seurat_phase, 
                     selection.method = "vst",
                     nfeatures = 2000, 
                     verbose = FALSE)
		     
# Scale the counts
seurat_phase <- ScaleData(seurat_phase)
```
```{r}
# Identify the 15 most highly variable genes
ranked_variable_genes <- VariableFeatures(seurat_phase)
top_genes <- ranked_variable_genes[1:15]

# Plot the average expression and variance of these genes
# With labels to indicate which genes are in the top 15
p <- VariableFeaturePlot(seurat_phase)
LabelPoints(plot = p, points = top_genes, repel = TRUE)
```
```{r}
# Perform PCA
seurat_phase <- RunPCA(seurat_phase)

# Plot the PCA colored by cell cycle phase
DimPlot(seurat_phase,
        reduction = "pca",
        group.by= "Phase",
        split.by = "Phase")
```
```{r}
# Check quartile values
summary(seurat_phase@meta.data$mitoRatio)

# Turn mitoRatio into categorical factor vector based on quartile values
seurat_phase@meta.data$mitoFr <- cut(seurat_phase@meta.data$mitoRatio, 
                   breaks=c(-Inf, 0.0144, 0.0199, 0.0267, Inf), 
                   labels=c("Low","Medium","Medium high", "High"))
```
```{r}
# Plot PCA colored by mitoFr
DimPlot(seurat_phase,
        reduction = "pca",
        group.by = "mitoFr",
        split.by = "mitoFr")

```
```{r}
# Split seurat object by condition to perform cell cycle scoring and SCT on all samples
split_seurat <- SplitObject(seurat_phase, split.by = "sample")
```
```{r}
ctrl_reps <- split_seurat[c("ctrl_1", "ctrl_2")]
```
```{r}
options(future.globals.maxSize = 4000 * 1024^2)
```
```{r}
for (i in 1:length(split_seurat)) {
    split_seurat[[i]] <- SCTransform(split_seurat[[i]], vars.to.regress = c("mitoRatio"), vst.flavor = "v2")
    }
```
```{r}
# Check which assays are stored in objects
split_seurat$ctrl@assays
```
```{r}
# Save the split seurat object
saveRDS(split_seurat, "data/split_seurat.rds")
```


```{r}
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")
```


```{r}
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")

# Check the structure of the loaded object
str(split_seurat)

# Access specific elements or properties
# For example, to check the assays:
# Access specific elements or properties
# For example, to check the assays for each Seurat object within the list:
for (i in seq_along(split_seurat)) {
  cat("Assays for Seurat object", i, ":", names(split_seurat[[i]]@assays), "\n")
}





```
```{r}



```

```{r}
library(Seurat)

# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")

# Check the structure of the loaded object
str(split_seurat)
```
```{r}
# Select the most variable features to use for integration
integ_features <- SelectIntegrationFeatures(object.list = split_seurat, 
                                            nfeatures = 3000) 
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)

```
```{r}
# Load the split seurat object into the environment
split_seurat <- readRDS("data/split_seurat.rds")

```
```{r}
## Don't run this during class
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
```
```{r}
## Don't run this during class
# Prepare the SCT list object for integration
split_seurat <- PrepSCTIntegration(object.list = split_seurat, 
                                   anchor.features = integ_features)
## Don't run this during class
# Find best buddies - can take a while to run
integ_anchors <- FindIntegrationAnchors(object.list = split_seurat, 
                                        normalization.method = "SCT", 
                                        anchor.features = integ_features)
 

```
```{r}


```





