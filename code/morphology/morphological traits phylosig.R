### Arthropod phylogenetic signal code

## PGLS
library(ape)
library(caper)
library(phytools)
library(openxlsx)

### Blattodea
# Read the Newick tree file
Blattodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK")

# Resolve polytomies (if any)
Blattodea_tree <- multi2di(Blattodea_tree)

### Coleoptera
# Read the Newick tree file
Coleoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK")

# Resolve polytomies (if any)
Coleoptera_tree <- multi2di(Coleoptera_tree)


### Hymenoptera
# Read the Newick tree file
Hymenoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK")

# Resolve polytomies (if any)
Hymenoptera_tree <- multi2di(Hymenoptera_tree)

### Mantodea
# Read the Newick tree file
Mantodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK")

# Resolve polytomies (if any)
Mantodea_tree <- multi2di(Mantodea_tree)

### Odonata
# Read the Newick tree file
Odonata_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK")

# Resolve polytomies (if any)
Odonata_tree <- multi2di(Odonata_tree)

### Orthoptera
# Read the Newick tree file
Orthoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK")

# Resolve polytomies (if any)
Orthoptera_tree <- multi2di(Orthoptera_tree)

### Phasmatodea
# Read the Newick tree file
Phasmatodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Phasmatodea/Phasmatodea.NWK")

# Resolve polytomies (if any)
Phasmatodea_tree <- multi2di(Phasmatodea_tree)

library(readxl)

# Define the file path
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"

# Get the names of all sheets in the Excel file
sheet_names <- excel_sheets(file_path)

# Read all sheets into a list of data frames
list_of_dfs <- lapply(sheet_names, function(sheet) {
  read_excel(file_path, sheet = sheet)
})

# Assign each sheet to a separate data frame in the global environment
names(list_of_dfs) <- sheet_names
list2env(list_of_dfs, envir = .GlobalEnv)

### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Blattodea data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Blattodea_data <- as.data.frame(read_excel(file_path, sheet = "Blattodea"))

# Create proper row names
rownames(Blattodea_data) <- paste(Blattodea_data$genus, Blattodea_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Blattodea_data <- Blattodea_data[, selected_traits]

## Step 2: Load and prepare Blattodea tree
Blattodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK")
Blattodea_tree <- multi2di(Blattodea_tree)
Blattodea_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Blattodea_data), Blattodea_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

Blattodea_data <- Blattodea_data[common_species, , drop = FALSE]
Blattodea_tree <- keep.tip(Blattodea_tree, common_species)

## Step 4: Create phylo4d object
Blattodea_phylotraits <- phylo4d(Blattodea_tree, Blattodea_data)

## Step 5: Moran's I test
Blattodea_moran <- abouheif.moran(Blattodea_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Blattodea_moran, main = "Blattodea Moran's I Test")
}

# Save plots
png("Blattodea_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Blattodea_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Blattodea_data)) {
  cat("\nTesting trait:", trait, "\n")
  current_trait <- setNames(Blattodea_data[,trait], rownames(Blattodea_data))
  
  # Pagel's lambda
  lambda <- phylosig(Blattodea_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Blattodea_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
}

# Save all results
saveRDS(results, "Blattodea_phylogenetic_signal_results.rds")

### Coleoptera
### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Coleoptera data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Coleoptera_data <- as.data.frame(read_excel(file_path, sheet = "Coleoptera"))

# Create proper row names
rownames(Coleoptera_data) <- paste(Coleoptera_data$genus, Coleoptera_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Coleoptera_data <- Coleoptera_data[, selected_traits]

## Step 2: Load and prepare Coleoptera tree
Coleoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK")
Coleoptera_tree <- multi2di(Coleoptera_tree)
Coleoptera_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Coleoptera_data), Coleoptera_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

Coleoptera_data <- Coleoptera_data[common_species, , drop = FALSE]
Coleoptera_tree <- keep.tip(Coleoptera_tree, common_species)

## Step 4: Create phylo4d object
Coleoptera_phylotraits <- phylo4d(Coleoptera_tree, Coleoptera_data)

## Step 5: Moran's I test
Coleoptera_moran <- abouheif.moran(Coleoptera_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Coleoptera_moran, main = "Coleoptera Moran's I Test")
}

# Save plots
png("Coleoptera_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Coleoptera_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Coleoptera_data)) {
  cat("\nTesting trait:", trait, "\n")
  current_trait <- setNames(Coleoptera_data[,trait], rownames(Coleoptera_data))
  
  # Pagel's lambda
  lambda <- phylosig(Coleoptera_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Coleoptera_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
}

# Save all results
saveRDS(results, "Coleoptera_phylogenetic_signal_results.rds")

# Print final confirmation
cat("\nColeoptera analysis completed successfully!\n")
cat("Results saved in:\n")
cat("- Coleoptera_moran_test.png/pdf\n")
cat("- Coleoptera_phylogenetic_signal_results.rds\n")

## Hymenoptera
### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Hymenoptera data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Hymenoptera_data <- as.data.frame(read_excel(file_path, sheet = "Hymenoptera"))

# Create proper row names
rownames(Hymenoptera_data) <- paste(Hymenoptera_data$genus, Hymenoptera_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Hymenoptera_data <- Hymenoptera_data[, selected_traits]

## Step 2: Load and prepare Hymenoptera tree
Hymenoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK")
Hymenoptera_tree <- multi2di(Hymenoptera_tree)
Hymenoptera_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Hymenoptera_data), Hymenoptera_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

if(length(common_species) < 5) {
  warning(paste("Only", length(common_species), "species matched between tree and data!"))
}

Hymenoptera_data <- Hymenoptera_data[common_species, , drop = FALSE]
Hymenoptera_tree <- keep.tip(Hymenoptera_tree, common_species)

## Step 4: Create phylo4d object
Hymenoptera_phylotraits <- phylo4d(Hymenoptera_tree, Hymenoptera_data)

## Step 5: Moran's I test
Hymenoptera_moran <- abouheif.moran(Hymenoptera_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Hymenoptera_moran, main = "Hymenoptera Moran's I Test")
}

# Save plots
png("Hymenoptera_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Hymenoptera_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Hymenoptera_data)) {
  cat("\n=== Testing trait:", trait, "===\n")
  current_trait <- setNames(Hymenoptera_data[,trait], rownames(Hymenoptera_data))
  
  # Pagel's lambda
  lambda <- phylosig(Hymenoptera_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Hymenoptera_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
  
  # Plot trait evolution
  png(paste0("Hymenoptera_", trait, "_trait_evolution.png"), width = 10, height = 6, units = "in", res = 300)
  phenogram(Hymenoptera_tree, current_trait, main = paste(trait, "trait evolution"))
  dev.off()
}

# Save all results
saveRDS(results, "Hymenoptera_phylogenetic_signal_results.rds")

# Print final summary
cat("\n\n===== Hymenoptera Analysis Summary =====\n")
cat("Analyzed", length(colnames(Hymenoptera_data)), "traits across", length(common_species), "species\n")
cat("Output files created:\n")
cat("- Hymenoptera_moran_test.png/pdf\n")
cat("- Hymenoptera_phylogenetic_signal_results.rds\n")
cat("- Trait-specific evolution plots (PNG)\n")

## Mantodea 
### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Mantodea data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Mantodea_data <- as.data.frame(read_excel(file_path, sheet = "Mantodea"))

# Create proper row names
rownames(Mantodea_data) <- paste(Mantodea_data$genus, Mantodea_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Mantodea_data <- Mantodea_data[, selected_traits]

## Step 2: Load and prepare Mantodea tree
Mantodea_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK")
Mantodea_tree <- multi2di(Mantodea_tree)
Mantodea_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Mantodea_data), Mantodea_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

if(length(common_species) < 5) {
  warning(paste("Only", length(common_species), "species matched between tree and data!"))
}

Mantodea_data <- Mantodea_data[common_species, , drop = FALSE]
Mantodea_tree <- keep.tip(Mantodea_tree, common_species)

## Step 4: Create phylo4d object
Mantodea_phylotraits <- phylo4d(Mantodea_tree, Mantodea_data)

## Step 5: Moran's I test
Mantodea_moran <- abouheif.moran(Mantodea_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Mantodea_moran, main = "Mantodea Moran's I Test")
}

# Save plots
png("Mantodea_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Mantodea_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Mantodea_data)) {
  cat("\n=== Testing trait:", trait, "===\n")
  current_trait <- setNames(Mantodea_data[,trait], rownames(Mantodea_data))
  
  # Pagel's lambda
  lambda <- phylosig(Mantodea_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Mantodea_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
  
  # Plot trait evolution
  png(paste0("Mantodea_", trait, "_trait_evolution.png"), width = 10, height = 6, units = "in", res = 300)
  phenogram(Mantodea_tree, current_trait, main = paste(trait, "trait evolution"))
  dev.off()
}

# Save all results
saveRDS(results, "Mantodea_phylogenetic_signal_results.rds")

# Print final summary
cat("\n\n===== Mantodea Analysis Summary =====\n")
cat("Analyzed", length(colnames(Mantodea_data)), "traits across", length(common_species), "species\n")
cat("Output files created:\n")
cat("- Mantodea_moran_test.png/pdf\n")
cat("- Mantodea_phylogenetic_signal_results.rds\n")
cat("- Trait-specific evolution plots (PNG)\n")

### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Odonata data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Odonata_data <- as.data.frame(read_excel(file_path, sheet = "Odonata"))

# Create proper row names
rownames(Odonata_data) <- paste(Odonata_data$genus, Odonata_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Odonata_data <- Odonata_data[, selected_traits]

## Step 2: Load and prepare Odonata tree
Odonata_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK")
Odonata_tree <- multi2di(Odonata_tree)
Odonata_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Odonata_data), Odonata_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

if(length(common_species) < 5) {
  warning(paste("Only", length(common_species), "species matched between tree and data!"))
}

Odonata_data <- Odonata_data[common_species, , drop = FALSE]
Odonata_tree <- keep.tip(Odonata_tree, common_species)

## Step 4: Create phylo4d object
Odonata_phylotraits <- phylo4d(Odonata_tree, Odonata_data)

## Step 5: Moran's I test
Odonata_moran <- abouheif.moran(Odonata_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Odonata_moran, main = "Odonata Moran's I Test")
}

# Save plots
png("Odonata_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Odonata_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Odonata_data)) {
  cat("\n=== Testing trait:", trait, "===\n")
  current_trait <- setNames(Odonata_data[,trait], rownames(Odonata_data))
  
  # Pagel's lambda
  lambda <- phylosig(Odonata_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Odonata_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
  
  # Plot trait evolution
  png(paste0("Odonata_", trait, "_trait_evolution.png"), width = 10, height = 6, units = "in", res = 300)
  phenogram(Odonata_tree, current_trait, main = paste(trait, "trait evolution"))
  dev.off()
}

# Save all results
saveRDS(results, "Odonata_phylogenetic_signal_results.rds")

# Print final summary
cat("\n\n===== Odonata Analysis Summary =====\n")
cat("Analyzed", length(colnames(Odonata_data)), "traits across", length(common_species), "species\n")
cat("Output files created:\n")
cat("- Odonata_moran_test.png/pdf\n")
cat("- Odonata_phylogenetic_signal_results.rds\n")
cat("- Trait-specific evolution plots (PNG)\n"


### Orthoptera
### Load required packages
library(picante)
library(ape)
library(adephylo)
library(ade4)
library(phylobase)
library(geiger)
library(phytools)
library(readxl)

## Step 1: Load and prepare Orthoptera data
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
Orthoptera_data <- as.data.frame(read_excel(file_path, sheet = "Orthoptera"))

# Create proper row names
rownames(Orthoptera_data) <- paste(Orthoptera_data$genus, Orthoptera_data$species, sep = "_")

# Select only the specified numeric traits
selected_traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", 
                     "wing.l", "body.l", "max.bf.iBite", "max.bf.specimen", 
                     "mean.bf.specimen")

Orthoptera_data <- Orthoptera_data[, selected_traits]

## Step 2: Load and prepare Orthoptera tree
Orthoptera_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK")
Orthoptera_tree <- multi2di(Orthoptera_tree)
Orthoptera_tree$node.label <- NULL

## Step 3: Match data and tree
common_species <- intersect(rownames(Orthoptera_data), Orthoptera_tree$tip.label)
cat("Number of species in analysis:", length(common_species), "\n")

if(length(common_species) < 5) {
  warning(paste("Only", length(common_species), "species matched between tree and data!"))
}

Orthoptera_data <- Orthoptera_data[common_species, , drop = FALSE]
Orthoptera_tree <- keep.tip(Orthoptera_tree, common_species)

## Step 4: Create phylo4d object
Orthoptera_phylotraits <- phylo4d(Orthoptera_tree, Orthoptera_data)

## Step 5: Moran's I test
Orthoptera_moran <- abouheif.moran(Orthoptera_phylotraits, method = "Abouheif")

# Plot function
plot_moran <- function() {
  plot(Orthoptera_moran, main = "Orthoptera Moran's I Test")
}

# Save plots
png("Orthoptera_moran_test.png", width = 8, height = 6, units = "in", res = 300)
plot_moran()
dev.off()

pdf("Orthoptera_moran_test.pdf", width = 8, height = 6)
plot_moran()
dev.off()

## Step 6: Phylogenetic signal tests for each trait
results <- list()

for(trait in colnames(Orthoptera_data)) {
  cat("\n=== Testing trait:", trait, "===\n")
  current_trait <- setNames(Orthoptera_data[,trait], rownames(Orthoptera_data))
  
  # Pagel's lambda
  lambda <- phylosig(Orthoptera_tree, current_trait, method = "lambda", test = TRUE, nsim = 999)
  cat("Lambda:", lambda$lambda, "p-value:", lambda$P, "\n")
  
  # Blomberg's K
  K <- phylosig(Orthoptera_tree, current_trait, method = "K", test = TRUE, nsim = 999)
  cat("K:", K$K, "p-value:", K$P, "\n")
  
  # Store results
  results[[trait]] <- list(lambda = lambda, K = K)
  
  # Plot trait evolution
  png(paste0("Orthoptera_", trait, "_trait_evolution.png"), width = 10, height = 6, units = "in", res = 300)
  phenogram(Orthoptera_tree, current_trait, main = paste(trait, "trait evolution"))
  dev.off()
}

# Save all results
saveRDS(results, "Orthoptera_phylogenetic_signal_results.rds")

# Print final summary
cat("\n\n===== Orthoptera Analysis Summary =====\n")
cat("Analyzed", length(colnames(Orthoptera_data)), "traits across", length(common_species), "species\n")
cat("Output files created:\n")
cat("- Orthoptera_moran_test.png/pdf\n")
cat("- Orthoptera_phylogenetic_signal_results.rds\n")
cat("- Trait-specific evolution plots (PNG)\n")
    