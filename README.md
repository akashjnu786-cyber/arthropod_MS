# arthropod_MS
This repository contains all code and data necessary to reproduce the analyses and figures from the manuscript on Insect Macroevolution submitted to Molecular Phylogenetics and Evolution.


## 📊 Data Availability

- **Supplementary Figures, Trees & Extended Data:** [10.5281/zenodo.17461834](https://doi.org/10.5281/zenodo.17461834)

---

## 📖 Project Overview

This project investigates the macroevolutionary dynamics of morphological and genomic traits across major insect orders. Using a phylogenetic comparative framework, we test for:

- Phylogenetic signal and trait evolution models
- Correlations between signal and lineage diversification
- PGLS and Phylogenetic mixed models
- Testing for pulsed evolution
- Pulsed vs. gradual evolutionary models

Key findings include a fundamental decoupling between genomic and morphological evolution, with pulsed evolution in genomic architecture linked to diversification in species-rich orders.

---

## 🗂️ Repository Structure

# For genomic trait analyses on the genome backbone tree:
source("02_code/backbone_analyses/genome_backbone_tree.R") # Uses: 01_data/phylogenies/backbone_trees/genome_orders.newick

# For morphological trait analyses on the morphology backbone tree:
source("02_code/backbone_analyses/morphological_backbone_tree.R") # Uses: 01_data/phylogenies/backbone_trees/morphology_orders.newick

### Note
All the analogous scripts work the same way for genomic traits as they work for morphological traits. Only the morphological code is shown here for brevity


## 🧬 Morphology Analysis Codes

This directory contains all R scripts for morphological trait evolution analysis across six insect orders.

### 📁 Directory Structure


### 🚀 How to Run

1. **Set up paths in R scripts:**
   - Replace all file paths in the R scripts with your local directory structure
   - Example: Change `"D:/Akash Ajay/insect bite force project/Morphology/..."` to your actual path

2. **Required R packages:**
   ```r
   install.packages(c("ape", "geiger", "readxl", "openxlsx", "dplyr", 
                     "tidyr", "ggplot2", "nortest", "purrr"))

### Phylogenetic signal -

README: Phylogenetic Signal Analysis
Overview
This R script performs phylogenetic signal analysis on morphological trait data across multiple insect orders using standard metrics (Pagel's λ, Blomberg's K, Moran's I).

Input Requirements
1. Data File: Insect_Orders_Averaged.xlsx

Required sheets: Blattodea, Coleoptera, Hymenoptera, Mantodea, Odonata, Orthoptera

Required columns in each sheet:

genus, species (combined to create species identifiers)

Morphological traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l, max.bf.iBite, max.bf.specimen, mean.bf.specimen

2. Phylogenetic Trees (Newick format)

Blattodea.NWK, Coleoptera.NWK, Hymenoptera.NWK, Mantodea.NWK, Odonata.NWK, Orthoptera.NWK

Trees should be placed in respective order directories

Analysis Performed
For each insect order, the script calculates:

Pagel's λ: Tests trait evolution against Brownian motion expectation

Blomberg's K: Measures phylogenetic signal strength

Moran's I: Spatial autocorrelation test for phylogenetic patterns

Trait evolution plots: Phenograms visualizing trait evolution on phylogenies

Output Files
For each order, the script generates:

[Order]_moran_test.png/pdf: Moran's I test plots

[Order]_phylogenetic_signal_results.rds: Complete results (λ, K values with p-values)

[Order]_[trait]_trait_evolution.png: Trait-specific evolution plots


### arthropod morphological trait with parameter output
Readme - arthropod morphological trait with parameter output


This R script fits and compares four standard evolutionary models to morphological trait data across multiple insect orders using the geiger package.

Input Requirements
1. Data File: Insect_Orders_Averaged.xlsx

Required sheets: Blattodea, Coleoptera, Hymenoptera, Mantodea, Odonata, Orthoptera

Required columns in each sheet:

genus, species (combined to create species identifiers)

Morphological traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l

2. Phylogenetic Trees (Newick format)

Blattodea.NWK, Coleoptera.NWK, Hymenoptera.NWK, Mantodea.NWK, Odonata.NWK, Orthoptera.NWK

Trees should be placed in respective order directories

Models Tested
For each trait and order, the script fits and compares:

BM (Brownian Motion): Neutral drift model

OU (Ornstein-Uhlenbeck): Constrained evolution with stabilizing selection

Delta: Time-dependent evolutionary rate variation

Kappa: Punctuated versus gradual evolution model

Output
For each order, the script generates:

[Order]_model_fits.txt: Detailed output file containing:

Model parameter estimates for each trait

AIC values for model comparison

Identification of the best-fitting model for each trait


### PGLS Code

README: Morphology Phylogenetic Generalized Least Squares (PGLS) Analysis
Overview
This R script performs Phylogenetic Generalized Least Squares (PGLS) analysis to test for evolutionary correlations between morphological traits across multiple insect orders, while accounting for phylogenetic relationships.

Input Requirements
1. Data File: Insect_Orders_Averaged.xlsx

Required sheets: Blattodea, Coleoptera, Hymenoptera, Mantodea, Odonata, Orthoptera

Required columns in each sheet:

genus, species (combined to create species identifiers)

Morphological traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l

2. Phylogenetic Trees (Newick format)

Blattodea.NWK, Coleoptera.NWK, Hymenoptera.NWK, Mantodea.NWK, Odonata.NWK, Orthoptera.NWK

Trees should be placed in respective order directories

Analysis Performed
For each insect order, the script:

Tests all pairwise trait combinations (response ~ predictor)

Fits PGLS models with maximum likelihood estimation of phylogenetic parameters (λ, κ, δ)

Includes robust error handling for singular matrices and convergence issues

Automatically selects best-fitting evolutionary model when full model fails

Output Files
all_pgls_results.xlsx: Contains two sheets:

All_PGLS_Results: Complete results for all trait pairs

Significant_PGLS_Results: Only significant correlations (P < 0.05)

<order>_pgls_results.xlsx: Order-specific results (e.g., Blattodea_pgls_results.xlsx)


### Phyr model
README: Phylogenetic Mixed Model Analysis for Bite Force
Overview
This R script performs phylogenetic mixed model analysis using the phyr package to identify morphological predictors of bite force evolution across multiple insect orders, accounting for phylogenetic relationships.

Input Requirements
1. Data File: Insect_Orders_Averaged_Renamed.xlsx

Required sheets: Blattodea, Coleoptera, Hymenoptera, Mantodea, Odonata, Orthoptera, Phasmatodea

Required columns in each sheet:

genus, species (combined to create species identifiers)

Morphological traits: iBite, head.w, head.h, head.l, th.w, body.l

2. Phylogenetic Trees (Newick format)

Trees for all seven orders in respective directories

Files: Blattodea.NWK, Coleoptera.NWK, etc.

Analysis Performed
For each insect order and predictor trait, the script:

Fits phylogenetic mixed models with bite force as response variable

Tests five predictor traits: Head_W, Head_H, Head_L, Thorax_W, Body_L

Includes phylogenetic random effects using (1 | sp__) syntax

Scales predictors and handles data quality issues automatically

Output Files
<predictor>_phyr_importance_per_order_bite_force.csv: Individual results for each predictor trait

merged_phyr_importance_bite_force.csv: Combined results across all predictors

phyr_log_bite_force.txt: Detailed analysis log with warnings

<predictor>_<order>_phyr_importance_bite_force.png/pdf: Visualization plots (800 DPI)


### phyr code wing length

README: Phylogenetic Mixed Model Analysis for Wing Length
Overview
This script performs the same phylogenetic mixed model analysis as the bite force code, but with Wing Length (Wing_L) as the response variable instead of Bite Force.

Key Differences from Bite Force Analysis
Aspect	Bite Force Analysis	Wing Length Analysis
Response Variable	Bite_F	Wing_L
Predictors	Head_W, Head_H, Head_L, Thorax_W, Body_L	Head_W, Head_H, Head_L, Thorax_W, Body_L
Output Files	*_bite_force.*	*_wing_length.*
Log File	phyr_log_bite_force.txt	phyr_log_wing_length.txt
Analysis Focus	Morphological predictors of bite force	Morphological predictors of wing length
Identical Components
Input data and tree files

Model structure: response ~ predictor + (1 | sp__)

Statistical methods and error handling

Output format and visualization style

Quality control and data processing steps

Output Files
merged_phyr_importance_wing_length.csv - Combined results

<predictor>_phyr_importance_per_order_wing_length.csv - Individual predictor results

<predictor>_<order>_phyr_importance_wing_length.png/pdf - Visualization plots


### random forest bite force code
README: Random Forest Analysis for Bite Force Predictors
Overview
This R script performs Random Forest analysis to identify the most important morphological predictors of bite force across multiple insect orders.

Input Requirements
1. Data File: Insect_Orders_Averaged.xlsx

Required sheets: Blattodea, Coleoptera, Hymenoptera, Mantodea, Odonata, Orthoptera

Required morphological traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l

Analysis Performed
For each insect order, the script:

Fits Random Forest models with bite force as response variable

Uses all other morphological traits as predictors

Calculates variable importance using Mean Decrease in Accuracy

Identifies the most influential predictors for bite force in each order

Output Files
rf_importance_per_order.csv: Importance scores for all predictors across orders

<order>_rf_importance.png/pdf: Visualization of predictor importance for each order

Methodological Note
This analysis uses a machine learning approach (Random Forest) that does not incorporate phylogenetic correction, providing a complementary perspective to the phylogenetically-informed analyses.

### random forest wing length
README: Random Forest Analysis for Wing Length
This script performs the same Random Forest analysis as the bite force code, but with Wing Length as the response variable instead of Bite Force.

Key Changes:

Response Variable: wing.l (Wing Length)

Output Files: *_wing_length.*

All other parameters, methods, and output formats remain identical to the bite force analysis.

### PIC vs BM vs OU Morphology code

README: Kolmogorov-Smirnov and Anderson-Darling Test for Evolutionary Model Fit
Overview
This R script performs Anderson-Darling (AD) tests to assess whether morphological trait evolution follows Brownian Motion (BM) or Ornstein-Uhlenbeck (OU) models by testing the normality of Phylogenetically Independent Contrasts (PICs).

Input Requirements
Data File: Insect_Orders_Averaged.xlsx with sheets for 6 insect orders

Phylogenetic Trees: Newick format trees for each order

Morphological Traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l

Analysis Performed
For each trait and order, the script:

Computes Phylogenetically Independent Contrasts (PICs)

Fits both BM and OU evolutionary models

Performs Anderson-Darling tests on PIC distributions

Includes Kolmogorov-Smirnov tests for comparison

Generates distribution plots with test results

Output Files
ad_test_results.xlsx: Complete test results with pivot tables

ad_test_output.txt: Detailed analysis log

plots_AD/: Directory with PIC distribution plots for all trait-order combinations


### best model plots simulation morphology 

README: Morphology Trait Evolution Simulation
Overview
This R script simulates morphological trait evolution across seven insect orders using the best-fitting evolutionary models identified in previous analyses.

Input Requirements
Data File: Insect_Orders_Averaged.xlsx with morphological trait data

Phylogenetic Trees: Newick format trees for all seven orders

Pre-defined Model Parameters: Best-fit evolutionary models (BM, OU, kappa, EB, delta, WN) for each trait and order

Analysis Performed
For each order and trait, the script:

Simulates trait evolution using the actual best-fit model parameters

Generates two types of visualizations:

Branching phenograms: Show trait evolution along phylogenetic branches

Simple phenograms: Plot trait values against evolutionary time

Creates both individual plots and combined PDFs for each order

Output Structure

morphology_all_orders_simulations/
├── [Order_Name]/
│   ├── branching_[Trait].png
│   ├── simple_[Trait].png
│   ├── [Order]_all_traits_branching.pdf
│   ├── [Order]_all_traits_simple.pdf
│   └── [Order]_model_parameters_summary.csv

### aic evaluation morphology

README: AIC Computation for Pulsed Evolution Models
Overview
This R script computes Akaike Information Criterion (AIC) values to compare Brownian Motion (BM) against various Pulsed Evolution (PE) models for morphological trait evolution across six insect orders.

Models Compared
BM: Standard Brownian Motion

PE1: Single-rate pulsed evolution

PE1_BM: Pulsed evolution + BM background

PE2: Two-rate pulsed evolution

PE2_BM: Two-rate pulses + BM background

PE3: Three-rate pulsed evolution

PE3_BM: Three-rate pulses + BM background

Input Requirements
Data File: Insect_Orders_Averaged.xlsx with morphological traits

Phylogenetic Trees: Newick format trees for all six orders

Traits: iBite, head.w, head.h, head.l, th.w, wing.l, body.l

Output Files
all_orders_bayou_results.xlsx: Complete AIC comparison results

all_orders_bayou_output.txt: Detailed optimization logs and diagnostics


### aic revaluation code BA

README: AIC Model Selection Interpretation
Overview
This R script applies the Burnham & Anderson ΔAIC framework to interpret model selection results, comparing Brownian Motion versus Pulsed Evolution models for morphological traits.

Analysis Performed
Calculates ΔAIC values (BM AIC - Best Pulsed AIC)

Applies Burnham & Anderson (2002) interpretation criteria:

|ΔAIC| < 2: Substantial support for both models

ΔAIC = 4-7: Considerably less support for higher-AIC model

|ΔAIC| > 10: Essentially no support for higher-AIC model

Provides recommended actions and final model choices

Input Requirements
AIC Results File: bayou_results.xlsx from previous analysis

Output
morphological_aic_BA_interpretation.xlsx: Complete ΔAIC analysis with:

ΔAIC values and support levels

Recommended actions for model selection

Final model choices based on statistical evidence


### GC skew compute
README: Parallel GC Skew Analysis
Overview
This R script performs parallel computation of GC skew from genomic FASTA files using the Biostrings package. GC skew is calculated as (G - C) / (G + C).

Key Features
Parallel Processing: Uses multiple CPU cores for faster computation

Robust Error Handling: Continues processing even if individual files fail

Automatic File Detection: Finds FASTA files (.fna, .fa, .fasta) with recursive search

Comprehensive Debugging: Detailed file discovery and validation output

Input Requirements
FASTA Files: Genomic sequences in .fna, .fa, or .fasta format (gzip compressed or uncompressed)

Directory: Path containing the FASTA files

Output
gc_skew_results.csv: Table with columns:

accession: File name

G_count, C_count: Nucleotide counts

GC_skew: Calculated GC skew value

