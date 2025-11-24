library(ape)
library(geiger)
library(readxl)
library(ggplot2)
library(phytools)
library(gridExtra)
library(dplyr)

# PATHS & SETTINGS
file_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/genomic_features_with_gc_skew1.xlsx"

orders <- c(
  "Blattodea","Coleoptera","Diptera","Hemiptera","Hymenoptera",
  "Lepidoptera","Odonata","Orthoptera","Phasmatodea","Psocodea","Trichoptera"
)

tree_paths <- list(
  Blattodea   = "D:/Akash Ajay/insect bite force project/genome features/Blattodea/Blattodea.NWK",
  Coleoptera  = "D:/Akash Ajay/insect bite force project/genome features/Coleoptera/Coleoptera.NWK",
  Diptera     = "D:/Akash Ajay/insect bite force project/genome features/Diptera/Diptera.NWK",
  Hemiptera   = "D:/Akash Ajay/insect bite force project/genome features/Hemiptera/Hemiptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/genome features/Hymenoptera/Hymenoptera.NWK",
  Lepidoptera = "D:/Akash Ajay/insect bite force project/genome features/Lepidoptera/Lepidoptera.NWK",
  Odonata     = "D:/Akash Ajay/insect bite force project/genome features/Odonata/Odonata.NWK",
  Orthoptera  = "D:/Akash Ajay/insect bite force project/genome features/Orthoptera/Orthoptera.NWK",
  Phasmatodea = "D:/Akash Ajay/insect bite force project/genome features/Phasmatodea/Phasmatodea.NWK",
  Psocodea    = "D:/Akash Ajay/insect bite force project/genome features/Psocodea/Psocodea.NWK",
  Trichoptera = "D:/Akash Ajay/insect bite force project/genome features/Trichoptera/Trichoptera.NWK"
)

# Create model parameters from your genomic features data
model_params_all <- list(
  Blattodea = list(
    "Coding_genes" = list(model = "delta", aic = 100.206, param1 = "delta", param1_value = 0.014696, param2 = "sig2", param2_value = 423893.97)
  ),
  Coleoptera = list(
    "Genome_Size" = list(model = "kappa", aic = 1675.125, param1 = "kappa", param1_value = 0.118765, param2 = "sig2", param2_value = 60918.839),
    "Genome_GC" = list(model = "kappa", aic = 470.477, param1 = "kappa", param1_value = 0.354741, param2 = "sig2", param2_value = 0.470599),
    "Chromosome_Number" = list(model = "OU", aic = 405.574, param1 = "alpha", param1_value = 0.794751, param2 = "sig2", param2_value = 58.57339),
    "Coding_genes" = list(model = "BM", aic = 336.079, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 334620.115),
    "Non_coding_genes" = list(model = "BM", aic = 232.508, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 13851.867)
  ),
  Diptera = list(
    "Genome_Size" = list(model = "kappa", aic = 4111.48, param1 = "kappa", param1_value = 0.354165, param2 = "sig2", param2_value = 2598.51),
    "Genome_GC" = list(model = "kappa", aic = 1287.89, param1 = "kappa", param1_value = 0.53183, param2 = "sig2", param2_value = 0.404031),
    "Chromosome_Number" = list(model = "kappa", aic = 290.533, param1 = "kappa", param1_value = 0.114369, param2 = "sig2", param2_value = 0.248949),
    "Coding_genes" = list(model = "kappa", aic = 1372.966, param1 = "kappa", param1_value = 0.595485, param2 = "sig2", param2_value = 290976.278),
    "Non_coding_genes" = list(model = "kappa", aic = 1361.137, param1 = "kappa", param1_value = 0.495959, param2 = "sig2", param2_value = 503833.798)
  ),
  Hemiptera = list(
    "Genome_Size" = list(model = "kappa", aic = 1216.317, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 234845.479),
    "Genome_GC" = list(model = "BM", aic = 318.527, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.037832),
    "Chromosome_Number" = list(model = "kappa", aic = 240.602, param1 = "kappa", param1_value = 0.089697, param2 = "sig2", param2_value = 4.501875),
    "Coding_genes" = list(model = "OU", aic = 489.38, param1 = "alpha", param1_value = 0.04298, param2 = "sig2", param2_value = 2821984.836),
    "Non_coding_genes" = list(model = "OU", aic = 364.271, param1 = "alpha", param1_value = 1.394913, param2 = "sig2", param2_value = 4191710.374)
  ),
  Hymenoptera = list(
    "Genome_Size" = list(model = "kappa", aic = 2348.784, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 22021.076),
    "Genome_GC" = list(model = "kappa", aic = 750.542, param1 = "kappa", param1_value = 0.128279, param2 = "sig2", param2_value = 1.369143),
    "Chromosome_Number" = list(model = "OU", aic = 387.13, param1 = "alpha", param1_value = 0.021713, param2 = "sig2", param2_value = 1.444489),
    "Coding_genes" = list(model = "OU", aic = 883.766, param1 = "alpha", param1_value = 0.026362, param2 = "sig2", param2_value = 286096.741),
    "Non_coding_genes" = list(model = "OU", aic = 802.93, param1 = "alpha", param1_value = 0.270549, param2 = "sig2", param2_value = 732286.085)
  ),
  Lepidoptera = list(
    "Genome_Size" = list(model = "kappa", aic = 7420.719, param1 = "kappa", param1_value = 0.540379, param2 = "sig2", param2_value = 2251.731),
    "Genome_GC" = list(model = "kappa", aic = 1708.471, param1 = "kappa", param1_value = 0.424805, param2 = "sig2", param2_value = 0.170868),
    "Chromosome_Number" = list(model = "kappa", aic = 1998.044, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 13.465618),
    "Coding_genes" = list(model = "OU", aic = 834.84, param1 = "alpha", param1_value = 1.82915, param2 = "sig2", param2_value = 54644458.5),
    "Non_coding_genes" = list(model = "kappa", aic = 671.186, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 942184.357)
  ),
  Odonata = list(
    "Genome_Size" = list(model = "BM", aic = 174.223, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 659.601),
    "Genome_GC" = list(model = "delta", aic = 37.222, param1 = "delta", param1_value = 0.258096, param2 = "sig2", param2_value = 0.013582),
    "Chromosome_Number" = list(model = "BM", aic = 23.993, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.010968)
  ),
  Orthoptera = list(
    "Genome_Size" = list(model = "OU", aic = 155.117, param1 = "alpha", param1_value = 0.024862, param2 = "sig2", param2_value = 382964.449),
    "Genome_GC" = list(model = "BM", aic = 20.825, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.003491)
  ),
  Phasmatodea = list(
    "Genome_Size" = list(model = "kappa", aic = 201.212, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 298253.008),
    "Genome_GC" = list(model = "OU", aic = 49.441, param1 = "alpha", param1_value = 2.718282, param2 = "sig2", param2_value = 11.885444),
    "Chromosome_Number" = list(model = "BM", aic = 42.626, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.129163),
    "Coding_genes" = list(model = "BM", aic = 84.338, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 379318.541),
    "Non_coding_genes" = list(model = "BM", aic = 72.965, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 22092.32)
  ),
  Psocodea = list(
    "Genome_Size" = list(model = "BM", aic = 41.72, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 3.476779),
    "Genome_GC" = list(model = "BM", aic = 22.776, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.030499)
  ),
  Trichoptera = list(
    "Genome_Size" = list(model = "BM", aic = 127.998, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 285.521),
    "Genome_GC" = list(model = "BM", aic = 40.105, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.016382)
  )
)

# Trait mapping for genomic features
trait_map <- c(
  "Genome_Size" = "Genome_Size",
  "Genome_GC" = "Genome_GC",
  "Chromosome Number" = "Chromosome_Number",
  "Coding genes" = "Coding_genes",
  "Non - coding genes" = "Non_coding_genes"
)

# Pretty display names
pretty_names <- c(
  "Genome_Size" = "Genome Size",
  "Genome_GC" = "Genome GC Content",
  "Chromosome_Number" = "Chromosome Number",
  "Coding_genes" = "Coding Genes",
  "Non_coding_genes" = "Non-coding Genes"
)

# IMPROVED function to prepare tree with better zero branch length handling
prepare_tree <- function(tree_path) {
  tryCatch({
    tree <- read.tree(tree_path)
    
    # Convert to bifurcating tree
    tree <- multi2di(tree)
    tree$node.label <- NULL
    
    # MORE ROBUST handling of zero branch lengths
    if (any(tree$edge.length == 0)) {
      cat("  Found", sum(tree$edge.length == 0), "zero branch lengths. Fixing...\n")
      # Replace zero branch lengths with a small fraction of the maximum branch length
      max_branch <- max(tree$edge.length, na.rm = TRUE)
      if (max_branch > 0) {
        tree$edge.length[tree$edge.length == 0] <- 1e-3 * max_branch
      } else {
        # If all branches are zero, set to a small constant
        tree$edge.length[tree$edge.length == 0] <- 1e-3
      }
    }
    
    # Additional check for very small branch lengths
    min_branch <- 1e-6 * max(node.depth.edgelength(tree))
    tree$edge.length[tree$edge.length < min_branch & tree$edge.length > 0] <- min_branch
    
    # Final validation
    if (any(is.na(tree$edge.length))) {
      tree$edge.length[is.na(tree$edge.length)] <- 1e-3
    }
    
    if (any(tree$edge.length <= 0)) {
      tree$edge.length[tree$edge.length <= 0] <- 1e-3
    }
    
    return(tree)
  }, error = function(e) {
    cat("  ERROR preparing tree:", e$message, "\n")
    return(NULL)
  })
}

# IMPROVED function to simulate trait with better error handling
simulate_trait_with_actual_params <- function(tree, trait_name, model_info, species_with_data) {
  tryCatch({
    set.seed(123) # For reproducibility
    
    model <- model_info$model
    sig2 <- model_info$param2_value
    
    # Prune tree to species that have data for this trait AND are in the tree
    trait_tree <- keep.tip(tree, species_with_data)
    
    # Additional tree validation
    if (any(trait_tree$edge.length <= 0)) {
      trait_tree$edge.length[trait_tree$edge.length <= 0] <- 1e-3
    }
    
    if (model == "BM") {
      # Brownian Motion
      sim <- fastBM(trait_tree, sig2 = sig2, internal = TRUE)
      
    } else if (model == "OU") {
      # Ornstein-Uhlenbeck
      alpha <- model_info$param1_value
      sim <- fastBM(trait_tree, model = "OU", alpha = alpha, sig2 = sig2, internal = TRUE)
      
    } else if (model == "kappa") {
      # Kappa - modify branch lengths
      kappa <- model_info$param1_value
      tree_kappa <- trait_tree
      tree_kappa$edge.length <- trait_tree$edge.length^kappa
      sim <- fastBM(tree_kappa, sig2 = sig2, internal = TRUE)
      
    } else if (model == "EB") {
      # Early Burst
      a <- model_info$param1_value
      sim <- fastBM(trait_tree, model = "EB", a = a, sig2 = sig2, internal = TRUE)
      
    } else if (model == "delta") {
      # Delta - modify branch lengths
      delta <- model_info$param1_value
      tree_delta <- trait_tree
      tree_delta$edge.length <- trait_tree$edge.length^delta
      sim <- fastBM(tree_delta, sig2 = sig2, internal = TRUE)
      
    } else if (model == "WN") {
      # White Noise - no phylogenetic structure
      trait_tips <- rnorm(length(trait_tree$tip.label), mean = 0, sd = 1)
      names(trait_tips) <- trait_tree$tip.label
      # For WN, ancestral states are meaningless, but we'll set them to mean for plotting
      trait_ancestral <- rep(mean(trait_tips), trait_tree$Nnode)
      names(trait_ancestral) <- paste0("node", 1:trait_tree$Nnode)
      return(list(tips = trait_tips, ancestral = trait_ancestral, tree = trait_tree))
    }
    
    # Extract tip and ancestral states
    trait_tips <- sim[1:length(trait_tree$tip.label)]
    trait_ancestral <- sim[(length(trait_tree$tip.label)+1):length(sim)]
    
    return(list(
      tips = trait_tips, 
      ancestral = trait_ancestral,
      tree = trait_tree
    ))
  }, error = function(e) {
    cat("    ERROR in simulation:", e$message, "\n")
    return(NULL)
  })
}

# IMPROVED function to create clean branching phenogram with error handling
make_clean_branching_phenogram <- function(tree, trait_data, trait_name, model_info, pretty_name, order_name) {
  tryCatch({
    # Extract tip and ancestral data
    trait_tips <- trait_data$tips
    trait_ancestral <- trait_data$ancestral
    
    # Combine all trait values (tips + ancestral)
    all_traits <- c(trait_tips, trait_ancestral)
    names(all_traits)[1:length(trait_tips)] <- tree$tip.label
    names(all_traits)[(length(trait_tips)+1):length(all_traits)] <- paste0("node", 1:tree$Nnode)
    
    # --- Node depths ---
    node_x <- node.depth.edgelength(tree)
    tip_x  <- node_x[1:length(tree$tip.label)]
    names(tip_x) <- tree$tip.label
    
    # --- Build edge lines ---
    edges <- tree$edge
    lines <- data.frame()
    for (i in 1:nrow(edges)) {
      a <- edges[i, 1]  # ancestor node
      d <- edges[i, 2]  # descendant node
      
      # Get trait values for this edge
      if (a <= length(tree$tip.label)) {
        y_anc <- all_traits[tree$tip.label[a]]
      } else {
        y_anc <- all_traits[paste0("node", a - length(tree$tip.label))]
      }
      
      if (d <= length(tree$tip.label)) {
        y_desc <- all_traits[tree$tip.label[d]]
      } else {
        y_desc <- all_traits[paste0("node", d - length(tree$tip.label))]
      }
      
      lines <- rbind(lines, data.frame(
        x = c(node_x[a], node_x[d]),
        y = c(y_anc, y_desc),
        group = i
      ))
    }
    
    # --- Tips ---
    tips <- data.frame(
      x = tip_x[names(trait_tips)],
      y = trait_tips
    )
    
    # --- Plot ---
    p <- ggplot() +
      # Plot branches
      geom_line(data = lines, aes(x, y, group = group), 
                color = "steelblue", size = 0.8, alpha = 0.8) +
      # Plot tip points (NO labels)
      geom_point(data = tips, aes(x, y), color = "red", size = 2) +
      # Labels and title
      labs(x = "Evolutionary Time", y = "Simulated Trait Value") +
      ggtitle(
        paste0(order_name, ": ", pretty_name, " - ", model_info$model, " Model"),
        subtitle = paste0("AIC: ", round(model_info$aic, 2), 
                          ifelse(!is.na(model_info$param2_value), 
                                 paste0(" | σ²: ", round(model_info$param2_value, 3)), ""),
                          ifelse(!is.na(model_info$param1), 
                                 paste0(" | ", model_info$param1, ": ", round(model_info$param1_value, 3)), 
                                 ""))
      ) +
      # Theme
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.margin = margin(10, 10, 10, 10),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.grid.minor = element_line(color = "grey95")
      )
    
    return(p)
  }, error = function(e) {
    cat("    ERROR creating branching plot:", e$message, "\n")
    return(NULL)
  })
}

# IMPROVED function to create simple phenogram with error handling
plot_clean_simple_phenogram <- function(tree, trait_data, trait_name, model_info, pretty_name, order_name) {
  tryCatch({
    trait_tips <- trait_data$tips
    
    # Calculate node depths (time)
    node_depths <- node.depth.edgelength(tree)
    tip_depths <- node_depths[1:length(tree$tip.label)]
    
    # Create data frame
    pheno_df <- data.frame(
      time = tip_depths,
      trait_value = trait_tips[tree$tip.label]
    )
    
    p <- ggplot(pheno_df, aes(x = time, y = trait_value)) +
      geom_point(color = "red", size = 2, alpha = 0.7) +
      geom_smooth(method = "loess", alpha = 0.3, color = "blue", se = FALSE) +
      theme_bw() +
      labs(title = paste(order_name, ":", pretty_name, "-", model_info$model, "Model"),
           subtitle = paste0("AIC: ", round(model_info$aic, 2)),
           x = "Time", y = "Simulated Trait Value") +
      theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
            plot.subtitle = element_text(hjust = 0.5, size = 10))
    
    return(p)
  }, error = function(e) {
    cat("    ERROR creating simple plot:", e$message, "\n")
    return(NULL)
  })
}

# Main processing loop for all orders
cat("Starting genomic features simulation for all orders...\n")

# Create main output directory
main_output_dir <- "D:/Akash Ajay/insect bite force project/genomic_features_simulations"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

# Check available sheets in the Excel file
cat("Checking available sheets in the Excel file...\n")
available_sheets <- excel_sheets(file_path)
cat("Available sheets:", paste(available_sheets, collapse = ", "), "\n")

# Process each order
for (order in orders) {
  cat("\n", rep("=", 50), "\n")
  cat("Processing order:", order, "\n")
  cat(rep("=", 50), "\n")
  
  # Create order-specific output directory
  order_output_dir <- file.path(main_output_dir, order)
  if (!dir.exists(order_output_dir)) {
    dir.create(order_output_dir, recursive = TRUE)
  }
  
  # Load and prepare tree with error handling
  tree_path <- tree_paths[[order]]
  if (!file.exists(tree_path)) {
    cat("✗ Tree file not found:", tree_path, "\n")
    next
  }
  
  tree <- prepare_tree(tree_path)
  if (is.null(tree)) {
    cat("✗ Failed to load/prepare tree for", order, "\n")
    next
  }
  
  cat("Tree loaded:", length(tree$tip.label), "tips\n")
  
  # Try to load data for this order from different possible sheet names
  order_data <- NULL
  possible_sheet_names <- c(
    order,
    tolower(order),
    toupper(order),
    gsub(" ", "_", order),
    gsub(" ", "", order)
  )
  
  for (sheet_name in possible_sheet_names) {
    if (sheet_name %in% available_sheets) {
      cat("Loading data from sheet:", sheet_name, "\n")
      tryCatch({
        order_data <- read_excel(file_path, sheet = sheet_name)
        order_data <- as.data.frame(order_data)
        break
      }, error = function(e) {
        cat("  Failed to load sheet:", sheet_name, "\n")
      })
    }
  }
  
  if (is.null(order_data)) {
    cat("✗ No data sheet found for", order, "\n")
    next
  }
  
  cat("Data for", order, ":", nrow(order_data), "species\n")
  
  # Create Species column - check which columns exist
  if ("Species" %in% colnames(order_data)) {
    # If Species column exists, clean it by replacing spaces with underscores
    order_data$Species <- gsub("\\s+", "_", trimws(as.character(order_data$Species)))
    cat("  Using existing Species column\n")
  } else if ("genus" %in% colnames(order_data) && "species" %in% colnames(order_data)) {
    order_data$Species <- paste(order_data$genus, order_data$species, sep = "_")
    order_data$Species <- gsub("\\s+", "_", trimws(as.character(order_data$Species)))
  } else if ("Genus" %in% colnames(order_data) && "Species" %in% colnames(order_data)) {
    order_data$Species <- paste(order_data$Genus, order_data$Species, sep = "_")
    order_data$Species <- gsub("\\s+", "_", trimws(as.character(order_data$Species)))
  } else {
    # If no genus/species columns, check if there's a species column directly
    species_col <- grep("^species$", colnames(order_data), value = TRUE, ignore.case = TRUE)
    if (length(species_col) > 0) {
      order_data$Species <- order_data[[species_col[1]]]
      order_data$Species <- gsub("\\s+", "_", trimws(as.character(order_data$Species)))
    } else {
      cat("✗ Could not find species information in the data\n")
      next
    }
  }
  
  # Additional cleaning: remove any extra whitespace and ensure proper format
  order_data$Species <- gsub("__+", "_", order_data$Species)  # Replace multiple underscores with single
  order_data$Species <- gsub("_$", "", order_data$Species)    # Remove trailing underscores
  order_data$Species <- gsub("^_", "", order_data$Species)    # Remove leading underscores
  
  rownames(order_data) <- make.unique(order_data$Species)
  
  # Also clean tree tip labels to ensure consistency
  tree$tip.label <- gsub("\\s+", "_", tree$tip.label)
  tree$tip.label <- gsub("__+", "_", tree$tip.label)
  tree$tip.label <- gsub("_$", "", tree$tip.label)
  tree$tip.label <- gsub("^_", "", tree$tip.label)
  
  # Find common species between tree and data
  common_species_all <- intersect(tree$tip.label, order_data$Species)
  cat("Common species between tree and data:", length(common_species_all), "\n")
  
  # Show some examples for debugging
  if (length(common_species_all) > 0) {
    cat("First few common species:", paste(head(common_species_all, 3), collapse = ", "), "\n")
  } else {
    cat("Tree tip examples:", paste(head(tree$tip.label, 3), collapse = ", "), "\n")
    cat("Data species examples:", paste(head(order_data$Species, 3), collapse = ", "), "\n")
  }
  
  if (length(common_species_all) < 4) {
    cat("✗ Insufficient common species for analysis\n")
    next
  }
  
  # Get model parameters for this order
  model_params <- model_params_all[[order]]
  
  # Process each trait
  branching_plots <- list()
  simple_plots <- list()
  successful_traits <- c()
  
  for (trait_code_name in names(model_params)) {
    cat("\n  Processing:", trait_code_name, "->", pretty_names[trait_code_name], "\n")
    
    # Find the original column name in the data
    original_trait_name <- names(trait_map)[trait_map == trait_code_name]
    
    if (length(original_trait_name) == 0) {
      cat("    ✗ Trait not found in mapping\n")
      next
    }
    
    # Check if the trait column exists in the data
    if (!original_trait_name %in% colnames(order_data)) {
      cat("    ✗ Trait column not found in data:", original_trait_name, "\n")
      # Try alternative column names
      alt_names <- grep(original_trait_name, colnames(order_data), value = TRUE, ignore.case = TRUE)
      if (length(alt_names) > 0) {
        cat("    Trying alternative column:", alt_names[1], "\n")
        original_trait_name <- alt_names[1]
      } else {
        next
      }
    }
    
    # Get species that have data for this trait (non-NA) AND are in the tree
    trait_data_raw <- order_data[[original_trait_name]]
    names(trait_data_raw) <- rownames(order_data)
    
    # Only use species that are BOTH in the tree AND have non-NA data
    species_with_data <- names(trait_data_raw)[!is.na(trait_data_raw)]
    species_with_data <- intersect(species_with_data, common_species_all)
    
    cat("    Species with data (and in tree):", length(species_with_data), "\n")
    
    if (length(species_with_data) < 4) {
      cat("    ✗ Insufficient data for simulation\n")
      next
    }
    
    # Get model info
    model_info <- model_params[[trait_code_name]]
    cat("    Model:", model_info$model, "| AIC:", model_info$aic, "\n")
    
    # Simulate trait evolution with error handling
    simulated_data <- simulate_trait_with_actual_params(tree, trait_code_name, model_info, species_with_data)
    
    if (is.null(simulated_data)) {
      cat("    ✗ Simulation failed for", trait_code_name, "\n")
      next
    }
    
    # Create plots with error handling
    branching_plot <- make_clean_branching_phenogram(
      simulated_data$tree, simulated_data, trait_code_name, model_info, pretty_names[trait_code_name], order
    )
    
    simple_plot <- plot_clean_simple_phenogram(
      simulated_data$tree, simulated_data, trait_code_name, model_info, pretty_names[trait_code_name], order
    )
    
    if (is.null(branching_plot) || is.null(simple_plot)) {
      cat("    ✗ Plot creation failed for", trait_code_name, "\n")
      next
    }
    
    branching_plots[[trait_code_name]] <- branching_plot
    simple_plots[[trait_code_name]] <- simple_plot
    successful_traits <- c(successful_traits, trait_code_name)
    
    cat("    ✓ Plots created successfully\n")
  }
  
  # Save individual plots for this order if we have any successful traits
  if (length(successful_traits) > 0) {
    cat("\n  Saving individual plots for", order, "...\n")
    cat("  Successful traits:", paste(successful_traits, collapse = ", "), "\n")
    
    for (trait_name in successful_traits) {
      clean_name <- gsub(" ", "_", trait_name)
      
      # Branching phenogram
      ggsave(file.path(order_output_dir, paste0("branching_", trait_name, ".png")), 
             branching_plots[[trait_name]], width = 10, height = 6, dpi = 300, bg = "white")
      
      # Simple phenogram
      ggsave(file.path(order_output_dir, paste0("simple_", trait_name, ".png")), 
             simple_plots[[trait_name]], width = 8, height = 5, dpi = 300, bg = "white")
    }
    
    # Create combined PDFs for this order
    cat("  Creating combined PDFs for", order, "...\n")
    
    # Combined branching plots
    pdf(file.path(order_output_dir, paste0(order, "_all_traits_branching.pdf")), width = 14, height = 18)
    grid.arrange(
      grobs = branching_plots,
      ncol = 2,
      top = paste(order, "Genomic Features Evolution\nSimulated with Actual Best-Fit Models")
    )
    dev.off()
    
    # Combined simple plots  
    pdf(file.path(order_output_dir, paste0(order, "_all_traits_simple.pdf")), width = 14, height = 18)
    grid.arrange(
      grobs = simple_plots,
      ncol = 2,
      top = paste(order, "Genomic Features Evolution\nSimulated with Actual Best-Fit Models")
    )
    dev.off()
    
    # Create model summary for this order
    model_summary <- data.frame(
      Trait = sapply(names(model_params), function(x) pretty_names[x]),
      Model = sapply(model_params, function(x) x$model),
      AIC = sapply(model_params, function(x) x$aic),
      sigma2 = sapply(model_params, function(x) ifelse(!is.na(x$param2_value), x$param2_value, NA)),
      Parameter1 = sapply(model_params, function(x) ifelse(!is.na(x$param1), paste0(x$param1, ": ", round(x$param1_value, 6)), "NA")),
      Species_Used = sapply(names(model_params), function(trait_code) {
        original_name <- names(trait_map)[trait_map == trait_code]
        if (length(original_name) > 0 && original_name %in% colnames(order_data)) {
          trait_data_raw <- order_data[[original_name]]
          names(trait_data_raw) <- rownames(order_data)
          species_with_data <- names(trait_data_raw)[!is.na(trait_data_raw)]
          length(intersect(species_with_data, common_species_all))
        } else {
          0
        }
      }),
      Status = ifelse(names(model_params) %in% successful_traits, "Success", "Failed")
    )
    
    write.csv(model_summary, file.path(order_output_dir, paste0(order, "_model_parameters_summary.csv")), row.names = FALSE)
    
    cat("  ✓", order, "processing complete\n")
    cat("  Successful traits:", length(successful_traits), "/", length(model_params), "\n")
    cat("  Files saved to:", order_output_dir, "\n")
  } else {
    cat("  ✗ No successful plots created for", order, "- all trait simulations failed\n")
  }
}

cat("\n", rep("=", 60), "\n")
cat("=== ALL ORDERS PROCESSING COMPLETE ===\n")
cat(rep("=", 60), "\n")
cat("Results saved to:", main_output_dir, "\n")
cat("Each order has its own subdirectory with:\n")
cat("- Individual branching and simple phenograms\n")
cat("- Combined PDFs for all traits\n")
cat("- Model parameter summary CSV\n")
cat("\nOrders processed:", paste(orders, collapse = ", "), "\n")