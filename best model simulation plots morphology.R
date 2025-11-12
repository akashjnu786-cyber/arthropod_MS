library(ape)
library(geiger)
library(readxl)
library(ggplot2)
library(phytools)
library(gridExtra)
library(dplyr)

# Define paths and parameters
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
orders <- c("Blattodea", "Coleoptera", "Hymenoptera", "Mantodea", "Odonata", "Orthoptera", "Phasmatodea")
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK",
  Phasmatodea = "D:/Akash Ajay/insect bite force project/Morphology/Phasmatodea/Phasmatodea.NWK"
)

# Create model parameters from your data
model_params_all <- list(
  Blattodea = list(
    "Bite_F" = list(model = "BM", aic = 296.555, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 1105.837),
    "Head_W" = list(model = "kappa", aic = 77.755, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 0.88),
    "Head_H" = list(model = "kappa", aic = 85.655, param1 = "kappa", param1_value = 0.110595, param2 = "sig2", param2_value = 0.863),
    "Head_L" = list(model = "WN", aic = 53.595, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Thorax_W" = list(model = "kappa", aic = 119.733, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 7.176),
    "Wing_L" = list(model = "BM", aic = 156.93, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 1.027),
    "Body_L" = list(model = "kappa", aic = 166.442, param1 = "kappa", param1_value = 0.107793, param2 = "sig2", param2_value = 49.551)
  ),
  Coleoptera = list(
    "Bite_F" = list(model = "WN", aic = 597.285, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Head_W" = list(model = "kappa", aic = 194.156, param1 = "kappa", param1_value = 0, param2 = "sig2", param2_value = 3.431),
    "Head_H" = list(model = "OU", aic = 148.254, param1 = "alpha", param1_value = 0.005739, param2 = "sig2", param2_value = 0.038),
    "Head_L" = list(model = "OU", aic = 176.947, param1 = "alpha", param1_value = 0.008042, param2 = "sig2", param2_value = 0.102),
    "Thorax_W" = list(model = "kappa", aic = 205.673, param1 = "kappa", param1_value = 0.032485, param2 = "sig2", param2_value = 4.102),
    "Wing_L" = list(model = "OU", aic = 240.726, param1 = "alpha", param1_value = 0.008929, param2 = "sig2", param2_value = 0.589),
    "Body_L" = list(model = "OU", aic = 279.574, param1 = "alpha", param1_value = 0.007667, param2 = "sig2", param2_value = 1.462)
  ),
  Hymenoptera = list(
    "Bite_F" = list(model = "BM", aic = 375.285, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 1367.183),
    "Head_W" = list(model = "BM", aic = 83.313, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.012),
    "Head_H" = list(model = "BM", aic = 90.867, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.016),
    "Head_L" = list(model = "WN", aic = 57.2, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.021),
    "Thorax_W" = list(model = "BM", aic = 98.235, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.021),
    "Wing_L" = list(model = "EB", aic = 141.063, param1 = "a", param1_value = -0.016406, param2 = "sig2", param2_value = 0.255),
    "Body_L" = list(model = "BM", aic = 157.357, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 0.224)
  ),
  Mantodea = list(
    "Bite_F" = list(model = "WN", aic = 389.173, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Head_W" = list(model = "kappa", aic = 106.851, param1 = "kappa", param1_value = 0.384445, param2 = "sig2", param2_value = 0.352),
    "Head_H" = list(model = "kappa", aic = 133.154, param1 = "kappa", param1_value = 0.297852, param2 = "sig2", param2_value = 1.381),
    "Head_L" = list(model = "kappa", aic = 63.374, param1 = "kappa", param1_value = 0.07184, param2 = "sig2", param2_value = 0.187),
    "Thorax_W" = list(model = "delta", aic = 103.84, param1 = "delta", param1_value = 2.730188, param2 = "sig2", param2_value = 0.021),
    "Wing_L" = list(model = "delta", aic = 211.594, param1 = "delta", param1_value = 2.918442, param2 = "sig2", param2_value = 1.546),
    "Body_L" = list(model = "kappa", aic = 231.15, param1 = "kappa", param1_value = 0.320638, param2 = "sig2", param2_value = 64.093)
  ),
  Odonata = list(
    "Bite_F" = list(model = "WN", aic = 272.635, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Head_W" = list(model = "kappa", aic = 76.846, param1 = "kappa", param1_value = 3.52e-155, param2 = "sig2", param2_value = 1.592),
    "Head_H" = list(model = "kappa", aic = 78.367, param1 = "kappa", param1_value = 0.028319, param2 = "sig2", param2_value = 1.706),
    "Head_L" = list(model = "kappa", aic = 66.395, param1 = "kappa", param1_value = 0.080775001, param2 = "sig2", param2_value = 0.805),
    "Thorax_W" = list(model = "kappa", aic = 69.067, param1 = "kappa", param1_value = 0.084148, param2 = "sig2", param2_value = 4.673),
    "Wing_L" = list(model = "kappa", aic = 136.084, param1 = "kappa", param1_value = 4.24e-217, param2 = "sig2", param2_value = 51.914),
    "Body_L" = list(model = "kappa", aic = 146.452, param1 = "kappa", param1_value = 6.15e-203, param2 = "sig2", param2_value = 95.531)
  ),
  Orthoptera = list(
    "Bite_F" = list(model = "WN", aic = 352.529, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Head_W" = list(model = "OU", aic = 79.933, param1 = "alpha", param1_value = 0.050805, param2 = "sig2", param2_value = 0.269),
    "Head_H" = list(model = "OU", aic = 93.03, param1 = "alpha", param1_value = 0.105251, param2 = "sig2", param2_value = 0.885),
    "Head_L" = list(model = "OU", aic = 59.501, param1 = "alpha", param1_value = 0.020751, param2 = "sig2", param2_value = 0.054),
    "Thorax_W" = list(model = "OU", aic = 74.348, param1 = "alpha", param1_value = 0.04277, param2 = "sig2", param2_value = 0.182),
    "Wing_L" = list(model = "OU", aic = 176.567, param1 = "alpha", param1_value = 0.047211, param2 = "sig2", param2_value = 20.502),
    "Body_L" = list(model = "OU", aic = 157.247, param1 = "alpha", param1_value = 0.075515, param2 = "sig2", param2_value = 12.434)
  ),
  Phasmatodea = list(
    "Bite_F" = list(model = "WN", aic = 233.683, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = NA),
    "Head_W" = list(model = "BM", aic = 69.957, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 5.923),
    "Head_H" = list(model = "BM", aic = 73.1, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 5.205),
    "Head_L" = list(model = "BM", aic = 81.674, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 7.766),
    "Thorax_W" = list(model = "BM", aic = 86.162, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 7.003),
    "Wing_L" = list(model = "BM", aic = 136.244, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 13.472),
    "Body_L" = list(model = "BM", aic = 163.406, param1 = NA, param1_value = NA, param2 = "sig2", param2_value = 90.301)
  )
)

# Trait mapping
trait_map <- c(
  "iBite" = "Bite_F",
  "head.w" = "Head_W", 
  "head.h" = "Head_H",
  "head.l" = "Head_L",
  "th.w" = "Thorax_W",
  "wing.l" = "Wing_L", 
  "body.l" = "Body_L"
)

# Pretty display names
pretty_names <- c(
  "Bite_F" = "Bite Force",
  "Head_W" = "Head Width",
  "Head_H" = "Head Height", 
  "Head_L" = "Head Length",
  "Thorax_W" = "Thorax Width",
  "Wing_L" = "Wing Length",
  "Body_L" = "Body Length"
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
    
    # Final validation - FIXED SYNTAX
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
cat("Starting morphology simulation for all orders...\n")

# Create main output directory
main_output_dir <- "D:/Akash Ajay/insect bite force project/morphology_all_orders_simulations"
if (!dir.exists(main_output_dir)) {
  dir.create(main_output_dir, recursive = TRUE)
}

# Load the main data file
all_data <- read_excel(file_path)
all_data <- as.data.frame(all_data)

# Clean species names in the main data
all_data$Species <- paste(all_data$genus, all_data$species, sep = "_")
all_data$Species <- trimws(as.character(all_data$Species))
all_data$Species <- gsub("\\s+", "_", all_data$Species)
rownames(all_data) <- make.unique(all_data$Species)

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
  
  # Filter data for current order
  order_data <- all_data[all_data$order == order, ]
  cat("Data for", order, ":", nrow(order_data), "species\n")
  
  # Find common species between tree and data
  common_species_all <- intersect(tree$tip.label, order_data$Species)
  cat("Common species between tree and data:", length(common_species_all), "\n")
  
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
      top = paste(order, "Morphology Trait Evolution\nSimulated with Actual Best-Fit Models")
    )
    dev.off()
    
    # Combined simple plots  
    pdf(file.path(order_output_dir, paste0(order, "_all_traits_simple.pdf")), width = 14, height = 18)
    grid.arrange(
      grobs = simple_plots,
      ncol = 2,
      top = paste(order, "Morphology Trait Evolution\nSimulated with Actual Best-Fit Models")
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
        if (length(original_name) > 0) {
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