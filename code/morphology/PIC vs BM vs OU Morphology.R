## PIC vs BM vs OU Morphology - WITH ANDERSON-DARLING TEST
## automates plotting

library(ape)
library(geiger)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nortest)  # For Anderson-Darling test
library(purrr)

# Create directory for plots
plot_dir <- "D:/Akash Ajay/insect bite force project/Morphology/plots_AD"
dir.create(plot_dir, showWarnings = FALSE)

# Start capturing console output
sink("D:/Akash Ajay/insect bite force project/Morphology/ad_test_output.txt", type = "output")

# Define paths and parameters
file_path <- "D:/Akash Ajay/insect bite force project/Morphology/Insect_Orders_Averaged.xlsx"
tree_paths <- list(
  Blattodea = "D:/Akash Ajay/insect bite force project/Morphology/Blattodea/Blattodea.NWK",
  Coleoptera = "D:/Akash Ajay/insect bite force project/Morphology/Coleoptera/Coleoptera.NWK",
  Hymenoptera = "D:/Akash Ajay/insect bite force project/Morphology/Hymenoptera/Hymenoptera.NWK",
  Mantodea = "D:/Akash Ajay/insect bite force project/Morphology/Mantodea/Mantodea.NWK",
  Odonata = "D:/Akash Ajay/insect bite force project/Morphology/Odonata/Odonata.NWK",
  Orthoptera = "D:/Akash Ajay/insect bite force project/Morphology/Orthoptera/Orthoptera.NWK"
)
traits <- c("iBite", "head.w", "head.h", "head.l", "th.w", "wing.l", "body.l")
trait_display <- list(
  "iBite" = "Bite_F",
  "head.w" = "Head_W",
  "head.h" = "Head_H",
  "head.l" = "Head_L",
  "th.w" = "Th_W",
  "wing.l" = "Wing_L",
  "body.l" = "Body_L"
)
orders <- names(tree_paths)

# Load trees
trees <- lapply(tree_paths, function(path) {
  tree <- tryCatch(read.tree(path), error = function(e) {
    cat("Error loading tree for", path, ":", e$message, "\n")
    return(NULL)
  })
  return(tree)
})

# Initialize results data frame
ad_results <- data.frame(
  Order = character(),
  Trait = character(),
  Sample_Size = numeric(),
  Dropped_Species = character(),
  AD_Statistic_BM = numeric(),
  AD_P_Value_BM = numeric(),
  AD_Statistic_OU = numeric(),
  AD_P_Value_OU = numeric(),
  KS_Statistic_BM = numeric(),  # Keep KS for comparison
  KS_P_Value_BM = numeric(),
  KS_Statistic_OU = numeric(),
  KS_P_Value_OU = numeric(),
  stringsAsFactors = FALSE
)

# Simulation functions (same as before)
sim_pics_bm <- function(tree, nsim = 1000, sigsq) {
  sim_pics_list <- list()
  successful_sims <- 0
  max_attempts <- nsim * 2
  
  for(i in 1:max_attempts) {
    if(successful_sims >= nsim) break
    
    sim_data <- tryCatch({
      rTraitCont(tree, model = "BM", sigma = sqrt(sigsq))
    }, error = function(e) NULL)
    
    if(!is.null(sim_data) && length(sim_data) == length(tree$tip.label)) {
      pics <- tryCatch({
        pic(sim_data, tree)
      }, error = function(e) NULL)
      
      if(!is.null(pics) && length(pics) > 0 && !any(is.na(pics))) {
        sim_pics_list[[successful_sims + 1]] <- pics
        successful_sims <- successful_sims + 1
      }
    }
  }
  
  if(length(sim_pics_list) > 0) {
    unlist(sim_pics_list)
  } else {
    numeric(0)
  }
}

sim_pics_ou <- function(tree, nsim = 1000, alpha, sigma, theta) {
  sim_pics_list <- list()
  successful_sims <- 0
  max_attempts <- nsim * 2
  
  for(i in 1:max_attempts) {
    if(successful_sims >= nsim) break
    
    sim_data <- tryCatch({
      rTraitCont(tree, model = "OU", 
                 alpha = alpha, 
                 sigma = sigma,
                 theta = theta,
                 root.value = theta)
    }, error = function(e) NULL)
    
    if(!is.null(sim_data) && length(sim_data) == length(tree$tip.label)) {
      pics <- tryCatch({
        pic(sim_data, tree)
      }, error = function(e) NULL)
      
      if(!is.null(pics) && length(pics) > 0 && !any(is.na(pics))) {
        sim_pics_list[[successful_sims + 1]] <- pics
        successful_sims <- successful_sims + 1
      }
    }
  }
  
  if(length(sim_pics_list) > 0) {
    unlist(sim_pics_list)
  } else {
    numeric(0)
  }
}

# Robust test function with AD test
robust_ad_test <- function(observed_pics, tree, trait_data, model = "BM", nsim = 1000) {
  if(length(observed_pics) < 3) {
    return(list(ad_statistic = NA, ad_p.value = NA, ks_statistic = NA, ks_p.value = NA))
  }
  
  if(model == "BM") {
    # Fit BM model
    fit <- tryCatch({
      fitContinuous(tree, trait_data, model = "BM")
    }, error = function(e) NULL)
    
    if(is.null(fit)) return(list(ad_statistic = NA, ad_p.value = NA, ks_statistic = NA, ks_p.value = NA))
    
    sigsq <- fit$opt$sigsq
    
    # Simulate under BM
    null_dist <- sim_pics_bm(tree, nsim, sigsq)
    
  } else if(model == "OU") {
    # Fit OU model
    fit <- tryCatch({
      fitContinuous(tree, trait_data, model = "OU")
    }, error = function(e) NULL)
    
    if(is.null(fit)) return(list(ad_statistic = NA, ad_p.value = NA, ks_statistic = NA, ks_p.value = NA))
    
    alpha <- fit$opt$alpha
    sigsq <- fit$opt$sigsq
    theta <- fit$opt$z0
    
    # Simulate under OU
    null_dist <- sim_pics_ou(tree, nsim, alpha, sqrt(sigsq), theta)
  } else {
    return(list(ad_statistic = NA, ad_p.value = NA, ks_statistic = NA, ks_p.value = NA))
  }
  
  # Ensure we have enough simulated data
  if(length(null_dist) < 100) {
    cat("Warning: Insufficient null distribution samples (", length(null_dist), ")\n")
    return(list(ad_statistic = NA, ad_p.value = NA, ks_statistic = NA, ks_p.value = NA))
  }
  
  # Perform both AD and KS tests for comparison
  results <- list()
  
  # Anderson-Darling test
  tryCatch({
    # For AD test against normal distribution
    ad_test <- ad.test(observed_pics)
    results$ad_statistic <- ad_test$statistic
    results$ad_p.value <- ad_test$p.value
  }, error = function(e) {
    results$ad_statistic <- NA
    results$ad_p.value <- NA
  })
  
  # KS test for comparison
  tryCatch({
    ks_test <- ks.test(observed_pics, null_dist)
    results$ks_statistic <- ks_test$statistic
    results$ks_p.value <- ks_test$p.value
  }, error = function(e) {
    results$ks_statistic <- NA
    results$ks_p.value <- NA
  })
  
  return(results)
}

# Main function to compute PICs and perform AD tests
perform_ad_test <- function(order, tree, data, trait) {
  # Ensure numeric data
  data[[trait]] <- suppressWarnings(as.numeric(data[[trait]]))
  if (all(is.na(data[[trait]]) | is.na(as.numeric(data[[trait]])))) {
    cat("Warning:", order, "for", trait, ": all data non-numeric or NA\n")
    return(data.frame(
      Order = order,
      Trait = trait,
      Sample_Size = 0,
      Dropped_Species = "All",
      AD_Statistic_BM = NA,
      AD_P_Value_BM = NA,
      AD_Statistic_OU = NA,
      AD_P_Value_OU = NA,
      KS_Statistic_BM = NA,
      KS_P_Value_BM = NA,
      KS_Statistic_OU = NA,
      KS_P_Value_OU = NA
    ))
  }
  
  trait_data <- data[[trait]]
  names(trait_data) <- data$species
  trait_data <- trait_data[tree$tip.label]
  
  # Remove NAs and non-numeric data
  valid <- !is.na(trait_data) & !is.na(suppressWarnings(as.numeric(trait_data)))
  dropped_species <- tree$tip.label[!valid]
  
  if (length(dropped_species) > 0) {
    cat("Dropped species for", order, "trait", trait, ":", paste(dropped_species, collapse = ", "), "\n")
  }
  
  if (sum(valid) < 10) {
    cat("Skipping", order, "for", trait, ": insufficient non-NA data (", sum(valid), "species)\n")
    return(data.frame(
      Order = order,
      Trait = trait,
      Sample_Size = sum(valid),
      Dropped_Species = paste(dropped_species, collapse = ", "),
      AD_Statistic_BM = NA,
      AD_P_Value_BM = NA,
      AD_Statistic_OU = NA,
      AD_P_Value_OU = NA,
      KS_Statistic_BM = NA,
      KS_P_Value_BM = NA,
      KS_Statistic_OU = NA,
      KS_P_Value_OU = NA
    ))
  }
  
  # Subset tree and data
  trait_data <- trait_data[valid]
  tree_sub <- drop.tip(tree, tree$tip.label[!valid])
  
  # Compute PICs
  pic_values <- tryCatch({
    pic(trait_data, tree_sub)
  }, error = function(e) {
    cat("Error computing PICs for", order, trait, ":", e$message, "\n")
    return(numeric(0))
  })
  
  pic_values <- pic_values[!is.na(pic_values)]
  
  # Diagnostic information
  cat("PIC summary for", order, trait, ":\n")
  cat("  N PICs:", length(pic_values), "\n")
  cat("  PIC mean:", mean(pic_values), "\n")
  cat("  PIC variance:", var(pic_values), "\n")
  cat("  PIC range:", range(pic_values), "\n")
  
  # Check for extreme outliers (safe version)
  if(length(pic_values) > 1 && !is.na(sd(pic_values)) && sd(pic_values) > 0) {
    if(any(abs(pic_values) > 10 * sd(pic_values))) {
      cat("Warning: Extreme PIC values detected\n")
    }
  } else {
    cat("Warning: No variation in PIC values (all identical)\n")
  }
  
  if (length(pic_values) < 3) {
    cat("Warning: Too few PICs for analysis\n")
    return(data.frame(
      Order = order,
      Trait = trait,
      Sample_Size = length(pic_values),
      Dropped_Species = paste(dropped_species, collapse = ", "),
      AD_Statistic_BM = NA,
      AD_P_Value_BM = NA,
      AD_Statistic_OU = NA,
      AD_P_Value_OU = NA,
      KS_Statistic_BM = NA,
      KS_P_Value_BM = NA,
      KS_Statistic_OU = NA,
      KS_P_Value_OU = NA
    ))
  }
  
  # BM tests using robust function
  bm_tests <- robust_ad_test(pic_values, tree_sub, trait_data, model = "BM", nsim = 1000)
  
  # OU tests using robust function
  ou_tests <- robust_ad_test(pic_values, tree_sub, trait_data, model = "OU", nsim = 1000)
  
  # Plot PIC distribution with both normal and simulated curves
  p <- ggplot(data.frame(PIC = pic_values), aes(x = PIC)) +
    geom_histogram(aes(y = ..density..), bins = 30, fill = "lightblue", alpha = 0.7) +
    stat_function(fun = dnorm, args = list(mean = 0, sd = sd(pic_values)), 
                  color = "red", size = 1, linetype = "dashed") +
    labs(title = paste(order, ":", trait_display[[trait]], "PIC Distribution"),
         subtitle = paste("AD p(BM) =", round(bm_tests$ad_p.value, 4), 
                          "KS p(BM) =", round(bm_tests$ks_p.value, 4)),
         x = "PIC Value", y = "Density") +
    theme(
      panel.background = element_rect(fill = "white"),
      plot.background = element_rect(fill = "white"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(plot_dir, paste0(order, "_", trait, "_pic_distribution.png")), plot = p)
  
  data.frame(
    Order = order,
    Trait = trait,
    Sample_Size = length(pic_values),
    Dropped_Species = paste(dropped_species, collapse = ", "),
    AD_Statistic_BM = bm_tests$ad_statistic,
    AD_P_Value_BM = bm_tests$ad_p.value,
    AD_Statistic_OU = ou_tests$ad_statistic,
    AD_P_Value_OU = ou_tests$ad_p.value,
    KS_Statistic_BM = bm_tests$ks_statistic,
    KS_P_Value_BM = bm_tests$ks_p.value,
    KS_Statistic_OU = ou_tests$ks_statistic,
    KS_P_Value_OU = ou_tests$ks_p.value
  )
}

# [Rest of the processing code remains the same as your original script...]
# Process each order
for (order in orders) {
  cat("\nProcessing", order, "\n")
  
  # Load data
  data <- tryCatch({
    as.data.frame(read_excel(file_path, sheet = order))
  }, error = function(e) {
    cat("Error loading data for", order, ":", e$message, "\n")
    for (trait in traits) {
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = 0,
        Dropped_Species = "Data load failed",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
    next
  })
  
  # [Data processing and tree matching code remains exactly the same...]
  # Use ID column for species matching, with cleaning
  if ("ID" %in% colnames(data)) {
    data$ID <- trimws(as.character(data$ID))
    # Standardize ID: replace spaces with underscores, remove special characters
    data$ID <- gsub("[[:space:]]+", "_", data$ID)
    data$ID <- gsub("[^[:alnum:]_]", "", data$ID)
    # Remove NA, empty, or invalid ID entries
    valid_rows <- !is.na(data$ID) & data$ID != ""
    if (sum(valid_rows) == 0) {
      cat("Error:", order, ": No valid ID entries after cleaning\n")
      for (trait in traits) {
        ad_results <- rbind(ad_results, data.frame(
          Order = order,
          Trait = trait,
          Sample_Size = 0,
          Dropped_Species = "No valid ID data",
          AD_Statistic_BM = NA,
          AD_P_Value_BM = NA,
          AD_Statistic_OU = NA,
          AD_P_Value_OU = NA,
          KS_Statistic_BM = NA,
          KS_P_Value_BM = NA,
          KS_Statistic_OU = NA,
          KS_P_Value_OU = NA
        ))
      }
      next
    }
    data <- data[valid_rows, ]
    data$species_lower <- tolower(data$ID)
    data$species <- data$ID  # Use ID as species for final use
  } else {
    cat("Error:", order, ": Missing 'ID' column in data\n")
    for (trait in traits) {
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = 0,
        Dropped_Species = "No valid ID data",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
    next
  }
  
  # Get tree
  tree <- trees[[order]]
  if (is.null(tree)) {
    cat("Error:", order, ": Tree is NULL or failed to load\n")
    for (trait in traits) {
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = 0,
        Dropped_Species = "Tree load failed",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
    next
  }
  tree <- multi2di(tree)  # Resolve polytomies
  tree$node.label <- NULL
  # Standardize tree tip labels: replace spaces with underscores, remove special characters
  tree$tip.label <- gsub("[[:space:]]+", "_", tree$tip.label)
  tree$tip.label <- gsub("[^[:alnum:]_]", "", tree$tip.label)
  tree_tip_lower <- tolower(tree$tip.label)
  
  # Match species (case-insensitive for finding matches, but use original ID names)
  data_tip_lower <- data$species_lower
  match_indices <- match(data_tip_lower, tree_tip_lower)
  common_indices_data <- !is.na(match_indices)
  common_species <- data$species[common_indices_data]
  
  # Ensure tree tips are subset to common species
  common_indices_tree <- tree$tip.label %in% common_species
  common_species <- intersect(common_species, tree$tip.label[common_indices_tree])
  
  cat("Data species sample (from ID):", paste(head(data$species, 5), collapse = ", "), "\n")
  cat("Tree tip sample:", paste(head(tree$tip.label, 5), collapse = ", "), "\n")
  cat("Number of exact matches:", length(common_species), "\n")
  
  if (length(common_species) < 10) {
    cat("Skipping", order, ": insufficient common species (", length(common_species), "). Potential mismatches detected.\n")
    cat("Non-matching data species:", paste(setdiff(data$species, tree$tip.label)[1:min(5, length(setdiff(data$species, tree$tip.label)))], collapse = ", "), "\n")
    cat("Non-matching tree tips:", paste(setdiff(tree$tip.label, data$species)[1:min(5, length(setdiff(tree$tip.label, data$species)))], collapse = ", "), "\n")
    for (trait in traits) {
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = length(common_species),
        Dropped_Species = "Insufficient matches",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
    next
  }
  
  # Subset data and tree to common species
  data <- data[common_indices_data, ]
  tree <- keep.tip(tree, common_species)
  
  # Verify alignment
  if (nrow(data) != length(tree$tip.label) || !all(data$species %in% tree$tip.label)) {
    cat("Error:", order, ": Mismatch in data and tree species after subsetting\n")
    for (trait in traits) {
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = 0,
        Dropped_Species = "Mismatch in data and tree species",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
    next
  }
  
  cat("Proceeding with", nrow(data), "matched species.\n")
  
  # Test each trait
  for (trait in traits) {
    if (trait %in% colnames(data)) {
      result <- perform_ad_test(order, tree, data, trait)
      ad_results <- rbind(ad_results, result)
    } else {
      cat("Warning:", order, "for", trait, ": trait column missing in data\n")
      ad_results <- rbind(ad_results, data.frame(
        Order = order,
        Trait = trait,
        Sample_Size = 0,
        Dropped_Species = "Trait column missing",
        AD_Statistic_BM = NA,
        AD_P_Value_BM = NA,
        AD_Statistic_OU = NA,
        AD_P_Value_OU = NA,
        KS_Statistic_BM = NA,
        KS_P_Value_BM = NA,
        KS_Statistic_OU = NA,
        KS_P_Value_OU = NA
      ))
    }
  }
}

# Create pivot tables for AD and KS P-values
pivot_table_ad_bm <- ad_results %>%
  select(Order, Trait, AD_P_Value_BM) %>%
  pivot_wider(names_from = Trait, values_from = AD_P_Value_BM, names_prefix = "AD_BM_")

pivot_table_ad_ou <- ad_results %>%
  select(Order, Trait, AD_P_Value_OU) %>%
  pivot_wider(names_from = Trait, values_from = AD_P_Value_OU, names_prefix = "AD_OU_")

pivot_table_ks_bm <- ad_results %>%
  select(Order, Trait, KS_P_Value_BM) %>%
  pivot_wider(names_from = Trait, values_from = KS_P_Value_BM, names_prefix = "KS_BM_")

pivot_table_ks_ou <- ad_results %>%
  select(Order, Trait, KS_P_Value_OU) %>%
  pivot_wider(names_from = Trait, values_from = KS_P_Value_OU, names_prefix = "KS_OU_")

# Combine pivot tables
pivot_table <- reduce(list(pivot_table_ad_bm, pivot_table_ad_ou, pivot_table_ks_bm, pivot_table_ks_ou), 
                      full_join, by = "Order")

# Format P-values
ad_results$AD_P_Value_BM <- format.pval(ad_results$AD_P_Value_BM, digits = 3, eps = 1e-100)
ad_results$AD_P_Value_OU <- format.pval(ad_results$AD_P_Value_OU, digits = 3, eps = 1e-100)
ad_results$KS_P_Value_BM <- format.pval(ad_results$KS_P_Value_BM, digits = 3, eps = 1e-100)
ad_results$KS_P_Value_OU <- format.pval(ad_results$KS_P_Value_OU, digits = 3, eps = 1e-100)

# Save results
write.xlsx(list(
  Raw_Results = ad_results,
  Pivot_Table = pivot_table
), "D:/Akash Ajay/insect bite force project/Morphology/ad_test_results.xlsx")

cat("\nResults saved to: D:/Akash Ajay/insect bite force project/Morphology/ad_test_results.xlsx\n")

# Print summary comparing AD vs KS
cat("\n===== AD vs KS Test Comparison =====\n")
for (order in orders) {
  order_data <- ad_results %>% filter(Order == order)
  if (nrow(order_data) > 0) {
    cat("\n", order, ":\n")
    cat("  AD significant BM deviations (P < 0.05):", sum(as.numeric(order_data$AD_P_Value_BM) < 0.05, na.rm = TRUE), "\n")
    cat("  KS significant BM deviations (P < 0.05):", sum(as.numeric(order_data$KS_P_Value_BM) < 0.05, na.rm = TRUE), "\n")
    cat("  AD significant OU deviations (P < 0.05):", sum(as.numeric(order_data$AD_P_Value_OU) < 0.05, na.rm = TRUE), "\n")
    cat("  KS significant OU deviations (P < 0.05):", sum(as.numeric(order_data$KS_P_Value_OU) < 0.05, na.rm = TRUE), "\n")
    
    # Compare sensitivity
    ad_more_sensitive <- sum(as.numeric(order_data$AD_P_Value_BM) < as.numeric(order_data$KS_P_Value_BM), na.rm = TRUE)
    cat("  AD more sensitive than KS in", ad_more_sensitive, "traits\n")
  }
}

sink()
cat("Console output saved to: D:/Akash Ajay/insect bite force project/Morphology/ad_test_output.txt\n")