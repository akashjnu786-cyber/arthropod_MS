## genome order raw tree taken from time tree
# Load libraries
library(ape)
library(phytools)
library(dplyr)
library(stringr)
library(ggtree)
library(taxize)# For taxonomic classification
library(ggplot2)

# Read your TimeTree Newick file
big_tree <- read.tree("D:/Akash Ajay/insect bite force project/genome features/genome_orders.NWK")

# Diagnostic: Show first 20 tip labels to understand the format
cat("=== Tree tip label format (first 20) ===\n")
print(head(big_tree$tip.label, 20))
cat("\nTotal tips in tree:", length(big_tree$tip.label), "\n\n")

# Define your vector of orders
orders_of_interest <- c("Archaeognatha", "Blattodea", "Coleoptera", "Dermaptera", "Dictyoptera",
                        "Diptera", "Ephemeroptera", "Lepidoptera", "Mantodea", "Mecoptera",
                        "Megaloptera", "Neolepidoptera", "Neuroptera", "Odonata", 
                        "Orthoptera", "Phasmatodea", "Plecoptera", "Siphonaptera",
                        "Strepsiptera", "Thysanoptera", "Trichoptera", "Zygentoma")

# First, try the direct approach (in case it works for some orders)
cat("=== Attempt 1: Direct string matching ===\n")
tips_to_keep <- big_tree$tip.label[
  sapply(orders_of_interest, 
         function(order) grep(order, big_tree$tip.label, ignore.case = TRUE)) %>%
    unlist() %>%
    unique()
]

cat("Found", length(tips_to_keep), "tips using direct matching\n")
if(length(tips_to_keep) > 0) {
  cat("Direct matches found:", paste(head(tips_to_keep, 5), collapse = ", "), "\n")
}

# Enhanced genus matching with more comprehensive genera list
cat("\n=== Attempt 2: Enhanced genus matching ===\n")

# More comprehensive representative genera for each order
order_representatives <- list(
  Archaeognatha = c("Machilis", "Pedetontus", "Meinertellus"),
  Blattodea = c("Blattella", "Periplaneta", "Blatta", "Ectobius", "Nauphoeta", "Coptotermes", "Loboptera", 
                "Supella", "Symploce", "Cryptocercus"),
  Coleoptera = c("Tribolium", "Carabus", "Tenebrio", "Dytiscus", "Coccinella", "Melolontha", "Leptinotarsa", 
                 "Coproporus", "Ophraella", "Anthonomus", "Phyllotreta", "Harmonia", "Pyrrhocoris"),
  Dermaptera = c("Forficula", "Labidura", "Anisolabis", "Proreus", "Chelidura"),
  Dictyoptera = c("Blattella", "Periplaneta", "Mantis", "Tenodera", "Blatta", "Cryptocercus"),
  Diptera = c("Drosophila", "Musca", "Anopheles", "Aedes", "Glossina", "Calliphora", "Lucilia", 
              "Sarcophaga", "Phlebotomus", "Culex"),
  Ephemeroptera = c("Ephemera", "Baetis", "Cloeon", "Caenis", "Hexagenia", "Paraleptophlebia"),
  Lepidoptera = c("Bombyx", "Manduca", "Helicoverpa", "Papilio", "Danaus", "Spodoptera", "Plutella", 
                  "Chilo", "Ostrinia", "Trichoplusia"),
  Mantodea = c("Mantis", "Tenodera", "Sphodromantis", "Iris", "Deroplatys"),
  Mecoptera = c("Panorpa", "Bittacus", "Boreus"),
  Megaloptera = c("Corydalus", "Sialis", "Chauliodes", "Protochauliodes"),
  Neolepidoptera = c("Bombyx", "Manduca", "Helicoverpa", "Spodoptera", "Plutella"),
  Neuroptera = c("Chrysopa", "Hemerobius", "Myrmeleon", "Cunctochrysa", "Nineta"),
  Odonata = c("Aeshna", "Sympetrum", "Libellula", "Calopteryx", "Coenagrion", "Enallagma", "Ischnura"),
  Orthoptera = c("Gryllus", "Acheta", "Locusta", "Schistocerca", "Tettigonia", "Romalea", "Dissosteira"),
  Phasmatodea = c("Carausius", "Extatosoma", "Phyllium", "Timema", "Diapheromera"),
  Plecoptera = c("Pteronarcys", "Perlodes", "Leuctra", "Nemoura", "Protonemura", "Isoperla"),
  Siphonaptera = c("Ctenocephalides", "Pulex", "Xenopsylla", "Tunga"),
  Strepsiptera = c("Xenos", "Stylops", "Halictophagus"),
  Thysanoptera = c("Frankliniella", "Thrips", "Scirtothrips", "Mycetothrips", "Haplothrips"),
  Trichoptera = c("Hydropsyche", "Rhyacophila", "Pycnopsyche", "Limnephilus", "Phryganea"),
  Zygentoma = c("Lepisma", "Thermobia", "Ctenolepisma", "Tricholepidion")
)

# Find tips that match these genera - more flexible matching
additional_tips <- c()
order_matches <- list()

for(order in names(order_representatives)) {
  order_matches[[order]] <- c()
  for(genus in order_representatives[[order]]) {
    # More flexible pattern matching
    pattern <- paste0("\\b", genus, "\\b")  # Word boundary matching
    matches <- grep(pattern, big_tree$tip.label, ignore.case = TRUE, value = TRUE)
    
    if(length(matches) == 0) {
      # Try substring matching if word boundary fails
      matches <- grep(genus, big_tree$tip.label, ignore.case = TRUE, value = TRUE)
    }
    
    if(length(matches) > 0) {
      cat("Found", length(matches), "tips for genus", genus, "(", order, "):", 
          paste(head(matches, 2), collapse = ", "), "\n")
      order_matches[[order]] <- c(order_matches[[order]], matches)
      additional_tips <- c(additional_tips, matches)
    }
  }
  if(length(order_matches[[order]]) > 0) {
    cat("Total for", order, ":", length(unique(order_matches[[order]])), "tips\n")
  }
}

# Combine with previously found tips
tips_to_keep <- unique(c(tips_to_keep, additional_tips))
cat("\nTotal unique tips found after enhanced matching:", length(tips_to_keep), "\n")

# Show which orders we found representatives for
found_orders <- sapply(order_matches, function(x) length(unique(x)) > 0)
cat("\nOrders with representatives found:\n")
cat(paste(names(found_orders)[found_orders], collapse = ", "), "\n")

# If still too few, try fuzzy matching
if(length(tips_to_keep) < 10) {
  cat("\n=== Attempt 3: Fuzzy matching with common species ===\n")
  
  # Try matching against known genomic species
  common_genomic_species <- c(
    # Drosophila
    "Drosophila_melanogaster", "Drosophila_simulans", "Drosophila_yakuba", "Drosophila_erecta",
    # Mosquitoes
    "Anopheles_gambiae", "Aedes_aegypti", "Culex_quinquefasciatus",
    # Beetles
    "Tribolium_castaneum", "Leptinotarsa_decemlineata",
    # Honey bee
    "Apis_mellifera",
    # Silkworm
    "Bombyx_mori",
    # Pea aphid (Hemiptera, but might be in tree)
    "Acyrthosiphon_pisum",
    # Red flour beetle
    "Tribolium_castaneum",
    # Pea aphid
    "Acyrthosiphon_pisum",
    # Body louse
    "Pediculus_humanus",
    # Monarch butterfly
    "Danaus_plexippus",
    # Flour beetle
    "Tribolium_castaneum",
    # Water flea (might be outgroup)
    "Daphnia_pulex"
  )
  
  fuzzy_matches <- c()
  for(species in common_genomic_species) {
    matches <- grep(species, big_tree$tip.label, ignore.case = TRUE, value = TRUE)
    if(length(matches) > 0) {
      cat("Found fuzzy match:", matches[1], "\n")
      fuzzy_matches <- c(fuzzy_matches, matches[1])
    }
  }
  
  tips_to_keep <- unique(c(tips_to_keep, fuzzy_matches))
  cat("Total tips after fuzzy matching:", length(tips_to_keep), "\n")
}

# Only proceed if we have any tips
if(length(tips_to_keep) > 0) {
  # Prune the tree to keep only the matched tips
  pruned_tree <- keep.tip(big_tree, tips_to_keep)
  cat("\nPruned tree has", length(pruned_tree$tip.label), "tips\n")
  
  # Classify tips into orders using the order_matches list first
  tree_data <- data.frame(
    tip_label = pruned_tree$tip.label,
    order = "Unknown"
  )
  
  # Assign orders based on our matching
  for(order in names(order_matches)) {
    if(length(order_matches[[order]]) > 0) {
      order_tips <- intersect(pruned_tree$tip.label, order_matches[[order]])
      if(length(order_tips) > 0) {
        tree_data$order[tree_data$tip_label %in% order_tips] <- order
      }
    }
  }
  
  # For remaining unclassified tips, try the original classification
  for(i in seq_along(orders_of_interest)) {
    order <- orders_of_interest[i]
    if(!order %in% names(order_matches) || sum(tree_data$order == order) == 0) {
      matches <- grep(order, pruned_tree$tip.label, ignore.case = TRUE)
      if(length(matches) > 0) {
        tree_data$order[matches] <- order
      }
    }
  }
  
  # Compute order counts
  order_counts <- table(tree_data$order)
  cat("\nOrder representation:\n")
  print(order_counts)
  
  # Only proceed if we have representatives from multiple orders
  if(length(unique(tree_data$order[tree_data$order != "Unknown"])) > 1) {
    # Get orders that actually have representatives (excluding Unknown)
    valid_orders <- names(order_counts)[order_counts > 0 & names(order_counts) != "Unknown"]
    
    cat("\nValid orders for backbone:", paste(valid_orders, collapse = ", "), "\n")
    
    # Create a list to store the MRCA node for each order
    order_mrca_list <- list()
    
    for(order_name in valid_orders) {
      tips_in_order <- tree_data$tip_label[tree_data$order == order_name]
      if(length(tips_in_order) > 1) {
        # For orders with multiple species, get their MRCA
        order_mrca_list[[order_name]] <- getMRCA(pruned_tree, tips_in_order)
      } else if(length(tips_in_order) == 1) {
        # For orders with single species, use that tip directly
        order_mrca_list[[order_name]] <- which(pruned_tree$tip.label == tips_in_order)
      }
    }
    
    # Create a vector of nodes to keep
    nodes_to_keep <- unlist(order_mrca_list)
    
    if(length(nodes_to_keep) > 1) {
      # Get the MRCA of all the order nodes
      root_node <- getMRCA(pruned_tree, nodes_to_keep)
      
      # Extract the clade containing all order nodes
      backbone_tree <- extract.clade(pruned_tree, root_node)
      
      # For each order, select one representative tip
      representative_tips <- c()
      for(order_name in valid_orders) {
        tips_in_order <- tree_data$tip_label[tree_data$order == order_name]
        # Take the first species as representative
        representative_tips <- c(representative_tips, tips_in_order[1])
      }
      
      # Create the final backbone tree with one tip per order
      final_backbone_tree <- keep.tip(backbone_tree, representative_tips)
      
      # Rename the tips to order names
      tip_mapping <- setNames(valid_orders, representative_tips)
      final_backbone_tree$tip.label <- tip_mapping[final_backbone_tree$tip.label]
      
      # Plot the clean backbone tree
      p <- ggtree(final_backbone_tree) + 
        geom_tiplab(align = TRUE, size = 5, fontface = "bold") +
        geom_nodelab(size = 3, hjust = -0.2) +
        xlim(0, max(final_backbone_tree$edge.length) * 1.2) +
        labs(title = "Evolutionary Relationships of Insect Orders (Genomic Tree)",
             subtitle = "Time-calibrated phylogeny from TimeTree") +
        theme_tree2()
      
      print(p)
      
      # Save the tree
      write.tree(final_backbone_tree, "D:/Akash Ajay/insect bite force project/genome features/genome_orders_backbone.tre")
      cat("Backbone tree saved with", length(final_backbone_tree$tip.label), "orders\n")
      
    } else {
      cat("Not enough orders found to create a backbone tree\n")
    }
  } else {
    cat("Not enough orders found to proceed (only 1 or 0 orders with representatives)\n")
    cat("Available tips for manual inspection:\n")
    print(head(pruned_tree$tip.label, 10))
  }
} else {
  cat("No tips found - check if the tree file path is correct and contains insect data\n")
}