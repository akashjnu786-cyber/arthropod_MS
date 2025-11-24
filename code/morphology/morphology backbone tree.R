## morphology backbone tree taken from time tree
# Load libraries
library(ape)
library(phytools)
library(dplyr)
library(stringr)
library(ggtree)
library(taxize)# For taxonomic classification
library(ggplot2)

# Read your TimeTree Newick file
big_tree <- read.tree("D:/Akash Ajay/insect bite force project/Morphology/morphology_orders.NWK")

# Define your vector of orders
orders_of_interest <- c("Blattodea", "Coleoptera", "Dermaptera", "Embioptera",
                        "Hymenoptera", "Mantodea", "Mantophasmatodea", 
                        "Megaloptera", "Neuroptera", "Odonata", 
                        "Orthoptera", "Phasmatodea", "Raphidioptera")

# First, try the direct approach (in case it works for some orders)
cat("=== Attempt 1: Direct string matching ===\n")
tips_to_keep <- big_tree$tip.label[
  sapply(orders_of_interest, 
         function(order) grep(order, big_tree$tip.label, ignore.case = TRUE)) %>%
    unlist() %>%
    unique()
]

cat("Found", length(tips_to_keep), "tips using direct matching\n")

# If direct matching found very few tips, use the representative genera approach
if(length(tips_to_keep) < 20) {
  cat("\n=== Attempt 2: Using representative genera ===\n")
  
  # Create a list of representative genera for each order
  order_representatives <- list(
    Blattodea = c("Blattella", "Periplaneta", "Blatta", "Ectobius", "Nauphoeta"),
    Coleoptera = c("Tribolium", "Carabus", "Tenebrio", "Dytiscus", "Coccinella", "Melolontha"),
    Dermaptera = c("Forficula", "Labidura", "Anisolabis"),
    Embioptera = c("Oligotoma", "Embia", "Haploembia"),
    Hymenoptera = c("Apis", "Formica", "Bombus", "Vespa", "Camponotus"),
    Mantodea = c("Mantis", "Tenodera", "Sphodromantis"),
    Mantophasmatodea = c("Mantophasma", "Karoophasma"),
    Megaloptera = c("Corydalus", "Sialis", "Chauliodes"),
    Neuroptera = c("Chrysopa", "Hemerobius", "Myrmeleon"),
    Odonata = c("Aeshna", "Sympetrum", "Libellula", "Calopteryx", "Coenagrion"),
    Orthoptera = c("Gryllus", "Acheta", "Locusta", "Schistocerca", "Tettigonia"),
    Phasmatodea = c("Carausius", "Extatosoma", "Phyllium"),
    Raphidioptera = c("Raphidia", "Inocellia")
  )
  
  # Find tips that match these genera
  additional_tips <- c()
  for(order in names(order_representatives)) {
    for(genus in order_representatives[[order]]) {
      matches <- grep(genus, big_tree$tip.label, ignore.case = TRUE, value = TRUE)
      if(length(matches) > 0) {
        cat("Found", length(matches), "tips for genus", genus, "(", order, ")\n")
        additional_tips <- c(additional_tips, matches)
      }
    }
  }
  
  # Combine with previously found tips
  tips_to_keep <- unique(c(tips_to_keep, additional_tips))
  cat("Total tips found after genus matching:", length(tips_to_keep), "\n")
} # <-- This closing brace was missing

# Only proceed if we have representatives from multiple orders
if(length(unique(tree_data$order)) > 1) {
  # Get orders that actually have representatives (excluding Unknown)
  valid_orders <- names(order_counts)[order_counts > 0 & names(order_counts) != "Unknown"]
  
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
    
    # Now we need to prune this tree to keep only one representative per order
    # Create a new tree with just the order representatives
    
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
      labs(title = "Evolutionary Relationships of Insect Orders",
           subtitle = "Time-calibrated phylogeny from TimeTree") +
      theme_tree2()
    
    print(p)
    
    # Save the tree
    write.tree(final_backbone_tree, "D:/Akash Ajay/insect bite force project/Morphology/insect_orders_backbone.tre")
    cat("Backbone tree saved with", length(final_backbone_tree$tip.label), "orders\n")
    
  } else {
    cat("Not enough orders found to create a backbone tree\n")
  }
} else {
  cat("Not enough orders found to proceed\n")
}