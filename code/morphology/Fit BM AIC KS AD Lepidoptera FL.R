#####################################################################
##  GENOMIC-TRAIT PIC DISTRIBUTION ANALYSIS (BM vs OU)             ##
##  – FULL CONSOLE OUTPUT FOR LOGGING                              ##
#####################################################################

library(ape)
library(geiger)
library(readxl)
library(openxlsx)
library(dplyr)
library(tidyr)
library(ggplot2)
library(nortest)
library(purrr)

## -------------------------------------------------
##  GLOBAL SETTINGS
## -------------------------------------------------
MIN_PICS_FOR_AD   <- 8L
MIN_PICS_FOR_PLOT <- 30L
MIN_SPECIES       <- 10L

## -------------------------------------------------
##  PATHS & TRAITS
## -------------------------------------------------
plot_dir <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/distribution_plots8"
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

sink("D:/Akash Ajay/insect bite force project/sravya work/GC skew/ks_ad_test_output8.txt",
     type = "output", split = TRUE)  # split = TRUE → also prints to R console

file_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/genomic_features_with_gc_skew1.xlsx"

tree_paths <- list(
  Lepidoptera = "D:/Akash Ajay/insect bite force project/genome features/Lepidoptera/Lepidoptera.NWK",
  Trichoptera = "D:/Akash Ajay/insect bite force project/genome features/Trichoptera/Trichoptera.NWK",
  Orthoptera  = "D:/Akash Ajay/insect bite force project/genome features/Orthoptera/Orthoptera.NWK"
)

traits <- c("Genome_Size", "Genome_GC", "GC_skew", "Chromosome_Number", "Coding_genes", "Non_coding_genes")
trait_display <- list(
  Genome_Size       = "Genome Size (Mb)",
  Genome_GC         = "Genome GC %",
  GC_skew           = "GC Skew",
  Chromosome_Number = "Chromosome Number",
  Coding_genes      = "Coding Genes",
  Non_coding_genes  = "Non-coding Genes"
)

orders <- names(tree_paths)

cat("=== STARTING ANALYSIS ===\n")
cat("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Working directory:", getwd(), "\n")
cat("Excel file:", file_path, "\n")
cat("Plot directory:", plot_dir, "\n\n")

## -------------------------------------------------
##  LOAD TREES
## -------------------------------------------------
make_ultrametric <- function(tree) {
  if (is.null(tree)) return(NULL)
  tree <- multi2di(tree)
  tree$node.label <- NULL
  zero <- which(tree$edge.length <= 0 | is.na(tree$edge.length))
  if (length(zero)) tree$edge.length[zero] <- 1e-6
  tree <- tryCatch(chronos(tree, lambda = 0, quiet = TRUE), error = function(e) tree)
  tree
}

trees <- list()
for (ord in orders) {
  p <- tree_paths[[ord]]
  if (!file.exists(p)) {
    cat("ERROR: Tree file NOT FOUND →", p, "\n")
    next
  }
  tr <- tryCatch(read.tree(p), error = function(e) NULL)
  if (is.null(tr)) {
    cat("ERROR: Failed to read tree →", p, "\n")
    next
  }
  tr <- make_ultrametric(tr)
  trees[[ord]] <- tr
  cat("Loaded tree for", ord, ":", length(tr$tip.label), "tips\n")
}

## -------------------------------------------------
##  RESULTS CONTAINER
## -------------------------------------------------
test_results <- data.frame(
  Order = character(), Trait = character(), Sample_Size = numeric(),
  AD_Statistic_BM = numeric(), AD_P_Value_BM = numeric(),
  AD_Statistic_OU = numeric(), AD_P_Value_OU = numeric(),
  KS_Statistic_BM = numeric(), KS_P_Value_BM = numeric(),
  KS_Statistic_OU = numeric(), KS_P_Value_OU = numeric(),
  stringsAsFactors = FALSE
)

## -------------------------------------------------
##  SIMULATION HELPERS
## -------------------------------------------------
sim_pics_bm <- function(tree, nsim = 1000, sigsq) {
  out <- numeric(0); ok <- 0
  for (i in seq_len(nsim * 2)) {
    if (ok >= nsim) break
    dat <- tryCatch(rTraitCont(tree, model = "BM", sigma = sqrt(sigsq)), error = function(e) NULL)
    if (is.null(dat) || length(dat) != length(tree$tip.label)) next
    pics <- tryCatch(pic(dat, tree), error = function(e) NULL)
    if (!is.null(pics) && !any(is.na(pics)) && var(pics) > 1e-10) {
      out <- c(out, pics); ok <- ok + 1
    }
  }
  out
}

sim_pics_ou <- function(tree, nsim = 1000, alpha, sigsq, theta) {
  out <- numeric(0); ok <- 0
  for (i in seq_len(nsim * 2)) {
    if (ok >= nsim) break
    dat <- tryCatch(rTraitCont(tree, model = "OU", alpha = alpha, sigma = sqrt(sigsq), theta = theta),
                    error = function(e) NULL)
    if (is.null(dat) || length(dat) != length(tree$tip.label)) next
    pics <- tryCatch(pic(dat, tree), error = function(e) NULL)
    if (!is.null(pics) && !any(is.na(pics)) && var(pics) > 1e-10) {
      out <- c(out, pics); ok <- ok + 1
    }
  }
  out
}

## -------------------------------------------------
##  ROBUST TEST
## -------------------------------------------------
robust_distribution_test <- function(obs_pics, tree, trait_vec, nsim = 1000) {
  n <- length(obs_pics)
  if (n < 3 || var(obs_pics) < 1e-10) {
    cat("   WARNING: Too few or zero-variance PICs (n =", n, ", var =", var(obs_pics), ")\n")
    return(list(bm_dist = numeric(0), ou_dist = numeric(0),
                ad_statistic_bm = NA, ad_p.value_bm = NA,
                ad_statistic_ou = NA, ad_p.value_ou = NA,
                ks_statistic_bm = NA, ks_p.value_bm = NA,
                ks_statistic_ou = NA, ks_p.value_ou = NA))
  }
  
  cat("   Fitting BM...\n")
  bm_fit <- tryCatch(fitContinuous(tree, trait_vec, model = "BM"), error = function(e) NULL)
  cat("   Fitting OU...\n")
  ou_fit <- tryCatch(fitContinuous(tree, trait_vec, model = "OU"), error = function(e) NULL)
  
  bm_dist <- ou_dist <- numeric(0)
  
  if (!is.null(bm_fit) && !is.na(bm_fit$opt$sigsq) && bm_fit$opt$sigsq > 0) {
    cat("   BM fit: sigsq =", round(bm_fit$opt$sigsq, 6), "\n")
    bm_dist <- sim_pics_bm(tree, nsim, bm_fit$opt$sigsq)
    cat("   Simulated", length(bm_dist), "BM PICs\n")
  } else {
    cat("   BM fit failed → using fallback\n")
  }
  
  if (!is.null(ou_fit) && !is.na(ou_fit$opt$alpha) && ou_fit$opt$alpha > 1e-6 && ou_fit$opt$sigsq > 0) {
    cat("   OU fit: alpha =", round(ou_fit$opt$alpha, 6),
        ", sigsq =", round(ou_fit$opt$sigsq, 6),
        ", theta =", round(ou_fit$opt$z0, 6), "\n")
    ou_dist <- sim_pics_ou(tree, nsim, ou_fit$opt$alpha, ou_fit$opt$sigsq, ou_fit$opt$z0)
    cat("   Simulated", length(ou_dist), "OU PICs\n")
  } else {
    cat("   OU fit failed → using weak OU fallback\n")
    ou_dist <- sim_pics_ou(tree, nsim, alpha = 0.1, sigsq = var(trait_vec) * 0.1, theta = mean(trait_vec))
  }
  
  std <- function(x) { if (length(x) == 0 || var(x) == 0) x else scale(x)[,1] }
  bm_dist_std <- std(bm_dist); ou_dist_std <- std(ou_dist)
  
  ad_res <- if (n >= MIN_PICS_FOR_AD) tryCatch(ad.test(obs_pics), error = function(e) list(statistic=NA, p.value=NA)) 
  else list(statistic=NA, p.value=NA)
  cat("   AD test on observed PICs: stat =", round(ad_res$statistic, 4), ", p =", round(ad_res$p.value, 4), "\n")
  
  ks_bm <- if (length(bm_dist_std) > 10) tryCatch(ks.test(obs_pics, bm_dist_std), error = function(e) list(statistic=NA, p.value=NA))
  else list(statistic=NA, p.value=NA)
  ks_ou <- if (length(ou_dist_std) > 10) tryCatch(ks.test(obs_pics, ou_dist_std), error = function(e) list(statistic=NA, p.value=NA))
  else list(statistic=NA, p.value=NA)
  
  cat("   KS BM: stat =", round(ks_bm$statistic, 4), ", p =", round(ks_bm$p.value, 4), "\n")
  cat("   KS OU: stat =", round(ks_ou$statistic, 4), ", p =", round(ks_ou$p.value, 4), "\n")
  
  list(
    bm_dist = bm_dist_std, ou_dist = ou_dist_std,
    ad_statistic_bm = ad_res$statistic, ad_p.value_bm = ad_res$p.value,
    ad_statistic_ou = ad_res$statistic, ad_p.value_ou = ad_res$p.value,
    ks_statistic_bm = ks_bm$statistic, ks_p.value_bm = ks_bm$p.value,
    ks_statistic_ou = ks_ou$statistic, ks_p.value_ou = ks_ou$p.value
  )
}

## -------------------------------------------------
##  PLOTTING (fixed ggplot2)
## -------------------------------------------------
create_distribution_plot <- function(order, trait, obs_pics, test_res, bm_dist, ou_dist) {
  if (length(obs_pics) < 3) return(NULL)
  
  df <- data.frame(
    Value = c(obs_pics, bm_dist, ou_dist),
    Type  = factor(c(rep("Observed PIC", length(obs_pics)),
                     rep("BM Simulated", length(bm_dist)),
                     rep("OU Simulated", length(ou_dist))),
                   levels = c("Observed PIC", "BM Simulated", "OU Simulated"))
  )
  
  p <- ggplot(df, aes(x = Value, fill = Type, colour = Type)) +
    geom_histogram(data = subset(df, Type == "Observed PIC"),
                   aes(y = after_stat(density)), bins = 30, alpha = 0.5) +
    geom_density(data = subset(df, Type == "Observed PIC"), linewidth = 1.2) +
    geom_density(data = subset(df, Type == "BM Simulated"), linewidth = 1, linetype = "dashed") +
    geom_density(data = subset(df, Type == "OU Simulated"), linewidth = 1, linetype = "dotted") +
    scale_fill_manual(values = c("Observed PIC" = "lightblue", "BM Simulated" = "red", "OU Simulated" = "darkgreen")) +
    scale_color_manual(values = c("Observed PIC" = "blue", "BM Simulated" = "red", "OU Simulated" = "darkgreen")) +
    labs(title = paste(order, ":", trait_display[[trait]]),
         subtitle = paste("AD p =", sprintf("%.3f", test_res$ad_p.value_bm),
                          "| KS BM p =", sprintf("%.3f", test_res$ks_p.value_bm),
                          "\nOU AD p =", sprintf("%.3f", test_res$ad_p.value_ou),
                          "| KS OU p =", sprintf("%.3f", test_res$ks_p.value_ou)),
         x = "Standardized PIC", y = "Density") +
    annotate("text", x = Inf, y = Inf, label = paste("n =", length(obs_pics)), hjust = 1.1, vjust = 1.1, size = 4) +
    theme_minimal(base_size = 12) +
    theme(legend.position = "top", panel.grid = element_blank(),
          panel.background = element_rect(fill = "white"),
          plot.background = element_rect(fill = "white"),
          panel.border = element_rect(colour = "black", fill = NA))
  
  # SAFE FILENAME: remove %, spaces, special chars
  safe_name <- gsub("%", "pct", trait_display[[trait]])  # % → pct
  safe_name <- gsub("[^[:alnum:]_]", "_", safe_name)     # only letters, numbers, _
  safe_name <- gsub("__+", "_", safe_name)               # collapse multiple _
  safe_name <- sub("_$", "", safe_name)                  # remove trailing _
  
  fn <- file.path(plot_dir, paste0(order, "_", safe_name, "_distribution.png"))
  
  tryCatch({
    ggsave(fn, p, width = 10, height = 8, dpi = 300, device = "png")
    cat("   PLOT SAVED →", fn, "\n")
  }, error = function(e) {
    cat("   FAILED TO SAVE PLOT:", e$message, "\n")
  })
  
  invisible(p)
}
## -------------------------------------------------
##  MAIN ANALYSIS
## -------------------------------------------------
perform_distribution_analysis <- function(order, tree, data, trait) {
  cat("\n--- ANALYZING TRAIT:", trait, "---\n")
  data[[trait]] <- suppressWarnings(as.numeric(data[[trait]]))
  valid_vals <- data[[trait]][!is.na(data[[trait]])]
  if (length(valid_vals) == 0) {
    cat("   All values NA → skip\n")
    return(NULL)
  }
  
  cat("   Non-NA values:", length(valid_vals), "\n")
  cat("   Min:", min(valid_vals), "| Max:", max(valid_vals), "| Mean:", mean(valid_vals), "| Var:", var(valid_vals), "\n")
  
  # Print first 5 values
  cat("   First 5 values:\n")
  print(head(data[!is.na(data[[trait]]), c("Species", trait)], 5))
  
  trait_vec <- setNames(data[[trait]], data$Species)[tree$tip.label]
  valid <- !is.na(trait_vec)
  n_valid <- sum(valid)
  cat("   Species with data on tree:", n_valid, "\n")
  
  if (n_valid < MIN_SPECIES) {
    cat("   Too few species (need >=", MIN_SPECIES, ") → skip\n")
    return(NULL)
  }
  
  trait_vec <- trait_vec[valid]
  tree_sub <- keep.tip(tree, names(trait_vec))
  
  if (var(trait_vec) < 1e-10) {
    cat("   Trait is constant (var < 1e-10) → skip\n")
    return(NULL)
  }
  
  pics <- tryCatch(pic(trait_vec, tree_sub), error = function(e) numeric(0))
  pics <- pics[!is.na(pics)]
  if (length(pics) < 3) {
    cat("   <3 PICs → skip\n")
    return(NULL)
  }
  
  cat("   PICs computed:", length(pics), "\n")
  cat("   First 10 PICs:", paste(head(pics, 10), collapse = ", "), "\n")
  
  pics_std <- scale(pics)[,1]
  test_out <- robust_distribution_test(pics_std, tree_sub, trait_vec)
  
  row <- data.frame(
    Order = order, Trait = trait, Sample_Size = length(pics),
    AD_Statistic_BM = test_out$ad_statistic_bm, AD_P_Value_BM = test_out$ad_p.value_bm,
    AD_Statistic_OU = test_out$ad_statistic_ou, AD_P_Value_OU = test_out$ad_p.value_ou,
    KS_Statistic_BM = test_out$ks_statistic_bm, KS_P_Value_BM = test_out$ks_p.value_bm,
    KS_Statistic_OU = test_out$ks_statistic_ou, KS_P_Value_OU = test_out$ks_p.value_ou,
    stringsAsFactors = FALSE
  )
  
  cat("   RESULT ROW ADDED:\n")
  print(row)
  
  if (length(test_out$bm_dist) > MIN_PICS_FOR_PLOT && length(test_out$ou_dist) > MIN_PICS_FOR_PLOT) {
    create_distribution_plot(order, trait, pics_std,
                             list(ad_p.value_bm = test_out$ad_p.value_bm,
                                  ad_p.value_ou = test_out$ad_p.value_ou,
                                  ks_p.value_bm = test_out$ks_p.value_bm,
                                  ks_p.value_ou = test_out$ks_p.value_ou),
                             test_out$bm_dist, test_out$ou_dist)
  } else {
    cat("   Not enough simulated PICs for plot\n")
  }
  
  row
}

## -------------------------------------------------
##  MAIN LOOP
## -------------------------------------------------
for (order in orders) {
  cat("\n", paste(rep("=", 60), collapse=""), "\n")
  cat("PROCESSING ORDER:", order, "\n")
  cat(paste(rep("=", 60), collapse=""), "\n")
  
  tree <- trees[[order]]
  if (is.null(tree)) {
    cat("   No valid tree → skip\n")
    next
  }
  
  if (!order %in% excel_sheets(file_path)) {
    cat("   Sheet '", order, "' not in Excel → skip\n")
    next
  }
  
  data <- tryCatch(read_excel(file_path, sheet = order), error = function(e) NULL)
  if (is.null(data) || nrow(data) == 0) {
    cat("   No data in sheet → skip\n")
    next
  }
  
  if (!"Species" %in% colnames(data)) {
    cat("   No 'Species' column → skip\n")
    next
  }
  
  data$Species <- trimws(gsub("[^[:alnum:]_]", "", gsub("\\s+", "_", data$Species)))
  data <- data[data$Species != "", ]
  data$species_lower <- tolower(data$Species)
  
  tree$tip.label <- gsub("[^[:alnum:]_]", "", gsub("\\s+", "_", tree$tip.label))
  tip_low <- tolower(tree$tip.label)
  
  match_idx <- match(data$species_lower, tip_low)
  keep <- !is.na(match_idx)
  common_sp <- data$Species[keep]
  
  if (length(common_sp) < MIN_SPECIES) {
    cat("   Only", length(common_sp), "matched species <", MIN_SPECIES, "→ skip\n")
    next
  }
  
  data <- data[keep, ]
  tree <- keep.tip(tree, common_sp)
  cat("   Matched", length(common_sp), "species. Starting trait analysis...\n")
  
  for (tr in traits) {
    if (!tr %in% colnames(data)) {
      cat("   ", tr, "not in data → skip\n")
      next
    }
    res <- perform_distribution_analysis(order, tree, data, tr)
    if (!is.null(res)) {
      test_results <- rbind(test_results, res)
    }
  }
}

## -------------------------------------------------
##  SAVE RESULTS
## -------------------------------------------------
excel_path <- "D:/Akash Ajay/insect bite force project/sravya work/GC skew/genomic_traits_distribution_test_results8.xlsx"
write.xlsx(test_results, excel_path, rowNames = FALSE)

cat("\n", paste(rep("=", 70), collapse=""), "\n")
cat("ANALYSIS COMPLETE!\n")
cat("Total result rows:", nrow(test_results), "\n")
cat("Excel saved →", excel_path, "\n")
cat("Plots saved →", plot_dir, "\n")
cat("Full log → ks_ad_test_output8.txt\n")
cat(paste(rep("=", 70), collapse=""), "\n")

sink()