## --------------------------------------------------------------
##  PGLS – GC-skew data (genomic_features_with_gc_skew1.xlsx)
##  – 100% working species matching + 1e-3 on zero branches
## --------------------------------------------------------------

library(ape)
library(geiger)
library(phytools)
library(caper)
library(readxl)
library(dplyr)
library(openxlsx)

# -----------------------------------------------------------------
# 1. PATHS & SETTINGS
# -----------------------------------------------------------------
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

raw_traits <- c(
  "Genome_Size","Genome_GC","Chromosome Number",
  "Coding genes","Non - coding genes","GC_skew"
)

selected_traits <- c(
  "Genome_Size","Genome_GC","Chromosome_Number",
  "Coding_genes","Non_coding_genes","GC_skew"
)

# -----------------------------------------------------------------
# 2. CLEAN DATA – same as working phylogenetic-signal script
# -----------------------------------------------------------------
clean_data <- function(df) {
  df <- as.data.frame(df)
  if (!"Species" %in% colnames(df)) stop("Column 'Species' missing")
  df$Species <- trimws(as.character(df$Species))
  df <- df[!is.na(df$Species), , drop = FALSE]
  
  df$species_clean <- gsub(" ", "_", df$Species)
  rownames(df) <- df$species_clean
  
  present <- intersect(raw_traits, colnames(df))
  if (length(present) == 0) stop("No trait columns found")
  df <- df[, c("species_clean", present), drop = FALSE]
  
  colnames(df) <- gsub("[ -]", "_", colnames(df))
  colnames(df) <- gsub("_{2,}", "_", colnames(df))
  
  rename_map <- setNames(selected_traits, gsub("[ -]", "_", raw_traits))
  idx <- match(colnames(df), names(rename_map))
  colnames(df)[!is.na(idx)] <- rename_map[idx[!is.na(idx)]]
  
  df
}

# -----------------------------------------------------------------
# 3. PREPARE TREE – add 1e-3 to zero branches
# -----------------------------------------------------------------
make_tree_ready <- function(tree) {
  if (is.null(tree)) return(NULL)
  tree <- multi2di(tree)
  tree$node.label <- NULL
  
  zero <- which(tree$edge.length <= 0 | is.na(tree$edge.length))
  if (length(zero)) tree$edge.length[zero] <- 1e-3
  
  if (!is.ultrametric(tree)) {
    tree <- tryCatch(force.ultrametric(tree, method = "extend"),
                     error = function(e) chronos(tree, lambda = 0, quiet = TRUE))
  }
  tree
}

# -----------------------------------------------------------------
# 4. ROBUST PGLS – FIXED: species_clean now in pair_df
# -----------------------------------------------------------------
run_pgls <- function(order, data, tree, traits) {
  cat("\n=== PGLS for", order, "===\n")
  
  common <- intersect(rownames(data), tree$tip.label)
  cat("  Common species:", length(common), "\n")
  if (length(common) < 5) {
    cat("  Skipping – <5 species matched\n")
    return(NULL)
  }
  
  data <- data[common, , drop = FALSE]
  tree <- keep.tip(tree, common)
  
  # Keep species_clean as a column for comparative.data
  data$species_clean <- rownames(data)
  
  comp_data <- comparative.data(phy = tree, data = data,
                                names.col = species_clean, vcv = TRUE)
  
  res <- data.frame(
    Order = character(), Response_Trait = character(),
    Predictor_Trait = character(), Slope = numeric(),
    P_value = numeric(), R_squared = numeric(), F_value = numeric(),
    df_num = numeric(), df_den = numeric(),
    logLik = numeric(), AIC = numeric(),
    lambda = numeric(), kappa = numeric(), delta = numeric(),
    Best_Model = character(), Sample_Size = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (resp in traits) {
    cat("  • Response:", resp, "\n")
    preds <- setdiff(traits, resp)
    for (pred in preds) {
      # ---- include species_clean in pair_df -----------------------
      pair_df <- data[, c("species_clean", resp, pred), drop = FALSE]
      pair_df <- pair_df[complete.cases(pair_df), , drop = FALSE]
      n <- nrow(pair_df)
      cat("      ", resp, "~", pred, " n =", n, "\n")
      
      if (n < 3) {
        res <- rbind(res, data.frame(Order=order, Response_Trait=resp,
                                     Predictor_Trait=pred, Slope=NA, P_value=NA,
                                     R_squared=NA, F_value=NA, df_num=NA, df_den=NA,
                                     logLik=NA, AIC=NA, lambda=NA, kappa=NA,
                                     delta=NA, Best_Model=NA, Sample_Size=n,
                                     stringsAsFactors=FALSE))
        next
      }
      
      comp_pair <- comparative.data(phy = tree, data = pair_df,
                                    names.col = species_clean, vcv = TRUE)
      formula <- as.formula(paste(resp, "~", pred))
      
      mod <- tryCatch(pgls(formula, data = comp_pair,
                           lambda="ML", kappa="ML", delta="ML"),
                      error = function(e) NULL)
      
      if (is.null(mod)) {
        m1 <- tryCatch(pgls(formula, data = comp_pair, lambda="ML", kappa=1, delta=1), error=function(e) NULL)
        m2 <- tryCatch(pgls(formula, data = comp_pair, lambda=1, kappa=1, delta="ML"), error=function(e) NULL)
        aics <- c(if(!is.null(m1)) AIC(m1) else Inf,
                  if(!is.null(m2)) AIC(m2) else Inf)
        best <- which.min(aics)
        mod <- if (aics[best] < Inf) get(paste0("m", best)) else NULL
        best_name <- if (aics[best] < Inf) c("Lambda","Delta")[best] else NA
      } else {
        best_name <- "Full"
      }
      
      if (!is.null(mod)) {
        s <- summary(mod)
        res <- rbind(res, data.frame(
          Order = order,
          Response_Trait = resp,
          Predictor_Trait = pred,
          Slope = round(coef(s)[2,1], 3),
          P_value = round(coef(s)[2,4], 4),
          R_squared = round(s$r.squared, 3),
          F_value = round(s$fstatistic["value"], 3),
          df_num = s$df[1], df_den = s$df[2],
          logLik = round(logLik(mod), 3),
          AIC = round(AIC(mod), 3),
          lambda = round(mod$param["lambda"], 3),
          kappa = round(mod$param["kappa"], 3),
          delta = round(mod$param["delta"], 3),
          Best_Model = best_name,
          Sample_Size = n,
          stringsAsFactors = FALSE
        ))
      } else {
        res <- rbind(res, data.frame(Order=order, Response_Trait=resp,
                                     Predictor_Trait=pred, Slope=NA, P_value=NA,
                                     R_squared=NA, F_value=NA, df_num=NA, df_den=NA,
                                     logLik=NA, AIC=NA, lambda=NA, kappa=NA,
                                     delta=NA, Best_Model=NA, Sample_Size=n,
                                     stringsAsFactors=FALSE))
      }
    }
  }
  res
}

# -----------------------------------------------------------------
# 5. MAIN LOOP
# -----------------------------------------------------------------
all_results <- data.frame()

for (ord in orders) {
  cat("\n=== PROCESSING", ord, "===\n")
  
  df_raw <- tryCatch(read_excel(file_path, sheet = ord),
                     error = function(e) { cat("  READ ERROR:", e$message, "\n"); NULL })
  if (is.null(df_raw)) next
  
  df <- tryCatch(clean_data(df_raw),
                 error = function(e) { cat("  CLEAN ERROR:", e$message, "\n"); NULL })
  if (is.null(df)) next
  
  tr <- tryCatch(read.tree(tree_paths[[ord]]), error = function(e) NULL)
  if (is.null(tr)) { cat("  TREE MISSING\n"); next }
  tr <- make_tree_ready(tr)
  
  res <- run_pgls(ord, df, tr, selected_traits)
  if (!is.null(res)) {
    all_results <- rbind(all_results, res)
    write.xlsx(res, paste0(ord, "_pgls_gc_skew_results.xlsx"))
  }
}

# -----------------------------------------------------------------
# 6. SAVE MASTER + SIGNIFICANT
# -----------------------------------------------------------------
sig <- all_results %>% filter(P_value < 0.05 & !is.na(P_value))

write.xlsx(
  list(All_Results = all_results, Significant_Results = sig),
  "PGLS_GC_Skew_All_Results.xlsx"
)

# -----------------------------------------------------------------
# 7. SUMMARY
# -----------------------------------------------------------------
cat("\n===== FINAL SUMMARY =====\n")
for (ord in orders) {
  sub <- all_results %>% filter(Order == ord)
  if (nrow(sub) == 0) next
  cat("\n", ord, ":\n")
  cat("  Pairs analysed :", nrow(sub), "\n")
  cat("  Significant (p<0.05):", sum(sub$P_value < .05, na.rm = TRUE), "\n")
  cat("  Min n :", min(sub$Sample_Size), " Max n :", max(sub$Sample_Size), "\n")
}

cat("\nDone! All PGLS results saved.\n")