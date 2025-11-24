# ===== Parallel GC Skew with Debugging =====

# Ensure Biostrings is available
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}
if (!requireNamespace("Biostrings", quietly = TRUE)) {
  BiocManager::install("Biostrings")
}
library(Biostrings)

# ---- Setup with Debugging ----
dir_path <- "C:/Users/LENOVO/OneDrive/Documents"

# First, let's see what files are actually in the directory
all_files <- list.files(dir_path, full.names = TRUE)
cat("Total files in directory:", length(all_files), "\n")
cat("First 10 files:\n")
print(head(all_files, 10))

# Now let's check specifically for FASTA files
files <- list.files(
  dir_path,
  pattern = "\\.(fna|fa|fasta)(\\.gz)?$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\nFASTA files found:", length(files), "\n")
if (length(files) > 0) {
  print(files)
}

# Let's also check with a simpler pattern
files_simple <- list.files(
  dir_path,
  pattern = "\\.fna$",
  full.names = TRUE,
  ignore.case = TRUE
)

cat("\n.fna files found:", length(files_simple), "\n")
if (length(files_simple) > 0) {
  print(files_simple)
}

# Check if your specific file exists
specific_file <- "C:/Users/LENOVO/Documents/GCA_000006945.2_ASM694v2_genomic.fna"
cat("\nDoes specific file exist?", file.exists(specific_file), "\n")

# If no files are found, let's stop here for debugging
if (length(files) == 0) {
  cat("No FASTA files found. Let's investigate the directory structure...\n")
  
  # Check subdirectories too
  all_files_recursive <- list.files(dir_path, pattern = "\\.(fna|fa|fasta)", 
                                    full.names = TRUE, recursive = TRUE, ignore.case = TRUE)
  cat("FASTA files found recursively:", length(all_files_recursive), "\n")
  if (length(all_files_recursive) > 0) {
    print(all_files_recursive)
  }
  
  stop("Please check the debug output above to understand why files aren't being found.")
}

# ---- If files are found, continue with the rest of the script ----
# ---- Functions ----
calculate_gc_skew <- function(sequence) {
  freqs <- alphabetFrequency(sequence)
  G_count <- sum(freqs[, "G"], na.rm = TRUE)
  C_count <- sum(freqs[, "C"], na.rm = TRUE)
  if ((G_count + C_count) == 0) NA_real_ else (G_count - C_count) / (G_count + C_count)
}

compute_gc_skew_for_file <- function(file_path) {
  # Wrapped in tryCatch for resilience per-file
  tryCatch({
    seqs <- readDNAStringSet(file_path)  # supports .gz
    freqs <- alphabetFrequency(seqs)
    G_total <- sum(freqs[, "G"], na.rm = TRUE)
    C_total <- sum(freqs[, "C"], na.rm = TRUE)
    gc_skew <- if ((G_total + C_total) == 0) NA_real_ else (G_total - C_total) / (G_total + C_total)
    
    data.frame(
      accession = basename(file_path),
      G_count   = G_total,
      C_count   = C_total,
      GC_skew   = gc_skew,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    # Return a row with NA if this file fails
    data.frame(
      accession = basename(file_path),
      G_count = NA_real_,
      C_count = NA_real_,
      GC_skew = NA_real_,
      stringsAsFactors = FALSE
    )
  })
}

# ---- Parallel execution (Windows-friendly) ----
# Use all but one core
cores_to_use <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(cores_to_use)

# Load packages & objects on workers
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(Biostrings))
  NULL
})
parallel::clusterExport(cl,
                        varlist = c("calculate_gc_skew", "compute_gc_skew_for_file"),
                        envir = environment()
)

# Process in parallel
results_list <- parallel::parLapply(cl, files, compute_gc_skew_for_file)

# Clean up cluster
parallel::stopCluster(cl)

# ---- Combine & save ----
results <- do.call(rbind, results_list)
output_file <- file.path(dir_path, "gc_skew_results.csv")
write.csv(results, output_file, row.names = FALSE)

cat("✅ Done. Results saved to:", output_file, "\n")
cat("Total files processed:", nrow(results), "\n")
cat("Failures (NA GC_skew):", sum(is.na(results$GC_skew)), "\n")