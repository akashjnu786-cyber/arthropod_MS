library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)

# Create the morphology dataset
morph_data <- data.frame(
  Order = c("Blattodea", "Blattodea", "Blattodea", "Blattodea", "Blattodea", "Blattodea", "Blattodea",
            "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera",
            "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera",
            "Mantodea", "Mantodea", "Mantodea", "Mantodea", "Mantodea", "Mantodea", "Mantodea",
            "Odonata", "Odonata", "Odonata", "Odonata", "Odonata", "Odonata", "Odonata",
            "Orthoptera", "Orthoptera", "Orthoptera", "Orthoptera", "Orthoptera", "Orthoptera", "Orthoptera",
            "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea"),
  Trait = rep(c("Bite Force", "Head Width", "Head Height", "Head Length", "Thorax Width", "Wing Length", "Body Length"), 7),
  Lambda = c(0.6591, 0.8699, 0.9218, 7.89e-05, 0.5687, 1.0195, 0.9048,
             7.55e-05, 0.7235, 0.7739, 0.6653, 0.7608, 0.6341, 0.6864,
             0.9999, 0.9999, 0.9999, 7.33e-05, 0.9999, 0.9999, 0.9999,
             5.60e-05, 0.6816, 0.5713, 0.8795, 0.6524, 0.5738, 0.6025,
             0.1852, 0.7125, 0.8631, 0.9509, 0.9546, 0.5781, 0.5744,
             7.33e-05, 0.9865, 0.9779, 0.9996, 0.9894, 0.9934, 0.9864,
             1.00e-04, 1.039, 0.8281, 0.9758, 1.0727, 1.0933, 1.0904),
  Lambda_P = c(0.0363, 0.0001, 0.0001, 1, 0.007, 0.0298, 0.0001,
               1, 0.0001, 0.0001, 0.0141, 0.0001, 0.008, 0.0039,
               0.0044, 0.0001, 0.0001, 1, 0.0001, 0.0001, 0.0012,
               1, 0.0322, 0.2823, 0.0135, 0.0742, 0.2305, 0.0685,
               0.4583, 0.0115, 0.0021, 0.0008, 0.0004, 0.0582, 0.0363,
               1, 0.0153, 0.0137, 0.1051, 0.009, 0.1338, 0.0286,
               1, 0.0053, 0.0283, 0.0012, 0.0476, 0.0541, 0.014),
  K = c(0.6788, 1.3175, 1.5666, 0.3062, 0.59, 0.9469, 1.3999,
        0.343, 0.6801, 0.5976, 0.5812, 0.8632, 0.469, 0.5726,
        1.3393, 1.2511, 1.0943, 0.4511, 1.2299, 2.6793, 1.1922,
        0.7672, 1.263, 0.9414, 1.376, 1.1259, 1.0172, 1.1029,
        0, 0, 0, 0, 0, 0, 0,
        0.0087, 0.1302, 0.0901, 0.2246, 0.1554, 0.124, 0.111,
        0.5475, 1.1175, 0.8625, 1.3241, 0.8446, 0.7435, 0.9269),
  K_P = c(0.0521, 0.001, 0.001, 0.7247, 0.2152, 0.013, 0.001,
          0.2513, 0.025, 0.017, 0.015, 0.003, 0.049, 0.009,
          0.001, 0.001, 0.007, 0.2553, 0.002, 0.001, 0.001,
          0.2272, 0.003, 0.0941, 0.005, 0.012, 0.018, 0.012,
          0.977, 0.4695, 0.2022, 0.0691, 0.033, 0.5776, 0.6937,
          0.5245, 0.015, 0.036, 0.002, 0.009, 0.007, 0.024,
          0.2843, 0.003, 0.026, 0.001, 0.021, 0.0671, 0.015)
)

# Create significance indicators
morph_data <- morph_data %>%
  mutate(
    Lambda_Sig = ifelse(Lambda_P < 0.05, "Significant", "Non-significant"),
    K_Sig = ifelse(K_P < 0.05, "Significant", "Non-significant"),
    # Create combined metrics for visualization
    Lambda_Value = ifelse(Lambda_Sig == "Significant", Lambda, NA),
    K_Value = ifelse(K_Sig == "Significant", K, NA)
  )

# Create the heatmap for Lambda
lambda_plot <- ggplot(morph_data, aes(x = Order, y = Trait)) +
  geom_tile(aes(fill = Lambda_Value), color = "white", size = 0.8) +
  geom_point(aes(shape = Lambda_Sig, size = Lambda_Sig), alpha = 0.8) +
  scale_fill_gradient2(
    name = "Lambda (λ) Value",
    low = "blue", 
    mid = "white",
    high = "red",
    midpoint = 0.5,
    na.value = "grey90",
    limits = c(0, 1.2)
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("Significant" = 16, "Non-significant" = 4),
    labels = c("Significant" = "P < 0.05", "Non-significant" = "P ≥ 0.05")
  ) +
  scale_size_manual(
    name = "Significance", 
    values = c("Significant" = 2, "Non-significant" = 1.5),
    guide = "none"
  ) +
  labs(
    title = "Phylogenetic Signal (Pagel's λ) for Morphological Traits",
    subtitle = "",
    x = "Insect Orders",
    y = "Morphological Traits"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),  # Add this
    axis.title.y = element_text(size = 14),  # Add this
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),  # Changed to 15
    plot.subtitle = element_text(hjust = 0.5)
  )

# Create the heatmap for K
k_plot <- ggplot(morph_data, aes(x = Order, y = Trait)) +
  geom_tile(aes(fill = K_Value), color = "white", size = 0.8) +
  geom_point(aes(shape = K_Sig, size = K_Sig), alpha = 0.8) +
  scale_fill_gradient2(
    name = "Blomberg's K Value",
    low = "lightblue", 
    mid = "white",
    high = "darkgreen",
    midpoint = 0.5,
    na.value = "grey90",
    limits = c(0, 1.6)
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("Significant" = 16, "Non-significant" = 4),
    labels = c("Significant" = "P < 0.05", "Non-significant" = "P ≥ 0.05")
  ) +
  scale_size_manual(
    name = "Significance", 
    values = c("Significant" = 2, "Non-significant" = 1.5),
    guide = "none"
  ) +
  labs(
    title = "Phylogenetic Signal (Blomberg's K) for Morphological Traits",
    subtitle = "",
    x = "Insect Orders",
    y = "Morphological Traits"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),  # Add this
    axis.title.y = element_text(size = 14),  # Add this
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),  # Changed to 15
    plot.subtitle = element_text(hjust = 0.5)
  )
# Display plots
print(lambda_plot)
print(k_plot)

# Save plots
ggsave("morphology_lambda_heatmap.png", lambda_plot, width = 12, height = 8, dpi = 800, bg = "white")
ggsave("morphology_k_heatmap.png", k_plot, width = 12, height = 8, dpi = 800, bg = "white")

# Create the genomic features dataset
genomic_data <- data.frame(
  Order = c("Blattodea", "Blattodea", "Blattodea", "Blattodea", "Blattodea",
            "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera", "Coleoptera", 
            "Diptera", "Diptera", "Diptera", "Diptera", "Diptera",
            "Hemiptera", "Hemiptera", "Hemiptera", "Hemiptera", "Hemiptera",
            "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera", "Hymenoptera",
            "Lepidoptera", "Lepidoptera", "Lepidoptera", "Lepidoptera", "Lepidoptera",
            "Odonata", "Odonata", "Odonata", "Odonata", "Odonata",
            "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea", "Phasmatodea"),
  Trait = c("Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes",
            "Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes"),
  Lambda = c(1.132948, 1.103162, NA, NA, NA,
             0.590418, 0.993005, 0.258911, 0.933583, 7.3e-05,
             0.996556, 0.995927, 0.885318, 0.937019, 0.945148,
             0.38761, 0.999927, 0.553997, 7.3e-05, 7.3e-05,
             0.699162, 0.976785, 0.652368, 0.684122, 0.205473,
             0.987357, 0.991442, 0.593737, 0.045627, 7.3e-05,
             0.530635, 1.171793, NA, NA, NA,
             0.662159, 0.100512, NA, NA, NA),
  Lambda_P = c(0.0023, 0.0024, NA, NA, NA,
               0.0008, 1e-04, 0.0155, 0.2389, 1,
               1e-04, 1e-04, 1e-04, 1e-04, 1e-04,
               1e-04, 1e-04, 0.0005, 1, 1,
               1e-04, 1e-04, 0.001, 0.4351, 0.4062,
               1e-04, 1e-04, 1e-04, 0.9428, 1,
               0.3026, 0.001, NA, NA, NA,
               0.0154, 0.703, NA, NA, NA),
  K = c(1.615199, 2.034449, NA, NA, NA,
        0.000874, 0.006149, 0.299505, 0.872627, 0.812008,
        0.004279, 0.14433, 0.086101, 0.206222, 0.096159,
        0.038918, 1.736331, 0.423634, 0.217049, 0.168052,
        6e-05, 0.000389, 0.23742, 0.342129, 0.256373,
        0.092202, 0.06001, 0.074821, 0.098289, 0.000564,
        0.753279, 2.348526, NA, NA, NA,
        0.075296, 0.019315, NA, NA, NA),
  K_P = c(0.004, 0.003, NA, NA, NA,
          0.3774, 0.0681, 0.8699, 0.1692, 0.2513,
          0.5976, 0.001, 0.004, 0.001, 0.1321,
          0.968, 0.001, 0.1031, 0.7658, 0.973,
          0.5856, 0.001, 0.001, 0.2002, 0.5556,
          0.002, 0.003, 0.013, 0.1862, 0.8929,
          0.1091, 0.001, NA, NA, NA,
          0.0621, 0.8819, NA, NA, NA)
)

# Reorder the trait levels
genomic_data <- genomic_data %>%
  mutate(Trait = factor(Trait, levels = c("Genome Size", "Genomic GC", "Chromosome Number", "Coding genes", "Non-coding genes")))

# Create significance indicators (treat NA as Non-significant)
genomic_data <- genomic_data %>%
  mutate(
    Lambda_Sig = case_when(
      is.na(Lambda_P) ~ "Non-significant",
      Lambda_P < 0.05 ~ "Significant",
      TRUE ~ "Non-significant"
    ),
    K_Sig = case_when(
      is.na(K_P) ~ "Non-significant",
      K_P < 0.05 ~ "Significant", 
      TRUE ~ "Non-significant"
    ),
    # Create combined metrics for visualization
    Lambda_Value = ifelse(Lambda_Sig == "Significant", Lambda, NA),
    K_Value = ifelse(K_Sig == "Significant", K, NA)
  )

# Create the heatmap for Lambda (genomic features)
lambda_genomic_plot <- ggplot(genomic_data, aes(x = Order, y = Trait)) +
  geom_tile(aes(fill = Lambda_Value), color = "white", size = 0.8) +
  geom_point(aes(shape = Lambda_Sig, size = Lambda_Sig), alpha = 0.8) +
  scale_fill_gradient2(
    name = "Lambda (λ) Value",
    low = "blue", 
    mid = "white",
    high = "red",
    midpoint = 0.5,
    na.value = "grey90",
    limits = c(0, 1.2)
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("Significant" = 16, "Non-significant" = 4),
    labels = c("Significant" = "P < 0.05", "Non-significant" = "P ≥ 0.05 / NA")
  ) +
  scale_size_manual(
    name = "Significance", 
    values = c("Significant" = 2, "Non-significant" = 1.5),
    guide = "none"
  ) +
  labs(
    title = "Phylogenetic Signal (Pagel's λ) for Genomic Features",
    subtitle = "",
    x = "Insect Orders",
    y = "Genomic Features"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Create the heatmap for K (genomic features)
k_genomic_plot <- ggplot(genomic_data, aes(x = Order, y = Trait)) +
  geom_tile(aes(fill = K_Value), color = "white", size = 0.8) +
  geom_point(aes(shape = K_Sig, size = K_Sig), alpha = 0.8) +
  scale_fill_gradient2(
    name = "Blomberg's K Value",
    low = "lightblue", 
    mid = "white",
    high = "darkgreen",
    midpoint = 0.5,
    na.value = "grey90",
    limits = c(0, 2.5)
  ) +
  scale_shape_manual(
    name = "Significance",
    values = c("Significant" = 16, "Non-significant" = 4),
    labels = c("Significant" = "P < 0.05", "Non-significant" = "P ≥ 0.05 / NA")
  ) +
  scale_size_manual(
    name = "Significance", 
    values = c("Significant" = 2, "Non-significant" = 1.5),
    guide = "none"
  ) +
  labs(
    title = "Phylogenetic Signal (Blomberg's K) for Genomic Features",
    subtitle = "",
    x = "Insect Orders",
    y = "Genomic Features"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 14),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
    plot.subtitle = element_text(hjust = 0.5)
  )

# Display plots
print(lambda_genomic_plot)
print(k_genomic_plot)

# Save plots
ggsave("genomic_lambda_heatmap.png", lambda_genomic_plot, width = 12, height = 8, dpi = 800, bg = "white")
ggsave("genomic_k_heatmap.png", k_genomic_plot, width = 12, height = 8, dpi = 800, bg = "white")