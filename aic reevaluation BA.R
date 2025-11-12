library(readxl)
library(dplyr)
library(openxlsx)

# Read morphological AIC data
aic_data_morphological <- read_excel("D:/Akash Ajay/insect bite force project/Morphology/bayou/bayou_results.xlsx")

# Apply Burnham & Anderson ΔAIC interpretation
aic_analysis_morphological <- aic_data_morphological %>%
  rowwise() %>%
  mutate(
    # Find best pulsed model AIC
    pulsed_aics = list(c(AIC_PE1, AIC_PE1_BM, AIC_PE2, AIC_PE2_BM, AIC_PE3, AIC_PE3_BM)),
    Best_Pulsed_AIC = min(unlist(pulsed_aics), na.rm = TRUE),
    
    # Calculate ΔAIC (BM - Best_Pulsed)
    Delta_AIC = AIC_BM - Best_Pulsed_AIC,  # POSITIVE means pulsed is better
    
    # Burnham & Anderson interpretation
    BA_Support_Level = case_when(
      is.na(Delta_AIC) ~ "Unable to calculate",
      Delta_AIC >= 10 ~ "ESSENTIALLY NO support for BM",
      Delta_AIC >= 4 ~ "CONSIDERABLY LESS support for BM",
      Delta_AIC >= 2 ~ "SUBSTANTIALLY LESS support for BM",
      Delta_AIC >= 0 ~ "SUBSTANTIAL support for both models",
      Delta_AIC >= -2 ~ "SUBSTANTIAL support for both models",
      Delta_AIC >= -4 ~ "SUBSTANTIALLY LESS support for pulsed",
      Delta_AIC >= -10 ~ "CONSIDERABLY LESS support for pulsed",
      TRUE ~ "ESSENTIALLY NO support for pulsed"
    ),
    
    # Recommended action
    Recommended_Action = case_when(
      Delta_AIC >= 4 ~ "PREFER pulsed model",
      Delta_AIC <= -4 ~ "PREFER BM model", 
      abs(Delta_AIC) < 2 ~ "RETAIN both for comparison",
      TRUE ~ "CAUTION: weak preference"
    ),
    
    # Final model choice
    Final_Model_Choice = case_when(
      Delta_AIC >= 2 ~ "Pulsed Evolution",
      Delta_AIC <= -2 ~ "Brownian Motion",
      TRUE ~ "Ambiguous (|ΔAIC| < 2)"
    )
  ) %>%
  select(Order, Trait, Sample_Size, AIC_BM, Best_Pulsed_AIC, Delta_AIC, 
         BA_Support_Level, Recommended_Action, Final_Model_Choice, Best_Model)

# Save results
write.xlsx(aic_analysis_morphological, 
           "D:/Akash Ajay/insect bite force project/Morphology/bayou/morphological_aic_BA_interpretation.xlsx")