# ==============================================================================
# IMPROVED GENE REGULATORY NETWORK PARAMETRIC ANALYSIS
# ==============================================================================
# Features:
# - Composite gene importance scoring with clear methodology
# - Gene categorization (Hub, Key Regulator, Moderate, Peripheral)
# - Directional classification (High>Mid>Low, etc.)
# - VitD-regulated gene emphasis
# - Three comparison approaches (Average, R1-only, Mid-normalized)
# - Parameter matrix sheets per replicate (Combined, High, Mid, Low)
# - VitD-weighted gene importance ranking
# ==============================================================================

# Install and load required packages
packages <- c("openxlsx", "stringr", "dplyr", "tidyr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
    install.packages(pkg, repos = "http://cran.r-project.org")
    library(pkg, character.only = TRUE)
  }
}

# ==============================================================================
# CONFIGURATION
# ==============================================================================
GENES <- c("ETS2", "FOS", "HIF1A", "EGR1", "ETS1", "JUN",
           "MYC", "NFKB1", "STAT1", "STAT3", "TP53")

# ==============================================================================
# PARSING FUNCTIONS
# ==============================================================================
extract_parameters <- function(filepath, replicate) {
  
  lines <- readLines(filepath, warn = FALSE)
  
  # Find replicate section
  rep_pattern <- paste0("REPLICATE ", replicate, " - OPTIMIZED")
  rep_start <- grep(rep_pattern, lines)
  if (length(rep_start) == 0) return(NULL)
  
  rep_end <- grep("^-{10,}$", lines)
  rep_end <- rep_end[rep_end > rep_start][1]
  if (is.na(rep_end)) rep_end <- length(lines)
  
  rep_text <- paste(lines[rep_start:rep_end], collapse = "\n")
  
  # Initialize results
  results <- data.frame(
    Interaction = character(),
    Target = character(),
    Source = character(),
    Type = character(),
    Gamma = numeric(),
    K = numeric(),
    n = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene equation
  for (gene in GENES) {
    eq_pattern <- paste0("d", gene, "/dt = (.+?)(?=\\n\\n|\\n[0-9]+\\.|Final)")
    eq_match <- str_match(rep_text, eq_pattern)
    
    if (is.na(eq_match[1,1])) next
    eq_text <- eq_match[1,2]
    
    # Extract activation terms
    for (source in c(GENES, "VitD")) {
      act_pattern <- paste0("([0-9.]+)\\s*\\*\\s*", source,
                            "\\^([0-9.]+)\\s*/\\s*\\(([0-9.]+)\\^")
      act_match <- str_match(eq_text, act_pattern)
      
      if (!is.na(act_match[1,1])) {
        act_part <- str_split(eq_text, "\\)\\s*\\*\\s*\\(")[[1]][1]
        if (grepl(paste0(source, "\\^"), act_part)) {
          results <- rbind(results, data.frame(
            Interaction = paste0(source, " -> ", gene),
            Target = gene,
            Source = source,
            Type = "Activation",
            Gamma = as.numeric(act_match[1,2]),
            K = as.numeric(act_match[1,4]),
            n = as.numeric(act_match[1,3]),
            stringsAsFactors = FALSE
          ))
        }
      }
    }
    
    # Extract repression terms
    for (source in c(GENES, "VitD")) {
      rep_pattern2 <- paste0("([0-9.]+)\\s*/\\s*\\(1\\s*\\+\\s*\\(",
                             source, "/([0-9.]+)\\)\\^([0-9.]+)\\)")
      rep_matches <- str_match_all(eq_text, rep_pattern2)[[1]]
      
      if (nrow(rep_matches) > 0) {
        results <- rbind(results, data.frame(
          Interaction = paste0(source, " -| ", gene),
          Target = gene,
          Source = source,
          Type = "Repression",
          Gamma = as.numeric(rep_matches[1,2]),
          K = as.numeric(rep_matches[1,3]),
          n = as.numeric(rep_matches[1,4]),
          stringsAsFactors = FALSE
        ))
      }
    }
  }
  
  results <- results[!duplicated(results$Interaction), ]
  return(results)
}

# ==============================================================================
# GENE IMPORTANCE - CORRECTED NORMALIZATION
# ==============================================================================
calculate_gene_importance <- function(general_file) {
  
  cat("\n================================================================================\n")
  cat("PART 1: GENE IMPORTANCE ANALYSIS (General Model)\n")
  cat("================================================================================\n\n")
  
  # Extract parameters from all replicates
  all_params <- list()
  for (rep in c("R1", "R2", "R3")) {
    params <- extract_parameters(general_file, rep)
    if (!is.null(params)) {
      params$Replicate <- rep
      all_params[[rep]] <- params
    }
  }
  
  all_params_df <- do.call(rbind, all_params)
  
  # Average gamma values across replicates
  avg_params <- all_params_df %>%
    group_by(Interaction, Target, Source, Type) %>%
    summarise(Avg_Gamma = mean(Gamma, na.rm = TRUE), .groups = "drop")
  
  # Initialize gene importance
  gene_importance <- data.frame(Gene = GENES, stringsAsFactors = FALSE)
  
  # ============================================================================
  # CONNECTIVITY METRICS (incoming + outgoing)
  # ============================================================================
  
  cat("Calculating connectivity...\n")
  
  gene_importance$Out_Activations <- sapply(GENES, function(g) {
    sum(avg_params$Source == g & avg_params$Type == "Activation")
  })
  
  gene_importance$Out_Repressions <- sapply(GENES, function(g) {
    sum(avg_params$Source == g & avg_params$Type == "Repression")
  })
  
  gene_importance$In_Activations <- sapply(GENES, function(g) {
    sum(avg_params$Target == g & avg_params$Type == "Activation")
  })
  
  gene_importance$In_Repressions <- sapply(GENES, function(g) {
    sum(avg_params$Target == g & avg_params$Type == "Repression")
  })
  
  gene_importance$Total_Connectivity <- 
    gene_importance$Out_Activations + 
    gene_importance$Out_Repressions + 
    gene_importance$In_Activations + 
    gene_importance$In_Repressions
  
  # ============================================================================
  # PARAMETER STRENGTH (incoming + outgoing gamma sums)
  # ============================================================================
  
  cat("Calculating parameter strength...\n")
  
  gene_importance$Out_Activation_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == g &
                               avg_params$Type == "Activation"], na.rm = TRUE)
  })
  
  gene_importance$Out_Repression_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == g &
                               avg_params$Type == "Repression"], na.rm = TRUE)
  })
  
  gene_importance$In_Activation_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Target == g &
                               avg_params$Type == "Activation"], na.rm = TRUE)
  })
  
  gene_importance$In_Repression_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Target == g &
                               avg_params$Type == "Repression"], na.rm = TRUE)
  })
  
  gene_importance$Total_Strength <- 
    gene_importance$Out_Activation_Strength + 
    gene_importance$Out_Repression_Strength + 
    gene_importance$In_Activation_Strength + 
    gene_importance$In_Repression_Strength
  
  # ============================================================================
  # NORMALIZE SEPARATELY, THEN COMBINE (CORRECT METHOD)
  # ============================================================================
  
  cat("Calculating composite scores...\n")
  
  # Normalize function
  normalize <- function(x) {
    if (max(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  # Normalize TotalDegree (connectivity) and TotalStrength separately
  gene_importance$NormDegree <- normalize(gene_importance$Total_Connectivity)
  gene_importance$NormStrength <- normalize(gene_importance$Total_Strength)
  
  # Composite score = average of normalized metrics
  gene_importance$CompositeScore <- 
    (gene_importance$NormDegree + gene_importance$NormStrength) / 2
  
  # Rank by CompositeScore (already the average of normalized metrics)
  gene_importance <- gene_importance[order(-gene_importance$CompositeScore), ]
  gene_importance$Rank <- 1:nrow(gene_importance)
  
  # ============================================================================
  # CATEGORIZATION
  # ============================================================================
  
  gene_importance$Role <- case_when(
    gene_importance$Rank <= 3 ~ "Hub Gene",
    gene_importance$Rank <= 6 ~ "Key Regulator",
    gene_importance$Rank <= 9 ~ "Moderate",
    TRUE ~ "Peripheral"
  )
  
  # VitD info (for reference)
  gene_importance$VitD_Regulated <- sapply(GENES, function(g) {
    vitd_gamma <- sum(avg_params$Avg_Gamma[avg_params$Source == "VitD" &
                                             avg_params$Target == g], na.rm = TRUE)
    total_incoming <- gene_importance$In_Activation_Strength[gene_importance$Gene == g] +
      gene_importance$In_Repression_Strength[gene_importance$Gene == g]
    
    if (total_incoming > 0 && vitd_gamma / total_incoming > 0.2) "Yes" else "No"
  })
  
  # Reorder columns to match expected output
  gene_importance <- gene_importance %>%
    select(Rank, Gene, Role,
           Out_Activation_Strength, Out_Repression_Strength,
           In_Activation_Strength, In_Repression_Strength,
           Out_Activations, Out_Repressions, In_Activations, In_Repressions,
           Total_Strength, Total_Connectivity,
           NormStrength, NormDegree, CompositeScore,
           VitD_Regulated)
  
  cat("\n✓ Gene importance analysis complete!\n")
  cat("\nTop 5 Most Important Genes:\n")
  print(gene_importance[1:5, c("Rank", "Gene", "Role", "CompositeScore")])
  
  return(gene_importance)
}

# ==============================================================================
# VITD-WEIGHTED GENE IMPORTANCE
# ==============================================================================
calculate_vitd_weighted_importance <- function(general_file) {
  
  cat("\n================================================================================\n")
  cat("PART 1b: VITD-WEIGHTED GENE IMPORTANCE ANALYSIS\n")
  cat("================================================================================\n\n")
  
  # Extract parameters from all replicates
  all_params <- list()
  for (rep in c("R1", "R2", "R3")) {
    params <- extract_parameters(general_file, rep)
    if (!is.null(params)) {
      params$Replicate <- rep
      all_params[[rep]] <- params
    }
  }
  
  all_params_df <- do.call(rbind, all_params)
  
  # Average gamma values across replicates
  avg_params <- all_params_df %>%
    group_by(Interaction, Target, Source, Type) %>%
    summarise(Avg_Gamma = mean(Gamma, na.rm = TRUE), .groups = "drop")
  
  # Initialize gene importance
  gene_importance <- data.frame(Gene = GENES, stringsAsFactors = FALSE)
  
  # ============================================================================
  # BASE CONNECTIVITY & STRENGTH 
  # ============================================================================
  
  gene_importance$Out_Activations <- sapply(GENES, function(g) {
    sum(avg_params$Source == g & avg_params$Type == "Activation")
  })
  gene_importance$Out_Repressions <- sapply(GENES, function(g) {
    sum(avg_params$Source == g & avg_params$Type == "Repression")
  })
  gene_importance$In_Activations <- sapply(GENES, function(g) {
    sum(avg_params$Target == g & avg_params$Type == "Activation")
  })
  gene_importance$In_Repressions <- sapply(GENES, function(g) {
    sum(avg_params$Target == g & avg_params$Type == "Repression")
  })
  gene_importance$Total_Connectivity <- 
    gene_importance$Out_Activations + gene_importance$Out_Repressions +
    gene_importance$In_Activations + gene_importance$In_Repressions
  
  gene_importance$Out_Activation_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == g & avg_params$Type == "Activation"], na.rm = TRUE)
  })
  gene_importance$Out_Repression_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == g & avg_params$Type == "Repression"], na.rm = TRUE)
  })
  gene_importance$In_Activation_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Target == g & avg_params$Type == "Activation"], na.rm = TRUE)
  })
  gene_importance$In_Repression_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Target == g & avg_params$Type == "Repression"], na.rm = TRUE)
  })
  gene_importance$Total_Strength <- 
    gene_importance$Out_Activation_Strength + gene_importance$Out_Repression_Strength +
    gene_importance$In_Activation_Strength + gene_importance$In_Repression_Strength
  
  # ============================================================================
  # VITD INTERACTION METRICS
  # ============================================================================
  
  cat("Calculating VitD interaction metrics...\n")
  
  # Direct VitD -> Gene activation/repression strength
  gene_importance$VitD_Activates_Gene_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == "VitD" &
                               avg_params$Target == g &
                               avg_params$Type == "Activation"], na.rm = TRUE)
  })
  
  gene_importance$VitD_Represses_Gene_Strength <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == "VitD" &
                               avg_params$Target == g &
                               avg_params$Type == "Repression"], na.rm = TRUE)
  })
  
  # Gene -> Other genes that are themselves VitD-regulated (indirect VitD coupling)
  # First identify VitD-regulated genes (VitD has any direct interaction)
  vitd_regulated_genes <- GENES[sapply(GENES, function(g) {
    vitd_strength <- sum(avg_params$Avg_Gamma[avg_params$Source == "VitD" &
                                                avg_params$Target == g], na.rm = TRUE)
    vitd_strength > 0
  })]
  
  cat("  VitD-regulated genes detected:", paste(vitd_regulated_genes, collapse = ", "), "\n")
  
  # How strongly does this gene regulate VitD-regulated genes (downstream VitD coupling)
  gene_importance$Downstream_VitD_Coupling <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Source == g &
                               avg_params$Target %in% vitd_regulated_genes], na.rm = TRUE)
  })
  
  # How strongly is this gene regulated BY VitD-regulated genes (upstream VitD coupling)
  gene_importance$Upstream_VitD_Coupling <- sapply(GENES, function(g) {
    sum(avg_params$Avg_Gamma[avg_params$Target == g &
                               avg_params$Source %in% vitd_regulated_genes], na.rm = TRUE)
  })
  
  # Total VitD interaction score: direct + indirect
  gene_importance$VitD_Direct_Total <- 
    gene_importance$VitD_Activates_Gene_Strength + 
    gene_importance$VitD_Represses_Gene_Strength
  
  gene_importance$VitD_Indirect_Total <- 
    gene_importance$Downstream_VitD_Coupling + 
    gene_importance$Upstream_VitD_Coupling
  
  gene_importance$VitD_Total_Interaction <- 
    gene_importance$VitD_Direct_Total + gene_importance$VitD_Indirect_Total
  
  # VitD fraction of total incoming strength (direct influence proportion)
  gene_importance$VitD_Direct_Fraction <- sapply(GENES, function(g) {
    total_in <- gene_importance$In_Activation_Strength[gene_importance$Gene == g] +
      gene_importance$In_Repression_Strength[gene_importance$Gene == g]
    vitd_direct <- gene_importance$VitD_Direct_Total[gene_importance$Gene == g]
    if (total_in > 0) vitd_direct / total_in else 0
  })
  
  # ============================================================================
  # NORMALIZE AND COMBINE WITH VITD BOOST
  # ============================================================================
  
  cat("Calculating VitD-weighted composite scores...\n")
  
  normalize <- function(x) {
    if (max(x, na.rm = TRUE) == 0) return(rep(0, length(x)))
    (x - min(x, na.rm = TRUE)) / (max(x, na.rm = TRUE) - min(x, na.rm = TRUE))
  }
  
  gene_importance$NormDegree    <- normalize(gene_importance$Total_Connectivity)
  gene_importance$NormStrength  <- normalize(gene_importance$Total_Strength)
  gene_importance$NormVitDTotal <- normalize(gene_importance$VitD_Total_Interaction)
  gene_importance$NormVitDFrac  <- normalize(gene_importance$VitD_Direct_Fraction)
  
  # Base score (same as original)
  gene_importance$BaseScore <- (gene_importance$NormDegree + gene_importance$NormStrength) / 2
  
  # VitD composite (direct + indirect normalized)
  gene_importance$VitD_Score <- (gene_importance$NormVitDTotal + gene_importance$NormVitDFrac) / 2
  
  # VitD-weighted composite:
  # 60% base network importance + 40% VitD interaction score
  gene_importance$VitD_Weighted_Score <- 
    0.60 * gene_importance$BaseScore + 
    0.40 * gene_importance$VitD_Score
  
  # Rank by VitD-weighted score
  gene_importance <- gene_importance[order(-gene_importance$VitD_Weighted_Score), ]
  gene_importance$VitD_Rank <- 1:nrow(gene_importance)
  
  # Original rank for comparison
  gene_importance_base_order <- gene_importance[order(-gene_importance$BaseScore), ]
  gene_importance$Original_Rank <- match(gene_importance$Gene, gene_importance_base_order$Gene)
  gene_importance$Rank_Change <- gene_importance$Original_Rank - gene_importance$VitD_Rank
  
  # Categorization
  gene_importance$Role <- case_when(
    gene_importance$VitD_Rank <= 3 ~ "Hub Gene",
    gene_importance$VitD_Rank <= 6 ~ "Key Regulator",
    gene_importance$VitD_Rank <= 9 ~ "Moderate",
    TRUE ~ "Peripheral"
  )
  
  gene_importance$VitD_Status <- case_when(
    gene_importance$VitD_Direct_Total > 0 ~ "Direct VitD Target",
    gene_importance$VitD_Indirect_Total > 0 ~ "Indirect VitD Coupling",
    TRUE ~ "No VitD Interaction"
  )
  
  # Select final columns
  gene_importance <- gene_importance %>%
    select(VitD_Rank, Original_Rank, Rank_Change, Gene, Role, VitD_Status,
           VitD_Activates_Gene_Strength, VitD_Represses_Gene_Strength,
           VitD_Direct_Total, VitD_Indirect_Total, VitD_Total_Interaction,
           VitD_Direct_Fraction,
           Total_Strength, Total_Connectivity,
           BaseScore, VitD_Score, VitD_Weighted_Score)
  
  cat("\n✓ VitD-weighted gene importance analysis complete!\n")
  cat("\nTop 5 Most Important Genes (VitD-weighted):\n")
  print(gene_importance[1:5, c("VitD_Rank", "Original_Rank", "Rank_Change", "Gene", "Role", 
                               "VitD_Status", "VitD_Weighted_Score")])
  
  return(gene_importance)
}

# ==============================================================================
# PARAMETER MATRIX (all 373 parameters x 4 models, per replicate)
# ==============================================================================
build_parameter_matrix <- function(general_file, high_file, mid_file, low_file, replicate) {
  
  cat(paste0("\n--- Building parameter matrix for ", replicate, " ---\n"))
  
  # Extract parameters for this replicate from each model
  combined_params <- extract_parameters(general_file, replicate)
  high_params     <- extract_parameters(high_file, replicate)
  mid_params      <- extract_parameters(mid_file, replicate)
  low_params      <- extract_parameters(low_file, replicate)
  
  # Get union of all interaction names
  all_interactions <- unique(c(
    if (!is.null(combined_params)) combined_params$Interaction,
    if (!is.null(high_params))     high_params$Interaction,
    if (!is.null(mid_params))      mid_params$Interaction,
    if (!is.null(low_params))      low_params$Interaction
  ))
  
  # Build helper to extract gamma by interaction name
  get_gamma <- function(params, interactions) {
    if (is.null(params)) return(rep(NA, length(interactions)))
    gamma_vals <- params$Gamma[match(interactions, params$Interaction)]
    return(gamma_vals)
  }
  
  get_meta <- function(params, interactions, col) {
    if (is.null(params)) return(rep(NA, length(interactions)))
    params[[col]][match(interactions, params$Interaction)]
  }
  
  # Reference params for metadata (use combined, fallback to any available)
  ref_params <- if (!is.null(combined_params)) combined_params else
    if (!is.null(high_params)) high_params else mid_params
  
  matrix_df <- data.frame(
    Interaction = all_interactions,
    Target      = get_meta(ref_params, all_interactions, "Target"),
    Source      = get_meta(ref_params, all_interactions, "Source"),
    Type        = get_meta(ref_params, all_interactions, "Type"),
    Gamma_Combined = get_gamma(combined_params, all_interactions),
    Gamma_High     = get_gamma(high_params,     all_interactions),
    Gamma_Mid      = get_gamma(mid_params,      all_interactions),
    Gamma_Low      = get_gamma(low_params,      all_interactions),
    stringsAsFactors = FALSE
  )
  
  # Add computed columns
  matrix_df$High_vs_Low_Diff <- matrix_df$Gamma_High - matrix_df$Gamma_Low
  matrix_df$High_vs_Low_FC   <- ifelse(!is.na(matrix_df$Gamma_Low) & matrix_df$Gamma_Low > 0,
                                       matrix_df$Gamma_High / matrix_df$Gamma_Low, NA)
  matrix_df$VitD_Interaction <- ifelse(matrix_df$Source == "VitD", "Yes", "No")
  
  # Sort: VitD interactions first, then by Target gene, then by magnitude
  matrix_df <- matrix_df %>%
    arrange(desc(VitD_Interaction == "Yes"), Target, Source, desc(abs(High_vs_Low_Diff)))
  
  cat(paste0("  ✓ Matrix built: ", nrow(matrix_df), " interactions found for ", replicate, "\n"))
  return(matrix_df)
}

# ==============================================================================
# DIFFERENTIAL ANALYSIS - ONLY REAL HIGH vs LOW DIFFERENCES
# ==============================================================================
approach_a_averaged <- function(high_file, mid_file, low_file) {
  
  cat("\n--- Approach A: Averaged Parameters ---\n")
  
  # Extract all parameters
  extract_all_reps <- function(filepath, condition) {
    all_reps <- list()
    for (rep in c("R1", "R2", "R3")) {
      params <- extract_parameters(filepath, rep)
      if (!is.null(params)) {
        params$Condition <- condition
        params$Replicate <- rep
        all_reps[[rep]] <- params
      }
    }
    do.call(rbind, all_reps)
  }
  
  high_all <- extract_all_reps(high_file, "High")
  mid_all  <- extract_all_reps(mid_file,  "Mid")
  low_all  <- extract_all_reps(low_file,  "Low")
  
  all_data <- rbind(high_all, mid_all, low_all)
  
  comparison <- all_data %>%
    group_by(Interaction, Target, Source, Type, Condition) %>%
    summarise(Avg_Gamma = mean(Gamma, na.rm = TRUE), .groups = "drop") %>%
    pivot_wider(names_from = Condition, values_from = Avg_Gamma)
  
  # Calculate differences
  comparison <- comparison %>%
    mutate(
      High_Low_Diff = High - Low,
      Abs_Diff = abs(High_Low_Diff),
      Fold_Change = ifelse(Low > 0, High / Low, NA),
      
      # Direction
      Direction = case_when(
        High > Mid & Mid > Low ~ "High > Mid > Low",
        High < Mid & Mid < Low ~ "High < Mid < Low",
        High > Low & Mid < Low ~ "High > Low > Mid",
        High > Low & Mid > High ~ "Mid > High > Low",
        Low > High & Mid > High ~ "Low > Mid > High",
        Low > High & Mid < High ~ "Low > High > Mid",
        High > Low ~ "High > Low",
        Low > High ~ "Low > High",
        TRUE ~ "Similar"
      ),
      
      VitD_Interaction = ifelse(Source == "VitD", "Yes", "No")
    ) %>%
    # FILTER OUT "Similar" - only keep real differences
    filter(Direction != "Similar") %>%
    arrange(desc(Abs_Diff))
  
  # Select columns
  comparison <- comparison %>%
    select(Interaction, Target, Source, Type, Direction,
           High, Mid, Low, High_Low_Diff, Fold_Change, VitD_Interaction)
  
  cat("✓ Found", nrow(comparison), "truly differential interactions\n")
  cat("\nTop 5 High vs Low Differentiators:\n")
  print(comparison[1:5, c("Interaction", "Direction", "High_Low_Diff")])
  
  return(comparison)
}

approach_b_r1_only <- function(high_file, mid_file, low_file) {
  
  cat("\n--- Approach B: R1 Only ---\n")
  
  high_r1 <- extract_parameters(high_file, "R1")
  mid_r1  <- extract_parameters(mid_file,  "R1")
  low_r1  <- extract_parameters(low_file,  "R1")
  
  if (is.null(high_r1) || is.null(mid_r1) || is.null(low_r1)) {
    cat("✗ R1 data not available\n")
    return(NULL)
  }
  
  high_r1$Condition <- "High"
  mid_r1$Condition  <- "Mid"
  low_r1$Condition  <- "Low"
  
  all_r1 <- rbind(high_r1, mid_r1, low_r1)
  
  comparison <- all_r1 %>%
    pivot_wider(names_from = Condition, values_from = Gamma)
  
  comparison <- comparison %>%
    mutate(
      High_Low_Diff = High - Low,
      Abs_Diff = abs(High_Low_Diff),
      Fold_Change = ifelse(Low > 0, High / Low, NA),
      
      Direction = case_when(
        High > Mid & Mid > Low ~ "High > Mid > Low",
        High < Mid & Mid < Low ~ "High < Mid < Low",
        High > Low & Mid < Low ~ "High > Low > Mid",
        High > Low & Mid > High ~ "Mid > High > Low",
        Low > High & Mid > High ~ "Low > Mid > High",
        Low > High & Mid < High ~ "Low > High > Mid",
        High > Low ~ "High > Low",
        Low > High ~ "Low > High",
        TRUE ~ "Similar"
      ),
      
      VitD_Interaction = ifelse(Source == "VitD", "Yes", "No")
    ) %>%
    filter(Direction != "Similar") %>%
    arrange(desc(Abs_Diff))
  
  comparison <- comparison %>%
    select(Interaction, Target, Source, Type, Direction,
           High, Mid, Low, High_Low_Diff, Fold_Change, VitD_Interaction)
  
  cat("✓ Found", nrow(comparison), "truly differential interactions\n")
  cat("\nTop 5 High vs Low Differentiators:\n")
  print(comparison[1:5, c("Interaction", "Direction", "High_Low_Diff")])
  
  return(comparison)
}

approach_c_normalized <- function(high_file, mid_file, low_file) {
  
  cat("\n--- Approach C: Normalized to Mid ---\n")
  
  extract_avg <- function(filepath) {
    all_reps <- list()
    for (rep in c("R1", "R2", "R3")) {
      params <- extract_parameters(filepath, rep)
      if (!is.null(params)) all_reps[[rep]] <- params
    }
    all_data <- do.call(rbind, all_reps)
    all_data %>%
      group_by(Interaction, Target, Source, Type) %>%
      summarise(Avg_Gamma = mean(Gamma, na.rm = TRUE), .groups = "drop")
  }
  
  high_avg <- extract_avg(high_file)
  mid_avg  <- extract_avg(mid_file)
  low_avg  <- extract_avg(low_file)
  
  comparison <- mid_avg %>%
    rename(Mid = Avg_Gamma) %>%
    left_join(high_avg %>% rename(High = Avg_Gamma),
              by = c("Interaction", "Target", "Source", "Type")) %>%
    left_join(low_avg %>% rename(Low = Avg_Gamma),
              by = c("Interaction", "Target", "Source", "Type"))
  
  comparison <- comparison %>%
    mutate(
      Norm_High = ifelse(Mid != 0, High / Mid, NA),
      Norm_Low  = ifelse(Mid != 0, Low / Mid, NA),
      Norm_Diff = Norm_High - Norm_Low,
      Abs_Norm_Diff = abs(Norm_Diff),
      Fold_Change = ifelse(Low > 0, High / Low, NA),
      
      Direction = case_when(
        Norm_High > 1.2 & Norm_Low < 0.8 ~ "High >> Mid >> Low",
        Norm_High > 1.1 & Norm_Low < 0.9 ~ "High > Mid > Low",
        Norm_High < 0.8 & Norm_Low > 1.2 ~ "Low >> Mid >> High",
        Norm_High < 0.9 & Norm_Low > 1.1 ~ "Low > Mid > High",
        abs(Norm_High - Norm_Low) < 0.2  ~ "Similar",
        Norm_High > Norm_Low ~ "High > Low",
        TRUE ~ "Low > High"
      ),
      
      VitD_Interaction = ifelse(Source == "VitD", "Yes", "No")
    ) %>%
    filter(Direction != "Similar") %>%
    arrange(desc(Abs_Norm_Diff))
  
  comparison <- comparison %>%
    select(Interaction, Target, Source, Type, Direction,
           Mid, High, Low, Norm_High, Norm_Low, Norm_Diff, 
           Fold_Change, VitD_Interaction)
  
  cat("✓ Found", nrow(comparison), "truly differential interactions\n")
  cat("\nTop 5 Normalized High vs Low Differentiators:\n")
  print(comparison[1:5, c("Interaction", "Direction", "Norm_Diff")])
  
  return(comparison)
}

# ==============================================================================
# MAIN FUNCTION
# ==============================================================================
run_parametric_analysis <- function(general_file, high_file, mid_file, low_file,
                                    output_file = "Parametric_Analysis_Results_recent.xlsx") {
  
  cat("================================================================================\n")
  cat("     PARAMETRIC ANALYSIS - CORRECTED VERSION\n")
  cat("================================================================================\n")
  
  if (!file.exists(general_file)) stop("General model file not found: ", general_file)
  
  # Part 1: Gene Importance
  gene_importance <- calculate_gene_importance(general_file)
  
  # Part 1b: VitD-Weighted Gene Importance
  vitd_gene_importance <- calculate_vitd_weighted_importance(general_file)
  
  # Part 2: Differential Analysis
  approach_a_results <- NULL
  approach_b_results <- NULL
  approach_c_results <- NULL
  
  if (file.exists(high_file) && file.exists(mid_file) && file.exists(low_file)) {
    cat("\n================================================================================\n")
    cat("PART 2: HIGH vs LOW DIFFERENTIATORS (filtering out similar)\n")
    cat("================================================================================\n")
    
    approach_a_results <- approach_a_averaged(high_file, mid_file, low_file)
    approach_b_results <- approach_b_r1_only(high_file, mid_file, low_file)
    approach_c_results <- approach_c_normalized(high_file, mid_file, low_file)
  }
  
  # Part 3: Parameter matrices per replicate
  param_matrix_r1 <- NULL
  param_matrix_r2 <- NULL
  param_matrix_r3 <- NULL
  
  if (file.exists(high_file) && file.exists(mid_file) && file.exists(low_file)) {
    cat("\n================================================================================\n")
    cat("PART 3: PARAMETER MATRICES PER REPLICATE\n")
    cat("================================================================================\n")
    
    param_matrix_r1 <- build_parameter_matrix(general_file, high_file, mid_file, low_file, "R1")
    param_matrix_r2 <- build_parameter_matrix(general_file, high_file, mid_file, low_file, "R2")
    param_matrix_r3 <- build_parameter_matrix(general_file, high_file, mid_file, low_file, "R3")
  }
  
  # Export to Excel
  cat("\n================================================================================\n")
  cat("EXPORTING TO EXCEL\n")
  cat("================================================================================\n")
  
  wb <- createWorkbook()
  
  header_style <- createStyle(
    fontSize = 11, fontColour = "#FFFFFF", halign = "center",
    fgFill = "#4472C4", textDecoration = "bold"
  )
  
  title_style <- createStyle(
    fontSize = 13, fontColour = "#FFFFFF", halign = "left",
    fgFill = "#203864", textDecoration = "bold"
  )
  
  hub_style     <- createStyle(fgFill = "#C6EFCE")
  key_style     <- createStyle(fgFill = "#FFEB9C")
  vitd_style    <- createStyle(fgFill = "#B4C7E7")
  vitd_hi_style <- createStyle(fgFill = "#9DC3E6")  # Stronger blue for direct VitD
  pos_diff_style <- createStyle(fgFill = "#E2EFDA")
  neg_diff_style <- createStyle(fgFill = "#FCE4D6")
  
  # ============================================================================
  # Sheet 1: Gene Importance
  # ============================================================================
  addWorksheet(wb, "1_Gene_Importance")
  writeData(wb, "1_Gene_Importance",
            "GENE IMPORTANCE: Normalized Connectivity + Normalized Strength",
            startRow = 1)
  mergeCells(wb, "1_Gene_Importance", rows = 1, cols = 1:15)
  addStyle(wb, "1_Gene_Importance", title_style, rows = 1, cols = 1)
  
  writeData(wb, "1_Gene_Importance", gene_importance, startRow = 3)
  addStyle(wb, "1_Gene_Importance", header_style, rows = 3,
           cols = 1:ncol(gene_importance), gridExpand = TRUE)
  
  hub_rows <- which(gene_importance$Role == "Hub Gene") + 3
  if (length(hub_rows) > 0) {
    addStyle(wb, "1_Gene_Importance", hub_style, rows = hub_rows,
             cols = 1:ncol(gene_importance), gridExpand = TRUE, stack = TRUE)
  }
  
  key_rows <- which(gene_importance$Role == "Key Regulator") + 3
  if (length(key_rows) > 0) {
    addStyle(wb, "1_Gene_Importance", key_style, rows = key_rows,
             cols = 1:ncol(gene_importance), gridExpand = TRUE, stack = TRUE)
  }
  
  setColWidths(wb, "1_Gene_Importance", cols = 1:ncol(gene_importance), widths = "auto")
  cat("  ✓ Sheet 1: Gene Importance\n")
  
  # ============================================================================
  # Sheet 2: VitD-Weighted Gene Importance (NEW)
  # ============================================================================
  addWorksheet(wb, "2_VitD_Weighted_Importance")
  writeData(wb, "2_VitD_Weighted_Importance",
            "GENE IMPORTANCE: VitD-Weighted (60% Network Score + 40% VitD Interaction Score)",
            startRow = 1)
  mergeCells(wb, "2_VitD_Weighted_Importance", rows = 1, cols = 1:17)
  addStyle(wb, "2_VitD_Weighted_Importance", title_style, rows = 1, cols = 1)
  
  # Legend row
  writeData(wb, "2_VitD_Weighted_Importance",
            "VitD Score = (Norm_VitD_Total_Interaction + Norm_VitD_Direct_Fraction) / 2  |  Rank_Change = Original_Rank - VitD_Rank (positive = rose in importance)",
            startRow = 2)
  
  writeData(wb, "2_VitD_Weighted_Importance", vitd_gene_importance, startRow = 4)
  addStyle(wb, "2_VitD_Weighted_Importance", header_style, rows = 4,
           cols = 1:ncol(vitd_gene_importance), gridExpand = TRUE)
  
  # Color hub/key rows
  hub_rows_v <- which(vitd_gene_importance$Role == "Hub Gene") + 4
  if (length(hub_rows_v) > 0) {
    addStyle(wb, "2_VitD_Weighted_Importance", hub_style, rows = hub_rows_v,
             cols = 1:ncol(vitd_gene_importance), gridExpand = TRUE, stack = TRUE)
  }
  key_rows_v <- which(vitd_gene_importance$Role == "Key Regulator") + 4
  if (length(key_rows_v) > 0) {
    addStyle(wb, "2_VitD_Weighted_Importance", key_style, rows = key_rows_v,
             cols = 1:ncol(vitd_gene_importance), gridExpand = TRUE, stack = TRUE)
  }
  
  # Highlight direct VitD targets
  direct_vitd_rows <- which(vitd_gene_importance$VitD_Status == "Direct VitD Target") + 4
  if (length(direct_vitd_rows) > 0) {
    addStyle(wb, "2_VitD_Weighted_Importance", vitd_hi_style, rows = direct_vitd_rows,
             cols = 6:12, gridExpand = TRUE, stack = TRUE)
  }
  
  setColWidths(wb, "2_VitD_Weighted_Importance", cols = 1:ncol(vitd_gene_importance), widths = "auto")
  cat("  ✓ Sheet 2: VitD-Weighted Gene Importance\n")
  
  # ============================================================================
  # Sheet 3: Approach A
  # ============================================================================
  if (!is.null(approach_a_results)) {
    addWorksheet(wb, "3_Approach_A_Averaged")
    writeData(wb, "3_Approach_A_Averaged",
              "APPROACH A: Top High vs Low Differentiators (Averaged, Filtered)",
              startRow = 1)
    mergeCells(wb, "3_Approach_A_Averaged", rows = 1, cols = 1:12)
    addStyle(wb, "3_Approach_A_Averaged", title_style, rows = 1, cols = 1)
    
    writeData(wb, "3_Approach_A_Averaged", approach_a_results, startRow = 3)
    addStyle(wb, "3_Approach_A_Averaged", header_style, rows = 3,
             cols = 1:ncol(approach_a_results), gridExpand = TRUE)
    
    addStyle(wb, "3_Approach_A_Averaged", hub_style, rows = 4:13,
             cols = 1:ncol(approach_a_results), gridExpand = TRUE, stack = TRUE)
    
    vitd_rows_a <- which(approach_a_results$VitD_Interaction == "Yes") + 3
    if (length(vitd_rows_a) > 0) {
      addStyle(wb, "3_Approach_A_Averaged", vitd_style, rows = vitd_rows_a,
               cols = 1:3, gridExpand = TRUE, stack = TRUE)
    }
    
    setColWidths(wb, "3_Approach_A_Averaged", cols = 1:ncol(approach_a_results), widths = "auto")
    cat("  ✓ Sheet 3: Approach A\n")
  }
  
  # ============================================================================
  # Sheet 4: Approach B
  # ============================================================================
  if (!is.null(approach_b_results)) {
    addWorksheet(wb, "4_Approach_B_R1_Only")
    writeData(wb, "4_Approach_B_R1_Only",
              "APPROACH B: Top High vs Low Differentiators (R1, Filtered)",
              startRow = 1)
    mergeCells(wb, "4_Approach_B_R1_Only", rows = 1, cols = 1:12)
    addStyle(wb, "4_Approach_B_R1_Only", title_style, rows = 1, cols = 1)
    
    writeData(wb, "4_Approach_B_R1_Only", approach_b_results, startRow = 3)
    addStyle(wb, "4_Approach_B_R1_Only", header_style, rows = 3,
             cols = 1:ncol(approach_b_results), gridExpand = TRUE)
    
    addStyle(wb, "4_Approach_B_R1_Only", hub_style, rows = 4:13,
             cols = 1:ncol(approach_b_results), gridExpand = TRUE, stack = TRUE)
    
    vitd_rows_b <- which(approach_b_results$VitD_Interaction == "Yes") + 3
    if (length(vitd_rows_b) > 0) {
      addStyle(wb, "4_Approach_B_R1_Only", vitd_style, rows = vitd_rows_b,
               cols = 1:3, gridExpand = TRUE, stack = TRUE)
    }
    
    setColWidths(wb, "4_Approach_B_R1_Only", cols = 1:ncol(approach_b_results), widths = "auto")
    cat("  ✓ Sheet 4: Approach B\n")
  }
  
  # ============================================================================
  # Sheet 5: Approach C
  # ============================================================================
  if (!is.null(approach_c_results)) {
    addWorksheet(wb, "5_Approach_C_Normalized")
    writeData(wb, "5_Approach_C_Normalized",
              "APPROACH C: Top High vs Low Differentiators (Normalized, Filtered)",
              startRow = 1)
    mergeCells(wb, "5_Approach_C_Normalized", rows = 1, cols = 1:14)
    addStyle(wb, "5_Approach_C_Normalized", title_style, rows = 1, cols = 1)
    
    writeData(wb, "5_Approach_C_Normalized", approach_c_results, startRow = 3)
    addStyle(wb, "5_Approach_C_Normalized", header_style, rows = 3,
             cols = 1:ncol(approach_c_results), gridExpand = TRUE)
    
    addStyle(wb, "5_Approach_C_Normalized", hub_style, rows = 4:13,
             cols = 1:ncol(approach_c_results), gridExpand = TRUE, stack = TRUE)
    
    vitd_rows_c <- which(approach_c_results$VitD_Interaction == "Yes") + 3
    if (length(vitd_rows_c) > 0) {
      addStyle(wb, "5_Approach_C_Normalized", vitd_style, rows = vitd_rows_c,
               cols = 1:3, gridExpand = TRUE, stack = TRUE)
    }
    
    setColWidths(wb, "5_Approach_C_Normalized", cols = 1:ncol(approach_c_results), widths = "auto")
    cat("  ✓ Sheet 5: Approach C\n")
  }
  
  # ============================================================================
  # Sheets 6-8: Parameter Matrices per Replicate (NEW)
  # ============================================================================
  write_matrix_sheet <- function(wb, sheet_name, matrix_data, replicate_label) {
    addWorksheet(wb, sheet_name)
    
    title_text <- paste0("PARAMETER MATRIX - ", replicate_label,
                         ": Gamma values for Combined, High, Mid, Low models",
                         " (", nrow(matrix_data), " interactions)")
    writeData(wb, sheet_name, title_text, startRow = 1)
    mergeCells(wb, sheet_name, rows = 1, cols = 1:12)
    addStyle(wb, sheet_name, title_style, rows = 1, cols = 1)
    
    writeData(wb, sheet_name, "Columns: Interaction | Target | Source | Type | Gamma_Combined | Gamma_High | Gamma_Mid | Gamma_Low | High_vs_Low_Diff | High_vs_Low_FC | VitD_Interaction",
              startRow = 2)
    
    writeData(wb, sheet_name, matrix_data, startRow = 4)
    addStyle(wb, sheet_name, header_style, rows = 4,
             cols = 1:ncol(matrix_data), gridExpand = TRUE)
    
    # Highlight VitD rows
    vitd_rows_m <- which(matrix_data$VitD_Interaction == "Yes") + 4
    if (length(vitd_rows_m) > 0) {
      addStyle(wb, sheet_name, vitd_style, rows = vitd_rows_m,
               cols = 1:ncol(matrix_data), gridExpand = TRUE, stack = TRUE)
    }
    
    # Color High > Low differences green, Low > High differences orange
    if ("High_vs_Low_Diff" %in% colnames(matrix_data)) {
      diff_col_idx <- which(colnames(matrix_data) == "High_vs_Low_Diff")
      pos_rows <- which(!is.na(matrix_data$High_vs_Low_Diff) & 
                          matrix_data$High_vs_Low_Diff > 0 &
                          matrix_data$VitD_Interaction != "Yes") + 4
      neg_rows <- which(!is.na(matrix_data$High_vs_Low_Diff) & 
                          matrix_data$High_vs_Low_Diff < 0 &
                          matrix_data$VitD_Interaction != "Yes") + 4
      if (length(pos_rows) > 0) {
        addStyle(wb, sheet_name, pos_diff_style, rows = pos_rows,
                 cols = diff_col_idx, gridExpand = TRUE, stack = TRUE)
      }
      if (length(neg_rows) > 0) {
        addStyle(wb, sheet_name, neg_diff_style, rows = neg_rows,
                 cols = diff_col_idx, gridExpand = TRUE, stack = TRUE)
      }
    }
    
    setColWidths(wb, sheet_name, cols = 1:ncol(matrix_data), widths = "auto")
  }
  
  if (!is.null(param_matrix_r1)) {
    write_matrix_sheet(wb, "6_Param_Matrix_R1", param_matrix_r1, "R1")
    cat("  ✓ Sheet 6: Parameter Matrix R1\n")
  }
  if (!is.null(param_matrix_r2)) {
    write_matrix_sheet(wb, "7_Param_Matrix_R2", param_matrix_r2, "R2")
    cat("  ✓ Sheet 7: Parameter Matrix R2\n")
  }
  if (!is.null(param_matrix_r3)) {
    write_matrix_sheet(wb, "8_Param_Matrix_R3", param_matrix_r3, "R3")
    cat("  ✓ Sheet 8: Parameter Matrix R3\n")
  }
  
  # ============================================================================
  # Sheet 9: Explanation
  # ============================================================================
  addWorksheet(wb, "9_Explanation")
  
  explanation <- data.frame(
    Content = c(
      "PARAMETRIC ANALYSIS - CORRECTED METHODOLOGY WITH EXTENSIONS",
      "",
      "PART 1: GENE IMPORTANCE (Sheet 1)",
      "==================================",
      "",
      "Step 1: Calculate RAW metrics",
      "   Total_Connectivity = Out_Activations + Out_Repressions + In_Activations + In_Repressions",
      "   Total_Strength = Out_Act_Strength + Out_Rep_Strength + In_Act_Strength + In_Rep_Strength",
      "   (Total_Strength is SUM of ALL gamma values)",
      "",
      "Step 2: Normalize SEPARATELY",
      "   NormDegree = (Total_Connectivity - min) / (max - min)",
      "   NormStrength = (Total_Strength - min) / (max - min)",
      "",
      "Step 3: CompositeScore = (NormDegree + NormStrength) / 2",
      "Step 4: Rank by CompositeScore",
      "",
      "Categories: Rank 1-3 = Hub Gene | 4-6 = Key Regulator | 7-9 = Moderate | 10-11 = Peripheral",
      "",
      "PART 1b: VITD-WEIGHTED GENE IMPORTANCE (Sheet 2)  [NEW]",
      "=========================================================",
      "",
      "Extends base scoring to reward genes that interact most with Vitamin D.",
      "",
      "VitD Metrics per gene:",
      "   VitD_Activates_Gene_Strength  = sum of gamma where VitD -> Gene (Activation)",
      "   VitD_Represses_Gene_Strength  = sum of gamma where VitD -> Gene (Repression)",
      "   VitD_Direct_Total             = VitD_Activates + VitD_Represses",
      "   VitD_Indirect_Total           = Downstream coupling (Gene regulates VitD-targets)",
      "                                 + Upstream coupling (VitD-targets regulate Gene)",
      "   VitD_Direct_Fraction          = VitD_Direct_Total / Total_Incoming_Strength",
      "   VitD_Score = (Norm_VitD_Total_Interaction + Norm_VitD_Direct_Fraction) / 2",
      "",
      "Final formula:",
      "   VitD_Weighted_Score = 0.60 * BaseScore + 0.40 * VitD_Score",
      "   (Genes with strong VitD interactions are boosted in ranking)",
      "",
      "Rank_Change = Original_Rank - VitD_Rank",
      "   Positive = gene rose in importance when VitD weighting applied",
      "   Negative = gene dropped when VitD weighting applied",
      "",
      "Color coding: Dark blue = Direct VitD Target | Yellow = Key Regulator | Green = Hub Gene",
      "",
      "PART 2: HIGH vs LOW DIFFERENTIATORS (Sheets 3-5)",
      "==================================================",
      "",
      "Only shows interactions where High and Low are actually different (filtered).",
      "Approach A: Averaged across R1, R2, R3 — ranked by |High - Low|",
      "Approach B: R1 only — same logic",
      "Approach C: Normalized to Mid — ranked by |Norm_High - Norm_Low|, threshold 0.2",
      "",
      "PART 3: PARAMETER MATRICES PER REPLICATE (Sheets 6-8)  [NEW]",
      "==============================================================",
      "",
      "Three sheets (R1, R2, R3), each showing ALL interactions found across all 4 models.",
      "Columns: Interaction | Target | Source | Type |",
      "         Gamma_Combined | Gamma_High | Gamma_Mid | Gamma_Low |",
      "         High_vs_Low_Diff | High_vs_Low_FC | VitD_Interaction",
      "",
      "Sorting: VitD interactions first, then by Target gene, then by |High-Low| magnitude.",
      "Color: Blue = VitD interaction | Green diff = High > Low | Orange diff = Low > High",
      "NA = interaction not found in that model for that replicate."
    )
  )
  
  writeData(wb, "9_Explanation", explanation, colNames = FALSE)
  setColWidths(wb, "9_Explanation", cols = 1, widths = 120)
  cat("  ✓ Sheet 9: Explanation\n")
  
  saveWorkbook(wb, output_file, overwrite = TRUE)
  
  cat("\n✓ ANALYSIS COMPLETE!\n")
  cat("Output:", output_file, "\n")
  cat("\nSheets produced:\n")
  cat("  1. Gene Importance (base)\n")
  cat("  2. VitD-Weighted Gene Importance  [NEW]\n")
  cat("  3. Approach A - Averaged differential\n")
  cat("  4. Approach B - R1-only differential\n")
  cat("  5. Approach C - Normalized differential\n")
  cat("  6. Parameter Matrix R1  [NEW]\n")
  cat("  7. Parameter Matrix R2  [NEW]\n")
  cat("  8. Parameter Matrix R3  [NEW]\n")
  cat("  9. Explanation\n\n")
  
  return(list(
    gene_importance      = gene_importance,
    vitd_gene_importance = vitd_gene_importance,
    approach_a           = approach_a_results,
    approach_b           = approach_b_results,
    approach_c           = approach_c_results,
    param_matrix_r1      = param_matrix_r1,
    param_matrix_r2      = param_matrix_r2,
    param_matrix_r3      = param_matrix_r3
  ))
}

# ==============================================================================
# USAGE
# ==============================================================================
cat("================================================================================
PARAMETRIC ANALYSIS - WITH VITD-WEIGHTED RANKING + PARAMETER MATRICES
================================================================================

results <- run_parametric_analysis(
  general_file = 'optimized_equations_combined.txt',
  high_file    = 'optimized_equations_High1.txt',
  mid_file     = 'optimized_equations_Mid1.txt',
  low_file     = 'optimized_equations_Low1.txt',
  output_file  = 'Parametric_Analysis_Results_recent1.xlsx'
)

NEW SHEETS:
-----------
Sheet 2 - VitD-Weighted Importance:
  VitD_Weighted_Score = 0.60 * BaseScore + 0.40 * VitD_Score
  VitD_Score captures both DIRECT VitD interactions and INDIRECT coupling
  Rank_Change shows how much each gene moved vs original ranking

Sheets 6-8 - Parameter Matrices (R1, R2, R3):
  All interactions as rows x 4 models (Combined/High/Mid/Low) as columns
  Shows Gamma values side by side + High-Low difference + fold change
================================================================================
")
