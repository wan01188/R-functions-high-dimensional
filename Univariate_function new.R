create_volcano_plot <- function(X, Y, x_id, y_id, outcome, method, control_var = NULL, already_normalized = FALSE, already_scaled = FALSE) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    install.packages("ggrepel")
  }
  library(ggplot2)
  library(ggrepel)
  
  if (method %in% c("Mann-Whitney", "Kruskal-Wallis", "ANOVA")) {
    
    if(method=="ANOVA"){
      colnames(X) <- make.names(colnames(X))
    }

    if (any(is.na(X))) {
      stop("Error: The input data X contains missing values. Please handle missing values before running the function.")
    }
    
    # Merge datasets
    proteomics_data <- merge(X, Y, by.x = x_id, by.y = y_id)
    
    if (!is.factor(proteomics_data[[outcome]])) {
      proteomics_data[[outcome]] <- as.factor(proteomics_data[[outcome]])
    }
    if (method == "Mann-Whitney" && length(levels(proteomics_data[[outcome]])) != 2) {
      stop("The outcome variable must have exactly two levels for Mann-Whitney U test.")
    }
    
    # Subset numeric columns and initialize results
    numeric_cols <- sapply(proteomics_data, is.numeric)
    protein_cols <- proteomics_data[, numeric_cols, drop = FALSE]
    results <- data.frame(Protein = colnames(protein_cols), p_value = NA, effect_size = NA)
    
    # Handle scaling and unscaling vectorized
    if (already_scaled) {
      means <- colMeans(protein_cols, na.rm = TRUE)
      sds <- apply(protein_cols, 2, sd, na.rm = TRUE)
      protein_cols <- sweep(protein_cols, 2, sds, `*`)
      protein_cols <- sweep(protein_cols, 2, means, `+`)
    }
    
    # Vectorized p-value and effect size calculation
    if (method == "Mann-Whitney") {
      group1 <- protein_cols[proteomics_data[[outcome]] == levels(proteomics_data[[outcome]])[2], , drop = FALSE]
      group0 <- protein_cols[proteomics_data[[outcome]] == levels(proteomics_data[[outcome]])[1], , drop = FALSE]
      
      effect_sizes <- if (already_normalized) {
        colMeans(group1, na.rm = TRUE) - colMeans(group0, na.rm = TRUE)
      } else {
        log2(colMeans(group1, na.rm = TRUE) / colMeans(group0, na.rm = TRUE))
      }
      
      p_values <- apply(protein_cols, 2, function(col) {
        test <- wilcox.test(col ~ proteomics_data[[outcome]], na.action = na.omit)
        test$p.value
      })
      results$p_value <- p_values
      results$effect_size <- effect_sizes
    }
    
    
    # Kruskal-Wallis test
    if (method == "Kruskal-Wallis") {
      kw_tests <- apply(protein_cols, 2, function(col) {
        test <- kruskal.test(col ~ proteomics_data[[outcome]])
        eta_squared <- (test$statistic - length(levels(proteomics_data[[outcome]])) + 1) / 
          (nrow(proteomics_data) - length(levels(proteomics_data[[outcome]])))
        c(test$p.value, eta_squared)
      })
      
      kw_tests <- t(kw_tests)
      results$p_value <- kw_tests[, 1]
      results$effect_size <- kw_tests[, 2]
    }
    
    # ANOVA
    if (method == "ANOVA") {
      formulas <- sapply(colnames(protein_cols), function(col) {
        if (!is.null(control_var)) {
          paste0(col, " ~ ", outcome, " + ", paste(control_var, collapse = " + "))
        } else {
          paste0(col, " ~ ", outcome)
        }
      })
      
      anova_results <- lapply(formulas, function(formula_str) {
        fit <- aov(as.formula(formula_str), data = proteomics_data)
        test_result <- summary(fit)[[1]]
        p_value <- test_result[["Pr(>F)"]][1]
        effect_size <- test_result[["Sum Sq"]][1] / sum(test_result[["Sum Sq"]])
        c(p_value, effect_size)
      })
      
      anova_results <- do.call(rbind, anova_results)
      results$p_value <- anova_results[, 1]
      results$effect_size <- anova_results[, 2]
    }
    
    # Calculate q-values
    results$q_value <- p.adjust(results$p_value, method = "BH")
    
    # Identify significance
    results$significance <- ifelse(results$p_value < 0.05, 
                                   ifelse(results$q_value < 0.05, "sig(BH-adj)", "sig(unadjusted)"), 
                                   "non-sig")
    
    results$log_pvalue <- -log10(results$p_value)
    results$log_qvalue <- -log10(results$q_value)
    results <- results[order(results$q_value), ]
    
    # Generate the plot
    if (method %in% c("Mann-Whitney", "Kruskal-Wallis")) {
      p_value_threshold <- 0.05
      q_value_threshold <- max(results$p_value[results$q_value < 0.05], na.rm = TRUE)
      log_pvalue_threshold <- -log10(p_value_threshold)
      log_qvalue_threshold <- -log10(q_value_threshold)
      effect_size_threshold <- 1  
      
      volcano_plot <- ggplot(results, aes(x = effect_size, y = log_pvalue, color = significance)) +
        geom_point(size = 3) +  
        scale_color_manual(values = c("non-sig" = "#999999", "sig(unadjusted)" = "#0072B2", "sig(BH-adj)" = "#56B4E9")) +
        geom_hline(yintercept = log_pvalue_threshold, linetype = "dashed", color = "darkred", size = 1) +
        geom_hline(yintercept = log_qvalue_threshold, linetype = "dashed", color = "darkgreen", size = 1) +
        geom_vline(xintercept = c(-effect_size_threshold, effect_size_threshold), linetype = "dashed", color = "blue", size = 1) +
        theme_minimal() +
        labs(title = "Volcano Plot", x = "Log fold change", y = "-log10(p-value)", color = "Significance") +
        ggrepel::geom_text_repel(aes(label = ifelse(significance %in% c("sig(BH-adj)", "sig(unadjusted)"), Protein, "")), 
                                 size = 3, max.overlaps = 25)
      
      return(list(results = results, volcano_plot = volcano_plot))
    } else {
      # Generate p-value plot for Kruskal-Wallis or ANOVA
      library(ggplot2)
      library(ggrepel)
      
      p_value_threshold <- 0.05
      q_value_threshold <- max(results$p_value[results$q_value < 0.05], na.rm = TRUE)
      log_pvalue_threshold <- -log10(p_value_threshold)
      log_qvalue_threshold <- -log10(q_value_threshold)
      
      ANOVA_pvalue_plot <- ggplot(results, aes(x = 1:nrow(results), y = log_pvalue, color = significance)) +
        geom_point(size = 3) +
        scale_color_manual(values = c("non-sig" = "#999999", "sig(unadjusted)" = "#0072B2", "sig(BH-adj)" = "#56B4E9")) +
        geom_hline(yintercept = log_pvalue_threshold, linetype = "dashed", color = "darkred", size = 1) +
        geom_hline(yintercept = log_qvalue_threshold, linetype = "dashed", color = "darkgreen", size = 1) +
        theme_minimal() +
        labs(title = paste(method, "Results"), x = "Compounds", y = "-log10(p-value)", color = "Significance") +
        ggrepel::geom_text_repel(aes(label = ifelse(significance == "sig(BH-adj)", Protein, "")), 
                                 size = 3, max.overlaps = 25)
      
      return(list(results = results, pvalue_plot = ANOVA_pvalue_plot))
    }
  }
  
  
  ####################################################
  if (method == "lm") {
    
    # Merge datasets
    proteomics_data <- merge(X, Y, by.x = x_id, by.y = y_id, all = TRUE)
    
    # Check if the outcome exists
    if (!(outcome %in% colnames(proteomics_data))) {
      stop("Error: Outcome variable is missing from the merged dataset.")
    }
    
    # Ensure valid column names
    colnames(proteomics_data) <- make.names(colnames(proteomics_data))
    
    # Select numeric columns and remove constant variables
    numeric_cols <- sapply(proteomics_data, is.numeric)
    protein_cols <- proteomics_data[, numeric_cols, drop = FALSE]
    protein_cols <- protein_cols[, sapply(protein_cols, function(x) length(unique(x)) > 1), drop = FALSE]
    
    # Check if there are valid predictors left
    if (ncol(protein_cols) == 0) {
      stop("Error: No valid numeric predictor columns found after filtering.")
    }
    
    # Initialize results
    results <- data.frame(Protein = colnames(protein_cols), p_value = NA, effect_size = NA)
    
    # Handle scaling if needed
    if (already_scaled) {
      means <- colMeans(protein_cols, na.rm = TRUE)
      sds <- apply(protein_cols, 2, sd, na.rm = TRUE)
      protein_cols <- sweep(protein_cols, 2, sds, `*`)
      protein_cols <- sweep(protein_cols, 2, means, `+`)
    }
    formulas <- lapply(colnames(protein_cols), function(col) {
      if (!is.null(control_var)) {
        as.formula(paste(col, "~", outcome, "+", paste(control_var, collapse = " + ")))
      } else {
        as.formula(paste(col, "~", outcome))
      }
    })
    
    lm_results <- lapply(formulas, function(formula) {
      fit <- tryCatch(lm(formula, data = proteomics_data, na.action = na.omit), error = function(e) NULL)
      
      if (is.null(fit) || length(coef(fit)) < 2) {
        return(c(NA, NA))
      }
      
      summary_fit <- summary(fit)
      return(c(coef(summary_fit)[2, "Estimate"], coef(summary_fit)[2, "Pr(>|t|)"]))
    })
    
    lm_results <- do.call(rbind, lm_results)
    results$effect_size <- lm_results[, 1]
    results$p_value <- lm_results[, 2]
    ### **Compute q-values**
    results$q_value <- p.adjust(results$p_value, method = "BH")
    results$significance <- ifelse(results$p_value < 0.05, 
                                   ifelse(results$q_value < 0.05, "sig(BH-adj)", "sig(unadjusted)"), 
                                   "non-sig")
    
    results$log_pvalue <- -log10(results$p_value)
    results$log_qvalue <- -log10(results$q_value)
    results <- results[order(results$q_value), ]
    
    ### **Generate the Correct Type of Plot**
      volcano_plot <- ggplot(results, aes(x = effect_size, y = log_pvalue, color = significance)) +
        geom_point(size = 3) +
        scale_color_manual(values = c("non-sig" = "#999999", "sig(unadjusted)" = "#0072B2", "sig(BH-adj)" = "#56B4E9")) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "darkred", size = 1) +
        theme_minimal() +
        labs(title = "Volcano Plot", x = "Regression Coefficient (β)", y = "-log10(p-value)", color = "Significance") +
        ggrepel::geom_text_repel(aes(label = ifelse(significance %in% c("sig(BH-adj)", "sig(unadjusted)"), Protein, "")), 
                                 size = 3, max.overlaps = 25)
      return(list(results = results, volcano_plot = volcano_plot))

  }
  
}
