############## violin/box plot
# create_violin_boxplot <- function(X, Y, id, outcome) {
#   library(tidyverse)
#   library(ggplot2)
#   # Merge data by the specified ID column
#   merged_data <- X %>%
#     inner_join(Y, by = id)
#   
#   # Reshape data to long format
#   long_data <- merged_data %>%
#     pivot_longer(
#       cols = -c(all_of(id), all_of(outcome)),  # Exclude ID and outcome columns
#       names_to = "Protein",
#       values_to = "Expression"
#     )
#   
#   # Generate the plot
#   plot <- ggplot(long_data, aes(
#     x = factor(.data[[outcome]], levels = c(1, 0), labels = c("alcohol", "Non-alcohol")), 
#     y = Expression, 
#     fill = factor(.data[[outcome]], levels = c(1, 0), labels = c("alcohol", "Non-alcohol"))
#   )) +
#     geom_violin(trim = FALSE, alpha = 0.8) +
#     geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.7) +
#     facet_wrap(~ Protein, scales = "free_y", ncol = 3) +  # Adjust columns
#     labs(
#       title = paste("Violin and Box Plots of Proteins by", outcome),
#       x = paste(outcome, "Status"),
#       y = "Protein Expression"
#     ) +
#     theme_minimal() +
#     theme(
#       plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
#       axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
#       axis.text.y = element_text(size = 12, face = "bold"),
#       axis.title.x = element_text(size = 14, face = "bold"),
#       axis.title.y = element_text(size = 14, face = "bold"),
#       strip.text = element_text(size = 14, face = "bold")
#     ) +
#     scale_fill_manual(
#       values = c("alcohol" = "#0072B2", "Non-alcohol" = "#56B4E9"), 
#       name = paste(outcome, "Status")
#     )
#   
#   return(plot)
# }






create_violin_boxplot <- function(X, Y, id, outcome) {
  library(tidyverse)
  library(ggplot2)
  # Merge data by the specified ID column
  merged_data <- X %>%
    inner_join(Y, by = id)
  
  # Reshape data to long format
  long_data <- merged_data %>%
    pivot_longer(
      cols = -c(all_of(id), all_of(outcome)),  # Exclude ID and outcome columns
      names_to = "Protein",
      values_to = "Expression"
    )
  
  # Generate dynamic labels for the legend and x-axis
  outcome_labels <- c("1" = outcome, "0" = paste0("Non-", outcome))
  
  # Generate the plot
  plot <- ggplot(long_data, aes(
    x = factor(.data[[outcome]], levels = c(1, 0), labels = outcome_labels), 
    y = Expression, 
    fill = factor(.data[[outcome]], levels = c(1, 0), labels = outcome_labels)
  )) +
    geom_violin(trim = FALSE, alpha = 0.8) +
    geom_boxplot(width = 0.2, position = position_dodge(width = 0.9), outlier.shape = NA, alpha = 0.7) +
    facet_wrap(~ Protein, scales = "free_y", ncol = 3) +  # Adjust columns
    labs(
      title = paste("Violin and Box Plots of Proteins by", outcome),
      x = paste(outcome, "Status"),
      y = "Protein Expression"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, face = "bold"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.title.x = element_text(size = 14, face = "bold"),
      axis.title.y = element_text(size = 14, face = "bold"),
      strip.text = element_text(size = 14, face = "bold")
    ) +
    scale_fill_manual(
      values = setNames(c("#0072B2", "#56B4E9"), outcome_labels), 
      name = paste(outcome, "Status")
    )
  
  return(plot)
}






# X <- X_1.0alcohol_withID_75iqr_original
# X[, -1] <- log2(X[, -1])
# X[,-1] <- scale(X[,-1])
# Y <- CSF_Y_updated
# colnames(X)[1] <- "Sample_ID"
# colnames(Y)[1] <- "Sample_ID"
# Y <- Y[,c(1,6)] # alcohol
# Y[,2] <- ifelse(Y[,2]=="Yes",1,0)
# violin_boxplot_alcohol_1.0 <- create_violin_boxplot(X, Y, "Sample_ID", "alcohol")
# ggsave("violin_boxplot_alcohol_1.0.png", plot = violin_boxplot_alcohol_1.0, width = 15, height = 20, dpi = 300)
# 
# 
# X <- X_1.0HAND_withID_75iqr_original
# X[, -1] <- log2(X[, -1])
# X[,-1] <- scale(X[,-1])
# Y <- CSF_Y_updated
# colnames(X)[1] <- "Sample_ID"
# colnames(Y)[1] <- "Sample_ID"
# Y <- Y[,c(1,4)] # HAND
# Y[,2] <- ifelse(Y[,2]=="Yes",1,0)
# violin_boxplot_HAND_1.0 <- create_violin_boxplot(X, Y, "Sample_ID", "HAND")
# ggsave("violin_boxplot_HAND_1.0.png", plot = violin_boxplot_HAND_1.0, width = 15, height = 20, dpi = 300)



