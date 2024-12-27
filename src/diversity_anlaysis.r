library(phyloseq)
library(ggplot2)
library(gridExtra)
library(dunn.test)  # For non-parametric post-hoc tests
library(rstatix)    # For effect sizes
library(multcomp)   # For multiple comparison procedures
library(FSA)        # For Dunn's test with Bonferroni adjustment

# Function to plot histogram of sample sizes
plot_sample_size_histogram <- function(physeq, breaks=50, main="Sample size distribution", xlab="Sample size (log10)", ylab="Frequency", col="#007c80") {
  hist(log10(sample_sums(physeq)), breaks=breaks, main=main, xlab=xlab, ylab=ylab, col=col)
}

# Function to filter samples based on log10 sample size
filter_samples_by_size <- function(physeq, lower_threshold=NULL, upper_threshold=NULL) {
  if (!is.null(lower_threshold)) {
    physeq <- subset_samples(physeq, lower_threshold < log10(sample_sums(physeq)))
  }
  if (!is.null(upper_threshold)) {
    physeq <- subset_samples(physeq, log10(sample_sums(physeq)) < upper_threshold)
  }
  return(physeq)
}

# Function to perform the entire workflow
analyze_sample_sizes <- function(physeq, lower_threshold=NULL, upper_threshold=NULL, breaks=50) {
  # Plot initial histogram
  plot_sample_size_histogram(physeq, breaks=breaks)
  
  # Filter samples
  ps_filtered <- filter_samples_by_size(physeq, lower_threshold, upper_threshold)
  
  # Plot histogram of filtered data
  plot_sample_size_histogram(ps_filtered, breaks=breaks, main="Filtered Sample size distribution")
  
  return(ps_filtered)
}





create_siamcat_plot <- function(psg, label_column, case_value, output_prefix) {
  # Create label
  sc_label <- create.label(meta=sample_data(psg), label=label_column, case=case_value)
  
  # Create SIAMCAT object
  siamcat_o <- siamcat(phyloseq=psg, label=sc_label)
  
  
  # Filter features
  siamcat_o <- filter.features(siamcat_o, filter.method='abundance', cutoff=1e-03)
  siamcat_o <- filter.features(siamcat_o, filter.method='prevalence', cutoff=0.05, feature.type='filtered')
  
  # Check associations
  siamcat_o <- check.associations(siamcat_o)
  
  association.plot(siamcat_o, panels=c("fc", "prevalence"), prompt = FALSE, verbose = 0)
  # Final association plot
  svg(paste0("../plots/",output_prefix, "_final_association.svg"), width=12, height=7)
  association.plot(siamcat_o, panels=c("fc", "prevalence"), prompt = FALSE, verbose = 0)
  dev.off()
  cat("plot saved to", paste0(output_prefix, "_final_association.svg"), "\n")
  
  return(siamcat_o)
}

# Function to run comprehensive post-hoc analysis
run_posthoc_analysis <- function(data, dependent, grouping) {
  # Print column names for debugging
  print(paste("Testing", grouping))
  print("Available columns:")
  print(colnames(data))
  
  # Ensure the grouping variable exists
  if (!(grouping %in% colnames(data))) {
    stop(paste("Grouping variable", grouping, "not found in data"))
  }
  
  # 1. Shapiro test for normality
  shapiro_results <- tapply(data[[dependent]], 
                            data[[grouping]], 
                            function(x) shapiro.test(x)$p.value)
  
  # 2. Levene's test for homogeneity of variance
  levene_formula <- as.formula(paste(dependent, "~", grouping))
  levene_result <- car::leveneTest(levene_formula, data = data)
  
  # 3. Choose appropriate post-hoc test based on assumptions
  if (all(shapiro_results > 0.05) && levene_result$`Pr(>F)`[1] > 0.05) {
    # Tukey HSD
    aov_result <- aov(levene_formula, data = data)
    tukey_result <- TukeyHSD(aov_result)
    posthoc_result <- list(
      test = "Tukey HSD",
      result = tukey_result
    )
  } else {
    # Dunn's test
    dunn_result <- FSA::dunnTest(levene_formula, 
                                 data = data,
                                 method = "bonferroni")
    posthoc_result <- list(
      test = "Dunn's test",
      result = dunn_result
    )
  }
  
  # 4. Calculate effect sizes using rstatix
  # Convert grouping variable to factor if it's not already
  data[[grouping]] <- as.factor(data[[grouping]])
  
  eff_size <- data %>%
    rstatix::cohens_d(as.formula(paste(dependent, "~", grouping))) %>%
    as.data.frame()
  
  return(list(
    shapiro = shapiro_results,
    levene = levene_result,
    posthoc = posthoc_result,
    effect_size = eff_size
  ))
}

# Visualization function
create_posthoc_visualization <- function(data, variable, posthoc_results) {
  # Basic boxplot
  p1 <- ggplot(data, aes_string(x=variable, y="Shannon")) +
    geom_boxplot(aes_string(fill=variable)) +
    geom_jitter(width=0.2, alpha=0.5) +
    theme_bw() +
    theme(axis.text.x = element_text(angle=45, hjust=1),
          legend.position="none") +
    labs(title=paste("Shannon Diversity by", variable),
         y="Shannon Diversity Index",
         x=variable) +
    scale_fill_viridis_d()
  
  # QQ plot
  p2 <- ggplot(data, aes(sample=Shannon)) +
    stat_qq() +
    stat_qq_line() +
    theme_bw() +
    labs(title="Normal Q-Q Plot")
  
  # Effect size plot
  if(!is.null(posthoc_results[[variable]]$effect_size)) {
    eff_data <- posthoc_results[[variable]]$effect_size
    p3 <- ggplot(eff_data, aes(x=group1, y=effsize)) +
      geom_col() +
      theme_bw() +
      theme(axis.text.x = element_text(angle=45, hjust=1)) +
      labs(title="Effect Sizes",
           x="Groups",
           y="Effect Size")
  } else {
    p3 <- NULL
  }
  
  # Arrange plots
  if(!is.null(p3)) {
    grid.arrange(p1, p2, p3, ncol=2)
  } else {
    grid.arrange(p1, p2, ncol=2)
  }
}

# Create summary table
create_summary_table <- function(posthoc_results, kruskal_results, wilcox_results) {
  summary_df <- data.frame(
    Variable = c(variables_to_test, binary_vars),
    Test_Used = c(
      sapply(posthoc_results, function(x) x$posthoc$test),
      rep("Wilcoxon", length(binary_vars))
    ),
    P_Value = c(
      sapply(kruskal_results, function(x) x$p.value),
      sapply(wilcox_results, function(x) x$p.value)
    ),
    Significant = NA
  )
  summary_df$Significant <- ifelse(summary_df$P_Value < 0.05, "*", "")
  return(summary_df)
}
