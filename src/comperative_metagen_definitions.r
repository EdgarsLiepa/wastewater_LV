## Function Definitions


# Function to read data and return a melted data frame
read_ARG_data <- function(file_location) {
  arg_t <- read.csv(file_location, header = TRUE, row.names = 1)
#  arg_t$ARG <- rownames(arg_t) # Add a new column "ARG" with row names
  return(arg_t)
}

# Function to plot antibiotic resistance genes (ARGs) from a given ARG table.
#
# Parameters:
#   ARG_table: A data frame or matrix with rows representing different ARGs and
#              columns potentially including counts, categories, or other relevant data.
#   plot_title: Optional. A character string specifying the title of the plot. Default
#               is "Antibiotic Resistance Genes".
#
# This function will generate a plot based on the provided ARG_table. The specific
# plotting details (e.g., type of plot, aesthetics) should be defined within the function
# body based on the structure and content of ARG_table.

# To use the plot with specific dimensions, wrap it in ggsave:
# plot <- plot_ARGs(amr_categories_counts, "AMR Categories")
# ggsave("amr_categories.pdf", plot, width = 15, height = 8)

# Or when printing to screen, wrap in a specified viewing device:
# dev.new(width = 15, height = 8)
# print(plot)

plot_ARGs <- function(ARG_table, plot_title = "ARGs counts", y_name = "Counts") {
  # Input validation
  if (!is.data.frame(ARG_table)) {
    stop("ARG_table must be a data frame")
  }
  
  # Create ARG column from rownames
  ARG_table$ARG <- rownames(ARG_table)
  
  # Melt the data
  data_melted <- reshape2::melt(ARG_table, 
                                id.vars = "ARG", 
                                variable.name = "Sample", 
                                value.name = "Counts")
  
  # Calculate totals and get top ARGs
  arg_totals <- data_melted %>%
    dplyr::group_by(ARG) %>%
    dplyr::summarise(Total = sum(Counts, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(Total))
  
  # Select top 10 ARGs
  top_args <- head(arg_totals$ARG, 10)
  
  # Create 'other' category
  data_melted <- data_melted %>%
    dplyr::mutate(ARG = ifelse(ARG %in% top_args, ARG, 'Other'))
  
  # Order samples by total counts
  sample_totals <- data_melted %>%
    dplyr::group_by(Sample) %>%
    dplyr::summarise(Total = sum(Counts, na.rm = TRUE)) %>%
    dplyr::arrange(dplyr::desc(Total))
  
  data_melted$Sample <- factor(data_melted$Sample, 
                               levels = sample_totals$Sample)
  
  # Create color palette with enough colors
  n_colors <- length(unique(data_melted$ARG))
  custom_colors <- colorRampPalette(RColorBrewer::brewer.pal(8, "Set2"))(n_colors)
  
  # Create the plot
  ggplot_obj <- ggplot(data_melted, aes(x = Sample, y = Counts, fill = ARG)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = custom_colors) +
    theme_bw() +
    theme(
      # Adjust plot dimensions
      aspect.ratio = 0.5,
      
      # Improve axis text
      axis.text.x = element_text(
        angle = 90,
        hjust = 1,
        vjust = 0.5,
        size = 10,
        face = "plain"
      ),
      axis.text.y = element_text(size = 10),
      
      # Improve axis titles
      axis.title = element_text(size = 12, face = "bold"),
      
      # Improve legend
      legend.position = "bottom",
      legend.box = "horizontal",
      legend.margin = margin(t = 10, b = 10),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 12, face = "bold"),
      legend.key.size = unit(0.8, "cm"),
      
      # Improve plot title
      plot.title = element_text(
        size = 14,
        face = "bold",
        hjust = 0.5,
        margin = margin(b = 20)
      ),
      
      # Improve grid lines
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      
      # Add some margins
      plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
    ) +
    # Improve legend layout
    guides(fill = guide_legend(
      ncol = 4,
      byrow = TRUE,
      title.position = "top",
      title.hjust = 0.5,
      label.position = "right",
      keywidth = 1,
      keyheight = 1
    )) +
    labs(
      x = "Sample",
      y = y_name,
      fill = "ARG",
      title = plot_title
    )
  
  return(ggplot_obj)
}

plotRareCurveForEachSample <- function(ARG) {
  # Number of samples
  num_samples <- nrow(ARG)
  
  # Plot a rarefaction curve for each sample
  for (i in 1:num_samples) {
    rarecurve(ARG[i, , drop = FALSE], step = 10, sample = TRUE, col = "blue", cex = 0.7,
              xlab = "Number of ARG sequences", ylab = "Number of ARG types", 
              main = rownames(ARG)[i], label = FALSE)
  }
}

plotRareCurve <- function(ARG_data) {
  # Transpose the data
  ARG <- ARG_data[, !(names(ARG_data) %in% c("Total", "Mean", "ARG"))]
  
  # Transpose the data
  # Note: If ARG_data is a data.table, you may need to convert it to a matrix for t() to work as expected
  ARG <- t(as.matrix(ARG))
  
  # Generate rarefaction curve for the overall dataset
  rarecurve(ARG, step = 10, sample = TRUE, col = "blue", cex = 0.7,
            xlab = "Number of ARG sequences", ylab = "Number of ARG types", label = FALSE)
  
  # Call the sub-function to plot rarefaction curves for each sample
  plotRareCurveForEachSample(ARG)
}


### Phyloseq object creation


# Helper function to replace special characters in row names
sanitizeRowNames <- function(data) {
  rownames(data) <- gsub("''", "_doubleprime", rownames(data))
  rownames(data) <- gsub("'", "_prime", rownames(data))
  return(data)
}

createPhyloseqObject <- function(amr_data, metadata_path) {
  
  # Read and process AMR genes data
  #  data <- tryCatch({
  #    read.csv(amr_genes_path, header = TRUE, row.names = 1)
  #  }, error = function(e) {
  #    stop("Failed to read AMR genes data: ", e$message)
  #  })
  
  amr_data <- sanitizeRowNames(amr_data)
  
  # Create taxonomy table from the row names of the data
  tax_t <- data.frame(Tax = rownames(amr_data), ARG = rownames(amr_data))
  # tax_t <- sanitizeRowNames(tax_t)
  rownames(tax_t) <- tax_t$Tax
  
  TAX <- tryCatch({
    tax_table(as.matrix(tax_t[, c('Tax', 'ARG')]))
  }, error = function(e) {
    stop("Failed to create taxonomy table: ", e$message)
  })
  
  # Read and process metadata
  metadata <- tryCatch({
    read.csv(metadata_path, header = TRUE, colClasses = c(Date = "character", Dairy_farming="character", Meat_production="character", Metal_processing="character", Washrooms="character"))
  }, error = function(e) {
    stop("Failed to read metadata: ", e$message)
  })
  metadata$Date <- as.Date(metadata$Date, format = "%d.%m")
  if (any(is.na(metadata$Date))) {
    stop("Failed to convert some dates. Please check the date format.")
  }
  
  # Filter samples based on metadata
  samples_to_keep <- colnames(amr_data)[colnames(amr_data) %in% rownames(metadata)]
  data_filtered <- amr_data[, samples_to_keep]
  
  # Create phyloseq object
  physeq <- tryCatch({
    phyloseq(otu_table(data_filtered, taxa_are_rows = TRUE), sample_data(metadata), TAX)
  }, error = function(e) {
    stop("Failed to create phyloseq object: ", e$message)
  })
  
  # Convert Date to factor
  sample_data(physeq)$Date <- as.factor(sample_data(physeq)$Date)
  
  return(physeq)
}

visualizePhyseqData <- function(physeq) {
  # Set global options and knitr options for plot dimensions
  options(width = 100)
  if(exists("knitr::opts_chunk")) {
    knitr::opts_chunk$set(fig.height = 6, fig.width = 9)
  }
  
  # Basic Composition Barplot
  p1 <- physeq %>%
    comp_barplot(tax_level = 'ARG') +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  print(p1)
  
  # Group Average Composition Barplot
  p2 <- physeq %>%
    comp_barplot(tax_level = 'ARG', n_taxa = 30, bar_width = 0.8, bar_outline_colour = NA) +
    coord_flip() + labs(x = NULL, y = NULL)
  print(p2)
  
  # Faceting by Taxonomic Composition by Samples
  if ("City" %in% colnames(sample_data(physeq))) {
    p3 <- physeq %>%
      ps_mutate(Group = factor(sample_data(physeq)$City)) %>%
      comp_barplot(tax_level = "ARG", n_taxa = 10, sample_order = "bray", bar_outline_colour = NA) +
      facet_grid(rows = vars(Group), scales = "free", space = "free") +
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "bottom")+
      coord_flip()
    print(p3)
  } else {
    warning("Column 'City' not found in sample_data. Skipping faceted plot.")
  }
}

visualizePhyseqRichnessAndBar <- function(physeq) {
  # Check if 'City' metadata is available
  if(!"City" %in% colnames(sample_data(physeq))) {
    stop("Column 'City' not found in sample_data. Ensure 'City' metadata is present.")
  }
  
  sample_data(physeq)$City <- factor(sample_data(physeq)$City, levels = sort(unique(sample_data(physeq)$City)))
  
  # Composition Bar Plot
  p1 <- plot_bar(physeq, fill = "City") + 
    geom_bar(aes(color = City, fill = City), stat = "identity", position = "stack") +
    geom_bar(stat = "identity", colour = "black") + 
    scale_y_continuous(expand = c(0, 0)) +
    theme_minimal() + coord_flip() # Use a minimal theme for better readability
  
  print(p1)
  
  # Richness Plot
  p2 <- plot_richness(physeq, measures = c("Observed", "Chao1", "Shannon", "Simpson")) +
    geom_boxplot() +
    theme_minimal() + # Consistent theme with the bar plot
    theme(legend.position = "bottom", axis.text.x = element_text(angle = 45, hjust = 1)) # Improve x-axis label readability 
  print(p2)
}


### Normalization



# Description:
#
# Function to process Kraken reports in base path.
# Sample ID is taken from metadata file and report file 
# is expected to be in the format: <sample_id>_k2_report.txt
#
# INPUTS: 
#
# metadata_path: Path to the metadata file
# reports_base_path: Base path to the Kraken reports
#
# OUTPUTS:
# 
# kraken2_table: A tabble containing the bacterial counts for each sample
#


process_kraken_reports <- function(metadata, reports_base_path) {
  # Initialize the results table as a data frame
  kraken2_table <- data.frame(Sample = character(), BacterialReads_Kraken = double(), stringsAsFactors = FALSE)
  
  # Iterate over each sample in the metadata
  for (i in rownames(metadata)) {
    report_path <- paste0(reports_base_path, i, "_k2_report.txt")
    
    # Check if the report file exists
    if (!file.exists(report_path)) {
      message("Report file does not exist for sample: ", i)
      next # Skip to the next iteration
    }
    
    # Attempt to read the Kraken report
    report <- tryCatch({
      read.table(report_path, header = FALSE, sep = "\t", col.names = c("percentage", "reads", "reads_in_clade", "rank", "taxid", "name"))
    }, error = function(e) {
      message("Error reading report for sample: ", i, "; Error: ", e$message)
      return(NULL)
    })
    
    # Skip this iteration if report couldn't be read
    if (is.null(report)) next
    
    # Filter for bacteria rows
    bacteria <- report[grepl("Bacteria", report$name), ]
    
    # Check if bacteria data is available
    if (nrow(bacteria) > 0) {
      # Add results to results table
      kraken2_table <- rbind(kraken2_table, data.frame(Sample = i, BacterialReads_Kraken = bacteria$reads[1], stringsAsFactors = FALSE))
    } else {
      message("No bacterial reads found for sample: ", i)
    }
  }
  
  # Return the results table
  return(kraken2_table)
}


# Example usage of the function
# Replace 'path/to/your/metadata.csv' and 'path/to/kraken/reports/' with actual paths
# kraken2_results <- process_kraken_reports("path/to/your/metadata.csv", "path/to/kraken/reports/")
# print(kraken2_results)


# g_length Function Description
#
# This function calculates the average gene length for a specified gene within a sample,
# based on RGI analysis results. It utilizes the Resistance Gene Identifier (RGI) output,
# which identifies antimicrobial resistance (AMR) genes.
#
# Parameters:
#   - gene: The gene of interest, specified by its ARO (Antibiotic Resistance Ontology) name.
#   - sample: The sample identifier. This parameter is included for consistency and future use
#             but is not currently utilized within the function.
#   - rgi_rez: A data frame containing the RGI analysis results. It is expected to have at least
#              'Best_Hit_ARO', 'Stop', and 'Start' columns, where 'Best_Hit_ARO' identifies the
#              gene, 'Stop' and 'Start' define the gene's position on the genome.
#
# The function filters 'rgi_rez' to find rows corresponding to the specified 'gene',
# then calculates the gene length as the difference between the 'Stop' and 'Start' positions.
# It returns the average length of all identified instances of the gene in the dataset.
# This average length is useful for normalizing gene expression levels in TPM calculations.
#
# Returns:
#   - The average length of the specified gene across all instances found in 'rgi_rez'.

g_length <- function(gene, sample, rgi_rez) {
  
  # find in CARD by aro name
  gene_length <- 0
  
  # find all rows with gene name
  gene_rows <- rgi_rez[rgi_rez$Best_Hit_ARO == gene,]
  
  # calculate length by Stop - Start
  gene_length <- gene_rows$Stop - gene_rows$Start
  
  # return mean length 
  return(mean(gene_length))
  
}

# calculate_TPM_RGI Function Description
#
# This function calculates the Transcripts Per Million (TPM) for genes identified in Resistance
# Gene Identifier (RGI) analysis across multiple samples. TPM is a normalization method for
# gene expression data that accounts for gene length and sequencing depth, allowing for the
# comparison of gene expression levels across different samples or conditions.
#
# Parameters:
#   - gene_t: A data frame or matrix where rows represent genes identified by RGI analysis
#             and columns represent samples. The values in this data frame should be raw
#             gene counts or another measure of gene expression prior to TPM normalization.
#   - basepath: A string that specifies the base directory path where RGI result files are
#               stored. These files should be named in a manner that matches the sample
#               names in 'gene_t' and suffixed with "_rgi.txt". Each file contains detailed
#               RGI analysis results, including gene lengths necessary for TPM calculation.
#
# The function iterates over each sample (column) in 'gene_t', loads the corresponding RGI
# result file, and calculates gene lengths. It then normalizes the raw gene counts by gene
# length to get Reads Per Kilobase (RPK) values, replaces any NaN values resulting from
# division by zero with 0, and sums these RPK values to adjust for sequencing depth. Finally,
# it calculates TPM values by dividing RPK values by the per-sample sum of RPKs and scaling
# by 1e6. The normalized TPM values replace the original values in 'gene_t'.
#
# Returns:
#   - A data frame or matrix equivalent to 'gene_t' but with values transformed to TPM, thus
#     facilitating the comparison of gene expression levels across samples.
#
# Note: This function requires the 'g_length' function to calculate gene lengths from RGI
# result files. Ensure 'g_length' is defined in the environment before calling this function.

calculate_TPM_RGI <- function(gene_t, basepath) {
  
  gene_t <- gene_t[, !(names(gene_t) %in% c("ARG"))]
  
  for (sample in colnames(gene_t)) {
    gene_lengths <- c()
    # Load sample table from RGI results
    print(paste(basepath, sample, "_rgi.txt", sep = ""))
    rgi_rez <- read.csv(paste(basepath, sample, "_rgi.txt", sep = ""), sep = "\t", row.names = "Contig", header = TRUE)
    
    for (gene in rownames(gene_t)) {
      # Get the length of the gene length per kilobase
      gene_len <- g_length(gene, sample, rgi_rez) / 1000
      
      # Add sample vector
      gene_lengths <- c(gene_lengths, gene_len)
    }
    
    # Divide sample column in gene_t by gene_lengths
    sample_row <- gene_t[,sample] / gene_lengths
    
    # Replace NaN values with 0
    sample_row[is.nan(sample_row)] <- 0
    
    # Calculate the per-sample sum of RPKs
    sum_RPK <- sum(sample_row)
    
    # Finally, divide each column of RPK by the sum of RPKs and multiply by 1e6 to get TPM
    gene_t[,sample] <- sample_row / sum_RPK * 1e6
  }
  
  return(gene_t)
}

### PcOA plot

plot_pcoa <- function(physeq, distance_method, color_variable) {
  # Calculate distances
  dist <- distance(physeq, method=distance_method)
  
  # PCoA
  pcoa <- ordinate(physeq, method="PCoA", distance=dist)
  
  # Plot PCoA
  p <- plot_ordination(physeq, pcoa, color=color_variable) + 
    geom_point(size=3) +
    theme_minimal()
  
  # connect points by city
  p <- p + stat_chull(ggplot2::aes(colour = City))
  
  # add descriptive labels
  p <- p + labs(title=paste("PCoA on", distance_method, "distances"))
  
  
  return(p)
}


### CARD DB

process_CARD_models <- function(CARD, target_ARO_name) {
  # Initialize an empty data frame to store results
  results_df <- data.frame(ShortName=character(), AntibioticCategory=character(), DrugClass=character(), stringsAsFactors=FALSE)
  
  # Define a function to process each model based on the ARO_name
  for (model in CARD) {
    # Check if model size is bigger than 1 and matches the target ARO_name
    if (length(model) > 2 && model$ARO_name == target_ARO_name) {
      antibioticCategory <- NA
      drugClass <- NA
      amrGeneFamily <-NA
      
      # Iterate over ARO_category to find categories
      for(aro_id in names(model$ARO_category)) {
        category <- model$ARO_category[[aro_id]]
        if(category$category_aro_class_name == "Antibiotic") {
          #antibioticCategory <- category$category_aro_name
          # apped string to antibioticCategory
          if(is.na(antibioticCategory)) {
            antibioticCategory <- category$category_aro_name
          } else {
            antibioticCategory <- paste(antibioticCategory, category$category_aro_name, sep = ";")
          }
        }
        if(category$category_aro_class_name == "Drug Class") {
          if(is.na(drugClass)) {
            drugClass <- category$category_aro_name
          } else {
            drugClass <- paste(drugClass, category$category_aro_name, sep = ";")
          }
        }
        if(category$category_aro_class_name == "AMR Gene Family") {
          if(is.na(amrGeneFamily)) {
            amrGeneFamily <- category$category_aro_name
          } else {
            amrGeneFamily <- paste(amrGeneFamily, category$category_aro_name, sep = ";")
          }
        }
      }
      
      # Append to results data frame  
      results_df <- rbind(results_df, data.frame(ShortName=model$CARD_short_name, AntibioticCategory=antibioticCategory, DrugClass=drugClass, AMRGeneFamily=amrGeneFamily, stringsAsFactors=FALSE))
      
      # Break the loop after finding the target ARO_name
      break
    }
  }
  
  # if data frame is empty, return NA
  if (nrow(results_df) == 0) {
    print(paste("No results found for", target_ARO_name))
  }
  
  return(results_df)
}

process_card_data <- function(arg_scaffolds_t, gene_list=NULL) {
  
  # Load the JSON data
  CARD <- fromJSON(CARDDB)
  
  # Initialize an empty data frame for storing results
  results_df <- data.frame(
    ShortName = character(),
    AntibioticCategory = character(),
    DrugClass = character(),
    AMRGeneFamily = character(),
    stringsAsFactors = FALSE
  )
  
  # Process each gene in the arg_scaffolds_t
  for (args in rownames(arg_scaffolds_t)) {
    results_df <- rbind(results_df, process_CARD_models(CARD, args))
  }
  
  # Handle specific genes that require manual adjustment
  if (!is.null(gene_list)) {
    for (gene in gene_list) {
      results_df <- rbind(results_df, process_CARD_models(CARD, gene))
    }
  }
  
  # Fix NA values
  results_df[is.na(results_df)] <- "NA"
  
  # Create taxonomy table
  card_tax_t <- data.frame(
    ShortName = results_df$ShortName,
    AntibioticCategory = results_df$AntibioticCategory,
    DrugClass = results_df$DrugClass,
    AMRGeneFamily = results_df$AMRGeneFamily
  )
  
  # remove duplicates
  card_tax_t <- card_tax_t[!duplicated(card_tax_t$ShortName), ]
  
  # Set row names
  rownames(card_tax_t) <- card_tax_t$ShortName
  
  return(card_tax_t)
}

create_drug_class_counts <- function(otu_table, tax_t, unique_drug_classes) {
  
  samples <- colnames(otu_table)
  # Initialize 'drug_class_counts' with correct structure
  drug_class_counts <- data.frame(matrix(nrow = length(unique_drug_classes), ncol = length(samples), dimnames = list(unique_drug_classes, samples)))
  
  # Iterate over unique drug classes
  for(drug_class in unique_drug_classes) {
    # Filter rows for the current drug class
    matched_rows <- grepl(drug_class, tax_t$DrugClass, fixed = TRUE)
    
    # Apply the filter to tax_t to get rows matching the condition
    filtered_tax_t <- otu_table[matched_rows, ]
    
    current_sum <- colSums(filtered_tax_t)
    
    # add current_sum to drug_class_counts
    df <- t(data.frame(current_sum))
    
    # Set row names separately
    rownames(df) <- drug_class
    
    # Assuming 'drug_class' is the row name and it exists in drug_class_counts
    row_to_update <- which(rownames(drug_class_counts) == drug_class)
    
    # Assuming df contains only one row and you want to update drug_class_counts with its values
    drug_class_counts[row_to_update, ] <- df
  }
  
  return(drug_class_counts)
}

create_amr_families_counts <- function(otu_table, tax_t, unique_amr_families) {
  
  samples <- colnames(otu_table)
  
  # Initialize 'unique_amr_families_counts' with correct structure
  unique_amr_families_counts <- data.frame(matrix(nrow = length(unique_amr_families), ncol = length(samples), dimnames = list(unique_amr_families, samples)))
  
  # Iterate over unique drug classes
  for(drug_class in unique_amr_families) {
    # Filter rows for the current drug class
    matched_rows <- grepl(drug_class, tax_t$AMRGeneFamily, fixed = TRUE)
    
    # Apply the filter to tax_t to get rows matching the condition
    filtered_tax_t <- otu_table[matched_rows, ]
    
    current_sum <- colSums(filtered_tax_t)
    
    # add current_sum to drug_class_counts
    df <- t(data.frame(current_sum))
    
    # Set row names separately
    rownames(df) <- drug_class
    
    # Assuming 'drug_class' is the row name and it exists in drug_class_counts
    row_to_update <- which(rownames(unique_amr_families_counts) == drug_class)
    
    # Assuming df contains only one row and you want to update drug_class_counts with its values
    unique_amr_families_counts[row_to_update, ] <- df
  }
  
  return(unique_amr_families_counts)
}

create_amr_categories_counts <- function(otu_table, tax_t, unique_amr_categories) {
  
  samples <- colnames(otu_table)
  
  # Initialize 'unique_amr_categories_counts' with correct structure
  unique_amr_categories_counts <- data.frame(matrix(nrow = length(unique_amr_categories), ncol = length(samples), dimnames = list(unique_amr_categories, samples)))
  
  # Iterate over unique drug classes
  for(drug_class in unique_amr_categories) {
    # Filter rows for the current drug class
    matched_rows <- grepl(drug_class, tax_t$AntibioticCategory, fixed = TRUE)
    
    # Apply the filter to tax_t to get rows matching the condition
    filtered_tax_t <- otu_table[matched_rows, ]
    
    current_sum <- colSums(filtered_tax_t)
    
    # add current_sum to drug_class_counts
    df <- t(data.frame(current_sum))
    
    # Set row names separately
    rownames(df) <- drug_class
    
    # Assuming 'drug_class' is the row name and it exists in drug_class_counts
    row_to_update <- which(rownames(unique_amr_categories_counts) == drug_class)
    
    # Assuming df contains only one row and you want to update drug_class_counts with its values
    unique_amr_categories_counts[row_to_update, ] <- df
  }
  
  return(unique_amr_categories_counts)
}

plot_stacked_bar <- function(drug_class_counts) {
  # Create a melted version of the data
  data_melted <- melt(drug_class_counts, id.vars = "Drugs", variable.name = "Sample", value.name = "Drug_Count")
  
  # Create a stacked bar plot
  ggplot(data_melted, aes(fill=Drugs, y=Drug_Count, x=Sample)) + 
    geom_bar(position="stack", stat="identity") + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

create_taxonomy_table <- function(tax_t) {
  # Create taxonomy table from the row names of the data
  tax_t_drugs <- data.frame(Tax = rownames(tax_t), ARG = rownames(tax_t))
  rownames(tax_t_drugs) <- tax_t_drugs$Tax
  
  return(tax_t_drugs)
}

create_phyloseq <- function(tax_t, sample_data, metadata) {
  # Create phyloseq object
  physeq_drug_class <- phyloseq(otu_table(select(tax_t,-Drugs), taxa_are_rows = TRUE), metadta = sample_data(metadata), TAX)
  
  return(physeq_drug_class)
}



# Create a function to compare normalization methods
compare_normalizations <- function(ps_object) {
  # Original counts
  raw_counts <- as.matrix(otu_table(ps_object))
  
  # CSS normalization
  ps_css <- microbiomeMarker::normalize(ps_object, method="CSS")
  css_counts <- as.matrix(otu_table(ps_css))
  
  # TSS normalization
  ps_tss <- transform_sample_counts(ps_object, function(x) x/sum(x))
  tss_counts <- as.matrix(otu_table(ps_tss))
  
  # CLR normalization
  ps_clr <- clr(otu_table(ps_object))
  clr_counts <- as.matrix(ps_clr)
  
  # Compare distributions
  par(mfrow=c(2,2))
  boxplot(raw_counts, main="Raw counts", las=2)
  boxplot(css_counts, main="CSS normalized", las=2)
  boxplot(tss_counts, main="TSS normalized", las=2)
  boxplot(clr_counts, main="CLR normalized", las=2)
  
  # Compare coefficient of variation
  cv <- function(x) sd(x)/mean(x)
  cv_stats <- data.frame(
    Raw = apply(raw_counts, 1, cv),
    CSS = apply(css_counts, 1, cv),
    TSS = apply(tss_counts, 1, cv),
    CLR = apply(clr_counts, 1, cv)
  )
  
  return(cv_stats)
}


library(DESeq2)
library(ALDEx2)

# DESeq2 analysis
run_deseq2 <- function(ps_object, variable) {
  # Convert to DESeq2 object
  dds <- phyloseq_to_deseq2(ps_object, design = as.formula(paste("~", variable)))
  
  # Run DESeq2
  dds <- DESeq(dds)
  res <- results(dds)
  
  # Get significant results
  sig_res <- res[which(res$padj < 0.05), ]
  return(sig_res)
}

# ALDEx2 analysis
run_aldex2 <- function(ps_object, variable) {
  # Prepare data
  otu_table <- as.matrix(otu_table(ps_object))
  meta <- sample_data(ps_object)
  
  # Run ALDEx2
  aldex_clr <- aldex.clr(otu_table, meta[[variable]], mc.samples=128)
  aldex_effect <- aldex.effect(aldex_clr)
  aldex_test <- aldex.ttest(aldex_clr)
  
  return(list(effect=aldex_effect, test=aldex_test))
}



library(igraph)
library(SpiecEasi)

# Create co-occurrence network
create_arg_network <- function(ps_object, correlation_threshold=0.6) {
  # Get abundance matrix
  abund_mat <- as.matrix(otu_table(ps_object))
  
  # Calculate correlations
  cor_mat <- cor(t(abund_mat), method="spearman")
  
  # Apply threshold
  cor_mat[abs(cor_mat) < correlation_threshold] <- 0
  
  # Create network
  net <- graph_from_adjacency_matrix(cor_mat, 
                                     mode="undirected", 
                                     weighted=TRUE, 
                                     diag=FALSE)
  
  # Calculate network properties
  centrality <- degree(net)
  betweenness <- betweenness(net)
  clusters <- cluster_fast_greedy(net)
  
  return(list(network=net, 
              centrality=centrality, 
              betweenness=betweenness, 
              clusters=clusters))
}
