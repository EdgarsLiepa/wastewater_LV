
library(dplyr)
library(readr)
library(stringr)

convert_funcscan_harmonized_to_gtf <- function(harmonized_file, output_gtf, 
                                               feature_type = "CDS",
                                               source = "funcscan",
                                               min_identity = 90) {

  
  # Read harmonized data
  data <- read_tsv(harmonized_file, 
                   show_col_types = FALSE,
                   quote = "\"",
                   locale = locale(encoding = "UTF-8")) %>%
    filter(sequence_identity >= min_identity)
  
  cat("Loaded", nrow(data), "annotations from funcscan harmonization\n")
  
  # Validate required columns for your format
  required_cols <- c("input_sequence_id", "input_gene_start", "input_gene_stop", 
                     "strand_orientation", "gene_symbol")
  missing_cols <- setdiff(required_cols, colnames(data))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Process the data
  gtf_processed <- data %>%
    # Filter valid entries
    filter(
      !is.na(input_sequence_id),
      !is.na(input_gene_start),
      !is.na(input_gene_stop),
      input_gene_start <= input_gene_stop
    ) %>%
    # Create GTF-compatible columns
    mutate(
      # Extract sequence ID (remove everything after | if present)
      #sequence_id = str_extract(input_sequence_id, "^[^|]+"),
      sequence_id = input_sequence_id,
      
      # Use input gene coordinates
      start = as.integer(input_gene_start),
      end = as.integer(input_gene_stop),
      
      # Normalize strand
      strand = case_when(
        strand_orientation == "+" ~ "+",
        strand_orientation == "-" ~ "-",
        TRUE ~ "."
      ),
      
      # Clean gene identifiers
      gene_id = str_replace_all(gene_symbol, '[";\\s]', "_"),
      
      # Create transcript ID
      transcript_id = paste0(gene_id, "_", row_number()),
      
      # Clean gene name/product
      product = if ("gene_name" %in% colnames(.)) {
        str_replace_all(gene_name, '[";]', "'")
      } else {
        gene_symbol
      },
      
      # Clean drug class
      resistance_class = if ("drug_class" %in% colnames(.)) {
        str_replace_all(drug_class, '[";]', "'")
      } else {
        NA_character_
      },
      
      # Detection tool
      detection_tool = if ("analysis_software_name" %in% colnames(.)) {
        analysis_software_name
      } else {
        "unknown"
      },
      
      # Add quality metrics
      identity = if ("sequence_identity" %in% colnames(.)) {
        round(as.numeric(sequence_identity), 2)
      } else {
        NA_real_
      },
      
      coverage = if ("coverage_percentage" %in% colnames(.)) {
        round(as.numeric(coverage_percentage), 2)
      } else {
        NA_real_
      }
    )
  
  # Build GTF attributes
    gtf_final <- gtf_processed %>%
    rowwise() %>%
    mutate(
      attributes = {
        # Core GTF attributes
        attrs <- c(
          paste0('gene_id "', gene_id, '"'),
          paste0('transcript_id "', transcript_id, '"'),
          paste0('gene_name "', gene_id, '"')
        )
        
        # Add product/gene name
        if (!is.na(product) && product != "") {
          attrs <- c(attrs, paste0('product "', product, '"'))
        }
        
        # Add resistance class
        if (!is.na(resistance_class) && resistance_class != "") {
          attrs <- c(attrs, paste0('resistance_class "', resistance_class, '"'))
        }
        
        # Add detection tool
        if (!is.na(detection_tool) && detection_tool != "") {
          attrs <- c(attrs, paste0('detection_tool "', detection_tool, '"'))
        }
        
        # Add quality metrics
        if (!is.na(identity)) {
          attrs <- c(attrs, paste0('sequence_identity "', identity, '"'))
        }
        
        if (!is.na(coverage)) {
          attrs <- c(attrs, paste0('coverage_percentage "', coverage, '"'))
        }
        
        # Add reference information if available
        if ("reference_database_name" %in% colnames(gtf_processed) && 
            !is.na(reference_database_name)) {
          attrs <- c(attrs, paste0('reference_db "', reference_database_name, '"'))
        }
        
        if ("reference_accession" %in% colnames(gtf_processed) && 
            !is.na(reference_accession)) {
          attrs <- c(attrs, paste0('reference_accession "', reference_accession, '"'))
        }
        
        # Add genetic variation type
        if ("genetic_variation_type" %in% colnames(gtf_processed) && 
            !is.na(genetic_variation_type)) {
          attrs <- c(attrs, paste0('variation_type "', genetic_variation_type, '"'))
        }
        
        paste(attrs, collapse = "; ") %>% paste0(";")
      }
    ) %>%
    ungroup() %>%
    # Create final GTF format
    transmute(
      seqname = sequence_id,
      source = detection_tool,
      feature = feature_type,
      start = start,
      end = end,
      score = ".",
      strand = strand,
      frame = "0",
      attributes = attributes
    )
  
  # Write GTF file
  cat("Writing GTF file...\n")
  
  # Write header with metadata
  header_lines <- c(
    "##gtf-version 2.2",
    paste0("##source ", source),
    paste0("##date ", Sys.Date()),
    "##feature-type antimicrobial_resistance_genes"
  )
  
  writeLines(header_lines, output_gtf)
  
  # Write data
  write.table(gtf_final, output_gtf, append = TRUE, col.names = FALSE, 
              row.names = FALSE, quote = FALSE, sep = "\t")
  
  # Summary statistics
  cat("Conversion summary:\n")
  cat("- Total AMR features:", nrow(gtf_final), "\n")
  cat("- Unique sequences:", length(unique(gtf_final$seqname)), "\n")
  cat("- Unique genes:", length(unique(gtf_processed$gene_id)), "\n")
  
  # Drug class summary
  if ("resistance_class" %in% colnames(gtf_processed)) {
    drug_summary <- gtf_processed %>%
      filter(!is.na(resistance_class), resistance_class != "") %>%
      count(resistance_class, sort = TRUE, name = "count")
    
    if (nrow(drug_summary) > 0) {
      cat("- Drug classes:\n")
      for (i in 1:min(10, nrow(drug_summary))) {
        cat("  ", drug_summary$resistance_class[i], ":", drug_summary$count[i], "\n")
      }
    }
  }
  
  # Detection tool summary
  if ("detection_tool" %in% colnames(gtf_processed)) {
    tools_summary <- gtf_processed %>%
      filter(!is.na(detection_tool), detection_tool != "") %>%
      count(detection_tool, sort = TRUE, name = "count")
    
    if (nrow(tools_summary) > 0) {
      cat("- Detection tools:\n")
      for (i in 1:nrow(tools_summary)) {
        cat("  ", tools_summary$detection_tool[i], ":", tools_summary$count[i], "\n")
      }
    }
  }
  
  return(gtf_final)
}



# Convert funcscan harmonization results
gtf_data <- convert_funcscan_harmonized_to_gtf(
  harmonized_file = "~/ARG_waste_water/results/hamronization_MAGs_report.tsv",
  output_gtf = "~/ARG_waste_water/results/mags_funcascan.gtf",
  source = "funcscan"
)


