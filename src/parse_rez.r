dat <- SCP::srt_to_adata(seurat_object)

dat$write_h5ad("CO0407_umap_obj.h5ad")



###############
library(sceasy)
library(reticulate)
loompy <- reticulate::import('loompy')

###############


  arg_t <- read.csv("ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_new.csv", header = TRUE, row.names = 1)
  summary(arg_t)
  # count all the reads
  
  
  
  total_reads <- rowSums(data_normalized_k2)
  total <- sum(data_normalized_k2)
  length(rownames(data_normalized_k2))
  data_normalized_k2$SUM <- total_reads # Add a new column "SUM" with total reads
  
  data_normalized_k2$ARG <- rownames(data_normalized_k2) # Add a new column "ARG" with row names
 
  
  # print top 10 args
  top10 <- rowSums(head(arrange(total_reads, desc(total_reads)), 10))
  
  # add to dataframe top10
  data <- data.frame(ARG = "Top 10", SUM = top10)
  
  
  # remove collumn by name
  arg_t <- arg_t[, !names(arg_t) %in% c("SUM")]
  data_normalized_k2 <- data_normalized_k2[, !names(data_normalized_k2) %in% c("ARG")]
  
  # Create the pie chart
  pie(data_normalized_k2$SUM, labels = paste(names(data_normalized_k2$SUM), "(", round(data_normalized_k2$SUM / sum(sample_data) * 100, 1), "%)", sep = ""), main = sample_id)
  
  #
  #
  #
  #
  
  amr_data <- arg_deep_t
  
  sanitizeRowNames <- function(data) {
    rownames(data) <- gsub("''", "_doubleprime", rownames(data))
    rownames(data) <- gsub("'", "_prime", rownames(data))
    return(data)
  }
  
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
  
  
  # Filter samples based on metadata
  samples_to_keep <- colnames(amr_data)[colnames(amr_data) %in% rownames(metadata)]
  data_filtered <- data[, samples_to_keep]
  
  # Create phyloseq object
  physeq <- tryCatch({
    phyloseq(otu_table(data_filtered, taxa_are_rows = TRUE), sample_data(metadata))
  }, error = function(e) {
    stop("Failed to create phyloseq object: ", e$message)
  })
  
  # Convert Date to factor
  sample_data(physeq)$Date <- as.factor(sample_data(physeq)$Date)
  
  
  