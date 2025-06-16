
#
#
# PROCESS CARD JSON FILE
# FROM CARD DB TO EXTRACT DRUG CLASSES
# AND ANTIBIOTIC CATEGORIES
# BY MATCHING ARO_NAME TO THE ARG NAME IN THE DATA SET
#
#

library(jsonlite)
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

# Load the JSON data
file_path <- "~/DB/CARD_02_2024/card.json"  # Update this to the path of your JSON file
CARD <- fromJSON(file_path)

data <- read.csv("~/ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_filtered.csv", header = TRUE, row.names = 1)

# Initialize an empty data frame to store results
results_df <- data.frame(ShortName=character(), AntibioticCategory=character(), DrugClass=character(), AMRGeneFamily=character(), stringsAsFactors=FALSE)
for(args in rownames(data)) {
  # combine results from process_CARD_models(CARD, args) in data frame
  results_df <- rbind(results_df, process_CARD_models(CARD, args))
}


# List o genes that for some reson have some small variantion in names and cant be found in CARD DB
{
  results_df <- rbind(results_df, process_CARD_models(CARD, "mef(J)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "mef(C)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(32)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(36)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(M)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(O)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(Q)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(W)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(X)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "Escherichia coli AcrAB-TolC with MarR mutations conferring resistance to ciprofloxacin and tetracycline"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "Escherichia coli AcrAB-TolC with AcrR mutation conferring resistance to ciprofloxacin, tetracycline, and ceftazidime"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "Escherichia coli ampC beta-lactamase"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(X6)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "Limosilactobacillus reuteri cat-TC"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "RSA1-1"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(X4)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "AAC(6')-Ie-APH(2'')-Ia bifunctional protein"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "ANT(3'')-II-AAC(6')-IId bifunctional protein"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "Plasmid-encoded cat (pp-cat)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "catA1"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(T)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "fosA5"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(37)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "EBR-1"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "tet(S)"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "vanR gene in vanO cluster"))
  results_df <- rbind(results_df, process_CARD_models(CARD, "AAC(6')-30/AAC(6')-Ib' bifunctional protein"))
}

#remove duplicates
results_df <- results_df[!duplicated(results_df), ]


# Fix NA values
results_df[is.na(results_df)] <- "NA"

# Save the results to a CSV file
write.csv(results_df, file = "~/ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_CARD_filtered.csv", row.names = FALSE)



# Create taxonomy table from the row names of the data

tax_t <- data.frame(ShortName = results_df$ShortName, AntibioticCategory = results_df$AntibioticCategory, DrugClass = results_df$DrugClass, AMRGeneFamily = results_df$AMRGeneFamily)
rownames(tax_t) <- tax_t$ShortName
tax_t <- sanitizeRowNames(tax_t)


#
# SEPERATE DRUG CLASSES INTO INDIVIDUAL NAMES
#

library(dplyr)
library(tidyr)

# Assuming tax_t is your data frame and DrugClass is the column of interest
unique_drug_classes <- tax_t %>%
  # Separate the DrugClass column into rows by splitting on ";"
  separate_rows(DrugClass, sep = ";") %>%
  # Select unique DrugClass names
  distinct(DrugClass) %>%
  # Pull the DrugClass column to get a vector of unique drug class names
  pull(DrugClass)

unique_drug_classes


unique_amr_families <- tax_t %>%
  # Separate the DrugClass column into rows by splitting on ";"
  separate_rows(AMRGeneFamily, sep = ";") %>%
  # Select unique DrugClass names
  distinct(AMRGeneFamily) %>%
  # Pull the DrugClass column to get a vector of unique drug class names
  pull(AMRGeneFamily)

unique_amr_families

unique_amr_categories <- tax_t %>%
  # Separate the DrugClass column into rows by splitting on ";"
  separate_rows(AntibioticCategory, sep = ";") %>%
  # Select unique DrugClass names
  distinct(AntibioticCategory) %>%
  # Pull the DrugClass column to get a vector of unique drug class names
  pull(AntibioticCategory)

unique_amr_categories

#
#
# CREATE A NEW WORKBOOK WITH 
# ONE SHEET FOR EACH DRUG CLASS
# AND WRITE THE FILTERED DATA TO THE SHEET
# 
#

library(openxlsx)


wb <- createWorkbook()
for(drug_class in unique_drug_classes) {
  # Filter rows for the current drug class
  matched_rows <- grepl(drug_class, tax_t$DrugClass)
  
  # Apply the filter to tax_t to get rows matching the condition
  filtered_tax_t <- otu_table(physeq_scaffolds)[matched_rows, ]
  colnames(filtered_tax_t) <- sample_data(physeq_scaffolds)$Sample
  
  # Create a new workbook
  sheet_name <- substr(drug_class, 1, 31)
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, filtered_tax_t, rowNames = TRUE)
  
}



# Specify file name, adjust path as necessary
file_name <- paste0("output_drugClasses.xlsx")
saveWorkbook(wb, file_name, overwrite = TRUE)



### 

## Create a new data frame where rows are drug classes and collumns are samples 
## and the values are the number of genes in the sample that are in the drug class

# Assuming  otu_table(physeq_contigs) is your filtered data frame

samples <- colnames(otu_table(physeq_scaffolds))
# Initialize 'drug_class_counts' with correct structure
drug_class_counts <- data.frame(matrix(nrow = length(unique_drug_classes), ncol = length(samples), dimnames = list(unique_drug_classes, samples)))

# Iterate over unique drug classes
for(drug_class in unique_drug_classes) {
  # Filter rows for the current drug class
  matched_rows <- grepl(drug_class, tax_t$DrugClass, fixed = TRUE)
  
  # Apply the filter to tax_t to get rows matching the condition
  filtered_tax_t <- otu_table(physeq_scaffolds)[matched_rows, ]
  
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

write.csv(drug_class_counts, file = "~/ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_CARD_drug_class_counts.csv", row.names = TRUE)

unique_amr_families_counts <- data.frame(matrix(nrow = length(unique_amr_families), ncol = length(samples), dimnames = list(unique_amr_families, samples)))

# Iterate over unique drug classes
for(drug_class in unique_amr_families) {
  # Filter rows for the current drug class
  matched_rows <- grepl(drug_class, tax_t$AMRGeneFamily, fixed = TRUE)
  
  # Apply the filter to tax_t to get rows matching the condition
  filtered_tax_t <- otu_table(physeq_scaffolds)[matched_rows, ]
  
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

write.csv(unique_amr_families_counts, file = "~/ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_CARD_amr_families_counts.csv", row.names = TRUE)



# Iterate over unique drug classes

unique_amr_categories_counts <- data.frame(matrix(nrow = length(unique_amr_categories), ncol = length(samples), dimnames = list(unique_amr_categories, samples)))

for(drug_class in unique_amr_categories) {
  # Filter rows for the current drug class
  matched_rows <- grepl(drug_class, tax_t$AntibioticCategory, fixed = TRUE)
  
  # Apply the filter to tax_t to get rows matching the condition
  filtered_tax_t <- otu_table(physeq_scaffolds)[matched_rows, ]
  
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

write.csv(unique_amr_categories_counts, file = "~/ARG_waste_water/Resistance_genes/AMR_genes_RGI_scaffolds_CARD_amr_categories_counts.csv", row.names = TRUE)


# stacked barplots with samples on X axis and drug classes on Y axis
library(ggplot2)
drug_class_counts$Drugs <- rownames(drug_class_counts)
data_melted <- melt(drug_class_counts, id.vars = "Drugs", variable.name = "Sample", value.name = "Drug_Count")
ggplot(data_melted)

ggplot(data_melted, aes(fill=Drugs, y=Drug_Count, x=Sample)) + 
  geom_bar(position="stack", stat="identity") + # verticle sample names
  theme(axis.text.x = element_text(angle = 90, hjust = 1))






# Create taxonomy table from the row names of the data
tax_t_drugs <- data.frame(Tax = rownames(drug_class_counts), ARG = rownames(drug_class_counts))
# tax_t <- sanitizeRowNames(tax_t)
rownames(tax_t_drugs) <- tax_t_drugs$Tax

TAX <- tryCatch({
  tax_table(as.matrix(tax_t_drugs[, c('Tax', 'ARG')]))
}, error = function(e) {
  stop("Failed to create taxonomy table: ", e$message)
})

tax_table(physeq_drug_class) <- tax_table(tax_t)

physeq_drug_class <- phyloseq(otu_table(select(drug_class_counts,-Drugs), taxa_are_rows = TRUE), metadta = sample_data(metadata), TAX)

head(tax_table(as.matrix(tax_t)))

