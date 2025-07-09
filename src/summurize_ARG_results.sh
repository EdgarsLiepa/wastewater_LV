#!/bin/bash

# Directory containing the .txt files
dirRGI="/home/edgars.liepa/NAS/bioinfo/edgars.liepa/Antibiotic_resistance/ww_publication/resistance_genes_RGI_contigs"


mkdir -p ${dirRGI}/hAMRonization


# hAMRonize the RGI results
for file in "$dirRGI"/*.txt; do   
    file_id=$(basename "$file" .txt);     
    file_id=${file_id%_rgi};      
    hamronize rgi --analysis_software_version 6.0.3 --reference_database_version 3.2.8 \
    --input_file_name "$file_id" "$file" > ${dirRGI}/hAMRonization/${file_id}_rgi.txt ; 
done


dirDeepARG="/home/edgars.liepa/NAS/bioinfo/edgars.liepa/Antibiotic_resistance/ww_publication/ARG_deepARG-ss"


mkdir -p ${dirDeepARG}/hAMRonization

# hAMRonize the DeepARG results
for file in "$dirDeepARG"/*.ARG; do   
    file_id=$(basename "$file" .ARG);     
    file_id=${file_id%.clean.deeparg.mapping};      
    hamronize deeparg --analysis_software_version 1.0.2 --reference_database_version 1.0.2 \
    --input_file_name "$file_id" "$file" > ${dirDeepARG}/hAMRonization/${file_id}.txt ; 
done

# Generate Summuray report
hamronize summarize -o arg_report.tsv -t tsv \
    ${dirRGI}/hAMRonization/*_rgi.txt \
    ${dirDeepARG}/hAMRonization/*.txt

hamronize summarize -o arg_report.html -t interactive \
    ${dirRGI}/hAMRonization/*_rgi.txt \
    ${dirDeepARG}/hAMRonization/*.txt


# endles foor loop to keep prints on the screen

while true; do
    echo "Press [CTRL+C] to stop.."
    sleep 1
done
