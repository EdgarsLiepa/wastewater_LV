# Latvian Waste water microbiome analysis

Microbiome and metagenome analysis of Latvian waste-water data.

In this study were analysed microbiome

-   Tax composition

-   Antibiotic resistance gene diversity

<!-- -->

-   AMR in Plasmid

-   AMR in MAGs

-   Compositional and Diversity,

-   r-code and input/output data

This Whole Metagenome Shotgun project has been deposited in ENA under the accession no. PRJEB79273 and metagenome assembled genome assemblies under the accession no. PRJEB80484.

## Pipeline

Input from

[*kraken-biom*](https://github.com/smdabdoub/kraken-biom)

[**mobilome-annotation-pipeline**](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline) [**miassembler**](https://github.com/EBI-Metagenomics/miassembler) [**genomes-generation**](https://github.com/EBI-Metagenomics/genomes-generation)

[**funcscan**](https://github.com/nf-core/funcscan)

Inputs:

-   Metadata: "City", "Sample", "Date", "Regional_Hospital", "Hospital_Type", "Population", "Industrial_wastewater_impact", "Industrial_wastewater_impact_from_food", "Dairy_farming", "Meat_production", "Metal_processing", "Washrooms"

-   RGI (*Resistance_genes/AMR_genes_RGI_scaffolds_filtered.csv*)

-   FuncScan (*results/hamronization_combined_report_fixed_names.tsv* and *results/hamronization_MAGs_report.tsv*)

-   MobileElements ( *results/all_amr_locations.tsv* )

-   Taxonomy data (*waste_water.biome*)

Notebooks:

-   [Publication_ARG.Rmd](Notebooks/Publication_ARG.Rmd)

-   [Publication.Rmd](Notebooks/Publication.Rmd)

-   [Publication_mobilome.RMD](Notebooks/Publication_mobilome.RMD)

## Results

plots.

Renders.
