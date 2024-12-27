#!/bin/bash

# ---
# --- Count reads with resistance qene
# ---

PROJECTDIR=~/NAS/bioinfo/edgars.liepa/Antibiotic_resistance/ww_publication/sequencingData/
RGIDIR=~/NAS/bioinfo/edgars.liepa/Antibiotic_resistance/ww_publication/ARG_RGI_scaffolds
MAPPEDFILEDIR=$PROJECTDIR/allignment/scaffolds/
ESCAPED_KEYWORD=$(printf '%s\n' "$MAPPEDFILEDIR" | sed -e 's/[]\/$*.^[]/\\&/g');


# echo parameters
echo "PROJECTDIR: $PROJECTDIR"
echo "RGIDIR: $RGIDIR"
echo "MAPPEDFILEDIR: $MAPPEDFILEDIR"

mkdir -p $RGIDIR/htseq_filtered

# IDS=("V300079975_L01_42" "V300079975_L01_16" "V300072773_L01_41" "V300072773_L01_44" "V300072773_L01_47" "V300072773_L01_125" "V300072773_L01_127" "V300072773_L01_67" "V300072773_L01_65" "V300072773_L01_68" "V300072773_L01_90" "V300079503_L02_63" "V300079975_L01_41" "V300079975_L01_70" "V300079503_L02_65" "V300079975_L01_45" "V300079503_L02_61" "V300079503_L02_67" "V300079975_L01_66" "V300079503_L02_58" "V300079503_L02_69" "V300079975_L01_68" "V300079975_L01_69" "V300079503_L02_62" "V300079503_L02_72" "V300079503_L02_64" "V300079503_L02_71" "V300079975_L01_71" "V300079975_L01_43" "V300079975_L01_44" "V300079975_L01_72" "V300079503_L02_57" "V300079503_L02_68" "V300079975_L01_46" "V300079503_L02_59" "V300079503_L02_70" "V300079975_L01_67" "V300079503_L02_60" "V300079503_L02_66" "V300079975_L01_65" "V300079975_L01_47" "V300079975_L01_13" "V300079975_L01_15" "V300079975_L01_48" "V300079975_L01_14")
IDS=("V300072773_L01_41")

for f in "${IDS[@]}"
do
    date
    echo "Start htseq-count: $f"
    OUT_File=$RGIDIR/htseq_filtered/${f}.csv
    # if out file exists then skip
    if [ -f "$OUT_File" ]; then
        echo "File $OUT_File exists. Skipping..."
        continue
    fi

    # PRINT command and files
    echo "htseq-count -r pos -m intersection-nonempty -t CDS --idattr=gene_id $MAPPEDFILEDIR/${f}_mapped_indexed_sorted.bam $RGIDIR/${f}_rgi_filtered.gtf -c $OUT_File"

    htseq-count -r pos -m intersection-nonempty -t CDS --idattr=gene_id $MAPPEDFILEDIR/${f}_mapped_indexed_sorted.bam $RGIDIR/${f}_rgi_filtered.gtf -c $OUT_File
done

date

