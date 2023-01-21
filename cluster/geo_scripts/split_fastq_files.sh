#/bin/bash
#pip install --upgrade pip
#pip install --upgrade niet

# YAML file
YAML_FILE="smn2_info.yaml"
# Original Illumina directory
ILLUMINA_DIR=$(niet illumina_dir ${YAML_FILE})
IFS=" "

libs=(smn2_lib smn2_dmso smn2_nvs smn2_rg)

for lib in "${libs[@]}"
do
    LID=$(niet -f ifs ${lib}.lid ${YAML_FILE})
    SAMPLES=( $(niet -f ifs ${lib}.samples ${YAML_FILE}) )
    BARCODES=( $(niet -f ifs ${lib}.barcodes ${YAML_FILE}) )
    READ1=$ILLUMINA_DIR$(niet ${lib}.read1 ${YAML_FILE})
    READ2=$ILLUMINA_DIR$(niet ${lib}.read2 ${YAML_FILE})
    # Make the SAMPLES and BARCODES as an array
    declare -a SAMPLES
    declare -a BARCODES
    
    echo " "
    echo "Split the $lib files with LID=$LID"
    echo "----------------------------------"
    
    # Keep only the sequences and corresponding lines of
    # reads starts with sample specific barcodes
    for i in "${!SAMPLES[@]}"
    do
        echo "Make read1 fastq file for" "${SAMPLES[$i]}" "with barcode" "${BARCODES[$i]}"
        zcat $READ1|grep -A 2 -B 1 "^${BARCODES[$i]}" --no-group-separator|gzip > "${SAMPLES[$i]}".fastq.gz
        if [ -z "$READ2" ]
        then
            echo "Make read2 fastq file for" "${SAMPLES[$i]}" "with barcode" "${BARCODES[$i]}"
            zcat $READ2|grep -A 2 -B 1 "^${BARCODES[$i]}" --no-group-separator|gzip > "${SAMPLES[$i]}"_read2.fastq.gz
        fi
    done
done
