#/bin/bash
pip install -q --upgrade pip
pip install -q --upgrade niet

# YAML file
YAML_FILE="smn2_info.yaml"
libs=(smn2_lib smn2_dmso smn2_nvs smn2_rg)
# Original Illumina directory
ILLUMINA_DIR=$(niet illumina_dir ${YAML_FILE})
IFS=" "


for lib in "${libs[@]}"
do
    LID=$(niet -f ifs ${lib}.lid ${YAML_FILE})
    SAMPLES=( $(niet -f ifs ${lib}.samples ${YAML_FILE}) )
    BARCODES=( $(niet -f ifs ${lib}.barcodes ${YAML_FILE}) )
    READ1=$(niet ${lib}.read1 ${YAML_FILE})
    READ2=$(niet ${lib}.read2 ${YAML_FILE})
    # Make the SAMPLES and BARCODES as an array
    declare -a SAMPLES
    declare -a BARCODES
    
    echo " "
    echo "Merge the fastq files for the LID=$LID"
    echo "----------------------------------"

    SAMPLES_FASTQ=( "${SAMPLES[@]/%/.fastq.gz}" )
    # Concatenate the FASTQ files in READ1
    pv "${SAMPLES_FASTQ[@]}" > $READ1

    if [ "$READ2" != "NONE" ]; then
        SAMPLES_FASTQ=( "${SAMPLES[@]/%/_read2.fastq.gz}" )
        pv "${SAMPLES_FASTQ[@]}" > $READ2
    fi

done
