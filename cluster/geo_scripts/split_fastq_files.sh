#/bin/bash
pip install --upgrade pip
pip install --upgrade niet

# YAML file
YAML_FILE="smn2_info.yaml"
# Original Illumina directory
ILLUMINA_DIR=$(niet illumina_dir ${YAML_FILE})

LID=$(niet smn2_lib.lid ${YAML_FILE})
SAMPLES=$(niet smn2_lib.samples ${YAML_FILE})
BARCODES=$(niet smn2_lib.barcodes ${YAML_FILE})
READ1=$ILLUMINA_DIR$(niet smn2_lib.read1 ${YAML_FILE})
READ2=$ILLUMINA_DIR$(niet smn2_lib.read2 ${YAML_FILE})
