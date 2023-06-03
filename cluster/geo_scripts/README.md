# Scripts for GEO [NOT COMPLETE]

There are two `yaml` files in this directory for each case:

```bash
ikbkap_info.yaml
smn2_info.yaml
```

The path to the illumina runs on Elzar is given on the top of these yaml files. The rest are the LID, Barcodes, Read names and Sample names. 

To read `yaml` files we will install the `niet` package in the beginning of each script.

## Split the Illumina FASTQ files:

There are two bash scripts which will read the `yaml` files and perform the task. The outputs of this step have the form `SAMPLE_NAME.fastq.gz`. If there is a second read, additionally we have `SAMPLE_NAME_read2.fatq.gz`.

To perform this task run

```bash
bash ikbkap_split_fastq_files.sh
bash smn2_split_fastq_files.sh
```

## Merge FASTQ files.

The files created in the previous step and potentially uploaded to the GEO server can be merged back with the following scripts:

```bash
bash ikbkap_merge_fastq_files.sh
bash smn2_merge_fastq_files.sh
```

