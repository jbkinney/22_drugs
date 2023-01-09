# Cluster analysis [UNDER CONSTRUCTION]

Directory contains scripts used process raw Illumina sequencing data (in the form of `.fastq.gz` files) from the MPSA experiments, thereby yielding the `.csv` in `../local/data/mpsa`. Note: these scripts are designed to run on Elzar (the High Performance Compute Cluster at Cold Spring Harbor Laboratory), and will likely have to be adapted to run in your environtment.

To process these data, execute
```
$ python run_pipeline analysis_smn2/
$ python run_pipeline analysis_elp1/
```
