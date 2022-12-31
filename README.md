# Ishigami, Wong et al. (2022): computational pipeline

**Reference:**
Yuma Ishigami\*, Mandy S. Wong\*, Carlos Martí-Gómez, Andalus Ayaz, Mahdi Kooshkbaghi, Sonya Hanson, David M. McCandlish, Adrian R. Krainer‡, Justin B. Kinney‡.  Specificity, synergy, and mechanisms of splice-modifying drugs. bioRxiv doi: https://doi.org/10.1101/2022.12.30.522303. (2022)
\*Equal contribution. 
‡Correspondence.

The computational pipeline used in this study was split into two parts: a "cluster" component that is designed to run on Elzar (the High Performance Compute Cluster  at Cold Spring Harbor Laboratory), and a "local" component suitable for execution on a standard laptop computer. 

<!-- Code for the cluster component is provided in the cluster/ directory. This code was configured specifically for the Elzar at CSHL and will likely have to be modified before it is run on a different cluster. To run this  pipeline, first download all Illumina sequencing data from SRA BioProject number PRJNA420342 and deposit into the directory cluster/data/illumina_runs/. These file names should match the file names listed in cluster/data/metadata.xlsx. The pipeline can then be run by executing

$ cd cluster

$ python2 run_pipeline.py

Note that this pipeline should be run under Python 2.7.11. 

Code for the local component is in the local/ directory. Unlike the cluster pipeline, this code should be run using Python >=3.6.3. First, copy the results from the cluster/saved directory to the local/from_pipeline directory. Next, download the human genome (both hg38.fa and hg19.fa) from the UCSC genome browser, decompress them, and place then in the data/ directory. Then do

$ cd local

$ python3 run_pipeline.py

This generates textual output (deposited in local/output/) as well as graphical output (deposited in local/plots). 
 -->
