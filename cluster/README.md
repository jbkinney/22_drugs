# Cluster analysis

Directory contains scripts used process raw Illumina sequencing data (in the form of `.fastq.gz` files) from the MPSA experiments, thereby yielding the `.csv` in `../local/data/mpsa`. Note: these scripts are designed to run on Elzar (the High Performance Compute Cluster at Cold Spring Harbor Laboratory), and will likely have to be adapted to run in your environment.

To process the SMN2 data, execute

```
$ python run_pipeline.py analysis_smn2
$ python run_postprocess.py analysis_smn2
```

Note that final step in the postprocess script, is calculating the psi values. Since the psi is normalized by consensus sequences, the `make_psi.py` is case dependent and should be replaced in the `postprocess` directory.

Similarly, to process the ELP1 data, execute
```
$ python run_pipeline.py analysis_elp1
$ python run_postprocess.py analysis_elp1
```
The barcode-splice site libraries for the IKBKAP (ELP1) is borrowed from 

Wong, Mandy S., Justin B. Kinney, and Adrian R. Krainer. "Quantitative activity profile and context dependence of all human 5′ splice sites." Molecular cell 71.6 (2018): 1012-1026.
