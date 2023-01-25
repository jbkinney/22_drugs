# Alternative splicing changes induced by risdiplam and branaplam in HeLa cells

This section of the repository contains the [nextflow](https://www.nextflow.io/) pipeline used to analyze the alternative splicing changes induced by risdiplam and branaplam across the whole transcriptome using RNA-seq data to study their 5'splice site sequence dependent effect. 

It employs a general pipeline to count reads mapping to alternative splicing events based on [rMATS](https://rnaseq-mats.sourceforge.net) and performs Quality Control on the sequenced samples using [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [rseqc](https://pythonhosted.org/RSeQC/). Next, we perform more adhoc quantitative analysis to infer exon and splice site sequence level drug effects using [AS_quant](https://bitbucket.org/cmartiga/as_quant/src). 

### Dependencies

Create a new python environment and install the package manually for [AS_quant](https://bitbucket.org/cmartiga/as_quant/src)

```bash
conda create -n as_quant python=3.7
git clone git@bitbucket.org:cmartiga/as_quant.git
cd as_quant
python setup.py install
```

Create other environments for the other tools. It may take some time to solve all the dependencies

```bash
conda create -f environment.yml
```

### Download raw data

Download the data prior to running the analysis

```bash
cd fastq
bash download.sh
```

### Run the pipeline

The pipeline downloads genome references directly and performs all analysis automatically by just running

```bash
conda activate nextflow
nextflow pipeline.nf
```

Producing as main output files

- exons.psi.csv: file containing the estimated PSI and their posterior credible intervals in each of the experimental conditions for each exon with sufficient read coverage
- allelic_manifolds.csv: file containin the estimated effect of each drug for every 5'splice site sequence that was observed at least 10 times in the data

